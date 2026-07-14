"""
multimodel_etccdi.py

Core pipeline for the multi-model ETCCDI extremes analysis of G6-1.5K-SAI and
G6-1.5K-HiLLA (deliverables D1/D2). Extracted from
multimodel_etccdi_indices_edited_0627.ipynb so demo notebooks import it instead
of carrying the definitions themselves.

Contents: run configuration (Period, ComparisonConfig, windows), per-model
bucket layouts and native-variable converters (pr, tasmax, tasmin, tas, psl,
ua850), E3SM regrid/merge cache, member loaders, ETCCDI index computation with
xclim, the Lee et al. warming/SAI/combined framework, per-step netCDF caching,
land masking, and the summary plotting kept for analysis notebooks:
plot_all_mod_comparisons, plot_contrast_summary, contrast_table,
plot_medit_contrast_bars.

Written by Francis Osei Tutu Afrifa, 2026. Reflective internship, Track A.
"""

# Written by Francis Osei Tutu Afrifa, 2026

# Section 1: Imports and Global configuration

import s3fs
import fsspec
import xarray as xr
import numpy as np
import pandas as pd
import os

import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.util import add_cyclic_point

from dataclasses import dataclass, field
from pathlib import Path
from scipy import stats
import gc
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning, module='scipy')

import xclim
from xclim.core.calendar import percentile_doy
from xclim.core.units import convert_units_to
from xclim.indices import (days_over_precip_thresh, tx90p, tx10p, tn90p, tn10p,
                           warm_spell_duration_index, cold_spell_duration_index)


import uuid
import shutil
import tempfile
import itertools
import subprocess

import resource

def mem_gb():
    """
    Current max resident set size in GB 
    Use this to check current storage peak so far at 
    any current process or stage when called with the log_mem(tag)
    
    """
    # ru_maxrss is in kB on Linux
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

def log_mem(tag):
    print(f"    [mem] {tag}: peak {mem_gb():.1f} GB")

    
s3 = s3fs.S3FileSystem()
REFL_BUCKET = 's3://reflective-persistent-prod-large/'

OUTDIR = Path('./Figures')
OUTDIR.mkdir(exist_ok=True)
(OUTDIR / 'cached_data').mkdir(exist_ok=True)


# Local cache for regridded E3SM files. Persists across kernel restarts.
E3SM_CACHE_DIR = Path('./e3sm_cache')
E3SM_CACHE_DIR.mkdir(exist_ok=True)
E3SM_MAP_S3_PATH = ('E3SMv3/G6-1.5K-HiLLA/maps/'
                    'map_ne30pg2_to_cmip6_180x360_aave.20200201.nc')
E3SM_MAP_LOCAL = E3SM_CACHE_DIR / 'map_ne30pg2_to_cmip6_180x360_aave.20200201.nc'

# Bucket area for the slim merged regridded files (one per case+variable+period).
# Set to my username area on the persistent bucket. These persist
# across kernel restarts and hub sessions, so regridding is a one-time cost.
E3SM_MERGED_S3_PREFIX = ('s3://reflective-persistent-prod-large/'
                         'fafrifa/E3SMv3_regridded/')

# Dataset (D2) output tree, etccdi-indices format. Defaults to your writable
# area so a test writes nothing into Alistair's namespace. To merge into the
# shared zarr, repoint this at his tree (with his OK) or sync the files across.
ETCCDI_OUTPUT_ROOT = ('s3://reflective-persistent-prod-large/'
                      'fafrifa/ETCCDI/ETCCDI_indices_annual')

print(f"xclim version: {xclim.__version__}")


# ============================================================================
# SECTION 2  -  Period and ComparisonConfig dataclasses
# ============================================================================
# Period and ComparisonConfig dataclasses define time windows and comparison framings.
# These are the same across all models so we keep them at module level.

@dataclass
class Period:
    """
    A named time window with start and end years, both inclusive.
    
    """
    name: str
    start_year: int
    end_year: int
    
    def slice(self):
        """
        xarray time slice covering the full year range.
        This is calendar-agnostic: works for noleap, 360-day,
        and Gregorian without tripping over month lengths.
        
        """
        return slice(str(self.start_year), str(self.end_year))


@dataclass
class ComparisonConfig:
    """
    Define a scenario-vs-reference comparison.
    
    scenario_period: the window the SAI scenario is evaluated over
    reference_period: the window the SSP245 reference is evaluated over
    baseline_for_percentile: the window used to define percentile thresholds
                            (always SSP245, generally the same as reference_period
                             in the 'same_temperature' framing)
    All three Lee et al. (2026) comparisons reference the same baseline percentile 
    threshold so their index counts are measured on a common 
    yardstick and the anomalies add up.
    
    """
    scenario_name: str
    scenario_period: Period
    reference_name: str
    reference_period: Period
    baseline_for_percentile: Period
    label: str = ''    # tag for plot titles, e.g. "Global warming (1->2)"
    # default keeps it optional
    
    @property
    def description(self):
        return (f"{self.scenario_name}({self.scenario_period.start_year}-"
                f"{self.scenario_period.end_year}) vs "
                f"{self.reference_name}({self.reference_period.start_year}-"
                f"{self.reference_period.end_year})")


"""
Standard periods used by all variables. These follow Wang et al. 2026 
for the baseline (2020-2039 = 1.5C reference) and 
Duffey et al. 2026 for the assessment window (last 20 years of the run).
All SSP245 runs extend past 2084.

"""
BASELINE   = Period('baseline',   2020, 2039)
ASSESSMENT = Period('assessment', 2065, 2084)

# Per-model assessment window. CESM is 2050-2069 because its SSP245 members
# r7-r10 end in 2070 (Alistair: use the last 20 years all members share). The
# others keep the 2065-2084 default. Note for later: E3SM G6-1.5K-SAI tasmax
# ends ~2070, so any SAI-involving comparison for E3SM will need its own window
# when we add the E3SM SAI reader.
ASSESSMENT_BY_MODEL = {
    'CESM':  Period('assessment', 2050, 2069),
    'UKESM': Period('assessment', 2065, 2084),
    'MIROC': Period('assessment', 2065, 2084),
    'E3SM':  Period('assessment', 2065, 2084),
}

SAI_WINDOW = Period('assessment', 2050, 2069)   # common to all four models' SAI runs

def assessment_for(model):
    """
    Assessment window for a model, defaulting to the 2065-2084 convention.
    
    """
    return ASSESSMENT_BY_MODEL.get(model, ASSESSMENT)


# ============================================================================
# SECTION 3  -  Three-way comparison configuration (Lee et al. 2026)
# ============================================================================
def make_comparison_config(comparison_type='combined', assessment=ASSESSMENT):
    """
    Return a ComparisonConfig for one of the three Lee et al. (2026) comparisons.
 
    Accepts the new names plus backward-compatible aliases so older call sites
    keep working:
      'warming'  / '1to2'
      'sai'      / 'same_time'        / '2to3'
      'combined' / 'same_temperature' / '1to3'
      
    """
    ct = comparison_type.lower()
 
    if ct in ('warming', '1to2'):
        # 1 -> 2: pure global warming. SSP245 on both sides.
        return ComparisonConfig(
            scenario_name='SSP245', scenario_period=assessment,
            reference_name='SSP245', reference_period=BASELINE,
            baseline_for_percentile=BASELINE,
            label='Global warming (1\u21922)',
        )
    elif ct in ('sai', 'same_time', '2to3'):
        # 2 -> 3: pure SAI effect at fixed date. Reference is the warmed world.
        return ComparisonConfig(
            scenario_name='HiLLA', scenario_period=assessment,
            reference_name='SSP245', reference_period=assessment,
            baseline_for_percentile=BASELINE,
            label='SAI (2\u21923)',
        )
    elif ct in ('combined', 'same_temperature', '1to3'):
        # 1 -> 3: warming + SAI. Residual at matched GMST against the baseline.
        return ComparisonConfig(
            scenario_name='HiLLA', scenario_period=assessment,
            reference_name='SSP245', reference_period=BASELINE,
            baseline_for_percentile=BASELINE,
            label='Warming + SAI (1\u21923)',
        )

    elif ct in ('sai_vs_hilla', 'injection_latitude'):
        # G6-1.5K-SAI vs G6-1.5K-HiLLA at matched time and target. Isolates the
        # injection-latitude effect on the index, both sides scored against the
        # same SSP245 baseline percentile. Anomaly = SAI - HiLLA: positive where
        # standard injection gives more of the index than high-latitude
        # injection. Swap the two names for the HiLLA - SAI sign.
        return ComparisonConfig(
            scenario_name='SAI', scenario_period=assessment,
            reference_name='HiLLA', reference_period=assessment,
            baseline_for_percentile=BASELINE,
            label='SAI vs HiLLA (injection latitude)',
        )
        
    else:
        raise ValueError(f"Unknown comparison_type: {comparison_type}")
 
 
# ============================================================================
# SECTION 4  -  Per-model bucket layout dispatch table
# ============================================================================
# Each (model, scenario) entry carries:
#   path_template : format string with {member}, {var}, {stream} placeholders
#   members       : tuple of member identifiers for that scenario
#   file_pattern  : '' if the variable lives in the path (UKESM/E3SM),
#                   '.{var}.' if the filename carries the var (CESM),
#                   '{var}_..._{member}.nc' if it carries var and member (MIROC)
#   streams       : optional dict mapping variable -> subdir (UKESM HiLLA only)
#   custom_reader : 'e3sm' to route to the regridding reader
MODEL_BUCKET_LAYOUT = {
    # ---- CESM2-WACCM ---------------------------------------------------------
    # HiLLA r1/r2/r3 under ADAY/ (note lowercase 'k' in G6-1.5k-HiLLA).
    # SSP245 r6-r10 under day/ (the clean members uploaded from Derecho).
    ('CESM', 'HiLLA'): {
        'path_template': 'CESM2-WACCM/G6-1.5k-HiLLA/{member}/ADAY',
        'members': ('r1', 'r2', 'r3'),
        'file_pattern': '.{var}.',
        'streams': None,
    },
    
    ('CESM', 'SAI'): {
        'path_template': 'CESM2-WACCM/G6-1.5k-SAI/{member}/ADAY',
        'members': ('r1', 'r2', 'r3'),
        'file_pattern': '.{var}.',
        'streams': None,
    },
    
    ('CESM', 'SSP245'): {
        'path_template': 'CESM2-WACCM/SSP2-4.5/{member}/day',
        'members': ('r6', 'r7', 'r8', 'r9', 'r10'),
        'file_pattern': '.{var}.',
        'streams': None,
    },
 
    # ---- UKESM1-1 ------------------------------------------------------------
    # HiLLA (physics f2) splits variables across streams.
    # SSP245 (physics f1) keeps the standard ap6/day/{var} layout; r2 excluded
    # because it was uploaded non-CMORized (precipitation_flux / air_temperature_*).
    ('UKESM', 'HiLLA'): {
        'path_template': 'UKESM1-1/G6-1p5K-HiLLA/{member}/{stream}/{var}',
        'members': ('r12i1p1f2', 'r2i1p1f2', 'r3i1p1f2'),
        'file_pattern': '',
        'streams': {
            'pr': 'apj/day',
            'tasmax': 'apd/day',
            'tasmin': 'apd/day',
        },
    },
    ('UKESM', 'SSP245'): {
        'path_template': 'UKESM1-1/SSP245/{member}/ap6/day/{var}',
        'members': ('r12i1p1f1', 'r3i1p1f1'),
        'file_pattern': '',
        'streams': None,
    },
 
    # ---- MIROC-ES2H ----------------------------------------------------------
    # Both scenarios in one flat directory, split by filename.
    ('MIROC', 'HiLLA'): {
        'path_template': 'MIROC-ES2H/G6-1.5K-HiLLA/day',
        'members': ('r01', 'r02', 'r03'),
        'file_pattern': '{var}_G6-1.5K-HiLLA_{member}.nc',
        'streams': None,
    },
    ('MIROC', 'SSP245'): {
        'path_template': 'MIROC-ES2H/G6-1.5K-HiLLA/day',
        'members': ('r01', 'r02', 'r03'),
        'file_pattern': '{var}_baseline_{member}.nc',
        'streams': None,
    },

    ('MIROC', 'SAI'): {
        'path_template': 'MIROC-ES2H/G6-1.5K-SAI/day',
        'members': ('r01', 'r02', 'r03'),
        'file_pattern': '{var}_G6-1.5K-SAI_{member}.nc',
        'streams': None,
    },
 
    # ---- E3SMv3 --------------------------------------------------------------
    # Both scenarios live under E3SMv3/G6-1.5K-HiLLA/. Native cubed-sphere grid;
    # the custom reader downloads, regrids with ncremap, and caches locally.
    ('E3SM', 'HiLLA'): {
        'path_template': '',
        'members': ('v3.LR.ssp245.g6_hilla.sai.0101',
                    'v3.LR.ssp245.g6_hilla.sai.0151',
                    'v3.LR.ssp245.g6_hilla.sai.0201'),
        'file_pattern': '',
        'streams': None,
        'custom_reader': 'e3sm',
    },
    ('E3SM', 'SSP245'): {
        'path_template': '',
        'members': ('v3.LR.ssp245_0101',
                    'v3.LR.ssp245_0151',
                    'v3.LR.ssp245_0201'),
        'file_pattern': '',
        'streams': None,
        'custom_reader': 'e3sm',
    },

    ('UKESM', 'SAI'): {
        'path_template': '', 'members': ('001', '002', '003'),
        'file_pattern': '', 'streams': None, 'custom_reader': 'ukesm_sai',
    },
    ('E3SM', 'SAI'): {
        'path_template': '', 'members': ('001', '002', '003'),
        'file_pattern': '', 'streams': None, 'custom_reader': 'e3sm_sai',
    },
}

def get_layout(model, scenario):
    """
    Return the layout dict for a (model, scenario) pair, or raise.
    
    """
    key = (model, scenario)
    if key not in MODEL_BUCKET_LAYOUT:
        raise KeyError(f"No bucket layout for ({model}, {scenario}). "
                       f"Available: {list(MODEL_BUCKET_LAYOUT.keys())}")
    return MODEL_BUCKET_LAYOUT[key]
 

# Labels that make my written paths match the etccdi-indices (GitHub) tree exactly.
MODEL_LABEL = {
    'CESM': 'CESM2-WACCM', 'UKESM': 'UKESM1-1',
    'MIROC': 'MIROC-ES2H', 'E3SM': 'E3SMv3',
}
SCENARIO_LABEL = {
    ('CESM', 'HiLLA'): 'G6-1.5K-HiLLA',  ('CESM', 'SAI'): 'G6-1.5K-SAI', ('CESM', 'SSP245'): 'SSP2-4.5',
    ('UKESM', 'HiLLA'): 'G6-1.5K-HiLLA', ('UKESM', 'SSP245'): 'SSP245',
    ('MIROC', 'HiLLA'): 'G6-1.5K-HiLLA', ('MIROC', 'SSP245'): 'baseline',
    ('MIROC', 'SAI'): 'G6-1.5K-SAI', 
    ('E3SM', 'HiLLA'): 'G6-1.5K-HiLLA',  ('E3SM', 'SSP245'): 'SSP245',
    ('UKESM', 'SAI'): 'G6-1.5K-SAI',   ('E3SM', 'SAI'): 'G6-1.5K-SAI',   # SCENARIO_LABEL
}
# layout member id -> the short tree label on GitHub
MEMBER_LABEL = {
    ('CESM', 'HiLLA'):   {'r1': 'r1', 'r2': 'r2', 'r3': 'r3'},
    ('CESM', 'SAI'): {'r1': 'r1', 'r2': 'r2', 'r3': 'r3'},
    ('CESM', 'SSP245'):  {'r6': 'r6', 'r7': 'r7', 'r8': 'r8', 'r9': 'r9', 'r10': 'r10'},
    ('UKESM', 'HiLLA'):  {'r12i1p1f2': 'r1', 'r2i1p1f2': 'r2', 'r3i1p1f2': 'r3'},
    ('UKESM', 'SSP245'): {'r12i1p1f1': 'r1', 'r3i1p1f1': 'r3'},   # provisional, see note
    ('MIROC', 'HiLLA'):  {'r01': 'r01', 'r02': 'r02', 'r03': 'r03'},
    ('MIROC', 'SSP245'): {'r01': 'r01', 'r02': 'r02', 'r03': 'r03'},
    ('MIROC', 'SAI'): {'r01': 'r01', 'r02': 'r02', 'r03': 'r03'},
    ('E3SM', 'HiLLA'):   {'v3.LR.ssp245.g6_hilla.sai.0101': 'r1',
                          'v3.LR.ssp245.g6_hilla.sai.0151': 'r2',
                          'v3.LR.ssp245.g6_hilla.sai.0201': 'r3'},
    ('E3SM', 'SSP245'):  {'v3.LR.ssp245_0101': 'r1',
                          'v3.LR.ssp245_0151': 'r2',
                          'v3.LR.ssp245_0201': 'r3'},
    ('UKESM', 'SAI'): {'001': 'r1', '002': 'r2', '003': 'r3'},           # MEMBER_LABEL
    ('E3SM',  'SAI'): {'001': 'r1', '002': 'r2', '003': 'r3'},           # MEMBER_LABEL
}


# ============================================================================
# SECTION 5  -  Unit converters and per-model variable metadata
# ============================================================================
def _convert_cesm_PRECT(ds):
    """
    CESM/E3SM PRECT is m/s. Multiply by 86400 (s/day) and 1000 (m->mm).
    
    """
    pr = (ds['PRECT'] * 86400 * 1000).astype('float32')
    pr.attrs['units'] = 'mm/day'
    pr.attrs['standard_name'] = 'precipitation_flux'
    return pr
 

def _convert_cmip_pr(ds, varname='pr'):
    """
    CMIP pr is kg m-2 s-1 == mm/s. Multiply by 86400 to get mm/day.
    Idempotent: if a custom reader already produced mm/day, pass it through.
    
    """
    da = ds[varname]
    if da.attrs.get('units') == 'mm/day':
        out = da.astype('float32')
        out.attrs['units'] = 'mm/day'
        out.attrs['standard_name'] = 'precipitation_flux'
        return out
    pr = (da * 86400).astype('float32')
    pr.attrs['units'] = 'mm/day'
    pr.attrs['standard_name'] = 'precipitation_flux'
    return pr
    
 
def _convert_cesm_TREFHTMX(ds):
    da = ds['TREFHTMX'].astype('float32')
    da.attrs['units'] = 'K'
    return da


def _convert_cesm_TREFHTMN(ds):
    da = ds['TREFHTMN'].astype('float32')
    da.attrs['units'] = 'K'
    return da
 
 
def _convert_cmip_tasmax(ds, varname='tasmax'):
    da = ds[varname].astype('float32')
    if 'units' not in da.attrs:
        da.attrs['units'] = 'K'
    return da
 
 
def _convert_cmip_tasmin(ds, varname='tasmin'):
    da = ds[varname].astype('float32')
    if 'units' not in da.attrs:
        da.attrs['units'] = 'K'
    return da


# E3SM uses CESM-lineage names (PRECT, TREFHTMX, TREFHTMN), so it reuses the
# CESM converters. UKESM and MIROC use CMOR-standard names.
VARIABLE_INFO_BY_MODEL = {
    'CESM': {
        'pr':     {'native_short': 'PRECT',    'cmor_short': 'pr',     'converter': _convert_cesm_PRECT},
        'tasmax': {'native_short': 'TREFHTMX', 'cmor_short': 'tasmax', 'converter': _convert_cesm_TREFHTMX},
        'tasmin': {'native_short': 'TREFHTMN', 'cmor_short': 'tasmin', 'converter': _convert_cesm_TREFHTMN},
    },
    'UKESM': {
        'pr':     {'native_short': 'pr',     'cmor_short': 'pr',     'converter': _convert_cmip_pr},
        'tasmax': {'native_short': 'tasmax', 'cmor_short': 'tasmax', 'converter': _convert_cmip_tasmax},
        'tasmin': {'native_short': 'tasmin', 'cmor_short': 'tasmin', 'converter': _convert_cmip_tasmin},
    },
    'MIROC': {
        'pr':     {'native_short': 'pr',     'cmor_short': 'pr',     'converter': _convert_cmip_pr},
        'tasmax': {'native_short': 'tasmax', 'cmor_short': 'tasmax', 'converter': _convert_cmip_tasmax},
        'tasmin': {'native_short': 'tasmin', 'cmor_short': 'tasmin', 'converter': _convert_cmip_tasmin},
    },
    'E3SM': {
        'pr':     {'native_short': 'PRECT',    'cmor_short': 'pr',     'converter': _convert_cesm_PRECT},
        'tasmax': {'native_short': 'TREFHTMX', 'cmor_short': 'tasmax', 'converter': _convert_cesm_TREFHTMX},
        'tasmin': {'native_short': 'TREFHTMN', 'cmor_short': 'tasmin', 'converter': _convert_cesm_TREFHTMN},
    },
}
 
 
def get_variable_info(model, variable):
    if model not in VARIABLE_INFO_BY_MODEL:
        raise KeyError(f"No variable info for model {model}")
    if variable not in VARIABLE_INFO_BY_MODEL[model]:
        raise KeyError(f"No variable info for {model}/{variable}")
    return VARIABLE_INFO_BY_MODEL[model][variable]
# ============================================================================
# SECTION 6  -  E3SM regridding machinery
# ============================================================================
# ---- SAI sources on Alistair's bucket (different layouts per model) ----
UKESM_SAI_S3_BASE = 's3://reflective-persistent-prod/alistairduffey/UKESM1-1/G6-1.5K-SAI/'
UKESM_SAI_VARMAP  = {'tasmax': 'temp', 'tasmin': 'temp_1', 'pr': 'precip'}

E3SM_SAI_S3_BASE = 's3://reflective-persistent-prod/alistairduffey/E3SMv3/G6-1.5K-SAI/'
E3SM_SAI_VARMAP  = {'tasmax': ('daily_T_max', 'TREFHTMX'),
                    'tasmin': ('daily_T_min', 'TREFHTMN'),
                    'pr':     ('daily_pr',    'PRECT')}

def _to_mm_day_repo(da, varname=None):
    """
    The etccdi-indices to_mm_day, so our SAI precip matches the values on GitHub.
    
    """
    units = (da.attrs.get('units') or '').strip().lower().replace(' ', '')
    is_kg = units in {'kgm-2s-1', 'kg/m^2/s', 'kg/m2/s', 'kgm^-2s^-1', 'kgm**-2s**-1'}
    is_ms = units in {'ms-1', 'm/s', 'm s-1', 'm/s-1'} or units == 'ms^-1'
    if is_kg:
        out = da * 86400
    elif is_ms:
        out = da * 86400 * 1000
    else:
        out = da * 86400 if (varname and varname.upper() in {'PRECT', 'PRECC', 'PRECL', 'PR'}) else da
    out = out.assign_attrs(da.attrs)
    out.attrs['units'] = 'mm/day'
    return out.astype('float32')


def open_netcdf_any(filepath, **kwargs):
    """
    Open a local netCDF trying h5netcdf then netcdf4.
 
    Compressed ncremap output (-7 -L 1) is HDF5 (h5netcdf reads it natively);
    uncompressed output is netCDF-5/CDF5 (only netcdf4 reads it). Trying both
    makes the reader robust regardless of the compression choice.
    
    """
    last_err = None
    for engine in ('h5netcdf', 'netcdf4'):
        try:
            return xr.open_dataset(filepath, engine=engine, **kwargs)
        except Exception as e:  # noqa: BLE001 - we want to try the next engine
            last_err = e
    raise last_err
 
 
def is_valid_netcdf_file(filepath):
    """
    True if the file exists, is non-trivial in size, and opens cleanly.
    
    """
    p = Path(filepath)
    if not p.exists() or p.stat().st_size < 1000:
        return False
    try:
        ds = open_netcdf_any(filepath)
        ds.close()
        return True
    except Exception:  # noqa: BLE001
        return False
 
 
def ensure_e3sm_map_file():
    """
    Download the ne30pg2 -> 180x360 mapping file from the bucket if absent.
    
    """
    if E3SM_MAP_LOCAL.exists():
        return str(E3SM_MAP_LOCAL)
    print("  Downloading E3SM map file...")
    src_uri = 's3://' + REFL_BUCKET.replace('s3://', '').rstrip('/') + '/' + E3SM_MAP_S3_PATH
    with fsspec.open(src_uri, mode='rb') as src, open(E3SM_MAP_LOCAL, 'wb') as dst:
        shutil.copyfileobj(src, dst)
    print(f"  Saved map file ({E3SM_MAP_LOCAL.stat().st_size / 1024 / 1024:.1f} MB)")
    return str(E3SM_MAP_LOCAL)
 

def run_ncremap(input_path, output_path, map_file):
    """
    Regrid one file with ncremap, single-threaded, with isolated scratch.
 
    Two hard-won settings:
      OMP_NUM_THREADS=1  -- the multi-threaded ncremap path segfaults on this
        hub (signal 11 inside libucs/UCX). Single-threaded avoids it.
      TMPDIR=<unique dir> -- ncremap writes intermediate temp files; if a prior
        call crashed it can leave a temp that makes the next call fail with
        'Permission denied'. A unique per-call scratch dir, removed afterwards,
        prevents one crash from poisoning later calls.
    The conda env bin is prepended to PATH so the call survives a kernel restart
    that dropped the interactive PATH edit. Output is verified non-empty.
    
    """
    output_path = Path(output_path)
    if output_path.exists():
        try:
            output_path.unlink()
        except OSError:
            pass
 
    scratch = output_path.parent / f".ncremap_tmp_{uuid.uuid4().hex[:8]}"
    scratch.mkdir(parents=True, exist_ok=True)
 
    env = os.environ.copy()
    conda_bin = '/srv/conda/envs/notebook/bin'
    if conda_bin not in env.get('PATH', ''):
        env['PATH'] = f"{conda_bin}:{env.get('PATH', '')}"
    env['TMPDIR'] = str(scratch)
    env['OMP_NUM_THREADS'] = '1'
 
    cmd = ['ncremap', '-i', str(input_path), '-m', str(map_file),
           '-o', str(output_path), '-7', '-L', '1', '--no_stdin']
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
        ok = output_path.exists() and output_path.stat().st_size > 1000
        if not ok:
            print(f"    ncremap produced no/empty output for {Path(input_path).name}")
        return ok
    except subprocess.CalledProcessError as e:
        print(f"    ncremap failed: {e.stderr[:300]}")
        return False
    except FileNotFoundError:
        raise RuntimeError("ncremap not found. Install: "
                           "conda install -c conda-forge nco -y")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
 
 
def _e3sm_merged_s3_key(case_name, var, period):
    """
    REFL_Bucket URI for the slim merged regridded file for one case+var+period.
    
    """
    return (E3SM_MERGED_S3_PREFIX.rstrip('/') +
            f"/{case_name}/merged_{var}_{period.start_year}_{period.end_year}.nc")
 
 
def read_e3sm_var(case_name, var, period, scenario_dir='G6-1.5K-HiLLA'):
    """
    Read one E3SM case+variable, backed by a SINGLE merged file on the bucket.
 
    Strategy (keeps local disk use to a few hundred MB no matter how much E3SM
    you process, and makes regridding a one-time cost that persists on S3):
 
      1. If the merged file exists on the bucket, stream and return it.
      2. Otherwise regrid each in-window 5-year chunk locally, one at a time:
         download raw -> ncremap -> keep ONLY the data variable -> slice to the
         window -> load -> delete both the raw and regridded chunk immediately.
      3. Concatenate the chunks in memory, write one slim compressed file to a
         local temp, upload it to the bucket, delete the local temp.
 
    The merged file holds only the data variable (no bounds/area/gw), so it is
    far smaller than the per-chunk files and there is one per case+var+period.
 
    Path layout of the source:
      E3SMv3/{scenario_dir}/{case}/day/{var}/gn/{YYYYMMDD}/{var}_YYYYMM_YYYYMM.nc
    The {YYYYMMDD} version subdir varies; the most recent is used.
    
    """
    merged_uri = _e3sm_merged_s3_key(case_name, var, period)
 
    # 1. Fast path: merged file already on the bucket.
    try:
        if s3.exists(merged_uri):
            print(f"  {case_name}/{var} {period.start_year}-{period.end_year}: "
                  f"streaming merged file from bucket")
            with fsspec.open(merged_uri, mode='rb') as fobj:
                return open_netcdf_any(fobj).load()
    except Exception as e:  # noqa: BLE001 - fall through to rebuild on any error
        print(f"  (bucket check/stream failed: {type(e).__name__}; rebuilding)")
 
    # Locate source chunks.
    gn_path = f'E3SMv3/{scenario_dir}/{case_name}/day/{var}/gn'
    full_gn = REFL_BUCKET + gn_path
    try:
        date_dirs = s3.ls(full_gn)
    except FileNotFoundError:
        print(f"  E3SM gn dir not found: {full_gn}")
        return None
    if not date_dirs:
        return None
    chosen = sorted(date_dirs, reverse=True)[0]
    if len(date_dirs) > 1:
        print(f"  Versions {[d.split('/')[-1] for d in date_dirs]}, "
              f"using {chosen.split('/')[-1]}")
    nc_files = sorted(f for f in s3.ls(chosen) if f.endswith('.nc'))
    if not nc_files:
        return None
 
    # Filter to chunks intersecting the window.
    in_window = []
    for f in nc_files:
        parts = f.split('/')[-1].replace('.nc', '').split('_')
        try:
            fs, fe = int(parts[-2][:4]), int(parts[-1][:4])
            if fe < period.start_year or fs > period.end_year:
                continue
        except (ValueError, IndexError):
            pass
        in_window.append(f)
    print(f"  {case_name}/{var} {period.start_year}-{period.end_year}: "
          f"regrid+merge {len(in_window)} chunks")
 
    map_file = ensure_e3sm_map_file()
    tmp_dir = E3SM_CACHE_DIR / f".chunks_{uuid.uuid4().hex[:8]}"
    tmp_dir.mkdir(parents=True, exist_ok=True)
 
    pieces = []
    time_span = (str(period.start_year), str(period.end_year))
    try:
        for s3_file in in_window:
            filename = s3_file.split('/')[-1]
            raw = tmp_dir / filename
            reg = tmp_dir / f"regridded_{filename}"
            print(f"    {filename}")
            with fsspec.open('s3://' + s3_file, mode='rb') as src, \
                 open(raw, 'wb') as dst:
                shutil.copyfileobj(src, dst)
            ok = run_ncremap(raw, reg, map_file)
            try:
                raw.unlink()
            except OSError:
                pass
            if not ok:
                print(f"      regrid failed, skipping {filename}")
                continue
            # Keep only the data variable; slice to the window; load eagerly.
            piece = open_netcdf_any(reg)
            keep = [v for v in piece.data_vars if v == var]
            piece = (piece[keep] if keep
                     else piece[[list(piece.data_vars)[0]]])
            piece = piece.sel(time=slice(*time_span)).load()
            pieces.append(piece)
            try:
                reg.unlink()
            except OSError:
                pass
 
        if not pieces:
            return None
        combined = xr.concat(pieces, dim='time')
        _, uniq = np.unique(combined.time.values, return_index=True)
        combined = combined.isel(time=sorted(uniq))
 
        # Write one slim compressed file locally, upload to the Reflective bucket, delete local.
        enc = {v: {'zlib': True, 'complevel': 1} for v in combined.data_vars}
        local_merged = tmp_dir / f"merged_{var}_{period.start_year}_{period.end_year}.nc"
        combined.to_netcdf(local_merged, encoding=enc)
        size_mb = local_merged.stat().st_size / 1024 / 1024
        print(f"    uploading merged {local_merged.name} ({size_mb:.0f} MB) -> bucket")
        with open(local_merged, 'rb') as src, fsspec.open(merged_uri, mode='wb') as dst:
            shutil.copyfileobj(src, dst)
        return combined
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
 

def read_ukesm_sai_var(member_num, variable, period):
    """
    UKESM1-1 G6-1.5K-SAI: Alistair's scipy single-file copies, one file per
    member holding temp (tasmax), temp_1 (tasmin), precip (pr, m/s) on dims t/ht
    or t/surface with latitude/longitude coords.

    Opened lazily with dask chunks so only the requested variable and window are
    read (the files hold both temperatures over the full run; loading that eagerly
    overruns memory). Normalised to a standard UKESM Dataset: t->time,
    latitude/longitude->lat/lon cast float64 to match the CMORized HiLLA grid,
    singleton squeezed, precip->mm/day, CMOR variable name.
    
    """
    src = UKESM_SAI_VARMAP[variable]
    sub = 'daily_pr' if variable == 'pr' else 'daily_Tmaxmin'
    prefix = 'PRECT' if variable == 'pr' else 'T'
    path = f'{UKESM_SAI_S3_BASE}{sub}/{prefix}_{member_num}.nc'
    span = (str(period.start_year), str(period.end_year))
    try:
        with fsspec.open(path, mode='rb') as fobj:
            ds = xr.open_dataset(fobj, engine='scipy', chunks={'t': 360})
            rename = {'t': 'time'}
            if 'latitude' in ds.coords or 'latitude' in ds.dims:
                rename['latitude'] = 'lat'
            if 'longitude' in ds.coords or 'longitude' in ds.dims:
                rename['longitude'] = 'lon'
            ds = ds.rename(rename)
            for singleton in ('ht', 'surface'):
                if singleton in ds.dims:
                    ds = ds.squeeze(singleton, drop=True)
            da = ds[src].sel(time=slice(*span)).load()   # only this var + window
    except FileNotFoundError:
        print(f"  UKESM SAI not found: {path}")
        return None

    if len(da.time) == 0:
        return None
    da = da.assign_coords(lat=da['lat'].astype('float64'),
                          lon=da['lon'].astype('float64'))
    if variable == 'pr':
        return _to_mm_day_repo(da, 'PRECT').rename('pr').to_dataset()
    da = da.rename(variable)
    da.attrs['units'] = 'K'
    return da.to_dataset()
    

def read_e3sm_sai_var(member_num, variable, period):
    """
    E3SMv3 G6-1.5K-SAI: Alistair's PRE-regridded 180x360 CDF-5 files (no
    ncremap). fsspec cannot stream CDF-5, so each chunk is downloaded, opened
    with netcdf4, sliced, then deleted; a slim merged file is cached on the
    bucket so later runs stream it. Returns a Dataset with the native E3SM name
    (the E3SM converter handles m/s -> mm/day).
    
    """
    sub, native = E3SM_SAI_VARMAP[variable]
    merged_uri = (E3SM_MERGED_S3_PREFIX.rstrip('/') +
                  f"/G6-1.5K-SAI_{member_num}/merged_{native}_"
                  f"{period.start_year}_{period.end_year}.nc")
    if s3.exists(merged_uri.replace('s3://', '')):
        with fsspec.open(merged_uri, mode='rb') as f:
            return xr.open_dataset(f, engine='h5netcdf').load()

    try:
        files = sorted(f for f in s3.ls(f'{E3SM_SAI_S3_BASE}{sub}/{member_num}')
                       if f.endswith('.nc'))
    except FileNotFoundError:
        print(f"  E3SM SAI dir not found: {E3SM_SAI_S3_BASE}{sub}/{member_num}")
        return None
    if not files:
        return None

    tmp_dir = E3SM_CACHE_DIR / f".sai_{uuid.uuid4().hex[:8]}"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    span = (str(period.start_year), str(period.end_year))
    pieces = []
    try:
        for s3_file in files:
            local = tmp_dir / s3_file.split('/')[-1]
            with fsspec.open('s3://' + s3_file, mode='rb') as src, open(local, 'wb') as dst:
                shutil.copyfileobj(src, dst)
            piece = open_netcdf_any(local)
            keep = [v for v in piece.data_vars if v == native] or [list(piece.data_vars)[0]]
            piece = piece[keep].sel(time=slice(*span)).load()
            if len(piece.time):
                pieces.append(piece)
            local.unlink(missing_ok=True)
        if not pieces:
            return None
        combined = xr.concat(pieces, dim='time')
        _, uniq = np.unique(combined.time.values, return_index=True)
        combined = combined.isel(time=sorted(uniq))
        enc = {v: {'zlib': True, 'complevel': 1} for v in combined.data_vars}
        local_merged = tmp_dir / f"merged_{native}_{period.start_year}_{period.end_year}.nc"
        combined.to_netcdf(local_merged, encoding=enc)
        with open(local_merged, 'rb') as src, fsspec.open(merged_uri, mode='wb') as dst:
            shutil.copyfileobj(src, dst)
        return combined
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        
    
# E3SM cases by scenario, for the standalone prebuild.
E3SM_CASES = {
    'SSP245': ('v3.LR.ssp245_0101', 'v3.LR.ssp245_0151', 'v3.LR.ssp245_0201'),
    'HiLLA':  ('v3.LR.ssp245.g6_hilla.sai.0101',
               'v3.LR.ssp245.g6_hilla.sai.0151',
               'v3.LR.ssp245.g6_hilla.sai.0201'),
}
 
 
def prebuild_e3sm_cache(variables=('PRECT', 'TREFHTMX', 'TREFHTMN'),
                        periods=(BASELINE, ASSESSMENT)):
    """
    Regrid and upload every merged E3SM file the analysis will need.
 
    Run this ONCE, ideally in its own kernel, before the analysis. It populates
    the bucket with one slim merged file per case+variable+period. Because the
    merged files live on S3, this is a one-time cost: later runs (this session or
    any future one) stream them and never call ncremap. HiLLA has no data before
    2035, so the baseline window is skipped for HiLLA. Already-present merged
    files are detected and skipped, so the prebuild is itself resumable.
    
    """
    for scenario, cases in E3SM_CASES.items():
        for case in cases:
            for var in variables:
                for period in periods:
                    if scenario == 'HiLLA' and period.start_year < 2035:
                        continue
                    print(f"\n### prebuild {scenario} {case} {var} "
                          f"{period.start_year}-{period.end_year}")
                    ds = read_e3sm_var(case, var, period)
                    if ds is not None:
                        print(f"   OK: {dict(ds.sizes)}")
                        del ds
                        gc.collect()
                    else:
                        print("   (no data returned)")
 
 
def list_e3sm_merged_on_bucket():
    """
    List the merged E3SM files currently on the bucket (for inspection).
    
    """
    try:
        files = s3.ls(E3SM_MERGED_S3_PREFIX.replace('s3://', '').rstrip('/'))
        for f in files:
            print(f"  {f}")
        return files
    except FileNotFoundError:
        print(f"  Nothing under {E3SM_MERGED_S3_PREFIX} yet.")
        return []


# ============================================================================
# SECTION 7  -  Generic per-variable reader (dispatches E3SM to its reader)
# ============================================================================
def read_model_var(model, scenario, member, variable, period):
    """
    Read one variable for one (model, scenario, member) over a time window.
 
    For E3SM the custom regridding reader is used. For the bucket-native models
    the file list is built from the layout's path_template and file_pattern,
    each chunk is sliced to the window, and overlapping boundaries are
    deduplicated after concatenation. Returns a Dataset (caller converts units),
    or None if nothing matched
    
    """
    layout = get_layout(model, scenario)
    var_info = get_variable_info(model, variable)
    native_var = var_info['native_short']
    
    cr = layout.get('custom_reader')
    if cr == 'e3sm':
        return read_e3sm_var(case_name=member, var=native_var,
                             period=period, scenario_dir='G6-1.5K-HiLLA')
    if cr == 'ukesm_sai':
        return read_ukesm_sai_var(member_num=member, variable=variable, period=period)
    if cr == 'e3sm_sai':
        return read_e3sm_sai_var(member_num=member, variable=variable, period=period)

        
    # Build the directory path (substitute stream for UKESM HiLLA).
    if layout.get('streams'):
        if variable not in layout['streams']:
            raise KeyError(f"No stream configured for {model}/{scenario}/{variable}.")
        stream = layout['streams'][variable]
        path = layout['path_template'].format(member=member, var=native_var, stream=stream)
    else:
        path = layout['path_template'].format(member=member, var=native_var)
    full_path = REFL_BUCKET + path
 
    try:
        all_files = s3.ls(full_path)
    except FileNotFoundError:
        print(f"  Path not found: {full_path}")
        return None
 
    # Select the variable's files. Three filename conventions are handled.
    overrides = layout.get('file_pattern_overrides', {})
    if variable in overrides:
        match = overrides[variable]
        files = sorted(f for f in all_files if match in f and f.endswith('.nc'))
    else:
        fp = layout['file_pattern']
        if fp == '':
            files = sorted(f for f in all_files if f.endswith('.nc'))
        elif '{member}' in fp:
            match = fp.format(var=native_var, member=member)
            files = sorted(f for f in all_files if match in f)
        else:
            match = fp.format(var=native_var)
            files = sorted(f for f in all_files if match in f)
 
    if not files:
        print(f"  No files for {model}/{scenario}/{member}/{variable} in {full_path}")
        return None
    print(f"  {model}/{scenario}/{member}/{variable}: {len(files)} files")
 
    # Eager .load() inside the with-block avoids fsspec I/O-on-closed-file.
    datasets = []
    time_span = (str(period.start_year), str(period.end_year))
    for f in files:
        with fsspec.open('s3://' + f, mode='rb') as fobj:
            ds_full = xr.open_dataset(fobj, engine='h5netcdf').load()
            ds = ds_full.sel(time=slice(*time_span)).load()
            if len(ds.time) > 0:
                datasets.append(ds)
            del ds_full, ds
 
    if not datasets:
        return None
    combined = xr.concat(datasets, dim='time')
    _, uniq = np.unique(combined.time.values, return_index=True)
    return combined.isel(time=sorted(uniq))
 
 
# ============================================================================
# SECTION 8  -  Load all members of one (model, scenario) for one variable
# ============================================================================
def load_scenario_members(model, variable, scenario, period, members=None):
    """
    
    Return {member: DataArray} in canonical units (mm/day or K).
 
    Members default to the layout list. The converter sets the units attribute
    so xclim accepts the array. Prints per-member timestep counts and means as
    a sanity trace.
    
    """
    layout = get_layout(model, scenario)
    var_info = get_variable_info(model, variable)
    converter = var_info['converter']
    if members is None:
        members = layout['members']
 
    print(f"Loading {model}/{scenario} {variable} ({var_info['native_short']}) "
          f"for {period.start_year}-{period.end_year}, members={members}")
 
    out = {}
    for m in members:
        ds = read_model_var(model, scenario, m, variable, period)
        if ds is not None:
            out[m] = converter(ds)
            del ds
            gc.collect()
 
    for m, da in out.items():
        print(f"  {m}: {len(da.time)} timesteps, "
              f"mean {float(da.mean()):.3f} {da.attrs.get('units', '')}")
    return out
# ============================================================================
# SECTION 9  -  Index registry
# ============================================================================
# NOTE on R95p/R99p: 
# R95p/R99p use the ETCCDI wet-day convention: the percentile is the Nth
# percentile of WET-day precipitation (days >= 1 mm/day), not of all days.
# This is enforced in compute_baseline_threshold via the 'wet_day_thresh' key:
# dry days are masked to NaN before percentile_doy, and percentile_doy is
# NaN-aware so the percentile is taken over wet days only. days_over_precip_thresh
# additionally applies its own wet-day mask at counting time (the 'thresh' kwarg),
# but it cannot correct a percentile that was biased by including dry days, which
# is why the masking must happen here at the threshold step. Temperature indices
# carry no 'wet_day_thresh' key and are computed over all days as usual.

INDEX_REGISTRY = {
    'R95p': {'variable': 'pr',     'percentile': 95, 'xclim_fn': days_over_precip_thresh,
             'threshold_kwarg': 'pr_per',     'extra_kwargs': {'thresh': '1 mm/day'}},
    'R99p': {'variable': 'pr',     'percentile': 99, 'xclim_fn': days_over_precip_thresh,
             'threshold_kwarg': 'pr_per',     'extra_kwargs': {'thresh': '1 mm/day'}},
    'TX90p': {'variable': 'tasmax', 'percentile': 90, 'xclim_fn': tx90p,
              'threshold_kwarg': 'tasmax_per', 'extra_kwargs': {}},
    'TX10p': {'variable': 'tasmax', 'percentile': 10, 'xclim_fn': tx10p,
              'threshold_kwarg': 'tasmax_per', 'extra_kwargs': {}},
    'TN90p': {'variable': 'tasmin', 'percentile': 90, 'xclim_fn': tn90p,
              'threshold_kwarg': 'tasmin_per', 'extra_kwargs': {}},
    'TN10p': {'variable': 'tasmin', 'percentile': 10, 'xclim_fn': tn10p,
              'threshold_kwarg': 'tasmin_per', 'extra_kwargs': {}},
    'WSDI': {'variable': 'tasmax', 'percentile': 90, 'xclim_fn': warm_spell_duration_index,
             'threshold_kwarg': 'tasmax_per', 'extra_kwargs': {'window': 6}},
    'CSDI': {'variable': 'tasmin', 'percentile': 10, 'xclim_fn': cold_spell_duration_index,
             'threshold_kwarg': 'tasmin_per', 'extra_kwargs': {'window': 6}},
}
 
 
# ============================================================================
# SECTION 10  -  Baseline percentile threshold
# ============================================================================
def compute_baseline_threshold(model, variable, percentile, baseline_period,
                               window=5, members=None, wet_day_thresh=None):
    """
    Day-of-year percentile threshold from SSP245 over the baseline window.
 
    percentile_doy with a 5-day centered window is the ETCCDI standard. The
    threshold is computed per member then averaged across members (Cindy Wang's
    convention; stacking members before percentile_doy would create duplicate
    (year, dayofyear) indices and crash). Returns (threshold, baseline_data) so
    the SSP245 baseline load can be reused when a comparison references it.

    wet_day_thresh : str or None
        For precipitation indices, a quantity like '1 mm/day'. When set, days
        below it are masked to NaN before percentile_doy so the percentile is
        taken over wet days only (ETCCDI convention; verified against the xclim
        docs for days_over_precip_thresh). The threshold is unit-converted to the
        data's units, so it is correct whether data is mm/day or kg/m2/s. Leave
        None for temperature, which uses all days. Note: in very arid cells with
        no wet days in a day-of-year window, the wet-day percentile is undefined
        (NaN); those cells correctly drop out of the precip index.
        
    """
    print(f"\nComputing baseline {model}/{variable} P{percentile} threshold from "
          f"SSP245 {baseline_period.start_year}-{baseline_period.end_year}...")
    baseline_data = load_scenario_members(model, variable, 'SSP245',
                                          baseline_period, members=members)
    if not baseline_data:
        raise RuntimeError(f"No baseline data for {model}/{variable}. "
                           f"Check MODEL_BUCKET_LAYOUT.")

    if wet_day_thresh is not None:
        print(f"  Wet-day masking before percentile (>= {wet_day_thresh})")
        
    print("  Computing percentile_doy per member...")
    per_member = {}
    for m, data in baseline_data.items():
        print(f"    {m}...")
        if wet_day_thresh is not None:
            # Express the wet-day threshold in the data's units, then mask dry
            # days to NaN. percentile_doy is NaN-aware, so the percentile is
            # computed over wet days only. Reassign attrs because .where() can
            # drop the units that downstream xclim checks rely on.
            thr = convert_units_to(wet_day_thresh, data)
            data_for_pct = data.where(data >= thr)
            data_for_pct.attrs.update(data.attrs)
        else:
            data_for_pct = data

        per_member[m] = percentile_doy(data_for_pct, per=percentile, window=window)
 
    stacked = xr.concat(list(per_member.values()),
                        dim=pd.Index(list(per_member.keys()), name='member'))
    threshold = stacked.mean('member')
    # Preserve percentile_doy attrs that downstream xclim functions check.
    threshold.attrs.update(per_member[next(iter(per_member))].attrs)
 
    print(f"  Threshold global mean: {float(threshold.mean()):.2f} "
          f"{threshold.attrs.get('units', '')}")
    print(f"  Range: {float(threshold.min()):.2f} to {float(threshold.max()):.2f}")
    return threshold, baseline_data
 
 
# ============================================================================
# SECTION 11  -  Index computation, ensemble mean, significance
# ============================================================================
def compute_index_for_members(data_dict, index_name, threshold, period):
    """
    Compute the time-mean index field for each member. {member: DataArray}.
 
    squeeze(drop=True) removes xclim's singleton 'percentiles' dimension, which
    otherwise breaks contourf downstream.
    
    """
    info = INDEX_REGISTRY[index_name]
    fn, thresh_kw, extra, var_kw = (info['xclim_fn'], info['threshold_kwarg'],
                                    info['extra_kwargs'], info['variable'])
    per_member = {}
    for m, data in data_dict.items():
        window = data.sel(time=period.slice())
        kwargs = {var_kw: window, thresh_kw: threshold, 'freq': 'YS', **extra}
        annual = fn(**kwargs).squeeze(drop=True)
        per_member[m] = annual.mean('time')
    return per_member
 
 
def ensemble_mean(per_member_dict):
    """
    Mean across members.
    
    """
    stacked = xr.concat(list(per_member_dict.values()),
                        dim=pd.Index(list(per_member_dict.keys()), name='member'))
    return stacked.mean('member')
 
 
def compute_significance(scenario_per_member, reference_per_member, equal_var=False):
    """
    Grid-point Welch t-test p-values between scenario and reference ensembles.
 
    p < 0.05 marks a difference unlikely to come from internal variability.
    Welch (equal_var=False) handles unequal member counts, which is the norm
    here (HiLLA n=3, SSP245 n up to 5).
    
    """
    s = xr.concat(list(scenario_per_member.values()),
                  dim=pd.Index(list(scenario_per_member.keys()), name='member'))
    r = xr.concat(list(reference_per_member.values()),
                  dim=pd.Index(list(reference_per_member.keys()), name='member'))
    _, p_val = stats.ttest_ind(s.values, r.values, axis=0, equal_var=equal_var)
    return xr.DataArray(p_val, coords={'lat': s.lat, 'lon': s.lon},
                        dims=['lat', 'lon'], name='p_value',
                        attrs={'test': 'Welch t-test',
                               'scenario_n': s.sizes['member'],
                               'reference_n': r.sizes['member']})
 
 
# ============================================================================
# SECTION 12  -  Top-level single comparison  (member resolution FIXED)
# ============================================================================
def run_comparison(model, index_name, config,
                   hilla_members=None, ssp245_members=None, sai_members=None,
                   precomputed_threshold=None, precomputed_baseline_data=None):
    """
    Run one anomaly comparison for one model and index.
 
    The member lists are resolved PER SIDE from each side's scenario name, so
    the warming comparison (SSP245 on both sides) sends SSP245 members to both
    loads instead of mismatching HiLLA members to an SSP245 path. For the
    HiLLA-vs-SSP245 comparisons the resolved lists are unchanged, so this is
    backward compatible and existing cached results stay valid.
 
    Returns a dict with scenario/reference ensemble means, anomaly, p-field,
    threshold, per-member arrays, and the config.
    
    """
    info = INDEX_REGISTRY[index_name]
    var, pct = info['variable'], info['percentile']
 
    hilla_layout = get_layout(model, 'HiLLA')
    ssp245_layout = get_layout(model, 'SSP245')
    if hilla_members is None:
        hilla_members = hilla_layout['members']
        
    if ssp245_members is None:
        ssp245_members = ssp245_layout['members']

    if sai_members is None and 'SAI' in (config.scenario_name, config.reference_name):
        sai_members = get_layout(model, 'SAI')['members']
 
    def members_for(scenario_name):
        if scenario_name == 'HiLLA':
            return hilla_members
        if scenario_name == 'SAI':
            return sai_members
        return ssp245_members
        
    scenario_members = members_for(config.scenario_name)
    reference_members = members_for(config.reference_name)
 
    print(f"\n{'='*70}")
    print(f"{model} {index_name}: {config.description}  [{config.label}]")
    print(f"  scenario ({config.scenario_name}) members:  {scenario_members}")
    print(f"  reference ({config.reference_name}) members: {reference_members}")
    print(f"{'='*70}")
 
    # 1. Threshold (always SSP245 baseline). Reused across the three Lee runs.
    #    info.get('wet_day_thresh') is '1 mm/day' for precip, None for temperature.
    if precomputed_threshold is not None:
        threshold = precomputed_threshold
        baseline_data = precomputed_baseline_data
        print("\nUsing precomputed baseline threshold.")
    else:
        threshold, baseline_data = compute_baseline_threshold(
            model=model, variable=var, percentile=pct,
            baseline_period=config.baseline_for_percentile,
            members=ssp245_members, wet_day_thresh=info.get('wet_day_thresh'))
 
    # 2. Scenario data: load, compute its index, then free the daily data so the
    #    reference load does not stack on top of it (this halves peak memory,
    #    which is what prevents the OOM on memory-limited hosts).
    scenario_data = load_scenario_members(
        model, var, config.scenario_name, config.scenario_period,
        members=scenario_members)
    
    log_mem("after scenario load")
    print(f"\nComputing {index_name} for scenario...")
    scenario_pm = compute_index_for_members(scenario_data, index_name,
                                            threshold, config.scenario_period)
    scenario_mean = ensemble_mean(scenario_pm)
    del scenario_data
    gc.collect()

    log_mem("after scenario freed")
    
    # 3. Reference data:reuse the baseline daily data only if the caller passed it
    #    and the periods match; otherwise load fresh. Free it after the index
    #    unless it is the caller-owned baseline.
    same_as_baseline = (
        config.reference_name == 'SSP245'
        and config.reference_period.start_year == config.baseline_for_percentile.start_year
        and config.reference_period.end_year == config.baseline_for_percentile.end_year
    )
    if same_as_baseline and baseline_data is not None:
        print("\nReusing baseline SSP245 data as reference (same period).")
        reference_data = baseline_data
        reused_baseline = True
    else:
        reference_data = load_scenario_members(
            model, var, config.reference_name, config.reference_period,
            members=reference_members)
        log_mem("after reference load")
        reused_baseline = False
 
    # 4. Reference index, ensemble means.
    print(f"Computing {index_name} for reference...")
    reference_pm = compute_index_for_members(reference_data, index_name,
                                             threshold, config.reference_period)
    reference_mean = ensemble_mean(reference_pm)
    if not reused_baseline:
        del reference_data
        gc.collect()


    reference_mean = reference_mean.assign_coords(lat=scenario_mean.lat, lon=scenario_mean.lon)
    reference_pm   = {m: v.assign_coords(lat=scenario_mean.lat, lon=scenario_mean.lon)
                      for m, v in reference_pm.items()}
    # 5. Anomaly. 6. Significance.
    anomaly = (scenario_mean - reference_mean).squeeze(drop=True)
    p_field = compute_significance(scenario_pm, reference_pm)
 
    print(f"\n{model} {index_name} [{config.label}] stats:")
    print(f"  scenario  global mean: {float(scenario_mean.mean()):.2f}")
    print(f"  reference global mean: {float(reference_mean.mean()):.2f}")
    print(f"  anomaly   global mean: {float(anomaly.mean()):.3f}  "
          f"range {float(anomaly.min()):.2f}..{float(anomaly.max()):.2f}")
    print(f"  significant @p<0.05:   {float((p_field < 0.05).mean()):.1%}")
 
    return {
        'model': model, 'index_name': index_name, 'config': config,
        'scenario': scenario_mean, 'reference': reference_mean,
        'anomaly': anomaly, 'p_field': p_field, 'threshold': threshold,
        'scenario_per_member': scenario_pm, 'reference_per_member': reference_pm,
    }
# ============================================================================
# SECTION 13  -  Lee et al. (2026) three-way framework + additivity check
# ============================================================================
def run_lee_framework(model, index_name, hilla_members=None, ssp245_members=None):
    """
    Run warming, sai, and combined for one model+index, sharing the threshold.
 
    The SSP245-baseline percentile threshold and baseline data are computed once
    and passed to all three comparisons (correct, since all three reference the
    same baseline yardstick, and fast, since percentile_doy is the costly step).
    Returns {'warming': res, 'sai': res, 'combined': res}.

    This version holds the baseline daily data in memory across all three
    comparisons. On a memory-limited host that can OOM; prefer
    run_lee_framework_cached, which caches each step to disk and frees the daily
    data between comparisons.
    
    """
    ssp245_layout = get_layout(model, 'SSP245')
    hilla_layout = get_layout(model, 'HiLLA')
    if hilla_members is None:
        hilla_members = hilla_layout['members']
    if ssp245_members is None:
        ssp245_members = ssp245_layout['members']
 
    info = INDEX_REGISTRY[index_name]
    threshold, baseline_data = compute_baseline_threshold(
        model=model, variable=info['variable'], percentile=info['percentile'],
        baseline_period=BASELINE, members=ssp245_members, wet_day_thresh=info.get('wet_day_thresh'))
 
    results = {}
    for key in ('warming', 'sai', 'combined'):
        cfg = make_comparison_config(key, assessment_for(model))
        print(f"\n----- {model} {index_name}: {cfg.label} -----")
        results[key] = run_comparison(
            model, index_name, cfg,
            hilla_members=hilla_members, ssp245_members=ssp245_members,
            precomputed_threshold=threshold, precomputed_baseline_data=baseline_data)
    return results
 
 
def check_additivity(results_dict):
    """
    Verify (warming + sai) == combined to floating-point precision.
    
    """
    w = results_dict['warming']['anomaly']
    s = results_dict['sai']['anomaly']
    c = results_dict['combined']['anomaly']
    resid = (w + s) - c
    print(f"  Additivity check (warming + sai - combined): "
          f"max|resid|={float(np.abs(resid).max()):.3e}, "
          f"mean|resid|={float(np.abs(resid).mean()):.3e}  (expect ~0)")
 

def check_member_distinctness(model, variable, scenario, period,
                              members=None, sample_lat=40.0, sample_lon=260.0):
    """
    Diagnose whether ensemble members are genuinely distinct or duplicates.
 
    Motivated by the E3SM SSP245 case where all three members reported an
    identical global-mean precip (2.510 mm/day). A shared global mean is a weak
    signal: distinct members share a climatology and will agree on long-window
    means and on the seasonal cycle, so neither matching means nor a high daily
    correlation proves duplication. The decisive signals are differences that are
    exactly zero ONLY when the underlying files are the same:
 
      1. Spatial RMS difference between member time-mean fields. Distinct members
         differ by a small but non-zero amount (internal variability not fully
         averaged out over 20 years). Exactly zero implies the same file.
      2. Daily series at one grid cell, pairwise. Distinct members have divergent
         day-to-day weather, so the max absolute daily difference is large (several
         K, or tens of mm/day) and the fraction of exactly-equal days is ~0.
         Duplicates give max|delta| = 0 and 100% equal days.
 
    sample_lon is in the data's longitude convention (0-360 here; 260 ~ central
    US). The exact cell does not matter for the verdict, only that it is the same
    cell across members, so method='nearest' is used and any convention is fine.
 
    Run this before trusting the significance test: duplicated members give zero
    inter-member variance and therefore invalid p-values.
    Returns the loaded {member: DataArray} dict so the load can be reused.
    
    """
    data = load_scenario_members(model, variable, scenario, period, members=members)
    keys = list(data.keys())
    print(f"\n=== Member distinctness: {model}/{scenario}/{variable} "
          f"{period.start_year}-{period.end_year} ===")
    if len(keys) < 2:
        print(f"  Only {len(keys)} member(s); nothing to compare.")
        return data
 
    for m in keys:
        da = data[m]
        print(f"  {m}: global mean={float(da.mean()):.6f}, "
              f"overall std={float(da.std()):.6f}")
 
    tmeans = {m: data[m].mean('time') for m in keys}
    series = {m: np.asarray(
        data[m].sel(lat=sample_lat, lon=sample_lon, method='nearest').values
    ).ravel() for m in keys}
 
    print(f"\n  Pairwise (sample cell ~lat {sample_lat}, lon {sample_lon}):")
    all_identical = True
    for a, b in itertools.combinations(keys, 2):
        rms = float(np.sqrt(float(((tmeans[a] - tmeans[b]) ** 2).mean())))
        sa, sb = series[a], series[b]
        n = min(len(sa), len(sb))
        sa, sb = sa[:n], sb[:n]
        max_abs = float(np.nanmax(np.abs(sa - sb))) if n else float('nan')
        frac_equal = float(np.mean(np.abs(sa - sb) < 1e-9)) if n else float('nan')
        identical = (rms < 1e-9) and (max_abs < 1e-9)
        all_identical = all_identical and identical
        verdict = "IDENTICAL (duplicate!)" if identical else "distinct"
        print(f"    {a} vs {b}: field RMS diff={rms:.3e}, "
              f"sample max|delta|={max_abs:.3e}, "
              f"equal-day fraction={frac_equal:.3f}  -> {verdict}")
 
    if all_identical:
        print("\n  VERDICT: members appear to be DUPLICATES. The significance "
              "test is invalid (zero inter-member variance). Trace the date-stamp "
              "version-directory selection in read_e3sm_var: all cases may be "
              "resolving to the same files.")
    else:
        print("\n  VERDICT: members are DISTINCT. Matching global means are "
              "coincidental, expected over a 20-year window.")
    return data


# ============================================================================
# SECTION 14  -  Caching: flatten a result dict to netCDF and back
# ============================================================================
def save_result(result, path):
    """
    Write one run_comparison result dict to a single netCDF file.
    
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    ds = xr.Dataset()
    ds['scenario'] = result['scenario']
    ds['reference'] = result['reference']
    ds['anomaly'] = result['anomaly']
    ds['threshold'] = result['threshold']
    ds['p_field'] = result['p_field']
    for tag in ('scenario_per_member', 'reference_per_member'):
        d = result[tag]
        ds[tag] = xr.concat(list(d.values()),
                            dim=pd.Index(list(d.keys()), name=f'{tag}_member'))
    ds.attrs['index_name'] = result['index_name']
    ds.attrs['model'] = result['model']
    ds.attrs['comparison_label'] = result['config'].label
    ds.to_netcdf(path)
 
 
def load_result(path, config, index_name, model):
    """
    Reconstruct a result dict from a cached netCDF file.
    
    """
    ds = xr.open_dataset(path)
    spm = {m: ds['scenario_per_member'].sel(scenario_per_member_member=m)
              .drop_vars('scenario_per_member_member')
           for m in ds['scenario_per_member_member'].values}
    rpm = {m: ds['reference_per_member'].sel(reference_per_member_member=m)
              .drop_vars('reference_per_member_member')
           for m in ds['reference_per_member_member'].values}
    return {
        'model': model, 'index_name': index_name, 'config': config,
        'scenario': ds['scenario'], 'reference': ds['reference'],
        'anomaly': ds['anomaly'], 'threshold': ds['threshold'],
        'p_field': ds['p_field'] if 'p_field' in ds else None,
        'scenario_per_member': spm, 'reference_per_member': rpm,
    }
 

def save_threshold(threshold, path):
    """
    Write the percentile threshold (a DataArray) to its own netCDF.
 
    Cached separately from the comparison results because it is the expensive
    step (percentile_doy over the baseline ensemble) and is identical across all
    three Lee comparisons. Caching it means a kernel restart never recomputes it.
    
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    threshold.to_dataset(name='threshold').to_netcdf(path)
 
 
def load_threshold(path):
    """
    Load a cached threshold DataArray, eagerly, and close the file handle.
    
    """
    ds = xr.open_dataset(path)
    thr = ds['threshold'].load()
    ds.close()
    return thr


def run_lee_framework_cached(model, index_name, hilla_members=None,
                             ssp245_members=None, free_memory=True):
    """
    Run warming/sai/combined for one model+index with per-step disk caching.
 
    Resumable at the finest useful granularity. The SSP245-baseline percentile
    threshold and each of the three comparisons are written to their own netCDF
    the instant they finish, so a kernel restart mid-run resumes from disk rather
    than starting the model over. On rerun, anything already cached is loaded and
    skipped.
 
    Memory-bounded by design, which is the fix for hosts that OOM:
      - The threshold is computed once (loading the baseline daily data), cached,
        and the baseline daily data is then freed.
      - Each comparison is run with precomputed_baseline_data=None, so it loads
        only the daily data it needs, frees the scenario data before loading the
        reference (see run_comparison), and frees the result before the next
        comparison. Peak memory is therefore one comparison's daily data, not all
        three plus the baseline held at once.
 
    The trade is extra I/O: the SSP245 baseline is reloaded as the reference for
    the warming and combined comparisons rather than kept in memory. That is the
    right trade when OOM is the blocker, and on rerun the cached steps are skipped
    so the reloads only happen for work not yet done.
    
    """
    keys = ('warming', 'sai', 'combined')
    cache = OUTDIR / 'cached_data'
    cache.mkdir(parents=True, exist_ok=True)
    
    a = assessment_for(model)
    win = f"{a.start_year}_{a.end_year}"
    comp_paths = {k: cache / f"{model.lower()}_{index_name.lower()}_{k}_{win}.nc"
                  for k in keys}
    thr_path = cache / f"{model.lower()}_{index_name.lower()}_threshold.nc"
 
    info = INDEX_REGISTRY[index_name]
    hilla_layout = get_layout(model, 'HiLLA')
    ssp245_layout = get_layout(model, 'SSP245')
    if hilla_members is None:
        hilla_members = hilla_layout['members']
    if ssp245_members is None:
        ssp245_members = ssp245_layout['members']
 
    # --- Threshold: load from cache, else compute once, cache, free baseline ---
    if thr_path.exists():
        print(f"Loading cached threshold {thr_path.name}")
        threshold = load_threshold(thr_path)
    else:
        threshold, baseline_data = compute_baseline_threshold(
            model=model, variable=info['variable'], percentile=info['percentile'],
            baseline_period=BASELINE, members=ssp245_members,
            wet_day_thresh=info.get('wet_day_thresh'))
        save_threshold(threshold, thr_path)
        print(f"  Saved {thr_path.name}")
        del baseline_data          # comparisons reload only what they need
        gc.collect()
 
    # --- Each comparison: skip if cached, else compute, save now, free ---
    for k in keys:
        if comp_paths[k].exists():
            print(f"Cached {comp_paths[k].name} exists, skipping compute")
            continue
        cfg = make_comparison_config(k, assessment_for(model))
        print(f"\n----- {model} {index_name}: {cfg.label} -----")
        res = run_comparison(
            model, index_name, cfg,
            hilla_members=hilla_members, ssp245_members=ssp245_members,
            precomputed_threshold=threshold, precomputed_baseline_data=None)
        save_result(res, comp_paths[k])
        print(f"  Saved {comp_paths[k].name}")
        if free_memory:
            del res
            gc.collect()
 
    # --- Load all three from cache for return / plotting ---
    return {k: load_result(comp_paths[k], make_comparison_config(k, assessment_for(model)),
                           index_name, model) for k in keys}


def run_comparison_cached(model, index_name, comparison_type,
                          hilla_members=None, ssp245_members=None,
                          sai_members=None, assessment=None, free_memory=True):
    """
    Run ONE comparison for a model+index with the same per-step disk caching
    as run_lee_framework_cached. For comparisons outside the warming/sai/combined
    trio, e.g. 'sai_vs_hilla'. The threshold cache file is shared with the trio
    ({model}_{index}_threshold.nc), so it is reused if present.
    
    """
    cache = OUTDIR / 'cached_data'
    cache.mkdir(parents=True, exist_ok=True)
    cfg = make_comparison_config(comparison_type, assessment or assessment_for(model))
    key = comparison_type.lower()
    comp_path = cache / f"{model.lower()}_{index_name.lower()}_{key}.nc"
    thr_path = cache / f"{model.lower()}_{index_name.lower()}_threshold.nc"

    info = INDEX_REGISTRY[index_name]
    if ssp245_members is None:
        ssp245_members = get_layout(model, 'SSP245')['members']

    if thr_path.exists():
        print(f"Loading cached threshold {thr_path.name}")
        threshold = load_threshold(thr_path)
    else:
        threshold, baseline_data = compute_baseline_threshold(
            model=model, variable=info['variable'], percentile=info['percentile'],
            baseline_period=BASELINE, members=ssp245_members,
            wet_day_thresh=info.get('wet_day_thresh'))
        save_threshold(threshold, thr_path)
        print(f"  Saved {thr_path.name}")
        del baseline_data
        gc.collect()

    if comp_path.exists():
        print(f"Cached {comp_path.name} exists, skipping compute")
    else:
        print(f"\n----- {model} {index_name}: {cfg.label} -----")
        res = run_comparison(
            model, index_name, cfg,
            hilla_members=hilla_members, ssp245_members=ssp245_members,
            sai_members=sai_members,
            precomputed_threshold=threshold, precomputed_baseline_data=None)
        save_result(res, comp_path)
        print(f"  Saved {comp_path.name}")
        if free_memory:
            del res
            gc.collect()

    return load_result(comp_path, cfg, index_name, model)
# ============================================================================
# SECTION 14B  -  Per-member annual index dataset writer (etccdi-indices format)
# ============================================================================
# Writes each percentile index as per-member ANNUAL fields into the Reflective
# etccdi-indices tree, so make_zarr_store.ipynb folds them into the shared zarr
# beside the 16 fixed indices. This will be the D2 dataset path. It does not touch the
# Lee comparison code, which stays the analysis path for the figures.

def save_nc_to_s3(obj, s3_path):
    """
    Write a DataArray/Dataset to S3 via a local temp file (same as save_nc pattern on GitHub).
    
    """
    with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
        tmp = f.name
    try:
        obj.to_netcdf(tmp, engine='h5netcdf')
        s3.put(tmp, s3_path.replace('s3://', ''))
    finally:
        os.remove(tmp)


def compute_annual_index_for_member(da, index_name, threshold):
    """
    Annual index field (time, lat, lon) for one member, against a fixed
    baseline percentile. This is compute_index_for_members without the time mean,
    so every year is kept and the file matches the fixed-index annual files. The
    singleton 'percentiles' dim xclim attaches is squeezed so the field stacks
    cleanly into (member, year, lat, lon).
    
    """
    info = INDEX_REGISTRY[index_name]
    kwargs = {info['variable']: da, info['threshold_kwarg']: threshold,
              'freq': 'YS', **info['extra_kwargs']}
    annual = info['xclim_fn'](**kwargs).squeeze(drop=True)
    annual.name = index_name.upper()
    return annual


def write_index_dataset(model, index_name,
                        scenarios=('HiLLA', 'SSP245'),
                        scenario_periods=None,
                        members_by_scenario=None,
                        ssp245_threshold_members=None,
                        output_root=ETCCDI_OUTPUT_ROOT,
                        overwrite=False):
    """
    Compute one percentile index as per-member annual fields and write them to
        {output_root}/{model_label}/{scenario_label}/{member_label}/{INDEX}.nc

    The percentile threshold is the SSP245 2020-2039 day-of-year percentile,
    computed once and applied to every scenario so counts are comparable across
    SSP245, HiLLA, and later SAI. Members are read one at a time and freed, so
    peak memory is one member's daily record.
    
    """
    info = INDEX_REGISTRY[index_name]
    variable = info['variable']
    model_label = MODEL_LABEL[model]
    idx_file = f'{index_name.upper()}.nc'

    default_periods = {
        'HiLLA':  Period('hilla',  2035, 2084),   # GeoMIP G6 window, edit if needed
        'SAI':    Period('sai',    2035, 2084),
        'SSP245': Period('ssp245', 2015, 2100),   # wide: reader keeps what exists
    }
    if scenario_periods:
        default_periods.update(scenario_periods)

    if ssp245_threshold_members is None:
        ssp245_threshold_members = get_layout(model, 'SSP245')['members']
    threshold, baseline_data = compute_baseline_threshold(
        model=model, variable=variable, percentile=info['percentile'],
        baseline_period=BASELINE, members=ssp245_threshold_members,
        wet_day_thresh=info.get('wet_day_thresh'))
    del baseline_data
    gc.collect()

    written = []
    for scenario in scenarios:
        scen_label = SCENARIO_LABEL[(model, scenario)]
        period = default_periods[scenario]
        layout = get_layout(model, scenario)
        members = (members_by_scenario or {}).get(scenario, layout['members'])
        label_map = MEMBER_LABEL[(model, scenario)]

        for member in members:
            member_label = label_map.get(member, member)
            s3_path = f'{output_root}/{model_label}/{scen_label}/{member_label}/{idx_file}'
            if not overwrite and s3.exists(s3_path.replace('s3://', '')):
                print(f'  exists, skip: {model_label}/{scen_label}/{member_label}/{idx_file}')
                written.append(s3_path)
                continue

            ds = read_model_var(model, scenario, member, variable, period)
            if ds is None:
                print(f'  [WARN] no data: {model}/{scenario}/{member}/{variable}')
                continue
            da = get_variable_info(model, variable)['converter'](ds)
            del ds

            annual = compute_annual_index_for_member(da, index_name, threshold).load()
            del da
            gc.collect()

            save_nc_to_s3(annual, s3_path)
            yrs = annual['time'].dt.year.values
            print(f'  wrote {model_label}/{scen_label}/{member_label}/{idx_file}  '
                  f'({len(yrs)} yrs {int(yrs.min())}-{int(yrs.max())})')
            written.append(s3_path)
            del annual
            gc.collect()

    print(f'Done {index_name.upper()} for {model_label}: {len(written)} files.')
    return written
# ============================================================================
# SECTION 15  -  Plotting
# ============================================================================
def convert_longitude_range(lon_min=-180.0, lon_max=180.0):
    if abs(lon_max - lon_min) >= 360:
        return 0.0, 360.0
    return lon_min % 360, lon_max % 360
 
 
def _basemap(ax, extent, proj_is_regional):
    if proj_is_regional:
        lon_min, lon_max = convert_longitude_range(extent[0], extent[1])
        ax.set_extent([lon_min, lon_max, extent[2], extent[3]], crs.PlateCarree())
    else:
        ax.set_global()
    ax.coastlines(linewidth=0.5)
    ax.add_feature(cfeature.BORDERS.with_scale('110m'), linewidth=0.3, alpha=0.5)
 

# Per-model land fraction. Fill 'path' from the bucket to use the model's own
# field; left None it falls back to regionmask. scale brings the field to 0-1
# (sftlf is 0-100, LANDFRAC is 0-1); thresh is the land cutoff.
LAND_FRACTION = {
    'CESM':  {'path': None, 'var': 'LANDFRAC', 'scale': 1.0,  'thresh': 0.5},
    'UKESM': {'path': None, 'var': 'sftlf',    'scale': 0.01, 'thresh': 0.5},
    'MIROC': {'path': None, 'var': 'sftlf',    'scale': 0.01, 'thresh': 0.5},
    'E3SM':  {'path': None, 'var': 'LANDFRAC', 'scale': 1.0,  'thresh': 0.5},
}
_land_mask_cache = {}

def get_land_mask(model, lat, lon):
    """
    Boolean DataArray, True over land, on the (lat, lon) grid passed in.
    Model land fraction if configured (Alistair's request), else regionmask
    Natural Earth, else all-True with a warning. Cached per (model, shape).
    
    """
    key = (model, len(lat), len(lon))
    if key in _land_mask_cache:
        return _land_mask_cache[key]
    cfg = LAND_FRACTION.get(model, {})
    mask = None
    if cfg.get('path'):
        try:
            with fsspec.open(cfg['path'], mode='rb') as f:
                lf = xr.open_dataset(f, engine='h5netcdf')[cfg['var']].squeeze(drop=True)
            ren = {}
            if 'latitude' in lf.coords:  ren['latitude'] = 'lat'
            if 'longitude' in lf.coords: ren['longitude'] = 'lon'
            if ren: lf = lf.rename(ren)
            lf = lf.assign_coords(lat=lf['lat'].astype('float64'),
                                  lon=lf['lon'].astype('float64'))
            mask = (lf * cfg['scale']) >= cfg['thresh']
        except Exception as e:
            print(f"  [land mask] {model}: land fraction read failed ({e}); falling back")
    if mask is None:
        try:
            import regionmask
            m = regionmask.defined_regions.natural_earth_v5_0_0.land_110.mask(lon, lat)
            mask = xr.DataArray(~np.isnan(m.values), dims=('lat', 'lon'),
                                coords={'lat': lat, 'lon': lon})
        except Exception:
            print(f"  [land mask] {model}: no land fraction and no regionmask; using all cells")
            mask = xr.DataArray(np.ones((len(lat), len(lon)), bool),
                                dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})
    _land_mask_cache[key] = mask
    return mask

def mask_to_land(field, model):
    """
    Field with ocean set to NaN, using the model's land mask.
    
    """
    return field.where(get_land_mask(model, field['lat'], field['lon']))


def _auto_sym_levels(arrays, n=21, q=98, fallback=60.0):
    """
    Symmetric contour levels from the robust max of |arrays| (NaN-aware),
    rounded up to a round number. Used to fit CESM's own colorbar to its signal.
    
    """
    import math
    vals = [np.nanpercentile(np.abs(np.asarray(a, float)), q)
            for a in arrays if np.isfinite(np.asarray(a, float)).any()]
    v = max(vals) if vals else fallback
    if not np.isfinite(v) or v <= 0:
        v = fallback
    step = 10 ** math.floor(math.log10(v))
    v = math.ceil(v / step) * step
    return np.linspace(-v, v, n)


def plot_all_mod_comparisons(all_results, index_name, models,
                             extent=(-125, 40, 10, 80),
                             cmap=None, levels=None, levels_cesm=None,
                             outlier_model='CESM', land_only=True,
                             units_label='days/year', p_threshold=0.05,
                             save_path=None,
                             cbar_height=0.03, cbar_y=0.05):
    """
    Lee Fig 5 analogue: models as rows, warming/SAI/combined as columns.

    Land-masked by default (variability differs land vs ocean). The equal-window
    models (UKESM/MIROC/E3SM, 2065-2084) share one colorbar; the short-window
    model (CESM-WACCM, 2050-2069) gets its own scale and colorbar so its weaker
    signal stays visible, per Alistair, without altering the data.
    
    """
    if cmap is None:
        cmap = 'BrBG' if index_name in PRECIP_INDICES else 'RdBu_r'
    if levels is None:
        levels = (np.arange(-12, 12.0001, 1) if index_name in PRECIP_INDICES
                  else np.arange(-60, 60.0001, 5))

    col_keys = ['warming', 'sai', 'combined']
    col_titles = ['Global warming (1\u21922)', 'SAI (2\u21923)', 'Warming + SAI (1\u21923)']

    def _anom(model, key):
        a = all_results.get(model, {}).get(index_name, {}).get(key, {}).get('anomaly')
        if a is None:
            return None
        return mask_to_land(a, model) if land_only else a

    if outlier_model in models and levels_cesm is None:
        arrs = [_anom(outlier_model, k).values for k in col_keys
                if _anom(outlier_model, k) is not None]
        fb = 12.0 if index_name in PRECIP_INDICES else 60.0
        levels_cesm = _auto_sym_levels(arrs, n=len(levels), fallback=fb) if arrs else levels

    ncol, nrow = 3, len(models)
    regional = extent is not None
    proj = crs.Mercator() if regional else crs.Robinson()
    if regional:
        lat_lo, lat_hi = sorted((extent[2] - 3, extent[3] + 3))

    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 10),
                             subplot_kw={'projection': proj}, squeeze=False,
                             sharey=True, sharex=True)
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(left=0.05, bottom=0.16, hspace=0.0, wspace=0.08)

    right_label_idx = {i * ncol + (ncol - 1) for i in range(nrow)}
    bottom_label_idx = {(nrow - 1) * ncol + j for j in range(ncol)}

    main_mappable = cesm_mappable = None
    for i, model in enumerate(models):
        rd = all_results.get(model, {}).get(index_name)
        lv = levels_cesm if model == outlier_model else levels
        for j, key in enumerate(col_keys):
            n = i * ncol + j
            ax = axes[i][j]
            if regional:
                lon_min, lon_max = convert_longitude_range(extent[0], extent[1])
                ax.set_extent([lon_min, lon_max, extent[2], extent[3]], crs.PlateCarree())
            else:
                ax.set_global()
            ax.coastlines(linewidth=0.8)
            ax.add_feature(cfeature.BORDERS.with_scale('110m'), linewidth=0.5, alpha=0.8)
            gl = ax.gridlines(draw_labels=True, x_inline=False, alpha=0.25)
            gl.top_labels = gl.bottom_labels = gl.left_labels = gl.right_labels = False
            if n in right_label_idx:
                gl.right_labels = True
                gl.yformatter = LATITUDE_FORMATTER
                gl.ylabel_style = {'color': 'green', 'weight': 'bold', 'fontsize': 14}
            if n in bottom_label_idx:
                gl.bottom_labels = True
                gl.xformatter = LONGITUDE_FORMATTER
                gl.xlabel_style = {'color': 'red', 'weight': 'bold', 'rotation': 45, 'fontsize': 14}

            if rd is None:
                ax.set_visible(False)
                continue

            anom = _anom(model, key)
            if regional:
                anom = anom.sortby('lat').sel(lat=slice(lat_lo, lat_hi))
            data, lons = add_cyclic_point(anom.values, coord=anom.lon, axis=-1)
            lon2d, lat2d = np.meshgrid(lons, anom.lat)
            cf = ax.contourf(lon2d, lat2d, data, levels=lv,
                             transform=crs.PlateCarree(), cmap=cmap, extend='both',
                             transform_first=True)
            if model == outlier_model:
                cesm_mappable = cf
            else:
                main_mappable = cf

            res = rd[key]
            if res.get('scenario_per_member'):
                p = compute_significance(res['scenario_per_member'], res['reference_per_member'])
                not_sig = (p >= p_threshold).where(p >= p_threshold)
                if land_only:
                    not_sig = mask_to_land(not_sig, model)
                if regional:
                    not_sig = not_sig.sortby('lat').sel(lat=slice(lat_lo, lat_hi))
                hlon2d, hlat2d = np.meshgrid(not_sig.lon, not_sig.lat)
                ax.contourf(hlon2d, hlat2d, not_sig, levels=[0.5, 1.5],
                            colors='none', hatches=['....'],
                            transform=crs.PlateCarree(), transform_first=True)
            if i == 0:
                ax.set_title(col_titles[j], fontsize=16)
            if j == 0:
                ax.text(-0.10, 0.5, model, transform=ax.transAxes, rotation=90,
                        va='center', ha='center', fontsize=16, fontweight='bold')

    def _add_cbar(mappable, x0, width, label):
        cax = fig.add_axes([x0, cbar_y, width, cbar_height])
        cb = fig.colorbar(mappable, cax=cax, orientation='horizontal', extend='both')
        cb.ax.tick_params(labelsize=13)
        cb.set_label(label, fontsize=13)

    main_models = [m for m in models if m != outlier_model]
    if cesm_mappable is not None and main_mappable is not None:
        _add_cbar(main_mappable, 0.09, 0.38, f"{' / '.join(main_models)}  ({units_label})")
        _add_cbar(cesm_mappable, 0.55, 0.38, f"{outlier_model} 2050-2069  ({units_label})")
    else:
        m = main_mappable or cesm_mappable
        if m is not None:
            _add_cbar(m, 0.31, 0.38, f"{index_name}  ({units_label})")

    fig.suptitle(f"{index_name}: warming / SAI / combined" + ("  (land)" if land_only else ""),
                 fontsize=16, y=0.92)
    if save_path:
        plt.savefig(save_path, dpi=100, bbox_inches='tight')
        print(f"  Saved: {save_path}")
    plt.show()


def plot_contrast_summary(contrast, models, indices=('R95p', 'TX90p', 'TN10p'),
                          extent=(-125, 40, 20, 70), land_only=True,
                          p_threshold=0.05, save_path=None,
                          cbar_height=0.03, cbar_y=0.12):
    """
    SAI minus HiLLA, rows = models, columns = indices, NH-land Lee style.
    One window (2050-2069) for all models, so no separate CESM colorbar; each
    index column gets its own. E3SM is drawn only for precipitation indices.
    
    """
    def valid(model, idx):
        return not (model == 'E3SM' and idx in TEMP_INDICES)

    def anom(idx, model):
        r = contrast.get(idx, {}).get(model)
        if r is None or not valid(model, idx):
            return None, None
        a = mask_to_land(r['anomaly'], model) if land_only else r['anomaly']
        return a, r

    ncol, nrow = len(indices), len(models)
    regional = extent is not None
    proj = crs.Mercator() if regional else crs.Robinson()
    if regional:
        lat_lo, lat_hi = sorted((extent[2] - 3, extent[3] + 3))

    def crop(da):
        return da.sortby('lat').sel(lat=slice(lat_lo, lat_hi)) if regional else da

    col_levels, col_cmap, col_map = {}, {}, {}
    for idx in indices:
        arrs = [crop(anom(idx, m)[0]).values for m in models if anom(idx, m)[0] is not None]
        fb = 12.0 if idx in PRECIP_INDICES else 30.0
        col_levels[idx] = _auto_sym_levels(arrs, n=21, fallback=fb)
        col_cmap[idx] = 'BrBG' if idx in PRECIP_INDICES else 'RdBu_r'

    low_vis = {j: max([i for i, m in enumerate(models) if anom(idx, m)[0] is not None] or [-1])
               for j, idx in enumerate(indices)}
    right_vis = {i: max([j for j, idx in enumerate(indices) if anom(idx, m)[0] is not None] or [-1])
                 for i, m in enumerate(models)}

    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 2.6 * nrow + 1),
                             subplot_kw={'projection': proj}, squeeze=False,
                             sharey=True, sharex=True)
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(left=0.06, bottom=0.14, hspace=-0.5, wspace=0.08)

    for i, model in enumerate(models):
        for j, idx in enumerate(indices):
            ax = axes[i][j]
            if regional:
                lon_min, lon_max = convert_longitude_range(extent[0], extent[1])
                ax.set_extent([lon_min, lon_max, extent[2], extent[3]], crs.PlateCarree())
            else:
                ax.set_global()
            ax.coastlines(linewidth=0.8)
            ax.add_feature(cfeature.BORDERS.with_scale('110m'), linewidth=0.5, alpha=0.8)
            gl = ax.gridlines(draw_labels=True, x_inline=False, alpha=0.25)
            gl.top_labels = gl.bottom_labels = gl.left_labels = gl.right_labels = False
            if j == right_vis[i]:
                gl.right_labels = True; gl.yformatter = LATITUDE_FORMATTER
                gl.ylabel_style = {'color': 'green', 'weight': 'bold', 'fontsize': 13}
            if i == low_vis[j]:
                gl.bottom_labels = True; gl.xformatter = LONGITUDE_FORMATTER
                gl.xlabel_style = {'color': 'red', 'weight': 'bold', 'rotation': 45, 'fontsize': 13}

            a, r = anom(idx, model)
            if a is None:
                ax.set_visible(False); continue
            a = crop(a)
            data, lons = add_cyclic_point(a.values, coord=a.lon, axis=-1)
            lon2d, lat2d = np.meshgrid(lons, a.lat)
            cf = ax.contourf(lon2d, lat2d, data, levels=col_levels[idx],
                             transform=crs.PlateCarree(), cmap=col_cmap[idx],
                             extend='both', transform_first=True)
            col_map[idx] = cf
            p = r.get('p_field')
            if p is not None:
                ns = (p >= p_threshold).where(p >= p_threshold)
                if land_only:
                    ns = mask_to_land(ns, model)
                ns = crop(ns)
                h_lon, h_lat = np.meshgrid(ns.lon, ns.lat)
                ax.contourf(h_lon, h_lat, ns, levels=[0.5, 1.5], colors='none',
                            hatches=['....'], transform=crs.PlateCarree(), transform_first=True)
            if i == 0:
                ax.set_title(idx, fontsize=15)
            if j == 0:
                ax.text(-0.12, 0.5, model, transform=ax.transAxes, rotation=90,
                        va='center', ha='center', fontsize=15, fontweight='bold')

    for j, idx in enumerate(indices):
        if idx not in col_map or low_vis[j] < 0:
            continue
        pos = axes[low_vis[j]][j].get_position()
        w = pos.width * 0.9
        cax = fig.add_axes([pos.x0 + (pos.width - w) / 2, cbar_y, w, cbar_height])
        cb = fig.colorbar(col_map[idx], cax=cax, orientation='horizontal', extend='both')
        cb.ax.tick_params(labelsize=12)
        cb.set_label(f"{idx}  SAI - HiLLA  (days/year)", fontsize=12)

    fig.suptitle("G6-1.5K-SAI minus G6-1.5K-HiLLA  (land, 2050-2069)", fontsize=15, y=0.85)
    if save_path:
        plt.savefig(save_path, dpi=100, bbox_inches='tight')
        print(f"  Saved: {save_path}")
    plt.show()


def land_stats(res, lat_band=None, lon_band=None):
    """
    Land-only mean anomaly and significant fraction for one comparison result,
    optionally within a lat band (lo, hi) and/or lon band (lo, hi). lon_band is
    matched on a -180..180 axis, so a Mediterranean box (-10, 40) keeps its
    western half on a 0-360 grid. Model is read from the result.
    
    """
    model = res['model']
    p = res['p_field']
    land = get_land_mask(model, p.lat, p.lon)
    if lat_band is not None:
        lo, hi = lat_band
        land = land & (p.lat >= lo) & (p.lat <= hi)
    if lon_band is not None:
        lo, hi = lon_band
        lon = ((p.lon + 180) % 360) - 180
        land = land & (lon >= lo) & (lon <= hi)
    return {
        'model': model, 'index': res['index_name'], 'comparison': res['config'].label,
        'land_mean': float(res['anomaly'].where(land).mean()),
        'land_sig':  float(((p < 0.05) & land).sum() / land.sum()),
    }

# Land-only table, E3SM temperature rows dropped
# When E3SM data is looked into, set drop_e3sm_temp=False and rerun

REGIONS = {
    'global_land':   {},
    'NH_midlat':     {'lat_band': (25, 75)},
    'Mediterranean': {'lat_band': (30, 45), 'lon_band': (-10, 40)},
}

def contrast_table(contrast, regions=REGIONS, drop_e3sm_temp=True):
    rows = []
    for idx, by_model in contrast.items():
        for model, res in by_model.items():
            if drop_e3sm_temp and model == 'E3SM' and idx in TEMP_INDICES:
                continue
            row = {'index': idx, 'model': model}
            for name, kw in regions.items():
                st = land_stats(res, **kw)
                row[f'{name} mean'] = round(st['land_mean'], 2)
                row[f'{name} %sig'] = round(100 * st['land_sig'], 1)
            rows.append(row)
    return pd.DataFrame(rows)


# ============================================================================
# SECTION 16  -  Index families and model eligibility
# ============================================================================
# E3SM daily temperature is unusable (TREFHTMX == TREFHTMN == TREFHT at source),
# so temperature indices run on three models and precipitation on four.
TEMP_INDICES = {'TX90p', 'TX10p', 'TN90p', 'TN10p'}
PRECIP_INDICES = {'R95p', 'R99p'}

def MODELS_FOR(idx):
    return (['CESM', 'UKESM', 'MIROC', 'E3SM'] if idx in PRECIP_INDICES
            else ['CESM', 'UKESM', 'MIROC'])


# ============================================================================
# SECTION 17  -  Extra daily variables for the circulation analyses
# ============================================================================
# tas: daily-mean temperature, the global-cooling denominator in the
# drying-per-cooling analysis. psl: sea-level pressure for the NAO index.
# ua850: 850 hPa zonal wind for the jet-latitude index (CESM U on levels,
# E3SM native U850, MIROC ua on pressure levels; UKESM has no ua on levels).

def _convert_cesm_TREFHT(ds):
    da = ds['TREFHT'].astype('float32')
    da.attrs['units'] = 'K'
    return da


def _convert_cmip_tas(ds, varname='tas'):
    da = ds[varname].astype('float32')
    if 'units' not in da.attrs:
        da.attrs['units'] = 'K'
    return da


def _convert_cesm_PSL(ds):
    da = ds['PSL'].astype('float32')
    da.attrs['units'] = 'Pa'
    return da


def _convert_cmip_psl(ds, varname='psl'):
    da = ds[varname].astype('float32')
    if 'units' not in da.attrs:
        da.attrs['units'] = 'Pa'
    return da


def _sel_850(da):
    # collapse a vertical dim to 850 (hPa or Pa); pass single-level fields through
    vdim = next((d for d in ('plev', 'lev', 'level', 'pfull', 'p', 'ilev')
                 if d in da.dims), None)
    if vdim is None:
        return da
    levs = da[vdim]
    target = 85000.0 if float(levs.max()) > 2000 else 850.0
    return da.sel({vdim: target}, method='nearest')


def _convert_cesm_U(ds):
    return _sel_850(ds['U']).astype('float32')


def _convert_e3sm_U850(ds):
    return ds['U850'].astype('float32')


def _convert_miroc_ua(ds):
    return _sel_850(ds['ua']).astype('float32')


for _m, _conv in [('CESM', _convert_cesm_TREFHT), ('E3SM', _convert_cesm_TREFHT),
                  ('UKESM', _convert_cmip_tas), ('MIROC', _convert_cmip_tas)]:
    VARIABLE_INFO_BY_MODEL[_m]['tas'] = {
        'native_short': 'TREFHT' if _m in ('CESM', 'E3SM') else 'tas',
        'cmor_short': 'tas', 'converter': _conv}

for _m in ('CESM', 'E3SM'):
    VARIABLE_INFO_BY_MODEL[_m]['psl'] = {
        'native_short': 'PSL', 'cmor_short': 'psl', 'converter': _convert_cesm_PSL}
for _m in ('UKESM', 'MIROC'):
    VARIABLE_INFO_BY_MODEL[_m]['psl'] = {
        'native_short': 'psl', 'cmor_short': 'psl', 'converter': _convert_cmip_psl}

VARIABLE_INFO_BY_MODEL['CESM']['ua850'] = {
    'native_short': 'U', 'cmor_short': 'ua', 'converter': _convert_cesm_U}
VARIABLE_INFO_BY_MODEL['E3SM']['ua850'] = {
    'native_short': 'U850', 'cmor_short': 'ua', 'converter': _convert_e3sm_U850}
VARIABLE_INFO_BY_MODEL['MIROC']['ua850'] = {
    'native_short': 'ua', 'cmor_short': 'ua', 'converter': _convert_miroc_ua}

# UKESM HiLLA splits variables across streams; tas and psl live in apd/day
MODEL_BUCKET_LAYOUT[('UKESM', 'HiLLA')]['streams']['tas'] = 'apd/day'
MODEL_BUCKET_LAYOUT[('UKESM', 'HiLLA')]['streams']['psl'] = 'apd/day'

E3SM_SAI_VARMAP['tas'] = ('daily_T', 'TREFHT')


# ============================================================================
# SECTION 18  -  Mediterranean SAI-vs-HiLLA contrast bars
# ============================================================================
def plot_medit_contrast_bars(tbl, save_path=None):
    """
    Mediterranean-mean contrast (SAI - HiLLA, land, 2050-2069), days/year.
    Warm-tail trio then cold-tail trio; E3SM excluded (no usable daily
    temperature). tbl is the frame from contrast_table(CONTRAST).
    """
    WARM = ['TX90p', 'TN90p', 'WSDI']      # fewer warm extremes under subtropical injection
    COLD = ['TX10p', 'TN10p', 'CSDI']      # more cold extremes
    idx = WARM + COLD
    models = ['CESM', 'UKESM', 'MIROC']
    COL = 'Mediterranean mean'
    colors = {'CESM': '#1b5e9c', 'UKESM': '#d2691e', 'MIROC': '#2e8b57'}

    def med(index_name, model):
        row = tbl[(tbl['index'] == index_name) & (tbl['model'] == model)]
        if row.empty:
            raise KeyError(f"no row for {model} {index_name} in the table")
        return float(row[COL].iloc[0])

    vals = {m: [med(i, m) for i in idx] for m in models}
    x = np.arange(len(idx))
    w = 0.26
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    fig.patch.set_facecolor('white')
    for i, m in enumerate(models):
        ax.bar(x + (i - 1) * w, vals[m], w, label=m, color=colors[m],
               edgecolor='black', linewidth=0.4)
    allv = [v for m in models for v in vals[m]]
    lo, hi = min(allv), max(allv)
    pad = 0.18 * (hi - lo)
    ax.set_ylim(lo - pad * 2.2, hi + pad * 1.6)
    ax.axhline(0, color='black', linewidth=0.9)
    ax.axvline(len(WARM) - 0.5, color='0.4', linewidth=1.0, linestyle='--')
    ax.set_xticks(x)
    ax.set_xticklabels(idx, fontsize=12)
    ax.set_ylabel(r'Mediterranean mean,  SAI $-$ HiLLA  (days/year)', fontsize=12)
    ax.set_title(r'Mediterranean: subtropical minus high-latitude injection  (land, 2050$-$2069)',
                 fontsize=13)
    trans = ax.get_xaxis_transform()
    ax.text((len(WARM) - 1) / 2, 0.05, 'warm-tail indices\n(fewer warm extremes)',
            transform=trans, ha='center', va='bottom', fontsize=10.5, color='#8c2d04')
    ax.text(len(WARM) + (len(COLD) - 1) / 2, 0.95, 'cold-tail indices\n(more cold extremes)',
            transform=trans, ha='center', va='top', fontsize=10.5, color='#08519c')
    ax.legend(frameon=False, fontsize=11, loc='lower right', ncol=3, bbox_to_anchor=(1.0, 0.02))
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(axis='y', alpha=0.25)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=100, bbox_inches='tight')
        print(f"saved {save_path}")
    plt.show()
