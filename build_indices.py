#!/usr/bin/env python
"""
build_indices.py  -  the 30 annual ETCCDI indices for the GeoMIP G6-1.5K
multi-model archive, on the Reflective Cloud Hub or on NCAR Derecho/Casper.

    python build_indices.py --site hub --model UKESM --out-root OUT --preflight
    python build_indices.py --site hub --model UKESM --lat-range 40 45 --out-root OUT
    python build_indices.py --site hub --model UKESM --full-grid --band-rows 20 \\
        --out-root s3://.../ETCCDI_indices_annual_bootstrap --output both
    python build_indices.py --site derecho --model MIROC --out-root /glade/.../tree

THE INDICES

Twelve percentile indices, built by etccdi_bootstrap.py:

    tasmax  TX90p TX10p WSDI           tasmin  TN90p TN10p CSDI
    pr      R95p R99p (mm)  R95pTOT R99pTOT (%)  N95 N99 (days)

Eighteen fixed-threshold indices, built by multimodel_etccdi.INDEX_REGISTRY:

    tasmax  SU ID TXX TXN              tasmin  FD TN TNX TNN
    pr      PRCPTOT SDII CDD CWD R10MM R20MM RX1D RX5D
    tas     GSL GDD

Every year each run has, 2015 to 2084, is computed. A year with fewer than 360
days is written as NaN. E3SMv3 gets no tasmax or tasmin index: its daily
maximum and minimum temperature equal the daily mean in the archive.

THE PERCENTILE METHOD

The threshold is taken over the SSP2-4.5 threshold members pooled into one
sample, over the years given by --baseline (default 2020 2039, the period
whose global mean temperature the G6-1.5K runs target). 

Temperature uses a 5-day day-of-year window; precipitation uses one threshold per cell from all wet days (pr >= 1 mm/day) of the baseline (Zhang et al. 2011). 

With --bootstrap (the default) the in-base years of the threshold members are
computed by the Zhang et al. (2005) bootstrap

--no-bootstrap uses the plain threshold in every year. 

SITES

--site hub reads the Reflective bucket through multimodel_etccdi

--site derecho reads GLADE through derecho_io and derecho_layout. The index code is
the same on both. The default is the ETCCDI_SITE environment variable, else
hub.

OUTPUT

    --output per-index   one file per index and member (default):
        <out-root>/<MODEL>/<SCENARIO>/<member>/<INDEX>_<MODEL>_<SCENARIO>_<member>_<y0>-<y1>.nc
    --output merged      two files per member, the percentile group and the
                         fixed group, as the team releases the data:
        <out-root>/<MODEL>/<SCENARIO>/<member>/ETCCDI_<group>_<MODEL>_<SCENARIO>_<member>_<y0>-<y1>.nc
                         The per-index files it is built from go to
                         <out-root>_per_index.
    --output both        per-index files under <out-root>, merged files under
                         <out-root>_merged.

The percentile group is called bootstrap when the bootstrap is on and
percentile when it is off. Thresholds are written under
<out-root>/_thresholds/<MODEL>/. 
A local --out-root writes to disk; on the Hub an s3:// one writes to the bucket.

MEMORY

Each member is read whole and cut to the latitude band before the next is
read. 
--lat-range computes one band

--full-grid computes the globe as bands of --band-rows rows and joins them, with a marker per finished band so a stopped run resumes at the next one.

"""

import os

# One thread per worker.
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
os.environ.setdefault('HDF5_USE_FILE_LOCKING', 'FALSE')

import argparse                                              # noqa: E402
import datetime                                              # noqa: E402
import gc                                                    # noqa: E402
import glob                                                  # noqa: E402
import json                                                  # noqa: E402
import multiprocessing as mp                                 # noqa: E402
import re                                                    # noqa: E402
import sys                                                   # noqa: E402
import tempfile                                              # noqa: E402
import time                                                  # noqa: E402
import warnings                                              # noqa: E402

import numpy as np                                           # noqa: E402
import xarray as xr                                          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 'src'))

import etccdi_bootstrap as eb                                # noqa: E402
SHIM = eb.ensure_xclim_compat()

import multimodel_etccdi as m                                # noqa: E402
from multimodel_etccdi import Period                         # noqa: E402
import naming                                                # noqa: E402

VARIABLES = ('tasmax', 'tasmin', 'pr', 'tas')
SCENARIOS = ('SSP245', 'HiLLA', 'SAI')
ALL_YEARS = Period('all', 2015, 2084)
DEFAULT_BASELINE = (2020, 2039)

PERCENTILE_BY_VAR = {'tasmax': tuple(eb.TEMPERATURE_INDICES['tasmax']),
                     'tasmin': tuple(eb.TEMPERATURE_INDICES['tasmin']),
                     'pr': tuple(eb.PRECIP_INDICES), 'tas': ()}
FIXED_BY_VAR = {v: tuple(n for n, i in m.INDEX_REGISTRY.items()
                         if i['percentile'] is None and i['variable'] == v)
                for v in VARIABLES}
PERCENTILE_INDICES = tuple(n for v in VARIABLES for n in PERCENTILE_BY_VAR[v])
FIXED_INDICES = tuple(n for v in VARIABLES for n in FIXED_BY_VAR[v])
ALL_INDICES = PERCENTILE_INDICES + FIXED_INDICES

# Models whose daily maximum and minimum temperature cannot carry an index.
NO_TEMPERATURE_EXTREMES = {
    'E3SM': ('E3SMv3 daily maximum and minimum near-surface air temperature '
             '(TREFHTMX, TREFHTMN) in this archive equal the daily mean '
             '(TREFHT), so no tasmax or tasmin index is computed for E3SMv3.'),
}


# --------------------------------------------------------------------------
# units and names of the percentile indices
# --------------------------------------------------------------------------

# ---- copied from build_indices_derecho_boot.py, 14 October 2026 -----------
# Day counts in 'd', not 'days', since 24 October 2026 on both sites
# (patch_boot_sources.py, repair_boot_attrs.py): xarray 2025.1 decodes units
# 'days' to timedelta64 on open.
UNITS = {'TX90p': 'd', 'TX10p': 'd', 'TN90p': 'd', 'TN10p': 'd',
         'WSDI': 'd', 'CSDI': 'd', 'N95': 'd', 'N99': 'd',
         'R95p': 'mm', 'R99p': 'mm', 'R95pTOT': '%', 'R99pTOT': '%'}
RANGE = {'d': (0.0, 366.0), 'mm': (0.0, 50000.0), '%': (0.0, 100.0)}

LONG_NAME = {
    'TX90p': 'Days with tasmax above the 90th percentile',
    'TX10p': 'Days with tasmax below the 10th percentile',
    'TN90p': 'Days with tasmin above the 90th percentile',
    'TN10p': 'Days with tasmin below the 10th percentile',
    'WSDI': 'Days in spells of at least 6 days with tasmax above the 90th '
            'percentile',
    'CSDI': 'Days in spells of at least 6 days with tasmin below the 10th '
            'percentile',
    'R95p': 'Annual total precipitation on wet days above the 95th percentile',
    'R99p': 'Annual total precipitation on wet days above the 99th percentile',
    'R95pTOT': 'Contribution of days above the 95th percentile to annual wet-'
               'day precipitation, 100 * R95p / PRCPTOT',
    'R99pTOT': 'Contribution of days above the 99th percentile to annual wet-'
               'day precipitation, 100 * R99p / PRCPTOT',
    'N95': 'Number of wet days above the 95th percentile',
    'N99': 'Number of wet days above the 99th percentile',
}

PRECIP_NOTE = (
    'In this tree R95p and R99p are annual totals in mm (Zhang et al. 2011, '
    'Table 1; climpact). In ETCCDI_indices_annual, R95p and R99p were counts '
    'of days against a day-of-year wet-day threshold; the nearest quantity '
    'here is N95 or N99, against an all-wet-day threshold.')
# ---------------------------------------------------------------------------



def say(msg=''):
    print(f'{time.strftime("%H:%M:%S")} {msg}', flush=True)


# --------------------------------------------------------------------------
# sites: where the daily fields come from, and what the members are called
# --------------------------------------------------------------------------

class HubSite:
    """
    Compute on the Reflective bucket through multimodel_etccdi.
    """
    name = 'hub'

    def members(self, model, scenario):
        return tuple(m.get_layout(model, scenario)['members'])

    def baseline_members(self, model):
        return tuple(m.baseline_members(model, 'SSP245'))

    def label(self, model, scenario, member):
        return m.MEMBER_LABEL.get((model, scenario), {}).get(member, member)

    def read(self, model, scenario, var, member, period, lat_range):
        """
        One member, one variable, cut to the band, on a cftime axis, in
        mm/day or K. load_scenario_members applies the module's converter;
        the readers cut the band before loading when READ_LAT_RANGE is set.
        E3SM is assembled from the regridded merged windows on the bucket
        (read_e3sm_band).
        """
        if model == 'E3SM':
            return read_e3sm_band(scenario, var, member, period, lat_range)
        got = m.load_scenario_members(model, var, scenario, period,
                                      members=[member])
        full = got.pop(member, None)
        del got
        if full is None:
            return None
        return m._to_cftime(cut_band(full, lat_range))

    def extra_attrs(self, model, scenario):
        return {}


class DerechoSite:
    """
    Compute on GLADE through derecho_io and derecho_layout.
    """
    name = 'derecho'

    def __init__(self):
        import derecho_io
        import derecho_layout
        self.dio, self.dl = derecho_io, derecho_layout

    def members(self, model, scenario):
        # Every SSP2-4.5 member on disk; for the G6 runs those whose record
        # reaches the end of the window, as derecho_layout records it.
        if scenario == 'SSP245':
            return tuple(self.dl.members(model, scenario))
        return tuple(self.dl.members(model, scenario,
                                     through_year=ALL_YEARS.end_year))

    def baseline_members(self, model):
        return tuple(self.dl.baseline_members(model, 'SSP245'))

    def label(self, model, scenario, member):
        lab = self.dl.out_label(model, scenario, member)
        return m.MEMBER_LABEL.get((model, scenario), {}).get(lab, lab)

    def read(self, model, scenario, var, member, period, lat_range):
        da = self.dio.read_member(model, scenario, var, member, period)
        if da is None:
            return None
        return m._to_cftime(cut_band(da, lat_range))

    def extra_attrs(self, model, scenario):
        out = {}
        note = getattr(self.dl, 'VERSION_NOTE', {}).get((model, scenario))
        if note:
            out['version_note'] = note
        note = getattr(self.dl, 'VALIDATION_NOTE', {}).get((model, scenario))
        if note:
            out['validation_note'] = note
        return out


SITE = None


def member_label(model, scenario, member):
    return SITE.label(model, scenario, member)


def scen_label(model, scenario):
    return naming.out_scenario(scenario)


def scenario_members(model, scenario, restrict=None):
    if restrict:
        return tuple(restrict)
    return SITE.members(model, scenario)


def cut_band(da, lat_range):
    """
    da on (time, lat, lon), restricted to the latitudes in lat_range by
    value and copied, so the full-grid read can be freed.
    """
    da = da.transpose('time', 'lat', 'lon')
    if lat_range is None:
        return da
    r0, r1 = lat_rows(da['lat'], *lat_range)
    out = da.isel(lat=slice(r0, r1)).copy(deep=True)
    del da
    gc.collect()
    return out


# --------------------------------------------------------------------------
# names and paths
# --------------------------------------------------------------------------

def member_dir(root, model, scenario, member):
    return (f'{root.rstrip("/")}/{naming.out_model(model)}/'
            f'{naming.out_scenario(scenario)}/'
            f'{member_label(model, scenario, member)}')


def out_path(root, model, scenario, member, name, span):
    return (f'{member_dir(root, model, scenario, member)}/'
            + naming.index_filename(name, model, scenario,
                                    member_label(model, scenario, member),
                                    years=span))


def threshold_path(root, model, var, convention):
    return (f'{root.rstrip("/")}/_thresholds/{naming.out_model(model)}/'
            f'{var}_threshold_{convention}.nc')


def output_roots(args):
    """
    (per-index root, merged root or None) for --output.
    """
    root = args.out_root.rstrip('/')
    if args.output == 'per-index':
        return root, None
    if args.output == 'merged':
        return root + '_per_index', root
    return root, root + '_merged'


# --------------------------------------------------------------------------
# files, local or on the bucket
# --------------------------------------------------------------------------

def is_s3(path):
    return str(path).startswith('s3://')


def find(pattern):
    """
    Glob locally or on the bucket; s3 results come back with s3://.
    """
    if is_s3(pattern):
        return sorted('s3://' + p for p in m.s3.glob(pattern[5:]))
    return sorted(glob.glob(pattern))


def exists(path):
    return m.s3.exists(path[5:]) if is_s3(path) else os.path.exists(path)


def open_any(path):
    """
    decode_timedelta=False: every day-count index here carries units
    'days', and xarray 2025.1 (the Hub's) decodes such a variable to
    timedelta64 on open without a warning, so 11 days reads back as
    950400000000000 ns. Found in the smoke test of 24 October 2026.
    """
    if is_s3(path):
        with m.s3.open(path[5:], 'rb') as f:
            return xr.open_dataset(f, engine='h5netcdf',
                                   decode_timedelta=False).load()
    return xr.open_dataset(path, decode_timedelta=False).load()


def write_any(ds, path):
    if is_s3(path):
        m.save_nc_to_s3(ds, path)
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ds.to_netcdf(path)


def lat_rows(lat, lo, hi):
    """
    First and stop row of the latitudes in [lo, hi]. By value, so a
    descending axis (MIROC-ES2H) gives the same band as an ascending one.
    """
    v = np.asarray(lat, dtype='float64')
    idx = np.flatnonzero((v >= lo) & (v <= hi))
    if idx.size == 0:
        raise SystemExit(f'no latitude in [{lo}, {hi}]; the grid runs '
                         f'{v.min():.2f} to {v.max():.2f}')
    return int(idx.min()), int(idx.max()) + 1



# --------------------------------------------------------------------------
# E3SM on the Hub
# --------------------------------------------------------------------------

# ---- E3SM: the regridded merged windows already on the bucket -------------
#
# The module's E3SM readers take one window per call and, when the merged file
# for that exact window is missing, regrid the raw cubed-sphere output on the
# pod. A full 2015-2084 window does not exist and would not fit (about 6.6 GB
# a member). So the driver assembles each member from the windows that do
# exist (fafrifa/E3SMv3_regridded/<case>/merged_PRECT_<y0>_<y1>.nc), day by
# day: each day from the earliest window holding it, band only.
#
# Time stamps: EAM, like CAM, stamps a daily mean at the end of its day:
# 2020-01-02 00:00 is the mean of 1 January. The module corrects CESM with the
# time_bnds midpoint (_fix_cesm_time); the merged E3SM files kept only the data
# variable, so no bounds survive and nothing corrected them. The merge also cut
# each window by stamp year, so every window lacks its first stamp (y0-01-01,
# in the previous raw chunk) and drops its last day (stamped y1+1-01-01).
# Found by the overlap check, 27 October 2026: 2050 had 365 steps in
# merged_PRECT_2040_2059 and 364 in merged_PRECT_2050_2069.
#
# So every stamp is moved back E3SM_STAMP_SHIFT, to the midpoint of its day as
# _fix_cesm_time places CESM, and the day a window drops at its end comes from
# the next window. --preflight reads time_bnds from one raw file to confirm the
# stamps mark the end of the day. Added 27 October 2026.

E3SM_STAMP_SHIFT = datetime.timedelta(hours=12)
E3SM_WINDOW_RE = re.compile(r'merged_(?P<var>[A-Z]+)_(?P<y0>\d{4})_(?P<y1>\d{4})'
                            r'\.nc$')
_E3SM_TIMES = {}      # file -> corrected time stamps, one read per process
_E3SM_SAID = set()    # coverage printed once per member and window


def e3sm_case(scenario, member):
    """
    Case directory under E3SM_MERGED_S3_PREFIX, as the module names it.
    """
    return f'G6-1.5K-SAI_{member}' if scenario == 'SAI' else member


def _e3sm_open(path):
    """
    (dataset, file object), lazy. The caller closes both.
    """
    f = m.s3.open(path[5:], 'rb') if is_s3(path) else open(path, 'rb')
    ds = xr.open_dataset(f, engine='h5netcdf', decode_timedelta=False)
    return ds, f


def _e3sm_times(path):
    """
    The file's time stamps moved to the middle of their day, as an
    object array of cftime dates.
    """
    if path not in _E3SM_TIMES:
        ds, f = _e3sm_open(path)
        try:
            t = ds['time'].values
        finally:
            ds.close()
            f.close()
        if t.size and isinstance(t[0], np.datetime64):
            raise RuntimeError(f'{path}: time decoded to datetime64; E3SM '
                               f'is noleap and should decode to cftime')
        _E3SM_TIMES[path] = np.array([x - E3SM_STAMP_SHIFT for x in t],
                                     dtype=object)
    return _E3SM_TIMES[path]


def _day(x):
    return (x.year, x.month, x.day)


def e3sm_windows(scenario, member, native):
    prefix = m.E3SM_MERGED_S3_PREFIX.rstrip('/')
    files = [p for p in find(f'{prefix}/{e3sm_case(scenario, member)}/'
                             f'merged_{native}_*.nc')
             if E3SM_WINDOW_RE.search(p)]
    files.sort(key=lambda p: (int(E3SM_WINDOW_RE.search(p)['y0']),
                              int(E3SM_WINDOW_RE.search(p)['y1'])))
    return files


def e3sm_plan(scenario, member, native, period):
    """
    [(file, positions)] with positions the indices into that file's time
    axis to read, and {year: days} of the assembled record. Each day of
    period comes from the earliest window holding it.
    """
    taken, plan, days = set(), [], {}
    for p in e3sm_windows(scenario, member, native):
        t = _e3sm_times(p)
        pos = []
        for i, x in enumerate(t):
            if not (period.start_year <= x.year <= period.end_year):
                continue
            k = _day(x)
            if k in taken:
                continue
            taken.add(k)
            pos.append(i)
            days[x.year] = days.get(x.year, 0) + 1
        if pos:
            plan.append((p, np.asarray(pos)))
    return plan, dict(sorted(days.items()))


def e3sm_coverage(scenario, member, plan, days):
    """
    One line per window used, then the years short of 365 days.
    """
    lab = member_label('E3SM', scenario, member)
    for p, pos in plan:
        t = _e3sm_times(p)[pos]
        say(f'  E3SM {scen_label("E3SM", scenario)} {lab}: {len(pos)} days, '
            f'{t[0].strftime("%Y-%m-%d")} to {t[-1].strftime("%Y-%m-%d")} '
            f'from {p.rsplit("/", 1)[-1]}')
    for y, n in days.items():
        if n < 365:
            say(f'  E3SM {lab}: {y} has {n} days, '
                + ('masked in the output' if n < 360 else
                   'kept (the rule is 360)'))


def _runs(pos):
    """
    Consecutive runs of sorted positions, as slices.
    """
    cut = np.flatnonzero(np.diff(pos) != 1) + 1
    return [slice(int(r[0]), int(r[-1]) + 1) for r in np.split(pos, cut)]


def read_e3sm_band(scenario, var, member, period, lat_range):
    info = m.get_variable_info('E3SM', var)
    native = info['native_short']
    plan, days = e3sm_plan(scenario, member, native, period)
    if not plan:
        return None
    key = (scenario, member, period.start_year, period.end_year)
    first = key not in _E3SM_SAID
    if first:
        _E3SM_SAID.add(key)
        e3sm_coverage(scenario, member, plan, days)
    parts = []
    for p, pos in plan:
        ds, f = _e3sm_open(p)
        try:
            sub = ds[[native]]
            if lat_range is not None:
                r0, r1 = lat_rows(sub['lat'], *lat_range)
                sub = sub.isel(lat=slice(r0, r1))
            for sl in _runs(pos):
                piece = sub.isel(time=sl).load()
                parts.append(piece.assign_coords(
                    time=('time', _e3sm_times(p)[sl], piece['time'].attrs)))
        finally:
            ds.close()
            f.close()
    ds = xr.concat(parts, dim='time') if len(parts) > 1 else parts[0]
    del parts
    t = ds['time'].values
    if not all(a < b for a, b in zip(t[:-1], t[1:])):
        raise RuntimeError(f'E3SM {scenario}/{member}: time axis not strictly '
                           f'increasing after assembly')
    if first:
        say(f'  E3SM {member_label("E3SM", scenario, member)}: first day '
            f'{t[0]}, last day {t[-1]} (stamps moved back '
            f'{E3SM_STAMP_SHIFT})')
    da = info['converter'](ds).transpose('time', 'lat', 'lon')
    da.attrs['src_e3sm_merged_files'] = ', '.join(p.rsplit('/', 1)[-1]
                                                  for p, _ in plan)
    da.attrs['src_time_correction'] = (
        f'time stamps moved back {E3SM_STAMP_SHIFT} to the middle of their '
        f'day: EAM stamps a daily mean at the end of the day and the merged '
        f'files carry no time_bnds')
    return m._to_cftime(da)



# --------------------------------------------------------------------------
# checks on the inputs
# --------------------------------------------------------------------------

def check_grid(base, targets):
    """
    Every series on the threshold members' grid. The block compute cuts by
    position and then labels everything with the first threshold member's
    coordinates, so a target on a shifted grid would be mislabelled rather
    than rejected. Added here; the Derecho driver assumes it.
    """
    ref = next(iter(base.values()))
    for k, da in list(base.items()) + list(targets.items()):
        for c in ('lat', 'lon'):
            a = np.asarray(ref[c].values, dtype='float64')
            b = np.asarray(da[c].values, dtype='float64')
            if a.shape != b.shape or not np.allclose(a, b, rtol=0,
                                                     atol=1e-4):
                raise RuntimeError(f'{k}: {c} differs from the threshold '
                                   f'members ({b.size} against {a.size} '
                                   f'points). Nothing computed.')


def _year_days(year, calendar):
    """
    Every day of one year on the calendar. Up to 1 January of the next
    year, not 31 December, which a 360_day calendar does not have.
    """
    return xr.date_range(f'{year}-01-01', f'{year + 1}-01-01', freq='D',
                         calendar=calendar, use_cftime=True,
                         inclusive='left')


def check_base(base, base_years, var):
    """
    Every threshold member must cover every base year completely.

    Complete means every day of the calendar year, not the 360 days the
    output mask asks for: the bootstrap (xclim build_bootstrap_year_da) puts
    one base year's values in place of another's, which fails on a year a
    day short. Tightened 27 October 2026, when E3SM baseline years came out
    of the merged windows one day short; until then the Derecho driver's
    360-day test was used, and every Derecho model passed the stricter one
    too, since their bootstraps ran.
    """
    want = set(range(base_years[0], base_years[1] + 1))
    for mem, da in base.items():
        t = da['time'].values
        cal = getattr(t[0], 'calendar', None) or 'standard'
        cnt = da['time'].dt.year.to_series().value_counts()
        short = {}
        for y in sorted(want):
            n, full = int(cnt.get(y, 0)), len(_year_days(y, cal))
            if n != full:
                short[y] = (n, full)
        if short:
            have = {(x.year, x.month, x.day) for x in t}
            lines = []
            for y, (n, full) in short.items():
                gone = [d.strftime('%m-%d') for d in _year_days(y, cal)
                        if (d.year, d.month, d.day) not in have]
                lines.append(f'{y}: {n} of {full} days, missing '
                             f'{", ".join(gone[:6])}'
                             + (' ...' if len(gone) > 6 else ''))
            raise RuntimeError(
                f'{var} baseline member {mem} is incomplete ({cal}): '
                + '; '.join(lines) + '. Nothing computed.')
        if var == 'pr' and da.attrs.get('units') != 'mm/day':
            raise RuntimeError(f'pr member {mem} arrived in '
                               f'{da.attrs.get("units")!r}; the wet-day test '
                               f'is 1 mm/day and needs mm/day')
        dims = tuple(da.dims)
        if dims != ('time', 'lat', 'lon'):
            raise RuntimeError(f'{var} member {mem} has dims {dims}; this '
                               f'driver expects (time, lat, lon)')



# --------------------------------------------------------------------------
# the percentile compute, in blocks of cells over forked workers
# --------------------------------------------------------------------------

_G = {}


def blocks(nlat, nlon, cells, rows=None):
    r0, r1 = rows if rows else (0, nlat)
    per_row = max(1, int(np.ceil(nlon / max(1, cells))))
    width = int(np.ceil(nlon / per_row))
    out = []
    for i in range(r0, r1):
        for j in range(0, nlon, width):
            out.append((i, i + 1, j, min(nlon, j + width)))
    return out


def _cut(da, blk):
    i0, i1, j0, j1 = blk
    return da.isel(lat=slice(i0, i1), lon=slice(j0, j1))


def _work(blk):
    warnings.simplefilter('ignore')
    t0 = time.time()
    g = _G
    base = {k: _cut(v, blk) for k, v in g['base'].items()}
    targets = {k: _cut(v, blk) for k, v in g['targets'].items()}
    if g['var'] == 'pr':
        res, thr = eb.precipitation_indices(
            targets, base, g['inbase'], g['base_years'],
            convention=g['convention'])
    else:
        res, thr = eb.temperature_indices(
            g['var'], targets, base, g['inbase'], g['base_years'],
            convention=g['convention'])
    out = {}
    for tkey, d in res.items():
        for name, da in d.items():
            extra = [x for x in da.dims if x not in ('time', 'lat', 'lon')]
            if extra:
                raise RuntimeError(f'{tkey}/{name} has extra dims {extra}')
            da = da.transpose('time', 'lat', 'lon')
            out[(tkey, name)] = (da['time'].values,
                                 np.asarray(da.values, dtype='float32'))
    thr = thr.transpose(..., 'lat', 'lon')
    return blk, out, np.asarray(thr.values, dtype='float32'), \
        {x: thr[x].values for x in thr.dims if x not in ('lat', 'lon')}, \
        time.time() - t0


def compute(var, base, targets, inbase, base_years, convention, workers,
            cells, rows=None):
    """
    ({(tkey, name): annual DataArray on the band}, threshold).
    """
    like = next(iter(base.values()))
    nlat, nlon = like.sizes['lat'], like.sizes['lon']
    blist = blocks(nlat, nlon, cells, rows)
    _G.clear()
    _G.update(var=var, base=base, targets=targets, inbase=inbase,
              base_years=base_years, convention=convention)

    say(f'  {len(blist)} block(s) of up to {blist[0][3] - blist[0][2]} cells, '
        f'{workers} worker(s), {len(targets)} target series, '
        f'{sum(1 for k in targets if k in inbase)} bootstrapped')
    fields, thr_full, thr_coords, times = {}, None, None, {}
    t0, done, slowest = time.time(), 0, 0.0
    ctx = mp.get_context('fork')
    with ctx.Pool(workers, maxtasksperchild=50) as pool:
        for blk, out, thr, tcoords, dt in pool.imap_unordered(_work, blist):
            i0, i1, j0, j1 = blk
            for key, (tt, vals) in out.items():
                if key not in fields:
                    times[key] = tt
                    fields[key] = np.full((len(tt), nlat, nlon), np.nan,
                                          dtype='float32')
                fields[key][:, i0:i1, j0:j1] = vals
            if thr_full is None:
                thr_coords = tcoords
                thr_full = np.full(thr.shape[:-2] + (nlat, nlon), np.nan,
                                   dtype='float32')
            thr_full[..., i0:i1, j0:j1] = thr
            done += 1
            slowest = max(slowest, dt)
            if done == 1 or done % max(1, len(blist) // 20) == 0 \
                    or done == len(blist):
                el = time.time() - t0
                say(f'    {done}/{len(blist)} blocks, {el / 60:.1f} min, '
                    f'slowest block {slowest:.0f} s, '
                    f'about {el / done * (len(blist) - done) / 60:.0f} min '
                    f'left')
    _G.clear()

    lat, lon = like['lat'], like['lon']
    das = {}
    for key, arr in fields.items():
        das[key] = xr.DataArray(arr, dims=('time', 'lat', 'lon'),
                                coords={'time': times[key], 'lat': lat,
                                        'lon': lon})
    tdims = tuple(thr_coords) + ('lat', 'lon')
    thr_da = xr.DataArray(thr_full, dims=tdims,
                          coords={**thr_coords, 'lat': lat, 'lon': lon})
    return das, thr_da


def check_values(da, name, where):
    """
    check_written, applied to the array before it is written, since the
    file may be on the bucket.
    """
    v = np.asarray(da.values, dtype='float64')
    fin = np.isfinite(v)
    if not fin.any():
        raise ValueError(f'{where}: every value is NaN')
    lo, hi = RANGE[UNITS[name]]
    vals = v[fin]
    tol = 1e-4 * max(1.0, abs(hi))
    if vals.min() < lo - tol or vals.max() > hi + tol:
        raise ValueError(f'{where}: range {vals.min():.3f} to '
                         f'{vals.max():.3f} outside {lo} to {hi} '
                         f'{UNITS[name]}')
    return f'{vals.min():.2f} to {vals.max():.2f}, mean {vals.mean():.2f}'



# --------------------------------------------------------------------------
# attributes
# --------------------------------------------------------------------------

def common_attrs(model, scenario, member, name, units, long_name, extra_src):
    return {'source_model': naming.out_model(model),
            'scenario': naming.out_scenario(scenario),
            'member': member_label(model, scenario, member),
            'source_member_id': str(member),
            'index': naming.canon_index(name),
            'units': units,
            'long_name': long_name,
            'created': time.strftime('%Y-%m-%d'),
            'producer': f'build_indices.py (site {SITE.name})',
            **SITE.extra_attrs(model, scenario),
            **extra_src}


def percentile_attrs(model, scenario, member, name, var, args, base_members,
                     inbase, extra_src):
    b = args.baseline
    a = {'percentile_baseline_scenario': 'SSP2-4.5',
         'percentile_baseline_period': f'{b.start_year}-{b.end_year}',
         'percentile_baseline_years': b.end_year - b.start_year + 1,
         'percentile_convention_changed': m.PERCENTILE_CONVENTION_CHANGED}
    if model in m.BASELINE_NOTE:
        a['baseline_caveat'] = m.BASELINE_NOTE[model]
    a.update(common_attrs(model, scenario, member, name, UNITS[name],
                          LONG_NAME[name], extra_src))
    a['baseline_members'] = ','.join(str(x) for x in base_members)
    a['threshold_convention'] = args.convention
    a['producer'] += ' with etccdi_bootstrap.py'
    in_threshold = (scenario == 'SSP245' and
                    str(member) in {str(x) for x in base_members})
    if not args.bootstrap and in_threshold:
        a['in_base_bootstrap'] = ('not applied (--no-bootstrap): every year '
                                  'uses the plain threshold')
    elif inbase:
        a['in_base_bootstrap'] = (
            f'yes, years {b.start_year}-{b.end_year} of this member are '
            f'bootstrapped (Zhang et al. 2005, J. Climate 18, 1641)')
    else:
        a['in_base_bootstrap'] = ('not applicable: this member did not enter '
                                  'the threshold, so every year uses the '
                                  'plain threshold')
    pooled = (args.convention != 'member_mean')
    if var == 'pr':
        a['percentile_method'] = (
            'one threshold per cell from all base-period days with pr >= 1 '
            'mm, no day-of-year window (Zhang et al. 2011, p. 855), type 8 '
            'quantile, ' + ('all members pooled' if pooled else
                            'per member then averaged over members'))
        a['wet_day_threshold'] = '1 mm/day'
        a['index_definition_note'] = PRECIP_NOTE
    else:
        a['percentile_method'] = (
            'xclim percentile_doy, 5-day centred window, type 8, ' +
            ('all members pooled into one sample, member j re-stamped '
             'j * 400 years (multimodel_etccdi.pooled_percentile_doy)'
             if pooled else 'per member then averaged over members'))
    return a


def fixed_attrs(model, scenario, member, name, da, extra_src):
    info = m.INDEX_REGISTRY[name]

    def fill(text):
        # Registry text may name its own threshold, e.g. 'above {thresh}'.
        try:
            return str(text).format(**info['extra_kwargs'])
        except (KeyError, IndexError, ValueError):
            return str(text)
    reg = {k: fill(v) for k, v in info.get('attrs', {}).items()}
    a = common_attrs(model, scenario, member, name,
                     da.attrs.get('units', naming.INDEX_UNITS.get(name, '')),
                     reg.get('long_name', da.attrs.get('long_name', name)),
                     extra_src)
    a['threshold_type'] = 'fixed'
    thresh = info['extra_kwargs'].get('thresh')
    if thresh is not None:
        a['fixed_threshold'] = str(thresh)
    if reg.get('comment'):
        a['comment'] = reg['comment']
    a['partial_years_masked'] = 'years with fewer than 360 days are NaN'
    import xclim
    a['xclim_version'] = xclim.__version__
    a['producer'] += ' with multimodel_etccdi.INDEX_REGISTRY'
    return a


def check_fixed(da, name, where):
    v = np.asarray(da.values, dtype='float64')
    fin = np.isfinite(v)
    if not fin.any():
        raise ValueError(f'{where}: every value is NaN')
    lo, hi = naming.index_range(name)
    vals = v[fin]
    if vals.min() < lo or vals.max() > hi:
        raise ValueError(f'{where}: range {vals.min():.3f} to {vals.max():.3f} '
                         f'outside {lo:g} to {hi:g}')
    return f'{vals.min():.2f} to {vals.max():.2f}, mean {vals.mean():.2f}'


def write_index(da, attrs, path, name, overwrite):
    if exists(path) and not overwrite:
        say(f'  exists, skipping: {path}')
        return False
    da = da.rename(naming.canon_index(name))
    da.attrs = {'units': attrs['units'], 'long_name': attrs['long_name']}
    ds = da.to_dataset()
    ds.attrs.update(attrs)
    write_any(ds, path)
    return True


# --------------------------------------------------------------------------
# one variable, one model
# --------------------------------------------------------------------------

def run_variable(args, model, var, base_members, out_root, written):
    """
    Every requested index on one variable, for every scenario. written
    collects (path, scenario, member, index) for each output.
    """
    pct = [n for n in PERCENTILE_BY_VAR[var] if n in args.indices]
    fix = [n for n in FIXED_BY_VAR[var] if n in args.indices]
    if model in NO_TEMPERATURE_EXTREMES and var in ('tasmax', 'tasmin'):
        if pct or fix:
            say(f'{model} {var}: skipped. {NO_TEMPERATURE_EXTREMES[model]}')
        return
    if not (pct or fix):
        return
    b = args.baseline
    base_years = (b.start_year, b.end_year)
    say(f'\n{"=" * 72}\n{model} {var}   percentile {pct or "none"}   '
        f'fixed {fix or "none"}\n{"=" * 72}')

    base = {}
    if pct:
        say(f'threshold: SSP2-4.5 {base_years[0]}-{base_years[1]}, members '
            f'{", ".join(base_members)}, '
            f'{"bootstrap on" if args.bootstrap else "no bootstrap"}')
        for mem in base_members:
            da = SITE.read(model, 'SSP245', var, mem, b, args.lat_range)
            if da is None:
                raise RuntimeError(f'threshold member {mem} unreadable for '
                                   f'{var}')
            base[mem] = da
        check_base(base, base_years, var)

    thr_written = False
    for scen in args.scenarios:
        t0 = time.time()
        mems = scenario_members(model, scen, args.members)
        say(f'\n--- {model} {scen} {var}, {ALL_YEARS.start_year}-'
            f'{ALL_YEARS.end_year}, members {", ".join(mems)}')
        targets = {}
        for mem in mems:
            da = SITE.read(model, scen, var, mem, ALL_YEARS, args.lat_range)
            if da is None:
                say(f'  SKIP {mem}: no data')
                continue
            targets[mem] = da
        if not targets:
            say(f'  no members read for {scen}; nothing to compute')
            continue
        check_grid(base or {next(iter(targets)): next(iter(targets.values()))},
                   targets)

        if pct:
            inbase = ({k: k for k in targets if k in base_members}
                      if scen == 'SSP245' and args.bootstrap else {})
            fields, thr = compute(var, base, targets, inbase, base_years,
                                  args.convention, args.workers,
                                  args.cells_per_task)
            if not thr_written:
                tp = threshold_path(out_root, model, var, args.convention)
                thr.attrs.update({
                    'baseline_members': ','.join(base_members),
                    'baseline_period': f'{base_years[0]}-{base_years[1]}',
                    'threshold_convention': args.convention,
                    'units': 'mm/day' if var == 'pr' else
                    next(iter(base.values())).attrs.get('units', '')})
                if args.lat_range:
                    thr.attrs['partial_grid'] = (
                        f'latitudes {args.lat_range[0]} to '
                        f'{args.lat_range[1]} only')
                write_any(thr.rename('threshold').to_dataset(), tp)
                say(f'  threshold -> {tp}')
                thr_written = True
            for mem, src in targets.items():
                complete = m._complete_bins(src, 'YS', min_days=360)
                extra_src = {k: v for k, v in src.attrs.items()
                             if k.startswith('src_')}
                for name in pct:
                    da = fields[(mem, name)]
                    complete_al = complete.reindex(time=da['time'])
                    if bool(complete_al.isnull().any()):
                        raise RuntimeError(f'{scen}/{mem}/{name}: annual time '
                                           f'axis does not match the '
                                           f'completeness mask')
                    da = da.where(complete_al.astype(bool))
                    yy = da['time'].dt.year.values
                    p = out_path(out_root, model, scen, mem, name,
                                 (int(yy.min()), int(yy.max())))
                    written.append((p, scen, mem, name))
                    a = percentile_attrs(model, scen, mem, name, var, args,
                                         base_members, mem in inbase,
                                         extra_src)
                    if args.lat_range:
                        a['partial_grid'] = (f'latitudes {args.lat_range[0]} '
                                             f'to {args.lat_range[1]} only')
                    rng = check_values(da, name, p)
                    if write_index(da, a, p, name, args.overwrite):
                        say(f'  wrote {p.rsplit("/", 1)[-1]}   {rng} '
                            f'{UNITS[name]}')
            del fields

        for mem, src in targets.items():
            extra_src = {k: v for k, v in src.attrs.items()
                         if k.startswith('src_')}
            for name in fix:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    da = m.compute_annual_index_for_member(
                        src, name, None, min_days=360).load()
                # Source-side scalar coordinates (UKESM's forecast_reference_time)
                # would stop the members of a model concatenating.
                da = da.reset_coords(drop=True)
                # The module squeezes singleton dimensions, so a band one
                # latitude row wide comes back without lat; put it back.
                for d in ('lat', 'lon'):
                    if d not in da.dims:
                        da = da.expand_dims({d: src[d].values})
                da = da.transpose('time', 'lat', 'lon')
                yy = da['time'].dt.year.values
                p = out_path(out_root, model, scen, mem, name,
                             (int(yy.min()), int(yy.max())))
                written.append((p, scen, mem, name))
                a = fixed_attrs(model, scen, mem, name, da, extra_src)
                if args.lat_range:
                    a['partial_grid'] = (f'latitudes {args.lat_range[0]} to '
                                         f'{args.lat_range[1]} only')
                rng = check_fixed(da, name, p)
                if write_index(da, a, p, name, args.overwrite):
                    say(f'  wrote {p.rsplit("/", 1)[-1]}   {rng} {a["units"]}')
                del da
                gc.collect()
        del targets
        gc.collect()
        say(f'--- {scen} {var} done in {(time.time() - t0) / 60:.1f} min')
    del base
    gc.collect()


# --------------------------------------------------------------------------
# --output merged / both: two files per member
# --------------------------------------------------------------------------
#
# The team releases the data as two files per member, one holding the
# percentile indices and one the fixed ones. The rule for attributes is the one
# merge_indices.py used for the released files: an attribute identical across
# every index in the file is global, one that differs or that only some carry
# is on the variable. Every merged file is read back and compared with the
# per-index files it came from before it is kept.

VAR_OF = {**{n: v for v in VARIABLES for n in PERCENTILE_BY_VAR[v]},
          **{n: m.INDEX_REGISTRY[n]['variable'] for n in FIXED_INDICES}}


def group_note(group, baseline):
    b = f'{baseline.start_year}-{baseline.end_year}'
    precip = ('Precipitation indices follow Zhang et al. (2011): one '
              'threshold per cell from all base-period wet days (>= 1 mm). '
              'R95p and R99p are annual totals in mm, R95pTOT and R99pTOT the '
              'percent of PRCPTOT, N95 and N99 day counts.')
    if group == 'bootstrap':
        return ('The twelve percentile-based indices with the Zhang et al. '
                '(2005) bootstrap applied to the in-base years of the '
                'SSP2-4.5 members that build the threshold. Threshold '
                f'reference: SSP2-4.5 {b}, pooled across members. Temperature '
                'indices use a 5-day day-of-year window. ' + precip +
                f' G6 years inside {b} are not bootstrapped, '
                'because no G6 data enter the threshold.')
    if group == 'percentile':
        return ('The twelve percentile-based indices, without the in-base '
                f'bootstrap. Threshold reference: SSP2-4.5 {b}, pooled across '
                'members. Temperature indices use a 5-day day-of-year window. '
                + precip)
    return ('The eighteen fixed-threshold indices. Most use an absolute '
            'threshold (0 C, 25 C, 20 C, 1 mm, 10 mm, 20 mm); TXX, TXN, TNX '
            'and TNN and the two precipitation depths RX1D and RX5D use no '
            'threshold at all, and PRCPTOT and SDII use a wet-day definition. '
            'Each variable carries what its own source file recorded.')


def merge_groups(args):
    pct = 'bootstrap' if args.bootstrap else 'percentile'
    return {pct: PERCENTILE_INDICES, 'fixed': FIXED_INDICES}


def _eq(a, b):
    try:
        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            return np.array_equal(np.asarray(a), np.asarray(b))
        return bool(a == b)
    except Exception:                                         # noqa: BLE001
        return False


def partition_attrs(per_index):
    names = list(per_index)
    keys = set()
    for a in per_index.values():
        keys |= set(a)
    glob_, var = {}, {n: {} for n in names}
    for k in sorted(keys):
        vals = [per_index[n].get(k, KeyError) for n in names]
        if all(v is not KeyError for v in vals) and \
                all(_eq(v, vals[0]) for v in vals):
            glob_[k] = vals[0]
        else:
            for n in names:
                if k in per_index[n]:
                    var[n][k] = per_index[n][k]
    return glob_, var


def merge_member(model, scen, mem, group, names, per_root, merged_root, args):
    label = member_label(model, scen, mem)
    tag = f'{scen:6s} {label:5s} {group}'
    found, missing, absent = {}, [], []
    for n in names:
        if model in NO_TEMPERATURE_EXTREMES and VAR_OF[n] in ('tasmax',
                                                              'tasmin'):
            absent.append(n)
            continue
        stem = naming.index_filename(n, model, scen, label, with_years=False)
        hits = find(f'{member_dir(per_root, model, scen, mem)}/'
                    f'{stem[:-3]}_*.nc')
        if len(hits) == 1:
            found[n] = hits[0]
        else:
            missing.append(n if not hits else f'{n} ({len(hits)} files)')
    if missing:
        say(f'  SKIP {tag}: no single per-index file for {", ".join(missing)}')
        return None
    arrays, src = {}, {}
    for n, p in found.items():
        ds = open_any(p)
        arrays[n] = ds[naming.canon_index(n)]
        src[n] = dict(ds.attrs)
    ref = next(iter(arrays.values()))
    for n, da in arrays.items():
        for c in ('lat', 'lon', 'time'):
            if not np.array_equal(np.asarray(da[c].values),
                                  np.asarray(ref[c].values)):
                raise RuntimeError(f'{tag}: {n} has a different {c} axis from '
                                   f'the other indices of this member')
    glob_, var = partition_attrs(src)
    out = xr.Dataset({n: arrays[n] for n in arrays})
    for n in arrays:
        out[n].attrs = {**arrays[n].attrs, **var[n],
                        'source_file': found[n].rsplit('/', 1)[-1]}
        out[n].encoding = {}
    out.attrs = glob_
    out.attrs.update({
        'index_group': group,
        'index_group_note': group_note(group, args.baseline),
        'indices': ','.join(sorted(arrays)),
        'n_indices': len(arrays),
        'merged_on': time.strftime('%Y-%m-%d'),
        'merged_from': f'per-index annual files under {per_root}',
        'merge_attribute_rule': (
            'An attribute identical across every index in this file is '
            'global; one that differs, or that only some indices carry, is on '
            'the variable.')})
    if absent:
        out.attrs['indices_not_computed'] = ','.join(sorted(absent))
        out.attrs['indices_not_computed_reason'] = \
            NO_TEMPERATURE_EXTREMES[model]
    yy = ref['time'].dt.year.values
    fname = (f'ETCCDI_{group}_{naming.out_model(model)}_'
             f'{naming.out_scenario(scen)}_{label}_{int(yy.min())}-'
             f'{int(yy.max())}.nc')
    path = f'{member_dir(merged_root, model, scen, mem)}/{fname}'
    if exists(path) and not args.overwrite:
        say(f'  exists, skipping: {path}')
        return path
    enc = {n: {'zlib': True, 'complevel': 4} for n in arrays}
    fd, tmp = tempfile.mkstemp(suffix='.nc')
    os.close(fd)
    try:
        out.to_netcdf(tmp, encoding=enc)
        with xr.open_dataset(tmp, decode_timedelta=False) as back:
            for n, da in arrays.items():
                if not (np.array_equal(back[n].values, da.values,
                                       equal_nan=True)
                        and np.array_equal(back['time'].values,
                                           da['time'].values)):
                    raise RuntimeError(f'{tag}: {n} reads back differently '
                                       f'from its per-index file; not written')
        if is_s3(path):
            m.s3.put(tmp, path[5:])
        else:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            import shutil
            shutil.move(tmp, path)
            tmp = None
    finally:
        if tmp and os.path.exists(tmp):
            os.remove(tmp)
    say(f'  merged {fname}   {len(arrays)} indices'
        + (f', {len(absent)} not computed by design' if absent else ''))
    return path


def merge_all(args, model, per_root, merged_root):
    say(f'\n{"=" * 72}\nmerging into {merged_root}\n{"=" * 72}')
    n = 0
    for scen in args.scenarios:
        for mem in scenario_members(model, scen, args.members):
            for group, names in merge_groups(args).items():
                if merge_member(model, scen, mem, group, names, per_root,
                                merged_root, args):
                    n += 1
    say(f'{n} merged file(s)')


# --------------------------------------------------------------------------
# --full-grid: a series of bands, then whole-grid files
# --------------------------------------------------------------------------

def grid_lat(model, var):
    """
    The grid's latitudes, from one year of the first threshold member read
    without a band.
    """
    keep = getattr(m, 'READ_LAT_RANGE', None)
    if hasattr(m, 'READ_LAT_RANGE'):
        m.READ_LAT_RANGE = None
    try:
        y = DEFAULT_BASELINE[0]
        da = SITE.read(model, 'SSP245', var, SITE.baseline_members(model)[0],
                       Period('probe', y, y), None)
    finally:
        if hasattr(m, 'READ_LAT_RANGE'):
            m.READ_LAT_RANGE = keep
    if da is None:
        raise RuntimeError(f'{model}: could not read {var} to find the grid')
    lat = np.asarray(da['lat'].values, dtype='float64')
    d = np.diff(lat)
    if not (np.all(d > 0) or np.all(d < 0)):
        raise RuntimeError(f'{model}: latitude axis is not monotonic')
    return lat


def set_band(args, band):
    args.lat_range = band
    if hasattr(m, 'READ_LAT_RANGE'):
        m.READ_LAT_RANGE = band


def run_full_grid(args, model, base_members, per_root):
    band_root = (args.band_root or f'{per_root}_bands').rstrip('/')
    variables = [v for v in VARIABLES if v in args.variables]
    probe = next(v for v in variables
                 if not (model in NO_TEMPERATURE_EXTREMES
                         and v in ('tasmax', 'tasmin')))
    lat = grid_lat(model, probe)
    bands = lat_bands(lat, args.band_rows)
    say(f'full grid: {lat.size} latitudes in {len(bands)} band(s) of up to '
        f'{args.band_rows} rows; bands under {band_root}')
    tag = f'{args.convention}_{"boot" if args.bootstrap else "noboot"}'
    for var in variables:
        per_band = []
        for k, (r0, r1, lo, hi) in enumerate(bands):
            broot = f'{band_root}/rows{r0:04d}-{r1 - 1:04d}'
            marker = f'{broot}/_DONE_{var}_{tag}.json'
            say(f'\n##### {var} band {k + 1}/{len(bands)}: rows {r0}-{r1 - 1}, '
                f'latitudes {lat[r0]:.2f} to {lat[r1 - 1]:.2f}')
            if exists(marker) and not args.overwrite:
                bw = [tuple(x) for x in json.loads(read_text(marker))]
                say(f'  done earlier ({len(bw)} files); skipping')
            else:
                set_band(args, (lo, hi))
                bw = []
                run_variable(args, model, var, base_members, broot, bw)
                write_text(marker, json.dumps(bw))
            per_band.append((broot, bw))
        set_band(args, None)
        if not per_band[0][1]:
            continue
        say(f'\n##### {var}: joining {len(bands)} band(s)')
        keys = [tuple(x[1:]) for x in per_band[0][1]]
        for broot, bw in per_band[1:]:
            if [tuple(x[1:]) for x in bw] != keys:
                raise RuntimeError(f'{broot} wrote a different set of files '
                                   f'from the first band; not joining')
        for i, (scen, mem, name) in enumerate(keys):
            ds = _join([bw[i][0] for _, bw in per_band], lat,
                       naming.canon_index(name))
            yy = ds['time'].dt.year.values
            p = out_path(per_root, model, scen, mem, name,
                         (int(yy.min()), int(yy.max())))
            if exists(p) and not args.overwrite:
                say(f'  exists, skipping: {p}')
                continue
            write_any(ds, p)
            say(f'  wrote {p.rsplit("/", 1)[-1]}')
        if any(n in args.indices for n in PERCENTILE_BY_VAR[var]):
            first = threshold_path(per_band[0][0], model, var,
                                   args.convention)
            if exists(first):
                ds = _join([threshold_path(br, model, var, args.convention)
                            for br, _ in per_band], lat, 'threshold')
                tp = threshold_path(per_root, model, var, args.convention)
                write_any(ds, tp)
                say(f'  threshold -> {tp}')

def write_text(path, text):
    if is_s3(path):
        with m.s3.open(path[5:], 'w') as f:
            f.write(text)
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(text)


def read_text(path):
    if is_s3(path):
        with m.s3.open(path[5:], 'r') as f:
            return f.read()
    with open(path) as f:
        return f.read()


def lat_bands(lat, rows):
    """
    [(r0, r1, lo, hi)]: row ranges of the grid and the value range that
    selects exactly those rows (the neighbouring rows are a grid step away).
    """
    out = []
    for r0 in range(0, lat.size, rows):
        r1 = min(lat.size, r0 + rows)
        lo, hi = sorted((float(lat[r0]), float(lat[r1 - 1])))
        out.append((r0, r1, lo - 1e-6, hi + 1e-6))
    return out


def _join(paths, lat, name=None):
    """
    Band files joined along latitude, checked against the grid.
    """
    parts = [open_any(p) for p in paths]
    ds = xr.concat(parts, dim='lat', data_vars='minimal', coords='minimal',
                   compat='override', join='exact', combine_attrs='override')
    got = np.asarray(ds['lat'].values, dtype='float64')
    if got.shape != lat.shape or not np.array_equal(got, lat):
        raise RuntimeError(f'joined latitudes do not match the grid '
                           f'({got.size} against {lat.size}): {paths[0]}')
    ds.attrs.pop('partial_grid', None)
    for v in ds.data_vars:
        ds[v].attrs.pop('partial_grid', None)
        ds[v].encoding = {}
    return ds



# --------------------------------------------------------------------------
# E3SM checks for --preflight
# --------------------------------------------------------------------------

def e3sm_raw_stamps(case, native):
    """
    'end', 'middle' or 'start': where one raw history file stamps its
    first day within time_bnds, or None if it cannot be read. Raw E3SM output
    may be HDF5 or CDF-5; the second is read with the module's netCDF3
    byte-range reader.
    """
    gn = f'{m.REFL_BUCKET.rstrip("/")}/E3SMv3/G6-1.5K-HiLLA/{case}/day/{native}/gn'
    try:
        dirs = sorted(m.s3.ls(gn[5:]), reverse=True)
        raw = sorted(f for f in m.s3.ls(dirs[0]) if f.endswith('.nc'))[0]
    except (FileNotFoundError, IndexError):
        say(f'  no raw file under {gn}')
        return None
    t = tb = units = None
    try:
        with m.s3.open(raw, 'rb') as f:
            with xr.open_dataset(f, engine='h5netcdf', decode_times=False,
                                 decode_timedelta=False) as ds:
                bn = ds['time'].attrs.get('bounds', 'time_bnds')
                t = ds['time'].values[:2]
                tb = ds[bn].values[:2]
                units = ds['time'].attrs.get('units')
    except Exception:                                   # noqa: BLE001
        if not hasattr(m, '_nc3_header'):
            say(f'  {raw.rsplit("/", 1)[-1]}: not HDF5, and the module has '
                f'no netCDF3 reader')
            return None
        with m.s3.open(raw, 'rb', block_size=2**20) as f:
            dims, _, V, _, recsize = m._nc3_header(f)
        bn = V['time']['attrs'].get('bounds', 'time_bnds')
        fs, rpath = m.s3, raw
        t = m._nc3_read(fs, rpath, V['time'], recsize, 0, 2)
        tb = m._nc3_read(fs, rpath, V[bn], recsize, 0, 2)
        units = V['time']['attrs'].get('units')
    t = [float(x) for x in t]
    tb = [[float(y) for y in x] for x in tb]
    say(f'  raw {raw.rsplit("/", 1)[-1]}: time {t}, {bn} {tb}, units {units}')
    i = 1 if tb[0][0] == tb[0][1] else 0     # skip a zero-width first bin
    if t[i] == tb[i][1]:
        return 'end'
    if t[i] == tb[i][0]:
        return 'start'
    if t[i] == (tb[i][0] + tb[i][1]) / 2:
        return 'middle'
    return None


def e3sm_preflight(args, bad):
    """
    Coverage of every case; where a raw file stamps its days; the chunk
    layout of the merged windows; and, for the first threshold member, the
    days two windows both hold, read from each: the windows were regridded
    separately and should agree exactly.
    """
    native = m.get_variable_info('E3SM', 'pr')['native_short']
    say('\nE3SM coverage, from the merged windows on the bucket, stamps '
        f'moved back {E3SM_STAMP_SHIFT}')
    for scen in args.scenarios:
        for mem in scenario_members('E3SM', scen, args.members):
            plan, days = e3sm_plan(scen, mem, native, ALL_YEARS)
            if not plan:
                bad(f'E3SM {scen} {mem}: no merged_{native} window found')
                continue
            full = [y for y, n in days.items() if n >= 360]
            say(f'  {scen_label("E3SM", scen):14s} '
                f'{member_label("E3SM", scen, mem):3s} {e3sm_case(scen, mem)}: '
                f'{len(full)} complete years, {full[0]}-{full[-1]}'
                + ('' if full == list(range(full[0], full[-1] + 1))
                   else '  GAPS'))
            e3sm_coverage(scen, mem, plan, days)

    bm = SITE.baseline_members('E3SM')
    say('\nwhere E3SM stamps a day, from one raw history file')
    where = e3sm_raw_stamps(bm[0], native)
    if where == 'end':
        say('  end of the day: the 12-hour shift puts each stamp mid-day')
    elif where is None:
        bad('could not confirm the E3SM stamp convention from a raw file')
    else:
        bad(f'raw stamps sit at the {where} of the day; the 12-hour shift '
            f'is wrong for these files')

    files = e3sm_windows('SSP245', bm[0], native)
    if files:
        ds, f = _e3sm_open(files[0])
        try:
            enc = ds[native].encoding
            say(f'\n  {files[0].rsplit("/", 1)[-1]}: {dict(ds[native].sizes)}, '
                f'chunks {enc.get("chunksizes")}, '
                f'compression {enc.get("compression")}. A band reads every '
                f'chunk it crosses, so full-latitude chunks mean each band '
                f'streams the whole file.')
        finally:
            ds.close()
            f.close()
    checked = 0
    for a, b_ in zip(files, files[1:]):
        ta, tb = _e3sm_times(a), _e3sm_times(b_)
        ka = {_day(x): i for i, x in enumerate(ta)}
        common = [(ka[_day(x)], j) for j, x in enumerate(tb) if _day(x) in ka]
        if not common:
            continue
        # The first and last 30 shared days: the windows' edges, where a
        # stamp or merge error would show, without streaming years of chunks.
        if len(common) > 60:
            common = common[:30] + common[-30:]
        ia = np.array([c[0] for c in common])
        ib = np.array([c[1] for c in common])
        vals = []
        for p, idx in ((a, ia), (b_, ib)):
            ds, f = _e3sm_open(p)
            try:
                sub = ds[native]
                if args.lat_range:
                    r0, r1 = lat_rows(sub['lat'], *args.lat_range)
                    sub = sub.isel(lat=slice(r0, r1))
                vals.append(np.concatenate(
                    [sub.isel(time=sl).values for sl in _runs(idx)]
                ).astype('float64'))
            finally:
                ds.close()
                f.close()
        d = float(np.nanmax(np.abs(vals[0] - vals[1])))
        t0, t1 = ta[ia[0]], ta[ia[-1]]
        say(f'  overlap {t0.strftime("%Y-%m-%d")} to '
            f'{t1.strftime("%Y-%m-%d")} ({len(common)} days): '
            f'{a.rsplit("/", 1)[-1]} against {b_.rsplit("/", 1)[-1]}, '
            f'max|d| {d:.2e} {native}')
        if d != 0:
            bad(f'E3SM windows disagree over shared days; the day source '
                f'matters')
        checked += 1
    if not checked:
        say('  no day held by two windows; overlap check skipped')



# --------------------------------------------------------------------------
# --preflight
# --------------------------------------------------------------------------

def preflight(args, model):
    import xclim
    ok = True

    def bad(msg):
        nonlocal ok
        ok = False
        say(f'  FAIL  {msg}')

    say(f'python {sys.version.split()[0]}, xclim {xclim.__version__}, '
        f'xarray {xr.__version__}, numpy {np.__version__}, shim: {SHIM}')
    say(f'site {SITE.name}; modules {m.__file__}, {eb.__file__}')
    b = args.baseline
    bm = SITE.baseline_members(model)
    say(f'\n{model}: baseline {b.start_year}-{b.end_year}, threshold members '
        f'{bm}, {"bootstrap on" if args.bootstrap else "no bootstrap"}, '
        f'convention {args.convention}')
    say(f'indices: {len(args.indices)} of {len(ALL_INDICES)}; variables '
        f'{args.variables}')
    if model in NO_TEMPERATURE_EXTREMES:
        say(f'  {NO_TEMPERATURE_EXTREMES[model]}')

    say('\nmembers this run would write')
    for scen in args.scenarios:
        mems = scenario_members(model, scen, args.members)
        say(f'  {scen:6s} {len(mems):2d} '
            f'{[member_label(model, scen, x) for x in mems]}')

    say('\none real read: one year of the first threshold member')
    var = next((v for v in args.variables
                if not (model in NO_TEMPERATURE_EXTREMES
                        and v in ('tasmax', 'tasmin'))), None)
    if var is None:
        bad('no variable this model can carry')
    else:
        da = SITE.read(model, 'SSP245', var, bm[0],
                       Period('probe', b.start_year, b.start_year),
                       args.lat_range)
        if da is None:
            bad('probe read returned nothing')
        else:
            t0 = da['time'].values[0]
            say(f'  {var}: dims {da.dims} shape {da.shape} units '
                f'{da.attrs.get("units")} calendar '
                f'{getattr(t0, "calendar", "?")} dtype {da.dtype}')
            say(f'  latitudes {float(da.lat.min()):.2f} to '
                f'{float(da.lat.max()):.2f}')
            if tuple(da.dims) != ('time', 'lat', 'lon'):
                bad(f'dims are {da.dims}, expected (time, lat, lon)')

    if model == 'E3SM' and SITE.name == 'hub':
        e3sm_preflight(args, bad)

    per_root, merged_root = output_roots(args)
    for what, root in (('per-index', per_root), ('merged', merged_root)):
        if root is None:
            continue
        say(f'\n{what} output {root}')
        if is_s3(root):
            say('  on the bucket')
            continue
        try:
            os.makedirs(root, exist_ok=True)
            probe = os.path.join(root, '.write_probe')
            open(probe, 'w').close()
            os.remove(probe)
            say('  writable')
        except OSError as e:
            bad(f'not writable: {e}')

    say(f'\nPREFLIGHT {"PASSED" if ok else "FAILED"}')
    return ok


# --------------------------------------------------------------------------

def parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('--site', choices=['hub', 'derecho'],
                    default=os.environ.get('ETCCDI_SITE', 'hub'))
    ap.add_argument('--model', required=True,
                    choices=['CESM', 'UKESM', 'MIROC', 'E3SM'])
    ap.add_argument('--indices', nargs='+', default=None, metavar='INDEX',
                    help=f'default: all {len(ALL_INDICES)}. Known: '
                         f'{" ".join(ALL_INDICES)}')
    ap.add_argument('--variables', nargs='+', default=None, choices=VARIABLES,
                    help='restrict to the indices of these daily variables')
    ap.add_argument('--scenarios', nargs='+', default=list(SCENARIOS),
                    choices=SCENARIOS)
    ap.add_argument('--members', nargs='+', default=None,
                    help='restrict the members written (testing); the '
                         'threshold always uses every threshold member')
    ap.add_argument('--baseline', nargs=2, type=int,
                    default=list(DEFAULT_BASELINE), metavar=('Y0', 'Y1'),
                    help='SSP2-4.5 years of the percentile threshold '
                         '(default 2020 2039)')
    ap.add_argument('--bootstrap', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='Zhang et al. (2005) bootstrap for the in-base years '
                         '(default on)')
    ap.add_argument('--convention', default=eb.CONVENTIONS[0],
                    choices=eb.CONVENTIONS,
                    help='pooled (default): one threshold from all threshold '
                         'members pooled; member_mean: per-member thresholds '
                         'averaged')
    ap.add_argument('--output', choices=['per-index', 'merged', 'both'],
                    default='per-index')
    ap.add_argument('--out-root', required=True,
                    help='local directory, or s3:// on the Hub')
    ap.add_argument('--lat-range', nargs=2, type=float, default=None,
                    metavar=('LO', 'HI'),
                    help='compute only latitudes LO..HI, in degrees (writes a '
                         'partial grid; local --out-root only)')
    ap.add_argument('--full-grid', action='store_true',
                    help='compute the globe as latitude bands and join them')
    ap.add_argument('--band-rows', type=int, default=20)
    ap.add_argument('--band-root', default=None,
                    help='where --full-grid keeps its bands (default '
                         '<per-index root>_bands)')
    ap.add_argument('--workers', type=int,
                    default=int(os.environ.get('NCPUS', 3)))
    ap.add_argument('--cells-per-task', type=int, default=96)
    ap.add_argument('--overwrite', action='store_true')
    ap.add_argument('--preflight', action='store_true')
    return ap


def build(model, out_root, site=None, indices=None, variables=None,
          scenarios=SCENARIOS, members=None, baseline=DEFAULT_BASELINE,
          bootstrap=True, convention=None, output='per-index',
          lat_range=None, full_grid=False, band_rows=20, band_root=None,
          workers=None, cells_per_task=96, overwrite=False, preflight=False):
    """
    The command line as a function, for notebooks. Arguments as the
    options of the same name; returns nothing, writes the files.
    """
    argv = ['--model', model, '--out-root', out_root,
            '--output', output, '--band-rows', str(band_rows),
            '--cells-per-task', str(cells_per_task),
            '--baseline', str(baseline[0]), str(baseline[1]),
            '--scenarios', *scenarios]
    if site:
        argv += ['--site', site]
    if indices:
        argv += ['--indices', *indices]
    if variables:
        argv += ['--variables', *variables]
    if members:
        argv += ['--members', *members]
    if not bootstrap:
        argv += ['--no-bootstrap']
    if convention:
        argv += ['--convention', convention]
    if lat_range:
        argv += ['--lat-range', str(lat_range[0]), str(lat_range[1])]
    if full_grid:
        argv += ['--full-grid']
    if band_root:
        argv += ['--band-root', band_root]
    if workers:
        argv += ['--workers', str(workers)]
    if overwrite:
        argv += ['--overwrite']
    if preflight:
        argv += ['--preflight']
    try:
        main(argv)
    except SystemExit as e:
        if e.code not in (0, None):
            raise


def main(argv=None):
    global SITE
    args = parser().parse_args(argv)

    SITE = HubSite() if args.site == 'hub' else DerechoSite()
    if args.site == 'derecho' and args.model == 'E3SM':
        sys.exit('E3SMv3 is not on GLADE; use --site hub')

    idx = list(args.indices or ALL_INDICES)
    canon = {n.lower(): n for n in ALL_INDICES}
    unknown = [n for n in idx if n.lower() not in canon]
    if unknown:
        sys.exit(f'unknown index {unknown}; known: {" ".join(ALL_INDICES)}')
    idx = [canon[n.lower()] for n in idx]
    if args.variables:
        idx = [n for n in idx if VAR_OF[n] in args.variables]
    args.indices = idx
    args.variables = [v for v in VARIABLES
                      if any(VAR_OF[n] == v for n in idx)]
    if not args.indices:
        sys.exit('no index selected')
    y0, y1 = args.baseline
    if y0 > y1:
        sys.exit('--baseline Y0 must not be after Y1')
    args.baseline = Period('baseline', y0, y1)

    if args.full_grid and args.lat_range:
        sys.exit('--full-grid and --lat-range exclude each other')
    if args.lat_range and is_s3(args.out_root):
        sys.exit('--lat-range writes a partial grid; use a local --out-root')
    if args.site == 'derecho' and is_s3(args.out_root):
        sys.exit('--site derecho writes to a local --out-root')
    if args.lat_range:
        if hasattr(m, 'READ_LAT_RANGE'):
            m.READ_LAT_RANGE = tuple(args.lat_range)
        args.lat_range = tuple(args.lat_range)

    if args.preflight:
        ok = preflight(args, args.model)
        if not ok:
            sys.exit(1)
        return

    t0 = time.time()
    per_root, merged_root = output_roots(args)
    bm = SITE.baseline_members(args.model)
    say(f'{args.model} on {SITE.name}: {len(args.indices)} indices, baseline '
        f'{y0}-{y1}, {"bootstrap" if args.bootstrap else "no bootstrap"}, '
        f'{args.convention}, workers {args.workers}, band '
        f'{"bands of %d rows" % args.band_rows if args.full_grid else args.lat_range or "globe"}'
        f'\nper-index output {per_root}'
        + (f'\nmerged output    {merged_root}' if merged_root else ''))
    if args.workers == 1 and args.bootstrap and \
            any(n in PERCENTILE_INDICES for n in args.indices):
        say('WARNING: one worker. The in-base bootstrap runs one block at a '
            'time and can take many hours; give --workers the number of CPUs '
            'in the job (and set OMP_NUM_THREADS=1).')
    if args.full_grid:
        run_full_grid(args, args.model, bm, per_root)
    else:
        written = []
        for var in args.variables:
            run_variable(args, args.model, var, bm, per_root, written)
    if merged_root:
        merge_all(args, args.model, per_root, merged_root)
    say(f'\nall done in {(time.time() - t0) / 60:.1f} min')


if __name__ == '__main__':
    main()
