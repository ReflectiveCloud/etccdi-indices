"""
derecho_io.py  -  read daily fields from GLADE, with the corrections the
pre-flight diagnostics turned up.

The hub reader (multimodel_etccdi.read_model_var) does the same job against S3.
This is its GLADE counterpart and carries three fixes that the hub reader does
not have, each traceable to a specific diagnostic output.


1. CESM timestamps are end-of-interval
   01_followup_diagnostics section B, on
   b.e21.BW.f09_g17.SSP245-G6-1p5K-SAI.001.cam.h1.TREFHTMX.20750101-20841231.nc:
   time[0] is 2075-01-01 while time_bnds[0] runs 2074-12-31 to 2075-01-01, and
   time equals the upper bound exactly. The day of 31 December 2074 is therefore
   labelled 1 January 2075, and every CESM record is shifted one day late. The
   visible symptom is the spurious one-day 2085 bin every member reported.

   Fix: replace time with the bounds midpoint before anything slices or
   resamples. Applied to CESM only; MIROC and UKESM were not end-stamped.

2. CESM HiLLA member 001 uses overlapping segments
   Section C: 001 is written as 2035-2039, then five-day-overlapping decades
   (2039-12-31 to 2050-01-05, 2050-01-01 to 2060-01-05, and so on), where 002
   and 003 use clean non-overlapping decades. open_mfdataset refuses it with a
   non-monotonic index. The overlap is real duplication, not corruption, and
   members 002 and 003 are unaffected.

   Fix: concatenate, then keep the first occurrence of each timestamp. The hub
   reader already does this, so the behaviour matches.

3. Two UKESM scenarios keep three temperature fields in one file
   Section A: G6-1.5K-SAI has pr/ and tas/ directories and no tasmax/ or
   tasmin/. The three fields are internal variables of the tas file. The
   UKESM1-1 SSP2-4.5 reference fetched from JASMIN on 19 September 2026 is
   written the same way, in tas_mean_min_max_day_*.nc.

   The mapping lives in derecho_layout.VAR_IN_FILE, keyed by (model, scenario).
   verify_multivar_temps() below checks it against both cell_methods and the
   ordering of the field means, and this module refuses to read tasmax or
   tasmin for such a scenario until it passes. The check is cached per
   (model, scenario): one scenario's pass is not evidence about another's file.

Not corrected here, because it does not need to be: MIROC latitude runs
descending where its land-fraction file runs ascending (section E; the
difference is pure ordering, 177.9 degrees at index 0, both float64, longitude
identical). Fields are left on their native grid and the ordering is handled at
masking time by mask_fixes.mask_after_reproject.
"""

import glob
import os

import numpy as np
import xarray as xr

import derecho_layout as dl

_CODER = None


def _coder():
    """xarray 2026.7 deprecates use_cftime as a keyword in favour of a coder."""
    global _CODER
    if _CODER is None:
        _CODER = xr.coders.CFDatetimeCoder(use_cftime=True)
    return _CODER


def _files(model, scenario, var, member):
    tmpl = dl.get(model, scenario, var)
    p = tmpl.format(ensemble_id=member)
    if p.endswith('/'):
        p = p + '*.nc'
    return sorted(glob.glob(p))


def _fix_cesm_time(ds):
    """Replace time with the bounds midpoint. No-op if bounds are absent."""
    bname = next((b for b in ['time_bnds', 'time_bounds'] if b in ds), None)
    if bname is None:
        print('    WARNING: CESM file has no time bounds; leaving time as found. '
              'Check whether it is end-stamped before trusting the window.')
        return ds
    b = ds[bname].values
    width = b[:, 1] - b[:, 0]
    mid = b[:, 0] + width / 2

    # Drop zero-width bins. The first record of a CAM run carries
    # time_bnds = [start, start], a bin covering no time at all: it is the model
    # state written at initialisation, not a day of output. Its midpoint lands
    # at the run start, 2015-01-01 00:00 for CESM SSP2-4.5, twelve hours before
    # the genuine 1 January record at 12:00. The two are distinct timestamps, so
    # de-duplication does not remove it, and the year gains a 366th day on a
    # noleap calendar. It shows up only in a window that includes the first year
    # of the run, which is why 2020-2039 is clean and 2015-2034, the new
    # baseline window, is not.
    keep = width != (width * 0)
    n_drop = int((~keep).sum())
    if n_drop:
        print(f'    dropping {n_drop} zero-width bin(s) at '
              f'{[str(x) for x in mid[~keep]]}')
        ds = ds.isel(time=np.flatnonzero(keep))
        mid = mid[keep]
    return ds.assign_coords(time=('time', mid))


def _dedupe_time(obj):
    """Keep the first occurrence of each timestamp. Works on a Dataset or a
    DataArray; CESM HiLLA member 001 is written with five-day overlaps between
    its decade segments and needs this."""
    _, uniq = np.unique(obj['time'].values, return_index=True)
    return obj.isel(time=sorted(uniq))


# cell_methods reduction -> the CMIP variable it corresponds to.
_METHOD_TO_VAR = {'time: maximum': 'tasmax', 'time: minimum': 'tasmin',
                  'time: mean': 'tas'}


def verify_multivar_temps(model, scenario, member):
    """
    Confirm which internal variable of a multi-variable temperature file is the
    maximum and which the minimum.

    Was verify_ukesm_vars(member), hardcoded to ('UKESM', 'SAI') in both the
    path it read and the mapping it checked. That was correct while UKESM SAI
    was the only scenario keeping three temperatures in one file. The UKESM1-1
    SSP2-4.5 reference from JASMIN does the same, so reading its tasmin sent a
    SSP245 request down this function and it went looking for a SAI file for
    member u-dp902, failing with 'no UKESM SAI tas file'. The scenario is now a
    parameter and the mapping is looked up under the pair actually being read.

    Two independent lines of evidence, because either alone can mislead.
    cell_methods states the reduction and is authoritative when present, but is
    absent from some non-CMORized output. The ordering of the means is always
    available but only says which field is warmest on average, which is weaker
    than it looks: two members of a triple could be close. Both are reported and
    a disagreement between them fails the check rather than picking a winner.

    Returns True if the derecho_layout.VAR_IN_FILE mapping holds.
    """
    paths = _files(model, scenario, 'tas', member)
    if not paths:
        raise FileNotFoundError(
            f'no {model} {scenario} tas file for member {member}; '
            f'template {dl.get(model, scenario, "tas")}')
    ds = xr.open_dataset(paths[0], decode_times=_coder())
    means, methods = {}, {}
    for name, v in ds.data_vars.items():
        if v.ndim >= 3 and 'time' in v.dims:
            # float64 before the reduction. A float32 accumulator over an array
            # this size stops advancing once the running sum outgrows the
            # spacing between representable values, which is the same failure
            # documented in read_member below and which biased these means low
            # by about 17 K. It biased all three the same way, so the ordering
            # survived, but a comparison between numbers that have all stalled
            # is not a comparison.
            means[name] = float(v.isel(time=slice(0, 365))
                                .astype('float64').mean().compute())
            cm = str(v.attrs.get('cell_methods', ''))
            methods[name] = next((var for key, var in _METHOD_TO_VAR.items()
                                  if key in cm), None)
            print(f'  {name:24s} mean {means[name]:8.3f} K  '
                  f'cell_methods={v.attrs.get("cell_methods", "-")!r}')
    ds.close()

    claimed = dl.VAR_IN_FILE[(model, scenario)]
    order = sorted(means, key=means.get, reverse=True)
    by_values = order[0] == claimed['tasmax'] and order[-1] == claimed['tasmin']
    print(f'  ordering high to low: {order}')
    print(f'  layout claims tasmax={claimed["tasmax"]}, '
          f'tasmin={claimed["tasmin"]}: '
          f'{"CONFIRMED" if by_values else "CONTRADICTED"} by the means')

    stated = {v: k for k, v in methods.items() if v is not None}
    if stated:
        by_meta = all(stated.get(k) == v for k, v in claimed.items()
                      if k in stated)
        print(f'  cell_methods says {stated}: '
              f'{"CONFIRMED" if by_meta else "CONTRADICTED"}')
        if by_meta != by_values:
            print('  cell_methods and the means DISAGREE. Neither is trusted '
                  'on its own here;\n  open the file and settle it before any '
                  'temperature index is computed.')
            return False
    else:
        print('  no cell_methods on any field; the means are the only evidence')

    return by_values


# Verification is per (model, scenario), not global. One flag for the whole
# module meant the first scenario checked stood in for every later one, so a
# second multi-variable scenario with a different internal layout would have
# inherited a pass it never earned.
_VARS_VERIFIED = {}

_BOUNDS_HINTS = ('bnds', 'bounds', '_bnds', 'time_bnds', 'lat_bnds', 'lon_bnds')


def _pick_var(ds, path):
    """
    The one data variable that is an actual field.

    The test is a time dimension, not the number of dimensions. Coordinate
    bounds (lat_bnds, lon_bnds) are two-dimensional in the file and carry no
    time axis, so requiring 'time' in dims separates them cleanly, and it stays
    correct whatever else a modelling centre writes alongside the field.
    """
    cands = [v for v in ds.data_vars
             if 'time' in ds[v].dims and ds[v].ndim >= 3
             and not any(h in v for h in _BOUNDS_HINTS)]
    if len(cands) != 1:
        raise ValueError(f'{os.path.basename(path)}: expected one field with a '
                         f'time dimension, found {cands} among '
                         f'{list(ds.data_vars)}. Name it in '
                         f'derecho_layout.VAR_IN_FILE.')
    return cands[0]



def read_member(model, scenario, var, member, period):
    """
    One variable, one member, sliced to a Period, in canonical units (mm/day for
    pr, K for temperature). Returns a DataArray, or None if no files match.
    """
    paths = _files(model, scenario, var, member)
    if not paths:
        print(f'  no files: {model}/{scenario}/{var}/{member}')
        return None

    pick = None
    if (model, scenario) in dl.VAR_IN_FILE:
        pick = dl.VAR_IN_FILE[(model, scenario)].get(var)
        if pick is not None and var in ('tasmax', 'tasmin'):
            key = (model, scenario)
            if key not in _VARS_VERIFIED:
                _VARS_VERIFIED[key] = verify_multivar_temps(
                    model, scenario, member)
            if not _VARS_VERIFIED[key]:
                raise RuntimeError(
                    f'{model}/{scenario}: the internal variable mapping in '
                    f'derecho_layout.VAR_IN_FILE is contradicted by the file. '
                    f'Fix the mapping before computing any temperature index; a '
                    f'swapped tasmax and tasmin produces plausible, wrong fields.')

    parts, src_attrs = [], {}
    for p in paths:
        ds = xr.open_dataset(p, decode_times=_coder())
        if not src_attrs:
            src_attrs = dict(ds.attrs)
        if model == 'CESM':
            ds = _fix_cesm_time(ds)
        name = pick if pick is not None else _pick_var(ds, p)
        if name not in ds:
            raise KeyError(f'{name!r} not in {os.path.basename(p)}; '
                           f'present: {list(ds.data_vars)}')
        # Extract the one variable before slicing and before concatenating.
        # Concatenating whole Datasets along time broadcasts every variable that
        # lacks a time dimension, so lat_bnds and lon_bnds acquire a length-7300
        # time axis: they then look like 3-D fields to any later selection, and
        # they cost real memory. Selecting first avoids both.
        sub = ds[name].sel(time=period.slice())
        if sub.sizes.get('time', 0):
            parts.append(sub.load())
        ds.close()

    if not parts:
        print(f'  nothing in {period.start_year}-{period.end_year}: '
              f'{model}/{scenario}/{var}/{member}')
        return None

    da = _dedupe_time(xr.concat(parts, dim='time'))

    da = da.rename({k: v for k, v in {'latitude': 'lat', 'longitude': 'lon'}.items()
                    if k in da.coords})

    # Carry the CMOR identity of the run onto the field. Without this the only
    # record of which model version and which realisation produced a value is
    # the directory it was read from, and directory names do not survive a copy.
    # It matters most for UKESM, where the GLADE SSP2-4.5 is UKESM1-0-LL while
    # the scenarios beside it are UKESM1-1-LL.
    for key in ('source_id', 'variant_label', 'realization_index',
                'physics_index', 'forcing_index'):
        if key in src_attrs:
            da.attrs[f'src_{key}'] = src_attrs[key]

    # When a source carries none of those, say so on the field rather than
    # leaving the absence to be inferred from missing attributes. The JASMIN
    # UKESM1-1 files have no CMOR identity at all: their only useful global
    # attribute is 'Data from Met Office Unified Model', so the member labels
    # r1, r2 and r3 come from the suite ids and nothing else. A reader three
    # months from now sees a directory called r1 and has no way to know that
    # unless the file says it.
    if not any(k.startswith('src_') for k in da.attrs):
        da.attrs['src_identity'] = (
            f'none in source; member label taken from the directory id '
            f'{member!r} in derecho_layout.OUT_LABEL')
        if 'source' in src_attrs:
            da.attrs['src_source'] = src_attrs['source']
    note = dl.VERSION_NOTE.get((model, scenario))
    if note:
        da.attrs['version_note'] = note

    file_units = da.attrs.get('units')
    if var == 'pr':
        da = da * dl.PR_SCALE[model]
        da.attrs['units'] = 'mm/day'
    else:
        # Read the units rather than asserting them. An earlier version set 'K'
        # unconditionally, which would have silently mislabelled a Celsius field
        # and handed xclim kelvin thresholds for Celsius data.
        if file_units in ('K', 'kelvin', 'Kelvin'):
            da.attrs['units'] = 'K'
        elif file_units in ('degC', 'C', 'celsius', 'Celsius'):
            da.attrs['units'] = 'degC'
        else:
            raise ValueError(
                f'{model}/{scenario}/{var}/{member}: unrecognised temperature '
                f'units {file_units!r}. Name them explicitly rather than '
                f'assuming; xclim thresholds are unit-sensitive.')

    n_years = len(np.unique(da['time'].dt.year.values))
    # float64 accumulator. A float32 mean over an array this size stagnates:
    # once the running sum reaches 2**33 the gap between representable values is
    # 1024, so adding a ~288 K value rounds to no change and the sum freezes.
    # Every member then returns exactly float32(2**33 / N), which for
    # 7300 x 192 x 288 is 21.280080795288086 whatever the data. That is what
    # produced six identical means across two scenarios and three members in the
    # 15 September diagnostic; the fields themselves were fine.
    mean = float(da.astype('float64').mean())
    print(f'  {model}/{scenario}/{var}/{member}: {da.sizes["time"]} steps, '
          f'{n_years} years, mean {mean:.3f} {da.attrs["units"]}')
    return da


def check_not_duplicate(model, scenario, member, period):
    """
    Confirm tas, tasmax and tasmin are three different fields.

    CESM SSP2-4.5 members r1-r5 on pangeo had tasmax and tasmin written as
    copies of the daily mean, which is why that ensemble was abandoned in favour
    of r6-r10. The GLADE CESM SAI and HiLLA ensembles are numbered 001 to 003,
    the same numbering as the members affected elsewhere, so this is worth
    running once per model before a production build rather than assuming raw
    CAM output is exempt. E3SMv3 fails this test by construction.

    A pass means the fields differ. A fail means the member cannot carry
    temperature indices.
    """
    fields = {}
    for var in ('tas', 'tasmax', 'tasmin'):
        try:
            fields[var] = read_member(model, scenario, var, member, period)
        except Exception as e:
            print(f'  {var}: unreadable ({type(e).__name__}: {e})')
    if len(fields) < 2:
        print('  too few fields to compare')
        return False
    ok = True
    names = sorted(fields)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            d = (fields[a] - fields[b]).astype('float64')
            mx = float(np.abs(d).max())
            rms = float((d ** 2).mean() ** 0.5)
            identical = mx == 0.0
            print(f'  {a} vs {b}: rms {rms:.4f}, max |diff| {mx:.4f}'
                  f'{"   IDENTICAL" if identical else ""}')
            ok &= not identical
    print(f'  {model}/{scenario}/{member}: '
          f'{"distinct fields" if ok else "DUPLICATE FIELDS, do not use"}')
    return ok


def load_members(model, scenario, var, period, members=None):
    """{member: DataArray}. Members default to those whose record reaches the
    end of the requested period, which excludes the two short MIROC HiLLA runs
    from the 2065-2084 window and keeps them for the baseline."""
    if members is None:
        members = dl.members(model, scenario, through_year=period.end_year)
    out = {}
    for m in members:
        da = read_member(model, scenario, var, m, period)
        if da is not None:
            out[m] = da
    if not out:
        raise RuntimeError(f'no members loaded for {model}/{scenario}/{var}')
    return out