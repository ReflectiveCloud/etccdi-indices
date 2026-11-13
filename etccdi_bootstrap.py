"""
etccdi_bootstrap.py  -  multi-member percentile thresholds with the Zhang et
al. (2005) bootstrap and the ETCCDI precipitation percentile indices as
defined by Zhang et al. (2011) and climpact.

This module is computation only. It reads nothing from disk and knows nothing
about models, paths or output trees; build_indices_derecho_boot.py supplies the
data and writes the files. 


WHY THIS EXISTS as of October 2026

The percentile indices showed a step at 2035, the first year after every
model's baseline. Zhang et al. (2005, J. Climate 18, 1641) explain it: a year
inside the base period is compared against a threshold its own data helped
estimate, so its exceedance rate is held near nominal, while years outside
meet a threshold that carries sampling error and exceed it more often. Zhang
et al. (2011, WIREs Clim. Change 2, 851, p. 857) give the fix: "This is done by
not using any data from a particular year within the base period to calculate
the index threshold applied to data for that year."

xclim implements that procedure (xclim.core.bootstrapping.bootstrap_func), but
only for a threshold computed by percentile_doy from the same single series the
index is computed on. Our thresholds come from several ensemble members and are
passed in from outside, so switching xclim's flag on would recompute a
single-member percentile for the in-base years while the out-of-base years kept
the multi-member one. This module does the bootstrap on the multi-member
threshold instead.


THE THRESHOLD, convention='pooled' (the default)

This is the live pipeline's threshold since 30 September 2026,
multimodel_etccdi.pooled_percentile_doy, reproduced line for line: member j's
base series is re-stamped j * 400 years later (400 keeps the Gregorian leap
cycle), the members are concatenated along time in the order given, and xclim's
percentile_doy is taken once over the result, 5-day centred window, type 8. The
window therefore crosses each join between members once, exactly as it does in
the live code. Out-of-base years computed here equal the existing archive bit
for bit; test T5 checks that against a verbatim copy of the live functions.

Bootstrap for in-base year x of threshold member m: xclim's own
build_bootstrap_year_da overwrites year x in m's base series with each of m's
other base years in turn, including its calendar-length handling. For each
altered series the pooled threshold is recomputed with the other members
unchanged, the index is computed on m's original year x against every altered
threshold with xclim's own index function, and the results are averaged. With
one member this is exactly xclim's bootstrap_func (test T2).

Drawing the replacement from the same member is the one choice the
multi-member case adds. Drawing from any pooled member-year would also remove
year x's data, at M times the cost; Zhang's replacement exists to keep the
sample size fixed, and either choice does that.

convention='member_mean' is the convention before 30 September: percentile_doy
per member, then the mean over members, NaN skipped. Kept because it is
tested and cheaper, not because it should be used: the live module documents
that it biases R99p about 20 percent low on CESM2-WACCM. It does not reproduce
the current archive.


PRECIPITATION

R95p, R99p and the derived indices use "the sample of all wet days in the base
period without consideration of the annual cycle" (Zhang 2011, p. 855): one
threshold per grid cell from every base-period day with at least 1 mm, every
member's wet days in one sample ('pooled') or per member then averaged
('member_mean'), type 8. No day-of-year window, so no join question arises.
The bootstrap is the same leave-one-year-out-with-replacement as above. This is
a change from the existing tree, whose R95p and R99p thresholds were day-of-year
percentiles of wet days in a 5-day window.


INDEX NAMES, following Mari Tye's clarification of 9 October 2026, which matches
Zhang 2011 Table 1 and climpact

    R95p, R99p        annual total precipitation on wet days above the 95th /
                      99th percentile                                   mm
    R95pTOT, R99pTOT  100 * R95p / PRCPTOT, the contribution of those days to
                      the annual total (Zhang 2011, Fig. 7b caption)    %
    N95, N99          number of wet days above the 95th / 99th percentile,
                      which is what the earlier tree called R95p and R99p  days

The temperature percentile indices are returned as counts of days, as in the
existing tree and as xclim returns them. ETCCDI defines them as a percentage of
days (Zhang 2011 Table 1); multiply by 100 / days-in-year to convert.


WHAT IT DOES NOT DO

No detrending. xclim's bootstrap does not detrend, Zhang (2005) does not, and
Mari's ARISE practice did not, for reasons she gave on 9 October 2026. Lala's
pooled-and-detrended approach for heatwave thresholds is a different method and
is not reproduced here.
"""

from __future__ import annotations

import operator
import sys
import types
import warnings

import numpy as np
import pandas as pd
import xarray as xr

CONVENTIONS = ('pooled', 'member_mean')     # first is the default
POOL_OFFSET_YEARS = 400           # as multimodel_etccdi.POOL_OFFSET_YEARS
BOOTSTRAP_DIM = '_bootstrap'      # the name xclim uses, so its index functions
                                  # average over it the way bootstrap_func does
WINDOW = 5                        # ETCCDI: five-day window centred on each day
ALPHA = BETA = 1.0 / 3.0          # Hyndman and Fan type 8, as ETCCDI specifies
WET_DAY_MM = 1.0                  # ETCCDI wet day, RR >= 1 mm


# --------------------------------------------------------------------------
# xclim compatibility
# --------------------------------------------------------------------------

def ensure_xclim_compat(verbose=False):
    """
    xclim 0.54 calls xarray.core.ops.get_op, which xarray 2026.7 removed.

    On Derecho the production patch, src/xclim_compat.py, is used so that this
    build and the existing one run xclim under the same shim. Elsewhere, for
    the test suite, a minimal equivalent is installed: get_op('__gt__') and its
    siblings resolve to the functions in the operator module, which is what the
    removed xarray function did.
    """
    try:
        import xclim_compat                                   # noqa: F401
        xclim_compat.patch(verbose=verbose)
        return 'xclim_compat'
    except ImportError:
        pass
    import xarray.core as xcore
    if not hasattr(xcore, 'ops'):
        mod = types.ModuleType('xarray.core.ops')
        mod.get_op = lambda name: getattr(operator, name)
        xcore.ops = mod
        sys.modules['xarray.core.ops'] = mod
    return 'local minimal shim'


# --------------------------------------------------------------------------
# Day-of-year thresholds, pooled across members
# --------------------------------------------------------------------------

def _shift_years(da, offset):
    """
    The same record stamped `offset` years later, calendar preserved, as
    multimodel_etccdi._shift_years does it. A numpy datetime64 axis cannot hold
    years past 2262, so it is first converted to cftime on the same
    (proleptic Gregorian) calendar; on Derecho derecho_io decodes every model
    to cftime already, so this branch is for the test suite.
    """
    if not offset:
        return da
    if np.issubdtype(da['time'].dtype, np.datetime64):
        da = da.convert_calendar('proleptic_gregorian', use_cftime=True)
    t = da['time'].values
    out = np.empty(t.shape, dtype=object)
    for i, x in enumerate(t):
        out[i] = x.replace(year=x.year + offset)
    return da.assign_coords(time=('time', out))


def pooled_series(series_list):
    """
    Members laid end to end, member j re-stamped j * 400 years, as
    multimodel_etccdi.pooled_percentile_doy builds its sample.
    """
    if len(series_list) == 1:
        big = series_list[0]
        big.attrs.update(series_list[0].attrs)
        return big
    # every member on cftime, or a datetime64 member 0 would be concatenated
    # with cftime members 1..n into an axis of mixed types
    series_list = [da.convert_calendar('proleptic_gregorian', use_cftime=True)
                   if np.issubdtype(da['time'].dtype, np.datetime64) else da
                   for da in series_list]
    parts = [_shift_years(da, j * POOL_OFFSET_YEARS)
             for j, da in enumerate(series_list)]
    big = xr.concat(parts, dim='time')
    big.attrs.update(series_list[0].attrs)
    return big


def pooled_percentile_doy(series_list, per, window=WINDOW, big=None):
    """
    Day-of-year percentile over several members pooled into one sample, the
    live pipeline's way: member j shifted j * 400 years, concatenated along
    time, one percentile_doy. With one member this is percentile_doy.

    series_list : base-period series, in the pipeline's member order
    per         : percentile or list of percentiles in [0, 100]
    big         : pooled_series(series_list), if already built
    """
    from xclim.core.calendar import percentile_doy
    if big is None:
        big = pooled_series(series_list)
    return percentile_doy(big, per=per, window=window).load()


def year_groups(series):
    """
    Calendar-year groups of a base series, as bootstrap_func forms them.
    """
    return series.resample(time='YS').groups


def _year_of(key):
    from xclim.core.bootstrapping import _get_year_label
    return int(_get_year_label(key))


def replacement_series(series, year):
    """
    xclim's build_bootstrap_year_da for one in-base year: dims (_bootstrap,
    time), one entry per other base year.
    """
    from xclim.core.bootstrapping import build_bootstrap_year_da
    groups = year_groups(series)
    key = next(k for k in groups if _year_of(k) == year)
    return build_bootstrap_year_da(series, groups, key)


def bootstrapped_pooled_doy(base, member, year, per, window=WINDOW, big=None):
    """
    Pooled thresholds with `year` of `member` replaced by each of that member's
    other base years. Returns dims (_bootstrap, ..., dayofyear, percentiles).

    The pooled sample is built once. For each replacement only the values of
    `year` in `member`'s stretch of it are overwritten, with the values xclim's
    build_bootstrap_year_da put there, and percentile_doy is taken again. The
    time axis never changes, since build_bootstrap_year_da keeps the replaced
    year's own time stamps.

    One percentile_doy per replacement rather than one call over the whole
    _bootstrap dimension: the vectorised call would broadcast every other
    member along _bootstrap and hold that many copies of the pooled sample at
    once. The numbers are the same either way.
    """
    from xclim.core.calendar import percentile_doy
    keys = list(base)
    j = keys.index(member)
    if big is None:
        big = pooled_series([base[k] for k in keys])
    series = base[member]
    start = sum(base[k].sizes['time'] for k in keys[:j])
    local = np.nonzero(series['time'].dt.year.values == year)[0]
    if local.size == 0:
        raise ValueError(f'{year} not in the base series of {member}')
    pos = start + local
    tax = big.get_axis_num('time')

    bda = replacement_series(series, year)
    out = []
    for b in range(bda.sizes[BOOTSTRAP_DIM]):
        rep = (bda.isel({BOOTSTRAP_DIM: b}, drop=True).isel(time=local)
               .transpose(*big.dims).values)
        vals = big.values.copy()
        idx = [slice(None)] * vals.ndim
        idx[tax] = pos
        vals[tuple(idx)] = rep
        alt = big.copy(data=vals)
        out.append(percentile_doy(alt, per=per, window=window).load())
    return xr.concat(out, dim=BOOTSTRAP_DIM)


def _member_mean(per_member, like_attrs=None):
    """
    Mean over members, NaN skipped, the way compute_baseline_threshold does
    it: concat along a 'member' index, then .mean('member'). Attributes of the
    first member are carried, since xclim's index functions read units and
    window from the threshold.
    """
    keys = list(per_member)
    stacked = xr.concat([per_member[k] for k in keys],
                        dim=pd.Index(keys, name='member'))
    out = stacked.mean('member')
    out.attrs.update(like_attrs if like_attrs is not None
                     else per_member[keys[0]].attrs)
    return out


def member_percentiles_doy(base, per, window=WINDOW):
    """
    {member: xclim percentile_doy of that member's base series}.
    """
    from xclim.core.calendar import percentile_doy
    return {k: percentile_doy(v, window=window, per=per).load()
            for k, v in base.items()}


def member_mean_bootstrapped_doy(base, per_member, member, year, per,
                                 window=WINDOW):
    """
    Member-mean thresholds with `year` of `member` replaced by each of that
    member's other base years. Returns dims (_bootstrap, ..., dayofyear,
    percentiles).

    The altered percentile is computed the way bootstrap_func computes it: one
    percentile_doy call on the series returned by build_bootstrap_year_da, with
    the _bootstrap dimension carried through. The other members' percentiles do
    not contain `year` and are reused unchanged.
    """
    from xclim.core.calendar import percentile_doy
    bda = replacement_series(base[member], year)
    altered = percentile_doy(bda, window=window, per=per).load()
    parts = {k: v for k, v in per_member.items() if k != member}
    parts[member] = altered
    # keep the member order of the plain threshold, so the floating-point sum
    # is formed in the same order
    ordered = {k: parts[k] for k in per_member}
    return _member_mean(ordered, per_member[member].attrs)


# --------------------------------------------------------------------------
# Temperature percentile indices
# --------------------------------------------------------------------------

def _xclim_temperature_indices():
    from xclim import indices as xi
    # name: (xclim function, data argument, threshold argument, percentile,
    #        extra keyword arguments)
    return {
        'TX90p': (xi.tx90p, 'tasmax', 'tasmax_per', 90, {}),
        'TX10p': (xi.tx10p, 'tasmax', 'tasmax_per', 10, {}),
        'WSDI': (xi.warm_spell_duration_index, 'tasmax', 'tasmax_per', 90,
                 {'window': 6}),
        'TN90p': (xi.tn90p, 'tasmin', 'tasmin_per', 90, {}),
        'TN10p': (xi.tn10p, 'tasmin', 'tasmin_per', 10, {}),
        'CSDI': (xi.cold_spell_duration_index, 'tasmin', 'tasmin_per', 10,
                 {'window': 6}),
    }


TEMPERATURE_INDICES = {
    'tasmax': ('TX90p', 'TX10p', 'WSDI'),
    'tasmin': ('TN90p', 'TN10p', 'CSDI'),
}


def temperature_index(name, data, per):
    """
    One xclim temperature percentile index, annual. `per` may carry a
    _bootstrap dimension, in which case the result is averaged over it, which
    is what bootstrap_func does.
    """
    fn, darg, parg, pct, extra = _xclim_temperature_indices()[name]
    thr = per.sel(percentiles=pct) if 'percentiles' in per.dims else per
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = fn(**{darg: data, parg: thr}, freq='YS', bootstrap=False,
                 **extra)
    if BOOTSTRAP_DIM in out.dims:
        out = out.mean(dim=BOOTSTRAP_DIM, keep_attrs=True)
    return out


def temperature_indices(var, series_by_target, base, inbase_members,
                        base_years, window=WINDOW, log=None,
                        convention=CONVENTIONS[0]):
    """
    Every temperature percentile index for one variable.

    var              : 'tasmax' or 'tasmin'
    series_by_target : {target_key: full series to compute indices on}
    base             : {threshold member: base-period series}, in the
                       pipeline's member order, which sets the pooling offsets
    inbase_members   : {target_key: threshold member id} for the targets whose
                       in-base years must be bootstrapped; any target not
                       listed uses the plain threshold for every year
    base_years       : (first, last) of the base period
    convention       : 'member_mean' or 'pooled', see the module docstring

    Returns ({target_key: {index name: annual DataArray}}, plain threshold).
    """
    if convention not in CONVENTIONS:
        raise ValueError(f'convention must be one of {CONVENTIONS}')
    names = TEMPERATURE_INDICES[var]
    pcts = sorted({_xclim_temperature_indices()[n][3] for n in names})

    if convention == 'member_mean':
        per_member = member_percentiles_doy(base, pcts, window)
        plain = _member_mean(per_member)

        def boot_thr(member, y):
            return member_mean_bootstrapped_doy(base, per_member, member, y,
                                                pcts, window)
    else:
        big = pooled_series(list(base.values()))
        plain = pooled_percentile_doy(None, pcts, window, big=big)

        def boot_thr(member, y):
            return bootstrapped_pooled_doy(base, member, y, pcts, window,
                                           big=big)

    out = {}
    for tkey, series in series_by_target.items():
        res = {n: temperature_index(n, series, plain) for n in names}
        member = inbase_members.get(tkey)
        if member is not None:
            years = [int(y) for y in np.unique(series.time.dt.year.values)
                     if base_years[0] <= int(y) <= base_years[1]]
            for y in years:
                if log:
                    log(f'      bootstrap {var} {tkey} {y}')
                boot = boot_thr(member, y)
                ydata = series.sel(time=str(y))
                for n in names:
                    val = temperature_index(n, ydata, boot)
                    res[n] = _replace_year(res[n], val, y)
        out[tkey] = res
    return out, plain


def _replace_year(annual, value_for_year, year):
    """
    Overwrite one year of an annual series with a recomputed value.
    """
    sel = annual.time.dt.year == year
    if int(sel.sum()) != 1:
        raise ValueError(f'expected exactly one {year} in the annual series')
    v = value_for_year
    if 'time' in v.dims:
        v = v.isel(time=0, drop=True)
    return annual.where(~sel, v.broadcast_like(annual.isel(time=0, drop=True)))


# --------------------------------------------------------------------------
# Precipitation percentile indices
# --------------------------------------------------------------------------

PRECIP_INDICES = ('R95p', 'R99p', 'R95pTOT', 'R99pTOT', 'N95', 'N99')


def wet_day_sample(series_list, wet=WET_DAY_MM):
    """
    All base-period values with pr >= wet, every member, other days NaN.
    Returns dims (..., stack_dim).
    """
    parts = []
    for s in series_list:
        parts.append(s.where(s >= wet).rename({'time': 'stack_dim'})
                     .drop_vars('stack_dim', errors='ignore'))
    return xr.concat(parts, dim='stack_dim')


def wet_day_threshold(series_list, per, wet=WET_DAY_MM, alpha=ALPHA, beta=BETA):
    """
    One threshold per grid cell from all base-period wet days, pooled across
    members, type 8, as Zhang (2011, p. 855) and climpact define it: "the sample
    of all wet days in the base period without consideration of the annual
    cycle". Returns dims (..., percentiles). A cell with no wet day in the base
    period gets NaN, and its indices are NaN rather than zero, since the
    threshold they would be measured against does not exist.
    """
    from xclim.core.utils import calc_perc
    if np.isscalar(per):
        per = [per]
    smp = wet_day_sample(series_list, wet)
    p = xr.apply_ufunc(
        calc_perc, smp,
        input_core_dims=[['stack_dim']],
        output_core_dims=[['percentiles']],
        kwargs={'percentiles': per, 'alpha': alpha, 'beta': beta,
                'copy': True},
        output_dtypes=[smp.dtype],
    )
    return p.assign_coords(percentiles=xr.DataArray(per, dims=('percentiles',)))


def precip_indices(pr, thr, wet=WET_DAY_MM):
    """
    Annual R95p, R99p (mm), R95pTOT, R99pTOT (%), N95, N99 (days).

    pr  : daily precipitation in mm/day
    thr : dims (..., percentiles) with percentiles 95 and 99; may also carry a
          _bootstrap dimension, averaged over at the end as bootstrap_func does

    Exceedance is strict, pr > threshold, as in xclim's days_over_precip_thresh
    and ETCCDI. A day must also be wet; with a threshold drawn from wet days
    that is already implied, and is kept explicit so it stays true if the
    threshold is ever NaN-filled.
    """
    wetd = pr >= wet
    prcptot = pr.where(wetd, 0.0).resample(time='YS').sum()
    out = {}
    for pct in (95, 99):
        t = thr.sel(percentiles=pct, drop=True)
        exceed = (pr > t) & wetd
        tot = pr.where(exceed, 0.0).resample(time='YS').sum()
        cnt = exceed.resample(time='YS').sum().astype('float32')
        # frac = 100.0 * tot / prcptot.where(prcptot > 0)
        # R95p <= PRCPTOT exactly; float32 division can still give 100.00001
        frac = (100.0 * tot / prcptot.where(prcptot > 0)).clip(0.0, 100.0)
        undefined = t.isnull()
        tot = tot.where(~undefined)
        cnt = cnt.where(~undefined)
        frac = frac.where(~undefined)
        if BOOTSTRAP_DIM in tot.dims:
            tot = tot.mean(BOOTSTRAP_DIM)
            cnt = cnt.mean(BOOTSTRAP_DIM)
            frac = frac.mean(BOOTSTRAP_DIM)
        out[f'R{pct}p'] = tot.astype('float32')
        out[f'R{pct}pTOT'] = frac.astype('float32')
        out[f'N{pct}'] = cnt.astype('float32')
    return out


def bootstrapped_wet_thresholds(base, member, year, per, wet=WET_DAY_MM):
    """
    Wet-day thresholds with `year` of `member` replaced by each of that
    member's other base years, using xclim's own build_bootstrap_year_da.
    Returns dims (_bootstrap, ..., percentiles).
    """
    bda = replacement_series(base[member], year)
    others = [v for k, v in base.items() if k != member]
    out = []
    for b in range(bda.sizes[BOOTSTRAP_DIM]):
        out.append(wet_day_threshold(others + [bda.isel({BOOTSTRAP_DIM: b})],
                                     per, wet))
    return xr.concat(out, dim=BOOTSTRAP_DIM)


def member_mean_wet_threshold(base, per, wet=WET_DAY_MM):
    """
    Per-member wet-day thresholds, then the mean over members, NaN skipped.
    Returns (mean, {member: threshold}).
    """
    pm = {k: wet_day_threshold([v], per, wet) for k, v in base.items()}
    return _member_mean(pm, {}), pm


def member_mean_bootstrapped_wet(base, per_member, member, year, per,
                                 wet=WET_DAY_MM):
    """
    Member-mean wet-day thresholds with `year` of `member` replaced by each
    of that member's other base years. Returns dims (_bootstrap, ...,
    percentiles).
    """
    bda = replacement_series(base[member], year)
    altered = wet_day_threshold([bda], per, wet)
    ordered = {k: (altered if k == member else v) for k, v in per_member.items()}
    return _member_mean(ordered, {})


def precipitation_indices(series_by_target, base, inbase_members, base_years,
                          bootstrap=True, wet=WET_DAY_MM, log=None,
                          convention=CONVENTIONS[0]):
    """
    Every precipitation percentile index. Arguments as temperature_indices.
    Returns ({target_key: {index name: annual DataArray}}, plain threshold).
    """
    if convention not in CONVENTIONS:
        raise ValueError(f'convention must be one of {CONVENTIONS}')
    if convention == 'member_mean':
        plain, per_member = member_mean_wet_threshold(base, [95, 99], wet)

        def boot_thr(member, y):
            return member_mean_bootstrapped_wet(base, per_member, member, y,
                                                [95, 99], wet)
    else:
        plain = wet_day_threshold(list(base.values()), [95, 99], wet)

        def boot_thr(member, y):
            return bootstrapped_wet_thresholds(base, member, y, [95, 99], wet)

    out = {}
    for tkey, series in series_by_target.items():
        res = precip_indices(series, plain, wet)
        member = inbase_members.get(tkey) if bootstrap else None
        if member is not None:
            years = [int(y) for y in np.unique(series.time.dt.year.values)
                     if base_years[0] <= int(y) <= base_years[1]]
            for y in years:
                if log:
                    log(f'      bootstrap pr {tkey} {y}')
                vals = precip_indices(series.sel(time=str(y)), boot_thr(member, y),
                                      wet)
                for n in PRECIP_INDICES:
                    res[n] = _replace_year(res[n], vals[n], y)
        out[tkey] = res
    return out, plain
