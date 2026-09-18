# Computing an ETCCDI index from this repository

Anyone with access to the Reflective Cloud Hub can compute any of the twelve
percentile-based and fixed-threshold indices without writing code, using
`compute_any_index.ipynb`.

## 1. Clone the repository on the hub

Open a terminal in JupyterLab and run:

    git clone https://github.com/ReflectiveCloud/etccdi-indices.git
    cd etccdi-indices

(While the percentile work is still on a branch, add `-b francis-percentile-indices`
to the clone command, or `git checkout francis-percentile-indices` afterwards.)

## 2. Open the notebook

In the JupyterLab file browser, navigate into `etccdi-indices/` and open
`compute_any_index.ipynb`. Choose the standard pangeo Python 3 kernel.

The notebook imports `multimodel_etccdi.py`, which sits beside it in the same
directory, so the import resolves with no path setup. Nothing to install: the
hub image already carries xclim, xarray, s3fs, and the rest.

## 3. Edit one cell and run

Only the USER CHOICES cell needs editing:

    INDEX     = 'TX90p'     # percentile: R95p R99p TX90p TX10p TN90p TN10p WSDI CSDI
                            # fixed-threshold: RX1D RX5D GSL GDD
    MODELS    = ['CESM']    # 'CESM', 'UKESM', 'MIROC', 'E3SM'
    SCENARIO  = 'HiLLA'     # 'SSP245', 'SAI', 'HiLLA'
    MEMBERS   = None        # None = all available
    FREQ      = 'YS'        # 'YS' annual, 'MS' monthly, 'QS-DEC' seasonal
    SEASON    = 'DJF'       # used with FREQ='QS-DEC': 'DJF', 'MAM', 'JJA' or 'SON'

    MAKE_PLOT = True        # quick-look map at the end
    SAVE_PLOTS= False       # True writes each figure as PNG

Then run the remaining subsequent cells. The notebook prints what it is loading, computes the index, and
plots the ensemble mean.

## 4. What you get back

    results[model]['per_member']   # {member: DataArray}
    results[model]['ens_mean']     # ensemble mean

At `FREQ='YS'` these are `(lat, lon)`: mean days per year.
At `FREQ='MS'` they are `(month, lat, lon)`: mean days per calendar month, the
seasonal cycle. Summing over `month` recovers the annual field.
At `FREQ='QS-DEC'` they are `(lat, lon)` for the season named in `SEASON`.

Each field carries `units`, `long_name` and a `comment` describing how it was
computed, so a file written to disk explains itself.

With `MAKE_PLOT = True` you also get a map with coastlines and national borders,
and for monthly runs a twelve-panel climatology plus the cosine-latitude
weighted seasonal cycle. `SAVE_PLOTS = True` writes each figure as a PNG named
for the index, model, scenario and window.

## 5. Note these costs before you pick a big configuration

The expensive step is the baseline percentile threshold, computed once per model
per variable and cached for the session. A single model with all members takes
roughly 30-40 minutes on the hub for a percentile index and four models scale
accordingly. Memory grows with the member count, since each member holds its
daily record.

Measured on the hub as of September 2026: one threshold from three MIROC baseline
members is about six minutes and peaks near 2 GB when members are loaded one at
a time. Loading every member before taking any percentile, which is what
`compute_baseline_threshold` does, raises the peak far higher, because
`percentile_doy` builds a rolling-window construct several times the size of its
input. A twelve-threshold run structured member-by-member took 73 minutes and
stayed under 2 GB. If a run dies without a traceback, memory is the first thing
to suspect.

Keep the kernel alive between runs on the same model and variable: the second
index is much faster because the threshold is reused.

The four fixed-threshold indices (`RX1D`, `RX5D`, `GSL`, `GDD`) skip the
percentile step entirely and finish faster. Monthly and seasonal frequency cost
no more than annual, since the threshold is shared.

E3SM's first access to a variable triggers cubed-sphere regridding and caching
to the bucket, which is slow once and fast afterwards.

## 6. Rules the code enforces to have it running smoothly

- E3SM is skipped for temperature indices built on daily temperature extremes
  (`TX90p`, `TX10p`, `TN90p`, `TN10p`, `WSDI`, `CSDI`): its archived daily max
  and min are identical to the daily mean at source, so no true daily extremes
  exist. Precipitation indices and the `tas`-based `GSL` and `GDD` keep all four
  models.

- `WSDI`, `CSDI` and `GSL` are refused at `FREQ='MS'` and `FREQ='QS-DEC'`: their
  spells and seasons cross month boundaries, so a sub-annual count is not the
  ETCCDI quantity.

- Incomplete bins are masked to NaN at every frequency, annual included. xclim
  returns a value for a partial bin rather than NaN, so a truncated final year
  would otherwise be written as a real year with a count low in proportion to
  its missing days. Two members in the archive needed this: MIROC-ES2H
  G6-1.5K-HiLLA r02 stops 62 days into 2084, and UKESM1-1 G6-1.5K-HiLLA
  r3i1p1f2 has 270 of 360 days of `pr` in the same year.

- Percentile thresholds come from SSP2-4.5 over the baseline window, whatever
  scenario you compute, so counts are comparable across SSP2-4.5, G6-1.5K-SAI,
  and G6-1.5K-HiLLA. The window defaults to 2020-2039, which is what the
  released dataset uses. `BASELINE_BY_MODEL` holds the newer 2015-2034 window
  (MIROC 2020-2034, because its SSP2-4.5 begins in 2020); pass it explicitly
  through `baseline_period` and write to a separate output root, since indices
  on two baselines are not comparable.

- `GSL` uses a hemisphere-aware wrapper, since a southern-hemisphere growing
  season straddles the calendar year: northern gridpoints use a 1 July mid-date
  on a calendar year, southern gridpoints a 1 January mid-date on a July-to-June
  year. Without this, southern values are pinned by the date convention rather
  than by temperature.

## 7. Sharing a result

Set `SAVE_TO_BUCKET = True` to write per-member fields into the shared S3 layout
via the pipeline's own writer:

    .../ETCCDI_indices_annual/{model}/{scenario}/{member}/{INDEX}.nc     # YS
    .../ETCCDI_indices_monthly/{model}/{scenario}/{member}/{INDEX}.nc    # MS
    .../ETCCDI_indices_djf/{model}/{scenario}/{member}/{INDEX}.nc        # QS-DEC, DJF
    .../ETCCDI_indices_mam|_jja|_son/...                                 # other seasons

Each season gets its own tree, so a DJF file cannot be mistaken for an annual
one. Only do this for configurations worth sharing.

A DJF value is labelled with the year its December falls in: a bin labelled 2070
covers December 2070 with January and February 2071. The convention is recorded
in the file attributes because it is invisible in the data.

## 8. Adding your own indicator

`INDEX_REGISTRY` in `multimodel_etccdi.py` maps an index name to its variable,
its percentile (or `None` for a fixed threshold), the xclim function and optional
`attrs` supplying `long_name` and `comment`. Add an entry and the notebook picks
it up with no other change: it appears in the printed index list and every cell
works as before.

`MODELS_FOR()` keys on the variable, so an index on `pr` or `tas` automatically
keeps all four models and one on `tasmax` or `tasmin` automatically drops E3SM.
If the index is spell-based, add it to `SPELL_INDICES` so the frequency guard
applies.

If your indicator needs a variable the loaders do not carry yet (wind, humidity),
that is a loader addition; open an issue or a pull request.

## 9. Versions

The pipeline pins `xclim==0.54.0`. Season-length partial-bin behaviour is version
dependent and `GSL` was validated on that version only. If the hub image moves
ahead, re-run the southern-hemisphere check in the repository before trusting
`GSL` output.

xclim 0.54.0 calls `xarray.core.ops.get_op`, which newer xarray has removed;
on xarray 2026.7.0 every percentile index fails with `AttributeError: module
'xarray.core' has no attribute 'ops'`. `xclim_compat.py` restores that one
function for the six comparison operators xclim actually asks for, and
self-tests them before returning. Import it before computing anything on an
environment where xarray is newer than the hub image:

    import xclim_compat; xclim_compat.patch()

Pinning xarray to the hub image's version is the cleaner fix where that is
possible.

## 10. Running somewhere other than the hub

The reading layer assumes the Reflective bucket. A GLADE port for NCAR's Derecho
and Casper lives in the analysis repository rather than here: it keeps this
module's index definitions, thresholds and registry, and replaces only the file
access. Its per-model baselines and its corrections for CAM end-of-interval
timestamps and overlapping file segments are described there.
