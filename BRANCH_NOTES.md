# Branch: francis-percentile-indices

Adds the percentile half of the ETCCDI set, four fixed-threshold indices, extra
output frequencies (monthly and all four seasons) and a user-facing "compute one
metric" notebook, covering two roadmap items in the README (extra metrics;
monthly-resolution indices).

## Contents

- `multimodel_etccdi.py` : one module unifying all four models (CESM2-WACCM6,
  UKESM1-1, MIROC-ES2H, E3SMv3) behind a single API: per-model loaders in native
  hub layouts, day-of-year percentile thresholds (SSP2-4.5, 5-day window,
  wet-day masking for precipitation), the Lee et al. (2026) three-way
  warming/SAI/combined decomposition, Welch significance, and a dataset writer
  targeting the shared S3 layout.

- `compute_any_index.ipynb` : press-go notebook where you choose an index,
  model(s), scenario, members, frequency and season in one cell. You get the
  field back with a cost-guidance table, a Cartopy quick-look map, an optional
  write to S3, and a recipe for adding new indicators via `INDEX_REGISTRY`.

- `xclim_compat.py` : restores `xarray.core.ops.get_op` for environments whose
  xarray is newer than the hub image. See Versions below.

## What this branch adds to the existing indices module

**Twelve (12) more indices:** The eight percentile indices (`R95p`, `R99p`, `TX90p`,
`TX10p`, `TN90p`, `TN10p`, `WSDI`, `CSDI`) plus four fixed-threshold ones:
`RX1D` and `RX5D` (1-day and 5-day precipitation maxima, following Cindy Wang's
definitions and short names so the two index sets merge without renaming), `GSL`
(growing season length) and `GDD` (growing degree days). Fixed-threshold indices
need no baseline percentile, so they run far faster and keep all four models.

**Output frequency-** `SUPPORTED_FREQ = ('YS', 'MS', 'QS-DEC')`: annual, monthly,
and seasonal. The frequency threads through `compute_index_for_members`,
`run_comparison`, `run_lee_framework` and both cached wrappers, and appears in
the comparison cache filenames so a seasonal run cannot silently load an annual
result. The tag is empty for `YS`, so caches written before the argument existed
are still found.

**All four DJF, MAM, JJA & SON seasons:** `QS-DEC` already bins every year into DJF, MAM, JJA and
SON, each stamped with the first month of its bin, so we've added a `season` argument
selecting the stamp month. DJF is binned identically however the option is used. One
completeness floor covers all four, thus every complete season is 90 to 92 days on the
calendars in this ensemble, and the largest partial stub a record end can produce
is January plus February at 59 days, so `min_days=85` accepts every complete
season and rejects every stub. Each season writes to its own output tree
(`ETCCDI_indices_djf`, `_mam`, `_jja`, `_son`). A DJF bin is labelled with the
year its December falls in, and that convention is written into the file
attributes because it is invisible in the data.

**Per-model baseline:** `BASELINE_BY_MODEL` and a `baseline_period` argument on
`write_index_dataset`, defaulting to the existing 2020-2039 window so nothing
already written changes basis by accident. `The newer window is 2015-2034, with
MIROC on 2020-2034 because its SSP2-4.5 begins in 2020` because 2015-2019 was never run, since the 2020 SSP2-4.5 was
initialised separately from combined prior runs plus a short spin-up. MIROC's
day-of-year percentiles therefore rest on 75 samples against 100, and this gap is
widest for `R99p`, where only wet days enter the estimate. Indices on two
baselines are not comparable and are written to separate output roots.

**Hemisphere-aware growing season-** `growing_season_length_global` wraps xclim's
`growing_season_length`: northern gridpoints use a 1 July mid-date on a calendar
year, southern gridpoints a 1 January mid-date on a July-to-June year, relabelled
to January stamps. Without it, southern values are pinned by the mid-date
constraint rather than by temperature (a validation on real data gave 181 days
for every year at a Patagonian point, which is exactly 1 January to 1 July on a
365-day calendar).

**Registry attributes:** Entries may carry an `attrs` dict supplying `long_name`
and `comment`, filled per computation from `extra_kwargs` so a field computed at
a non-default threshold is labelled with the threshold actually used. xclim's
`indices` layer sets `units` and nothing else, so without this a released file
carries no human-readable name.

**Subtropical comparison configs:** `subtropical` and `subtropical_combined` give
the 2-to-3 and 1-to-3 Lee contrasts for G6-1.5K-SAI. The existing `sai` key
carries G6-1.5K-HiLLA despite its name, so before this the subtropical scenario
was covered only by the direct SAI-versus-HiLLA contrast. Its labels are now
scenario-explicit.

**Assessment window argument:** `run_lee_framework_cached` takes `assessment`, so
all four models can be placed on the shared 2050-2069 window rather than each
using its own longest usable period. The window is already in the cache name, so
the two periods coexist.

## Data-quality rules the code enforces

- E3SMv3 daily max/min temperature are byte-identical to the daily mean at the
  source, so E3SM is excluded from temperature extremes. `MODELS_FOR()` keys on
  the variable, so precipitation and `tas`-based indices keep all four models
  automatically.

- `WSDI`, `CSDI` and `GSL` are annual-only by construction, since spells and
  seasons cross month boundaries. The code refuses them at `MS` and `QS-DEC`.

- Percentile thresholds are identical across scenarios for a given baseline, so
  counts are comparable between SSP2-4.5, G6-1.5K-SAI, and G6-1.5K-HiLLA.

- `_complete_bins` masks partial bins at every frequency, annual included,
  because xclim returns a value rather than NaN for a partial bin and a truncated
  year or season would otherwise be counted as a real low value. The annual path
  originally had no such mask, since no partial year existed in the archive when
  it was written. Two later arrivals did have one: MIROC-ES2H G6-1.5K-HiLLA r02
  stops 62 days into 2084 and UKESM1-1 G6-1.5K-HiLLA r3i1p1f2 has 270 of 360
  days of `pr` in the same year. Both fall inside the 2065-2084 assessment window
  for those models; the MIROC case biased that member's window mean low by about
  1.4 days per year, measured against an independent build of the same index.
  Sixteen files were regenerated on 15 September 2026, with the pre-fix copies
  kept under `ETCCDI_indices_annual_prebugfix`.

## Known limitations (currently not fixed here)

`run_comparison_cached` builds its cache filename from model, index and
comparison type, with no assessment window. That is safe while every call site
passes the same window, which is currently true, but two windows would collide.
Adding the window would rename existing caches and force a recompute, so it is
left with a comment on the line.

CAM h1 output labels each day with the end of its averaging interval, so 31
December carries a 1 January stamp of the following year. Every CESM member on
the bucket therefore shows a one-day final year, at 2070, 2085 or 2100 depending
on the run. Those years fall outside both the baseline and CESM's 2050-2069
assessment window, so no released number is affected, and the complete-bin mask
now drops the stray bin whenever those files are rebuilt. The related one-day
shift of the whole CESM record is uncorrected here; the GLADE port in the
analysis repository corrects it from `time_bnds`.

## Versions

`xclim==0.54.0` is pinned in `requirements.txt`. Season-length partial-bin
behaviour is version dependent and `GSL` was validated on that version only.

xclim 0.54.0 calls `xarray.core.ops.get_op`, which newer xarray has removed. On
xarray 2026.7.0 every percentile index fails at `percentile_doy` with
`AttributeError: module 'xarray.core' has no attribute 'ops'`. `xclim_compat.py`
restores that function for the six comparison operators xclim asks for and
self-tests them before returning; pinning xarray to the hub image's version is
the cleaner fix where that is possible.
