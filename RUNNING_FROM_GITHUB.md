# Running the ETCCDI indices from this repository

There are two ways in. `compute_any_index.ipynb` reads or computes one index and plots
it, with no code to write. `build_indices.py` builds all 30 indices for a
model and writes them as files, which is how the released dataset was made.
Both run on the Reflective Cloud Hub and on NCAR's Derecho and Casper, and
both use the same index code.

## 1. Get the code

On the Hub, open a terminal in JupyterLab:

    git clone https://github.com/ReflectiveCloud/etccdi-indices.git
    cd etccdi-indices

Until the percentile work is merged, clone the fork's branch instead:

    git clone -b francis-percentile-indices https://github.com/OT-Afrifa/etccdi-indices.git

The Hub image already has xclim, xarray, s3fs and the rest. Elsewhere, install
`requirements.txt` into Python 3.11 or later. On Derecho or Casper, also set
`export ETCCDI_SITE=derecho`, so the code reads GLADE instead of the bucket.

## 2. The notebook: look at or compute one index

Open `compute_any_index.ipynb` with the standard Python 3 kernel, edit the
USER CHOICES cell and run all cells.

    INDEX     = 'TX90p'                     # any of the 30
    MODELS    = ['CESM', 'UKESM', 'MIROC']  # 'CESM', 'UKESM', 'MIROC', 'E3SM'
    SCENARIOS = ['SSP245', 'SAI', 'HiLLA']
    MODE      = 'released'                  # or 'compute'
    REGION    = None                        # or (lat0, lat1, lon0, lon1)
    WINDOW    = (2050, 2069)                # years the change maps average over

`MODE='released'` reads the released annual files and takes seconds per
member. `MODE='compute'` computes the index over `LAT_BAND`: percentile indices
through `build_indices.build` (so `BASELINE` and `BOOTSTRAP` can be changed),
fixed-threshold indices through the registry call the build script uses (so
`FREQ` can be monthly or seasonal).

You get `series[model][scenario][member]`, a `(time, lat, lon)` DataArray per
member, then a time series over `REGION` with every member drawn and maps of
the change over `WINDOW`: SSP2-4.5 against 2020-2039, and each G6 experiment
against SSP2-4.5 in the same years. The maps carry no significance test.

## 3. The build script: all indices as files

Check first. The preflight reads one year, lists what would be written and
writes nothing:

    python build_indices.py --site hub --model UKESM --out-root OUT --preflight

Then one latitude band for a quick look, or the globe in bands:

    python build_indices.py --site hub --model UKESM --lat-range 40 45 --out-root OUT
    python build_indices.py --site hub --model UKESM --full-grid --workers 4 \
        --out-root OUT --output both

| Option | Default | |
|---|---|---|
| `--site hub\|derecho` | `$ETCCDI_SITE`, else hub | where the daily data are read from |
| `--indices`, `--variables` | all 30 | a subset |
| `--baseline Y0 Y1` | 2020 2039 | SSP2-4.5 years of the percentile threshold |
| `--bootstrap`, `--no-bootstrap` | bootstrap | in-base years of the threshold members |
| `--output per-index\|merged\|both` | per-index | one file per index, or two per member |
| `--workers` | `$NCPUS`, else 3 | processes for the percentile indices |

From Python, `build_indices.build(model, out_root, ...)` takes the same options
as keyword arguments. The README has the full list and the output layout.

## 4. What it costs

The in-base bootstrap sets the run time of the percentile indices: each block
of grid cells is resampled once per reference year for every threshold
member. One MIROC-ES2H block of 86 cells with ten members takes about 50
minutes on one CPU for each of tasmax and tasmin. Give `--workers` the number
of CPUs in the job, choose `--cells-per-task` so the number of blocks is close
to a multiple of the workers, and set `OMP_NUM_THREADS=1`. The script warns
when it would bootstrap on a single worker.

G6 scenarios need no bootstrap, so they take minutes, as do the
fixed-threshold indices and `--no-bootstrap` runs. A full-grid run resumes at
the next unfinished band if it stops.

## 5. Rules the code enforces

- **E3SMv3 has no tasmax or tasmin indices.** Its archived daily maximum and
  minimum temperature equal the daily mean, so there are no daily extremes to
  count. Its precipitation indices and the daily-mean GSL and GDD are computed.
  Its time stamps, written at the end of each day, are moved back 12 hours to
  the middle of the day before any index is taken.
- **Short years are NaN.** A year with fewer than 360 days is masked at every
  frequency. xclim returns a value for a partial year rather than NaN, which
  would otherwise enter the files as a real year with counts low in proportion
  to the missing days. MIROC-ES2H G6-1.5K-HiLLA r02, r05 and r10 end part-way
  through 2084, so that year is NaN for them.
- **Thresholds come from SSP2-4.5 only.** Whatever scenario is computed, its
  percentile indices are counted against the SSP2-4.5 threshold, so values are
  comparable across SSP2-4.5, G6-1.5K-SAI and G6-1.5K-HiLLA. Precipitation
  thresholds use wet days only (at least 1 mm/day), following Zhang et al.
  (2011).
- **Spell and season indices are annual.** WSDI, CSDI and GSL are refused at
  monthly and seasonal frequency, because their spells and seasons cross the
  bin boundaries.
- **GSL is hemisphere-aware.** Northern grid points use a 1 July mid-date on a
  calendar year and southern grid points a 1 January mid-date on a July-to-June
  year, so a southern growing season is not cut at the new year.
- **Seasonal labels.** At `FREQ='QS-DEC'`, a DJF value is labelled with the
  year its December falls in: 2070 covers December 2070 with January and
  February 2071.

## 6. Adding your own indicator

A fixed-threshold indicator is one entry in `multimodel_etccdi.INDEX_REGISTRY`:
the daily variable, `percentile: None`, the xclim function and its keyword
arguments, and optional `attrs` with the long name and comment. The build
script and the notebook pick it up with no other change. If it is spell-based,
add it to `SPELL_INDICES` so the frequency guard applies.

A percentile indicator belongs in `etccdi_bootstrap.py`, so that it gets the
pooled threshold and the bootstrap. The registry still holds older entries
named `R95p` and `R99p` that count wet days; the released `R95p` and `R99p`
are totals in mm, and the counts are `N95` and `N99`.

A variable the readers do not carry yet (wind, humidity) needs a reader change
in `multimodel_etccdi.py` for the Hub and `derecho_io.py` for GLADE; open an
issue or a pull request.

## 7. Versions

xclim is pinned at 0.54.0, the version every released file was built with.
Partial-season behaviour in GSL depends on the xclim version, and GSL was
checked on 0.54.0 only.

xclim 0.54.0 calls `xarray.core.ops.get_op`, which xarray 2026.7 removed.
`xclim_compat.py` restores it for the comparison operators xclim uses and tests
them before returning. `build_indices.py` loads it automatically, so nothing
needs importing by hand.
