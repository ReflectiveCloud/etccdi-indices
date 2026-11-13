# etccdi-indices

Annual ETCCDI climate extremes indices for the GeoMIP G6-1.5K experiments in
four models, computed with [xclim](https://xclim.readthedocs.io). The code runs
on the Reflective Cloud Hub, reading the Reflective bucket, and on NCAR
Derecho/Casper, reading GLADE. The index code is the same at both sites; only
the reading differs.

| Model | SSP2-4.5 | G6-1.5K-SAI | G6-1.5K-HiLLA |
|---|---|---|---|
| CESM2-WACCM | r6-r10 | r1-r3 | r1-r3 |
| UKESM1-1 | r2, r3, r12 | r1-r3 | r1-r3 |
| MIROC-ES2H | r01-r10 | r01-r10 | r01-r10 |
| E3SMv3 (Hub only, no tasmax or tasmin indices) | r1-r3 | r1-r3 | r1-r3 |

## The indices

Thirty annual indices, every year each run has from 2015 to 2084.

Twelve percentile indices, relative to an SSP2-4.5 reference period:

| Variable | Indices |
|---|---|
| tasmax | TX90p, TX10p, WSDI |
| tasmin | TN90p, TN10p, CSDI |
| pr | R95p, R99p (mm), R95pTOT, R99pTOT (% of PRCPTOT), N95, N99 (days) |

Eighteen fixed-threshold indices:

| Variable | Indices |
|---|---|
| tasmax | SU, ID, TXX, TXN |
| tasmin | FD, TN (tropical nights), TNX, TNN |
| pr | PRCPTOT, SDII, CDD, CWD, R10MM, R20MM, RX1D, RX5D |
| tas | GSL, GDD |

TXX, TXN, TNX and TNN are in kelvin. A year with fewer than 360 days is written
as NaN.

## Method

The percentile threshold is taken over the SSP2-4.5 threshold members pooled
into one sample. The default reference period is 2020-2039, the period whose
global mean temperature the G6-1.5K runs target (Lee et al. 2026); any other
can be given with `--baseline`. Temperature thresholds use a 5-day centred
day-of-year window. Precipitation follows Zhang et al. (2011): one threshold
per grid cell from all wet days (at least 1 mm/day) of the reference period.

In-base years of the threshold members are computed with the Zhang et al.
(2005) bootstrap by default (`--no-bootstrap` turns it off). G6 years are never
bootstrapped, since no G6 data enter the threshold.

The fixed-threshold indices are defined once, in
`multimodel_etccdi.INDEX_REGISTRY`, and computed by the same call from the
script and from the notebooks.

## Files

| File | Purpose |
|---|---|
| `build_indices.py` | Builds all 30 indices for one model, at either site |
| `multimodel_etccdi.py` | Index registry, bucket readers, pooled thresholds |
| `etccdi_bootstrap.py` | Percentile indices and the in-base bootstrap |
| `naming.py` | Output directory and file names |
| `derecho_layout.py`, `derecho_io.py` | GLADE paths, members and readers (Derecho only) |
| `xclim_compat.py` | Lets xclim 0.54.0 run on xarray 2026.7 |
| `calculate_yearly_indices_{CESM,UKESM,MIROC,E3SM}.ipynb` | The eighteen fixed-threshold indices for one model, through `build_indices` |

## Running

Check first. The preflight reads one year of one member, lists the members to
be written and checks the output root, and writes nothing:

```bash
python build_indices.py --site hub --model UKESM --out-root OUT --preflight
```

One latitude band, for a quick look or a test:

```bash
python build_indices.py --site hub --model UKESM --lat-range 40 45 --out-root OUT
```

The globe, in bands of 20 latitude rows, joined at the end. A stopped run
resumes at the next unfinished band:

```bash
python build_indices.py --site hub --model UKESM --full-grid --band-rows 20 \
    --workers 4 --out-root OUT --output both
```

On Derecho or Casper, inside a batch or interactive job with enough memory
for the members of one scenario:

```bash
export ETCCDI_SITE=derecho OMP_NUM_THREADS=1
python build_indices.py --model MIROC --full-grid --workers 16 --out-root /glade/...
```

The in-base bootstrap sets the run time. Each band is split into blocks of
`--cells-per-task` grid cells (default 96), and the blocks run on `--workers`
processes. Give `--workers` the number of CPUs in the job, and choose
`--cells-per-task` so the number of blocks is close to a multiple of the
workers. One MIROC block of 86 cells (ten members) takes about 50 minutes on
one CPU for each of tasmax and tasmin. `OMP_NUM_THREADS=1` stops numpy in each
worker from starting threads of its own.

Main options:

| Option | Default | |
|---|---|---|
| `--site hub\|derecho` | `$ETCCDI_SITE`, else hub | where the daily data are read from |
| `--model` | | CESM, UKESM, MIROC or E3SM |
| `--indices`, `--variables` | all 30 | a subset |
| `--scenarios` | SSP245 HiLLA SAI | |
| `--baseline Y0 Y1` | 2020 2039 | SSP2-4.5 years of the percentile threshold |
| `--bootstrap`, `--no-bootstrap` | bootstrap | in-base years |
| `--output per-index\|merged\|both` | per-index | see below |
| `--lat-range LO HI`, `--full-grid` | | one band, or the globe in bands |
| `--workers` | `$NCPUS`, else 3 | processes for the percentile indices; set it to the job's CPUs |

From Python or a notebook, `build_indices.build(model, out_root, ...)` takes
the same options as keyword arguments.

## Output

```
<out-root>/<MODEL>/<SCENARIO>/<member>/<INDEX>_<MODEL>_<SCENARIO>_<member>_<y0>-<y1>.nc
<out-root>/_thresholds/<MODEL>/<var>_threshold_pooled.nc
```

`--output merged` writes two files per member instead, the twelve percentile
indices (`ETCCDI_bootstrap_...`, or `ETCCDI_percentile_...` without the
bootstrap) and the eighteen fixed ones (`ETCCDI_fixed_...`). The per-index
files they are built from go to `<out-root>_per_index`. `--output both`
writes the per-index tree under `<out-root>` and the merged files under
`<out-root>_merged`. Every merged file is read back and checked against its
per-index files before it is kept.

The released data are on the Reflective bucket:

```
s3://reflective-persistent-prod-large/fafrifa/ETCCDI/ETCCDI_indices_annual_bootstrap/
s3://reflective-persistent-prod-large/fafrifa/ETCCDI/ETCCDI_indices_annual_bootstrap_merged/
```

## Notes on the data

- **CESM2-WACCM.** CAM stamps a daily mean at the end of its day; the readers
  move each stamp to the midpoint of its `time_bnds`. SSP2-4.5 members 1-5
  have a defect in daily maximum and minimum temperature, so r6-r10 are used,
  and r6-r8 build the threshold.
- **UKESM1-1.** A 360-day calendar. All runs have stratospheric heterogeneous
  chemistry on. The G6 daily records begin in 2034.
- **MIROC-ES2H.** SSP2-4.5 begins in 2020. HiLLA r02, r05 and r10 end
  part-way through 2084.
- **E3SMv3.** Hub only, read from regridded windows on a 180 x 360 grid, with
  each time stamp moved back 12 hours to the middle of its day. Its daily
  maximum and minimum temperature equal the daily mean in this archive, so no
  tasmax or tasmin index is computed.

## Requirements

Python 3.11 or later and the packages in `requirements.txt`. xclim is pinned
at 0.54.0, the version every released file was built with.

## References

Lee, W. R., and coauthors, 2026: Atmos. Chem. Phys., 26, 7463-7483,
https://doi.org/10.5194/acp-26-7463-2026.

Zhang, X., G. Hegerl, F. W. Zwiers, and J. Kenyon, 2005: Avoiding
inhomogeneity in percentile-based indices of temperature extremes. J. Climate,
18, 1641-1651.

Zhang, X., L. Alexander, G. C. Hegerl, P. Jones, A. Klein Tank, T. C.
Peterson, B. Trewin, and F. W. Zwiers, 2011: Indices for monitoring changes in
extremes based on daily temperature and precipitation data. WIREs Climate
Change, 2, 851-870.

## Credits

The original notebooks and the zarr store are by Cindy Wang and Alistair
Duffey. The percentile pipeline, the bootstrap and the multi-site build are by
Francis Osei Tutu Afrifa (Reflective summer research internship, 2026).
