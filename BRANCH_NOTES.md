# Branch notes: francis-percentile-indices

What this branch changes relative to `ReflectiveCloud/etccdi-indices`, as of
30 October 2026. The README describes how to run everything.

## New

- `build_indices.py` builds all 30 annual indices for one model, at either
  site (`--site hub` reads the Reflective bucket, `--site derecho` reads
  GLADE). Options set the SSP2-4.5 reference period (`--baseline`), the
  in-base bootstrap (`--bootstrap`, `--no-bootstrap`) and the output layout
  (`--output per-index|merged|both`).
- `etccdi_bootstrap.py` computes the twelve percentile indices: pooled
  thresholds over the SSP2-4.5 reference members, the Zhang et al. (2005)
  in-base bootstrap, and precipitation thresholds from wet days only (Zhang et
  al. 2011).
- `naming.py` holds the output directory and file names.
- `derecho_layout.py` and `derecho_io.py` hold the GLADE paths, members and
  readers used with `--site derecho`.

## Changed

- `multimodel_etccdi.py`: index registry for the eighteen fixed-threshold
  indices, bucket layouts for all four models (UKESM1-1 SSP2-4.5 het-chem-on
  suites, UKESM1-1 G6-1.5K-SAI suite files, ten MIROC-ES2H members, a complete
  copy of UKESM1-1 G6-1.5K-HiLLA r3i1p1f2), and masking of years with fewer
  than 360 days.
- The four `calculate_yearly_indices_*.ipynb` notebooks compute the eighteen
  fixed-threshold indices through `build_indices.build`, so the notebooks and
  the script run the same code.

## Checked

The script reproduces the released trees on both systems. On the Hub, UKESM1-1
at 40-45N agrees with `ETCCDI_indices_annual_bootstrap_merged` to within
float32 rounding. On Casper, MIROC-ES2H at 40-45N gives 843 of 903 files
identical and the remaining 60 (PRCPTOT and SDII) within float32 rounding.

## Released data

    s3://reflective-persistent-prod-large/fafrifa/ETCCDI/ETCCDI_indices_annual_bootstrap/
    s3://reflective-persistent-prod-large/fafrifa/ETCCDI/ETCCDI_indices_annual_bootstrap_merged/
