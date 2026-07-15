# Computing an ETCCDI index from this repository

Anyone with access to the Reflective Cloud Hub can compute any of the eight
percentile indices without writing code, using `compute_any_index.ipynb`.

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

    INDEX     = 'TX90p'     # R95p R99p TX90p TX10p TN90p TN10p WSDI CSDI
    MODELS    = ['CESM']    # 'CESM', 'UKESM', 'MIROC', 'E3SM'
    SCENARIO  = 'HiLLA'     # 'SSP245', 'SAI', 'HiLLA'
    MEMBERS   = None        # None = all available
    FREQ      = 'YS'        # 'YS' annual, 'MS' monthly (seasonal cycle)

Then Run All. The notebook prints what it is loading, computes the index, and
plots the ensemble mean.

## 4. What you get back

    results[model]['per_member']   # {member: DataArray}
    results[model]['ens_mean']     # ensemble mean

At `FREQ='YS'` these are `(lat, lon)`: mean days per year.
At `FREQ='MS'` they are `(month, lat, lon)`: mean days per calendar month, the
seasonal cycle. Summing over `month` recovers the annual field.

## 5. Cost, before you pick a big configuration

The expensive step is the baseline percentile threshold, computed once per
model per variable and cached for the session. A single model with all members
takes roughly 30-40 minutes on the hub; four models scale accordingly. Memory
grows with the member count, since each member holds its daily record. Keep the
kernel alive between runs on the same model and variable: the second index is
much faster because the threshold is reused.

E3SM's first access to a variable triggers cubed-sphere regridding and caching
to the bucket, which is slow once and fast afterwards.

## 6. Rules the code enforces, so you cannot get them wrong

- E3SM is skipped for temperature indices: its daily max and min are identical
  to the daily mean at source, so no true daily extremes exist. Precipitation
  is unaffected.
- WSDI and CSDI are refused at `FREQ='MS'`: their spells cross month
  boundaries, so a monthly count is not the ETCCDI quantity.
- Percentile thresholds always come from SSP2-4.5 2020-2039, whatever scenario
  you compute, so counts are comparable across SSP2-4.5, G6-1.5K-SAI, and
  G6-1.5K-HiLLA.

## 7. Sharing a result

Set `SAVE_TO_BUCKET = True` to write per-member fields into the shared S3
layout via the pipeline's own writer:

    .../ETCCDI_indices_annual/{model}/{scenario}/{member}/{INDEX}.nc     # YS
    .../ETCCDI_indices_monthly/{model}/{scenario}/{member}/{INDEX}.nc    # MS

Only do this for configurations worth sharing.

## 8. Adding your own indicator

`INDEX_REGISTRY` in `multimodel_etccdi.py` maps an index name to its variable,
percentile, and xclim function. Add an entry and the notebook picks it up with
no other change: it appears in the printed index list and every cell works as
before. If your indicator needs a variable the loaders do not carry yet (wind,
humidity), that is a loader addition; open an issue or a pull request.
