# Branch: francis-percentile-indices

Adds the percentile half of the ETCCDI set and a user-facing "compute one
metric" notebook, matching two roadmap items in the README (extra metrics;
monthly-resolution indices).

## Contents

- `multimodel_etccdi.py` - one module unifying all four models (CESM2-WACCM,
  UKESM1-1, MIROC-ES2H, E3SMv3) behind a single API, computing the eight
  percentile indices (R95p, R99p, TX90p, TX10p, TN90p, TN10p, WSDI, CSDI):
  per-model loaders in native hub layouts, day-of-year percentile thresholds
  (SSP2-4.5 2020-2039, 5-day window, wet-day masking for precipitation), the
  Lee et al. three-way warming/SAI/combined decomposition, Welch significance,
  and a dataset writer targeting the shared S3 layout. Supports annual (YS)
  and monthly (MS) frequency.
- `compute_any_index.ipynb` - press-go notebook: choose an index, model(s),
  scenario, members, and frequency in one cell; get the field back, with a
  cost-guidance table, a quick-look map, an optional write to S3, and a recipe
  for adding new (e.g. ecological) indicators via INDEX_REGISTRY.

## Data-quality rules the code enforces

- E3SMv3 daily max/min temperature are byte-identical at the source, so E3SM
  is excluded from temperature indices (precipitation is unaffected).
- WSDI/CSDI are annual-only by construction (spells cross month boundaries);
  the code refuses them at monthly frequency.

Percentile thresholds are identical across scenarios, so counts are comparable
between SSP2-4.5, G6-1.5K-SAI, and G6-1.5K-HiLLA.
