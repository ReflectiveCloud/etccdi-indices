# etccdi-indices
Code for creating ETCCDI indices from G6-1.5K simulations on the Reflective Cloud Hub. 

*This repo is currently in active development. We offer no guarantees that it is currently working properly, please check output data for bugs!*

At present, we calculate a set of annual ETCCDI indices for four models (UKESM1, E3SMv3, CESm2-WACCM and MIROC-ES2H), over the two scenarios G6-1.5K-SAI and G6-1.5K-HiLLA (and the SSP2-4.5 baseline). We are working to extend this, by:
- adding G6-1.5K-MCB
- adding extra metrics
- including monthly resolution indices
Contributions to the code-base are welcome! 

## Source
This repo is built off of work by Cindy Wang's Quantifying the effectiveness of multiple SAI strategies across different dimensions - ETCCDI indices. Please cite DOI 10.5281/zenodo.18880765

## Usage policy
These indices are processed G6-1.5K GeoMIP data, so any analyses using the data should follow [GeoMIP policy](https://climate.envsci.rutgers.edu/geomip/data.html) and invite modelling teams which produced the simulations to contribute as co-authors, if using within around the first year. We also ask that any users extend this offer to Cindy Wang, Alistair Duffey and John Orcutt as dataset producers. Contributions to the code 
