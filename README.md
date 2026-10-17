
# Data and Code Availability

This repository contains the datasets and analysis codes associated with our study on landscape interactions and ecosystem resilience in global drylands.

## Datasets

The datasets contain GeoTIFF raster layers describing landscape interaction strength, vegetation recovery-rate proxies, their temporal trends, and the local co-occurrence of high- and low-vegetation states.

- `Landscape_interaction_strength.tif`: Landscape interaction strength.
- `Landscape_interaction_strength_kendall_tau.tif`: Kendall’s tau for the temporal trend in landscape interaction strength.
- `Landscape_interaction_strength_kendall_p.tif`: P-value for the corresponding trend test.
- `RecoverRate_ac1.tif`: Recovery-rate proxy based on lag-1 autocorrelation (AC1).
- `RecoverRate_ac1_kendall_tau.tif`: Kendall’s tau for the temporal trend in the AC1-based proxy.
- `RecoverRate_ac1_kendall_p.tif`: P-value for the corresponding trend test.
- `RecoverRate_var.tif`: Recovery-rate proxy based on variance.
- `RecoverRate_var_kendall_tau.tif`: Kendall’s tau for the temporal trend in the variance-based proxy.
- `RecoverRate_var_kendall_p.tif`: P-value for the corresponding trend test.
- `Local_co-occurrence_states.tif`: Binary indicator of local state co-occurrence (0 = non-co-occurrence; 1 = co-occurrence).

## Analysis Codes

All analysis codes are provided in the `Code_upload` folder. All analyses were conducted using R version 4.1.3.

## Usage and Citation

Please consult the associated study for detailed definitions, temporal coverage, spatial resolution, and statistical procedures. Please cite the associated study when using these datasets or codes.
