# zoonotic-C
Calculations and figures to support Collins et al. review (in prep) of zoogeochemical carbon fluxes in the ocean. 

Notes on provenance of specific files in [data](data/):

* Files in [data/raw/sala_et_al_2021](data/raw/sala_et_al_2021) are from [Sala et al., 2021. Protecting the global ocean for biodiversity, food and climate. *Nature* **592**, 397–402](https://doi.org/10.1038/s41586-021-03371-z). Data files were retrieved from https://doi.org/10.25349/D9N89M on June 6, 2022.

* The files "fseq_OCIM2_48L.mat" and "fseq_OCIM2_48L.nc" in [data/raw/siegel_et_al_2021_v2](data/raw/siegel_et_al_2021_v2) are from [Siegel et al., 2021. Assessing the sequestration time scales of some ocean-based carbon dioxide reduction strategies. *Environ. Res. Lett.* **16** 104003](https://iopscience.iop.org/article/10.1088/1748-9326/ac0be0#erlac0be0s5). The files contain the same data in MATLAB and NetCDF formats and were retrieved from https://doi.org/10.6084/m9.figshare.15228690.v2 on June 6, 2022 **Update: Thsese files are not stored on GitHub due to size. User replicating this analysis will have to download the datasets from the figshare link and then move the files to the correct location.**

* The file "plot_sequestration_fraction.m" in [data/raw/siegel_et_al_2021_v2](data/raw/siegel_et_al_2021_v2) is the original, unmodified MATLAB script that Siegel et al. provided along with the data files containing the sequestration fractions. A modified version of this file, [gen_fracs_to_constrain_trawlCO2.m](gen_fracs_to_constrain_trawlCO2.m), was used to generate the necessary data for the analysis in [ConstrainCO2Flux.R](ConstrainCO2Flux.R)