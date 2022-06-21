# zoonotic-C
Calculations and figures to support Collins et al. review (in prep) of zoogeochemical carbon fluxes in the ocean. 

Notes on provenance of specific files in [data](data/):

* Files in [data/raw/sala_et_al_2021](data/raw/sala_et_al_2021) are from [Sala et al., 2021. Protecting the global ocean for biodiversity, food and climate. *Nature* **592**, 397–402](https://doi.org/10.1038/s41586-021-03371-z). Data files were retrieved from https://doi.org/10.25349/D9N89M on June 6, 2022. Resolution of the CO2 flux dataset is 1 km<sup>2</sup>.

* The files "fseq_OCIM2_48L.mat" and "fseq_OCIM2_48L.nc" in [data/raw/siegel_et_al_2021_v2](data/raw/siegel_et_al_2021_v2) are from [Siegel et al., 2021. Assessing the sequestration time scales of some ocean-based carbon dioxide reduction strategies. *Environ. Res. Lett.* **16** 104003](https://iopscience.iop.org/article/10.1088/1748-9326/ac0be0#erlac0be0s5). The files contain the same data in MATLAB and NetCDF formats and were retrieved from https://doi.org/10.6084/m9.figshare.15228690.v2 on June 6, 2022. Per the published paper, model resolution is 2 degrees with 48 vertical levels. **Update: These files are not stored on GitHub due to size. User replicating this analysis will have to download the datasets from the figshare link and then move the files to the correct location.**

* The file "plot_sequestration_fraction.m" in [data/raw/siegel_et_al_2021_v2](data/raw/siegel_et_al_2021_v2) is the original, unmodified MATLAB script that Siegel et al. provided along with the data files containing the sequestration fractions. A modified version of this file, [gen_fracs_to_constrain_trawlCO2.m](gen_fracs_to_constrain_trawlCO2.m), was used to pull out the [necessary benthic sequestration fractions](data/derived/benthic_seqfractions) for the analysis in [ConstrainCO2Flux.R](ConstrainCO2Flux.R).

* The files in [data/derived/benthic_seqfractions](data/derived/benthic_seqfractions) are 25, 50 and 100 year sequestration fractions for bottom depths, plus necessary metadata, pulled from the [original Siegel et al. model output](data/raw/siegel_et_al_2021_v2) using the script [gen_fracs_to_constrain_trawlCO2.m](gen_fracs_to_constrain_trawlCO2.m).
