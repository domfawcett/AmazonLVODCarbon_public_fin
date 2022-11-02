# AmazonLVODCarbon_public_fin

Author contact: d.fawcett@exeter.ac.uk

This repository includes the R codes used for generating the main results of the manuscript "Declining Amazon biomass due to deforestation and subsequent degradation losses exceeding gains". For required datasets please refer to the information below.
The graphics presented in the main manuscript were generated using:

**GenerateFig1.R**

**GenerateFig2and3.R**

**GenerateFig4.R**

## Datasets

Processed datasets for generation of Figures 1-4 are available at Zenodo: https://doi.org/10.5281/zenodo.7245450

Original input data is available from their respective sources: 

Mapbiomas Amazonia C2 (MapBiomas, 2021): <https://amazonia.mapbiomas.org>

Secondary forest maps (Silva Junior, Heinrich, et al., 2020):
<https://github.com/celsohlsj/gee_brazil_sv>

TMF dataset (Vancutsem et al., 2021): <https://forobs.jrc.ec.europa.eu/TMF/download/>

ESA CCI Biomass (Santoro & Cartus, 2021): 
<http://dx.doi.org/10.5285/84403d09cef3485883158f4df2989b0c> 

SMOS-IC V2 L-VOD daily data was provided by Jean-Pierre Wigneron (jean-pierre.wigneron@inrae.fr) and processed as described in the manuscript and SI.

## GEE Scripts repository

These GEE Scripts pre-process and export land cover data for further analysis. Access to private GEE assets used must be requested from their authors as indicated.

<https://code.earthengine.google.com/?accept_repo=users/dfawcett/RECCAP_pub>

## R Scripts

These R Scripts generate the results presented in Fawcett et al. 

### analyseLVOD_ccgfilt_diff_Amazon_bycountry_excl2011_V18_fin.R
This code generates Amazon country specific changes in AGC attributed to different processes 
Excluding 2011 which was a la Nina year and saw a large increase in AGC

Outputs:

- statistics of gross and net biomass changes and trends for each country and the entire Amazon (SI table S4), Figure S11


### analyseLVOD_ccgfilt_diff_Amazon_Main_V18_fin.R
This code 1. generates modelled values of AGC stocks and change and 2. compares them to annual L-VOD AGC

Outputs:

- Figure 2 a-e and 3 a and b, Supplementary Information Figure S5, S8, S15, S22, AGC change datasets for use by other codes

### analyseLVOD_ccgfilt_diff_Amazon_stats_V17_fin.R
This code calculates correlations between the L-VOD AGC and the modelled AGC change including different processes (deforestation, degradation, SF growth, IF change)

Outputs:

- Correlation statistics and errors for <90% old growth forest fraction grid cells, statistics reported in Suppl. Inf. Table S6


### getStaticIntactForestReferenceV16_fin.R
This code uses the ESA CCI biomass map of 2017 and precomputed values of intact forest cover to generate a reference biomass map of intact forest AGC for each 1 km cell used in further processing

Outputs:

- intact forest reference AGC maps, plus and minus SD of the original maps

### GenerateFig1.R

Outputs:

- Fig. 1 b and data for Fig. 1 a finalised in ArcMap, plus supplementary figures

### GenerateFig2and3.R

Outputs:

- Panels for Fig. 2 and 3, supplementary figures

### GenerateFig4.R

Outputs:

- Fig. 4 and statistics, supplementary figures


