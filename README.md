Hello!
Welcome to the Oxyria_digyna_phenology repository! If you use any of these scripts please cite Stunz et al. 2025. 
All code run in Stunz et al. 2025 can be found within the Temperature_2025.R and Oxyria_analysis_2025.R files. 
Briefly, these scripts detail how temperature and elevation data were obtained and how linear regression models and spatial autocorrelation models were built for multiple O. digyna groups and variables. All files except for CRU.7z* can be found in the 'main' branch of this repository. To find CRU.7z (large file with all temperature data) see the 'master' branch.

#---------------------FILE SPECIFICATIONS----------------------#

Oxyria_data.csv contains geographical and temporal data of
preserved Oxyria digyna GBIF occurrences, as well as whether
flowers were present or not.

O_digyna_all_occurrences.csv contains geographical and temporal
data of all available O. digyna GBIF occurrences.

The temperature script (Temperature_2025.R) unpacks the
temperature data (CRU.nc)* and creates .csv files, which is
required in order to run sections of the Oxyria analysis script
(Oxyria_analysis_2025.R). This script is adapted from:
Bartlein, P. (2024). netCDF in R. Github.
https://pjbartlein.github.io/REarthSysSci/netCDF.html

The Oxyria analysis script loads, reformats, and performs
analyses on all data available in the repository. More info on
the process is available in the script. The code for Spatial 
Autoregressive Models under "SAR" is adapted from:
Faurby, S. (2022). Spatial analyses on a round(ish) planet.
https://github.com/sorenfaurby/Spatial-R-PhD-course/blob/a3af750e6ddd568a6931f35f65117ae23a9fa4db/Tutorials/tutorial_4.pdf

The ne_110m_land files are world map projections required for
the SAR-models, and are all loaded through the shapefile
(ne_110m_land.shp)

ChatGPT was used as a tool in development of R code.

*CRU.7z needs to be unzipped before running any lines in
"Temperature_2025.R".
