# Overview

The encroachment of woody plants into grasslands is a global phenomenon with implications for biodiversity and ecosystem function. Understanding and predicting the pace of expansion and the underlying processes that control it are key challenges in the study and management of woody encroachment. Theory from spatial population biology predicts that the occurrence and speed of expansion should depend sensitively on the nature of conspecific density dependence. If fitness is maximized at the low-density encroachment edge, then shrub expansion should be “pulled” forward. However, encroaching shrubs have been shown to exhibit positive feedbacks, whereby shrub establishment modifies the environment in ways that facilitate further shrub recruitment and survival. In this case there may be a fitness cost to shrubs at low density causing expansion to be “pushed” from behind the leading edge. We studied the spatial dynamics of creosotebush (*Larrea tridentata*), which has a history of encroachment into Chihuahuan Desert grasslands over the past century. We used demographic data from observational censuses and seedling transplant experiments to test the strength and direction of density dependence in shrub fitness along a gradient of shrub density at the grass–shrub ecotone. We also used seed-drop experiments and wind data to construct a mechanistic seed-dispersal kernel, then connected demography and dispersal data within a spatial integral projection model (SIPM) to predict the dynamics of shrub expansion. Contrary to expectations based on potential for positive feedbacks, the shrub encroachment wave is “pulled” by maximum fitness at the low-density front. However, the predicted pace of expansion was strikingly slow (ca. 8 cm/year), and this prediction was supported by independent resurveys of the ecotone showing little to no change in the spatial extent of shrub cover over 12 years. Encroachment speed was acutely sensitive to seedling recruitment, suggesting that this population may be primed for pulses of expansion under conditions that are favorable for recruitment. Our integration of observations, experiments, and modeling reveals not only that this ecotone is effectively stalled under current conditions but also why that is so and how that may change as the environment changes.

*The corresponding publication for this repository can be found [here](https://doi.org/10.1002/ecm.1574). Note that the raw manuscripts and appendices in this repository may differ slightly from the published version.*

<br/>

# Files

*Note: Due to the large number of files in the repository, not all files are described below.*

## Data

**LT_DemographyData** *(.csv)* - Spreadsheet containing 2013-2017 census data.

**LT_TransectData** *(.csv)* - Spreadsheet containing shrub size and location data at the time of setup.

**LT_TransectLengths** *(.csv)* - Spreadsheet containing transect lengths.

**LT_TransectResurvey** *(.csv)* - Spreadsheet containing shrub cover data from the 2001 and 2013 resurveys.

**LT_TransplantExp** *(.csv)* - Spreadsheet containing the transplant experiment data.

**LT_XLSX** *(.xlsx)* - Spreadsheet containing all data in the "LT" data files, except for the transect resurvey.

**SD_Summary** *(.csv)* - Spreadsheet containing seed drop summary data.

**SD_Trials** *(.csv)* - Spreadsheet containing position versus time data for the seed drop experiment.

**SD_XLSX** *(.xlsx)* - Spreadsheet containing all data in the "SD" data files.

**Metadata** *(.txt)* - Plain text file describing fields for the various datasheets.

**Derived** *(folder)* - Folder containing data files derived from simulations and model runs, kept here so simulation outputs don't have to be generated again every time the scripts are run.

**EDI** *(folder)* - Folder containing data files prepared as submissions to the Environmental Data Initiative (EDI).

**Weather** *(folder)* - Folder containing SEV weather data from one of six time periods: 1988-1994, 1995-1999, 2000-2004, 2005-2009, 2010-2014, and 2015-2019. The data files can be found [here](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.1.14), and relevant metadata can be found [here](https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sev.1.14). Also contains data on monsoon precipitation.

## Scripts

**00_RunAll** *(.R)* - Runs scripts 01-08; this is computationally intensive and may take a while.

**01_SeedVelocities** *(.R)* - Finds the distribution of seed terminal velocities using seed drop data.

**02_WindSpeeds** *(.R)* - Finds the distribution of wind speeds using wind speed data from the SEV.

**03_Dispersal** *(.R)* - Sets up dispersal kernel functions.

**04_CDataPrep** *(.R)* - Cleans up the demographic data and calculates variables that will be used in size- and density-dependent demographic analyses.

**05_CDataAnalysis** *(.R)* - Generates models of demographic rates as a function of size and density.

**06_BootRes** *(.R)* - Sets up resampling for bootstraps.

**07_SIPM** *(.R)* - Constructs the transition matrix that projects population growth, and combines this with dispersal data to find the speed at which the encroaching shrub wave travels.

**08_MainFigures** *(.R)* - Generates figures on size- and density-dependence, wave speeds and population growth, and dispersal kernels.

**Extras** *(folder)* - Various extra scripts not launched from the master file, but are still used in other functions such as figure generation.

## Others

**Manuscript** *(folder)* - A folder with code and various other objects used to generate the manuscript, including figures.
