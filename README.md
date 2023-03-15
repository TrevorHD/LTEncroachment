# Overview

A combination of statistical and mechanistic models used to quantify creosotebush (*Larrea tridentata*) encroachment at the Sevilleta National Wildlife Refuge in central New Mexico. Here, we use dispersal kernels as well as size- and density-dependent demography to model encroachment as a moving wave pulled forward by individuals at the low-density vanguard. Data used here are collected via annual demographic surveys over the course of several years, a transplant experiment to assess the survival rates of seedlings, and seed laboratory experiments to determine dispersal-related traits of seeds.

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
