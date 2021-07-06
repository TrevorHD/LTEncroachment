# Overview

A combination of statistical and mechanistic models used to quantify creosotebush (Larrea tridentata) encroachment at the Sevilleta National Wildlife Refuge in central New Mexico. Here, we use dispersal kernels as well as size- and density-dependent demography to model encroachment as a moving wave pulled forward by individuals at the low-density vanguard. Data used here are collected annual demographic surveys over the course of several years, a transplant experiment to assess the survival rates of seedlings, and seed velocity data gathered from experiments in a laboratory setting.

<br/>

# Files

## Data

**LT_XLSX** *(.xlsx)* - Spreadsheet containing all of the demographic data on the shrubs that were a part of the study. Also contains data on transect length and the transplant experiment.

**SD_XLSX** *(.xlsx)* - Spreadsheet containing position and time for each seed drop trial.

**Metadata** *(.txt)* - Plain text file describing fields for the various datasheets.

**WeatherX** *(.csv)* - SEV weather data from one of six time periods: 1988-1994, 1995-1999, 2000-2004, 2005-2009, 2010-2014, and 2015-2019. The data files can be found [here](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.1.14), and relevant metadata can be found [here](https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sev.1.14).

**Creosote_transect_resurvey** *(.csv)* - Spreadsheet containing percent shrub cover from the transect resurveys; this file name will likely be changed later.

*Note: The data folder also contains CSV files with prefixes LT or SD; these extra files just CSV replicas of individual sheets in the XLSX files, and exist to increase accessibility to those who may not have the software to read XLSX files.*

## Scripts

**00_RunAll** *(.R)* - Runs scripts 01-07. This process takes a few minutes.

**01_SeedVelocities** *(.R)* - Finds the distribution of seed terminal velocities using data from **Seed_Drop**.

**02_WindSpeeds** *(.R)* - Finds the distribution of wind speeds using wind speed data from the SEV.

**03_Dispersal** *(.R)* - Creates dispersal kernels using the methods from Skarpaas and Shea (2007).

**04_CDataPrep** *(.R)* - Cleans up the demographic data from **LT_Data** and calculates variables that will be used in size- and density-dependent demographic analyses.

**05_CDataAnalysis_NS** *(.R)* - Generates models of demographic rates as a function of size and density.

**06_BootRes.R** *(.R)* - Resamples data for each bootstrap replicate.

**07_SIPM** *(.R)* - Constructs the transition matrix that projects population growth, and combines this with dispersal data to find the speed at which the encroaching shrub wave travels.

**08_MainFigures** *(.R)* - Generates figures on size- and density-dependence, wave speeds and population growth, and dispersal kernels.

## Others

**Manuscript** *(folder)* - A folder with code and various other objects used to generate the manuscript.

**05B_nonnormal_growth** *(.R)* - Models ofnon-normal growth rates as a function of size and density; not yet fully implemented.

**wavespeed_sensitivities** *(.R)* - Estimates wavespeeds for a range of per-seed recruitment rates; not yet fully implemented.

## Extras

**S1_SupportingMaterial** *(.R)* - Code not directly used in the analyses but contributes to our understanding of the data.

**S2_DeprecatedCode.R** *(.R)* - Bits of old code that is no longer in use; they are not organised in any particular order.

**S3_TransplantAnalysis.R** *(.R)* - Code to analyse data from the transplant experiment.
