# Overview

A combination of statistical and mechanistic models used to quantify creosotebush (Larrea tridentata) encroachment at the Sevilleta National Wildlife Refuge in central New Mexico. Here, we use dispersal kernels as well as size- and density-dependent demography to model encroachment as a moving wave pulled forward by individuals at the low-density vanguard. Data used here are collected annual demographic surveys over the course of several years, a transplant experiment to assess the survival rates of seedlings, and seed velocity data gathered from experiments in a laboratory setting.

# Files

## Data

**LT_Data** *(.xlsx)* - Spreadsheet containing all of the demographic data on the shrubs that were a part of the study. Also contains data on transect length and the transplant experiment. See the file metadata for more information.

**Seed_Drop** *(.xlsx)* - Spreadsheet containing position and time for each seed drop trial. See the file metadata for more information.

## Scripts

**00_RunAll** *(.R)* - Runs scripts 01-07. This process takes a few minutes.

**01_SeedVelocities** *(.R)* - Finds the distribution of seed terminal velocities using data from **Seed_Drop**.

**02_WindSpeeds** *(.R)* - Finds the distribution of wind speeds using wind speed data from the Sevilleta NWR; this data can be found [here](http://sevlter.unm.edu/content/meteorology-data-sevilleta-national-wildlife-refuge-new-mexico-1988-present). The 2011-2015 data was not used because it did not play nice, likely due to an apparent mismatch between the header row and the rest of the data.

**03_Dispersal** *(.R)* - Creates dispersal kernels using the methods from Skarpaas and Shea (2007).

**04_CDataPrep** *(.R)* - Cleans up the demographic data from **LT_Data** and calculates variables that will be used in size- and density-dependent demographic analyses.

**05_BootRes.R** *(.R)* - Resamples data for each bootstrap replicate.

**05_CDataAnalysis_BS** *(.R)* - Generates models of demographic rates under better survival (BS) conditions; these conditions were observed in a census that occurred after higher than average rainfall.

**05_CDataAnalysis_NS** *(.R)* - Generates models of demographic rates under normal survival (NS) conditions.

**06_SIPM** *(.R)* - Constructs the transition matrix that projects population growth, and combines this with dispersal data to find the speed at which the encroaching shrub wave travels.

**07_MainFigures** *(.R)* - Generates figures on size- and density-dependence, wave speeds and population growth, and dispersal kernels.

## Supporting Information

**S1_SupportingMaterial** *(.R)* - Code not directly used in the analyses but contributes to our understanding of the data.

**S2_DeprecatedCode.R** *(.R)* - Bits of old code that is no longer in use; they are not organised in any particular order.

**S3_TransplantAnalysis.R** *(.R)* - Code to analyse data from the transplant experiment.
