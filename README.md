# Overview

This project seeks to create a mathematical model of creosotebush (Larrea tridentata) encroachment at the Sevilleta National Wildlife Refuge in central New Mexico. Here, we use dispersal kernels and density dependent demography to model encroachment as a moving wave, and then combine the two to generate this model.

# Files

**LT_Data** *(.xlsx)* - Spreadsheet containing all of the demographic data on the shrubs that were a part of the study. Also contains data on transect length and the transplant experiment. See the file metadata for more information.

**Seed_Drop** *(.xlsx)* - Spreadsheet containing position and time for each seed drop trial. See the file metadata for more information.

**00_RunAll** *(.R)* - Runs scripts 01-07. This process takes a few minutes.

**01_SeedVelocities** *(.R)* - Finds the distribution of seed terminal velocities using data from **Seed_Drop**.

**02_WindSpeeds** *(.R)* - Finds the distribution of wind speeds using wind speed data from the Sevilleta NWR; this data can be found [here](http://sevlter.unm.edu/content/meteorology-data-sevilleta-national-wildlife-refuge-new-mexico-1988-present). The 2011-2015 data was not used because it did not play nice, likely due to an apparent mismatch between the header row and the rest of the data.

**03_Dispersal** *(.R)* - Creates dispersal kernels using the methods from Skarpaas and Shea (2007).

**04_CDataPrep** *(.R)* - Cleans up the demographic data from **LT_Data** and calculates variables that will be used in size- and density-dependent demographic analyses.

**05_CDataAnalysis_BS** *(.R)* - Generates models of demographic rates under better survival (BS) conditions; these conditions were observed in a census that occurred after higher than average rainfall.

**05_CDataAnalysis_NS** *(.R)* - Generates models of demographic rates under normal survival (NS) conditions.

**06_SIPM** *(.R)* - Constructs the transition matrix that projects population growth, and combines this with dispersal data to find the speed at which the encroaching shrub wave travels.

**07_MainFigures** *(.R)* - Generates figures on size- and density-dependence, wave speeds and population growth, and dispersal kernels.

**S1_SupportingMaterial** *(.R)* - Code not directly used in the analyses but contributes to our understanding of the data.
