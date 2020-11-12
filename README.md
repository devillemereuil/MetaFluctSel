Code and data for "Fluctuating optimum and temporally variable selection on breeding date in birds and mammals"
===============================================================================================================

Content
-------

- The data (distributed under CC BY-NC-SA, see `COPYING`) are available in the `Data` folder
- The R code is available in the `Analyses/Scripts R` folder
- The STAN models are available in the `Analyses/STAN models` folder

(Both the R and STAN code is distributed under GPL-3, see `COPYING` and `GPL`)

Data description
----------------

A description of the data is attached as an attribute to the R serialised data object `Data/alldata_final.rds` and contain the following information:

- **DataID** Internal ID for each dataset received for the meta-analysis, can correspond to several species and populations at once
- **Species** Scientific name of the species
- **Population**  Loose name for the location of the studied population
- **Year** Breeding year, rebased as integer for some datasets
- **ID** Breeding female ID, randomly anonymised for each dataset
- **Pheno** Phenological trait (parturition or lay date, depending on the dataset), provided as a date rebased on the current year
- **Fitness** Fitness measurement (number of offspring surviving to fledge or wean, or to the year after, depending on the dataset)

Contacts
--------

In the best interests of everyone, please contact the copyright holders (see COPYING) if you want to re-analyse the data for any other purpose than replication of the study, so as to avoid redundant, wasted work and to have a guided understanding of the biological sense behind the data.
