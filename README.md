# TacCA_sourcedata
- Source data for TarCA (`Fig1-5`)

Each folder represents the source data needed to draw the corresponding main figure. 
Each file in the folder contains the source data needed to draw the corresponding panel. 
For example, folder `Fig2` contains the source data needed to draw main figure 2, and file `Fig2/2bc.RData` contains the source data needed to draw main figure 2 panel **b** and **c**.  

Use command `readRDS("fileName")` in R to load RData files. For example, use command `readRDS("Fig2/2bc.RData")` to load `Fig2/2bc.RData`.

- QFM

`QFM` folder contains scripts used for processing QFM data.

- CoSpar

`CoSpar` folder contains scripts used for processing CoSpar data.
