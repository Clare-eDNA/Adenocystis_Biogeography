# Adenocystis Biogeography
Scripts for Adenocystis Biogeography
Please note: Some of the scripts are specific to the New Zealand e-Science Infrastructure (NeSI). This includes my username (adacl879) and project (uoo03666) and if you are using these scripts, they should be changed to whatever your account number / HPC is. 

Please feel free to reach out to clare.im.adams@gmail.com if you are having troubles.

Please also note that 

#Folder: Metadata
This folder contains the sampling location information. 

#Folder: Stacks
This folder contains scripts for Stacks analyses. First, you have a demultiplexing script (and the barcode file is included as barcodes_all.txt). Then you have a script for trimming off adapters, and then after a FastQC/MultiQC look, we decided to trim the sequences to 70 base pairs. This is because Stacks prefers to have sequences all of a uniform length. We then rename our trimmed files to something that Stacks will recognize. We tested the large M Stacks parameter on the BCS population, so there is a script for that, and then also tried out various values for the small m Stacks parameter. We ended up with choosing a Large M of 3, n of 3, and small m of 2 because these gave us the most SNPs without too many diminishing returns. 

#Folder: R Analyses

