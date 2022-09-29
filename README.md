# Adenocystis Biogeography
Scripts for Adenocystis Biogeography

The premise of this study was to examine the st

Please note: Some of the scripts are specific to the New Zealand e-Science Infrastructure (NeSI). This includes my username (adacl879) and project (uoo03666) and if you are using these scripts, they should be changed to whatever your account number / HPC is.

Please feel free to reach out to clare.im.adams@gmail.com if you are having troubles.

Please also note that this GBS data was not super great when it came back from the sequencing service, as we only got 12M reads instead of 80+M reads. Something went wrong with sequencing, hence why we only have ~500 SNPs. However, we believe that these are enough to paint broad brush strokes around what ecological patterns and processes may be going on. 

#Folder: Metadata
This folder contains the sampling location information. 

#Folder: Stacks
This folder contains scripts for Stacks analyses. First, you have a demultiplexing script (and the barcode file is included as barcodes_all.txt). Then you have a script for trimming off adapters, and then after a FastQC/MultiQC look, we decided to trim the sequences to 70 base pairs. This is because Stacks prefers to have sequences all of a uniform length. We then rename our trimmed files to something that Stacks will recognize. We tested the large M Stacks parameter on the BCS population, so there is a script for that, and then also tried out various values for the small m Stacks parameter. We ended up with choosing a Large M of 3, n of 3, and small m of 2 because these gave us the most SNPs without too many diminishing returns. Then, SNPs are filtered via VCFtools.

#Folder: R Analyses
This folder describes how to make maps in R, and includes some scratch paper/test scripts for data munging and testing out analyses in R. Probably the file you want is ScriptsForStatAnalyses.R, which produces some graphs and statistical analyses (including NJ trees and such). This is not complete, as a LEA analysis was not completed.

And then we have a bunch of random graphs in the main folder, perhaps of note would be the PCA plots. 

Mitochondrial data can be found in Jon Waters' High-Capacity Storage at the University of Otago HCS under sci-zoology-jonwaterslab/Clare-data/Clare-Postdoc-Work.

