# ===============
# scRNA_workflow
# ===============
This R script is written to provide easy to flow template using SEURAT package.
The script integrates all samples, preprocess the samples, split the samples
into seurate object lists, and implements DoubletFinder package to estimate pK
value for each sample separately. Here, I emphasize that parameter estimate
for DoubletFinder should be implemented separately, unfortunately, the output
of DoubletFinder does not allow us to combine all outputs and push them into
seurat object. This is a bit glitchy because the DoubletFinder adds a separate
column heading for each sample. I added a block of code to implement the 
DoubletFinder in iteration and override the seurat output format so that we can
combine them as a single seurate object after DoubletFinder implementation. 
The script implements SingleR and annotates cell type using several databases
like DICE and MONACO. 

The cell type annotation is the most challenging but extremely important part of
scRNA workflow. I am examining performance of celltyping using different database
and will document the results in this repo. 


