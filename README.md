# Code-for-Single-Cell-Analysis
Analysis code used in Zhang et al,Single cell transcriptome profiling of immune tissues from germ-free and specific pathogen-free piglet.2025
The code to identify the cell composition and functional characteristics of different immune tissues in pigs.
This workflow is used to identify the cell composition and functional characteristics of different immune tissues in pigs. The raw data can be accessed through the  Genome Sequence Archive (GSA) with accession number CRA018161.
The cell annotation of this project are provided  `Cells_annotation.zip` ,  and can be quickly integrated into new Seurat object using the `AddMetaData` function from the Seurat package.

##  Part1: Single-cell data preprocessing
1. Refer to the 10X CellRange workflow. First, extract the protein-coding gene information from the GTF file, and then construct the CellRange reference file.`01.get_PCG_ref.sh`
2.Quantification of single-cell data.`02.cellranger.count.sh`

##  Part1: Cell type identification
1. According to Seurat's workflow, we first filter the gene expression matrix from CellRanger, retaining valid cells and genes.
2. Cells are categorized through dimensionality reduction and clustering, and the identity of each cluster is characterized based on known gene markers.`03.Merge_and_remove_batch_effects.R`
3. Cell differentiation potential analysis is performed using the CytoTRACE package, running with the default parameters.`04.CytoTRACE.R`
