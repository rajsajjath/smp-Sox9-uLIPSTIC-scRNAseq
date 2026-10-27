These RScripts contains code to run the published scRNA-seq/CITE-seq analyses from Parigi et al., 
Nature Immunology (2026) on skin immune cells in contact with Sox9-expressing cells, as assessed
by the uLIPSTIC technology. 

"SMP_scRNAseq_object_processing.R" contains the analyses used for figures, and draws from the 
functions written in "SMP_scRNAseq_object_processing_functions.R".

To prepare the raw counts matrix used as input, raw fastq files were aligned to the mouse genome 
(mm10, v1.20 from 10X Genomics) and quantified using the ‘count’ function in CellRanger (v6.0.0) 
with default parameters (see 'Methods').
