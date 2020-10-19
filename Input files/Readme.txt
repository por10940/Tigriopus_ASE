This folder contains input files for the model

9679UniqueReads.csv = filtered raw unique mapped reads 
  each Column name = Cross_Treatment_Replicate_Read
    Cross: Parents = SD and SCN (SC in the manuscript) F1 hybrid = SCfxSDm and SDfxSCm
    Treatment: 20 = control, 35 = heat stress
    Replicate: 2015 = replicate#1, 2016 = replicate#2 (each biological replicate was made independently on that year)
    Read: number of reads mapped to SD or SC allele
    
9679ContigName.csv = file contain matching reference contig names and gene names from blast
  The gene numbers (first column) in the input file are orders of the contig names in the 9679ContigNames.csv file
  
Model Input folder contains files formatted for running ASE_Bayesian_model.R
  Model requires 1 parent file and 1 hybrid file input files of the same treatment to run. Each temperature treatment was run independently.
  parent20_9679.csv = parental reads from control treatment
  parent35_9679.csv = parental reads from heat stress streatment
  hybrid20_9679.csv = hybrid reads from control treatment
  hybrid35_9679.csv = hybrid reads from heat stress treatment

Parent file (parent*.csv) column names
  Gene: contig order in 9769ContigNames.csv
  Rep: replicate # (1 or 2)
  Temp.hot: control treatment = 0, heat stress treatment = 1
  N.Correct: # of parental reads mapped to the correct reference (#SD parental reads that mapped to SD reference or #SC parental reads that mapped to SC reference)
  N.total: # of parental reads mapped to SD refercne + # of parental reads mapped to SC reference
  N.is.SD: SD population 
  N.is.SC:
  
Hybrid File (hybrid*.csv) column names
  Gene: contig order in 9769ContigNames.csv
  Rep: replicate # (1 or 2)
  Temp.hot: control treatment = 0, heat stress treatment = 1
  Mom.is.SD: SCfxSDm = 0, SDfxSCm = 1 (direction of reciprocal cross)
  N.off.SD: # of F1 hybrid reads mapped to SD reference
  N.off.total: # of F1 hybrid reads mapped to SD reference + # of F1 hybrid reads mapped to SC reference
