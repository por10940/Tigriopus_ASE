This folder contains input files for the model

9679UniqueReads.csv = filtered raw unique mapped reads 
  each Column name = Cross_Treatment_Replicate_Read
    Cross: Parents = SD and SCN (SC in the manuscript) F1 hybrid = SCfxSDm and SDfxSCm
    Treatment: 20 = control, 35 = heat stress
    Replicate: 2015 = replicate#1, 2016 = replicate#2 (each biological replicate was made independently on that year)
    Read: number of reads mapped to SD or SC allele
    
9679ContigName.csv = file contain matching reference contig names and gene names from blast
  The gene number in the input file is corresponding with the contig names in the 9679ContigNames.csv file
  
Model Input folder contains files formatted for running ASE_Bayesian_model.R
