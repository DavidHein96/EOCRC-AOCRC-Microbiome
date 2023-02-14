ASV Workflow
================
<david.hein@utsouthwestern.edu>
2023-02-14

- <a href="#1-load-packages" id="toc-1-load-packages">1 Load Packages</a>
- <a href="#2-initial-instructions" id="toc-2-initial-instructions">2
  Initial instructions</a>
- <a href="#3-set-up-inital-things" id="toc-3-set-up-inital-things">3 Set
  up inital things</a>
- <a href="#4-run-dada2" id="toc-4-run-dada2">4 Run dada2</a>
  - <a href="#41-run-these-chunks-to-prepare-data"
    id="toc-41-run-these-chunks-to-prepare-data">4.1 Run these chunks to
    prepare data</a>
  - <a href="#42-run-core-program" id="toc-42-run-core-program">4.2 Run core
    program</a>
- <a href="#5-making-the-files-for-analysis"
  id="toc-5-making-the-files-for-analysis">5 Making the files for
  analysis</a>
- <a href="#6-session-info" id="toc-6-session-info">6 Session info</a>

<br>

# 1 Load Packages

``` r
library(dada2)
library(dplyr)
library(plyr)
library(ggplot2)
library(stringr)
library(DECIPHER)
library(phangorn)
# If these packages are not present, uncomment and run the following code 

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("dada2")

#install.packages("dplyr")
#install.packages("plyr")
#install.packages("ggplot2")
```

<br>

# 2 Initial instructions

1.  Retrieve the directory path of the file containing the FASTQ files
2.  Create a sub directory in that directory called “filtered”
3.  Create another sub directory somewhere else where you want to run
    all of your analysis from, set this as the working directory for R.
    This R markdown file should be in that directory
4.  In the working directory make sure there is the
    silva_nr99_v138.1_train_set file and the silva_species_assignments
    file <br>

# 3 Set up inital things

``` r
# Set this as the path of where your FASTQ files are 
path<-"nov2022FASTQ"

# Set this to the group name of who the files belong to, and if you want you can throw a date in there
group_name <- "Sanford_Nov_2022"
```

<br>

# 4 Run dada2

## 4.1 Run these chunks to prepare data

``` r
list.files(path)

#Here make sure the FASTQs follow the typical naming convention of SampleID_randomstuff_R1_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

Check out quality plots by changing the number, pick a good cut off to
trim the end of the reads, the reverse read will probably need to be
trimmed more

``` r
plotQualityProfile(fnFs[130])
plotQualityProfile(fnRs[130])
```

Change the trunc length and run this chunk, look at the output, some
files will have so few reads that they may need to be removed

``` r
# Filtering files and placing in the directory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Here you set the trim at the truncLen parameter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)

# Learning the error rate for DADA2 algorithm
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

print(out)
```

<br>

## 4.2 Run core program

This chunk takes the longest

``` r
# Running the core algo
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merging fwd and rev reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Making table of reads and removing chimeras
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Assigning taxa with a pretrained classifier, up to genus level
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
# Adds species only with an exact match to the ASV 
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

ASVs <- DNAStringSet(colnames(seqtab.nochim))
names(ASVs)<- paste0("ASV",1:ncol(seqtab.nochim))
alignment <- AlignSeqs(ASVs,anchor=NA,processors = NULL)
```

<br>

# 5 Making the files for analysis

Produces a file showing counts and taxonomy for all ASVs

``` r
## 1. file with all ASVs ##
# transpose and set as df with the dna sequences as a variable
seqtab_nochim_t <- t(seqtab.nochim)
seqtab_nochim_t_df <- data.frame(seqtab_nochim_t)

# set as df with dna sequence variable
taxa_df <- data.frame(taxa)
df_final <- cbind(taxa_df,seqtab_nochim_t_df)
write.csv(df_final,file = paste0(group_name,"_full_taxonomy_with_ASV.csv"),row.names = FALSE,quote=FALSE)
```

<br><br>

# 6 Session info

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] phangorn_2.10.0     ape_5.6-2           DECIPHER_2.26.0    
    ##  [4] RSQLite_2.2.20      Biostrings_2.66.0   GenomeInfoDb_1.34.6
    ##  [7] XVector_0.38.0      IRanges_2.32.0      S4Vectors_0.36.1   
    ## [10] BiocGenerics_0.44.0 stringr_1.5.0       ggplot2_3.4.0      
    ## [13] plyr_1.8.8          dplyr_1.0.10        dada2_1.26.0       
    ## [16] Rcpp_1.0.9         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] MatrixGenerics_1.10.0       Biobase_2.58.0             
    ##  [3] bit64_4.0.5                 RcppParallel_5.1.5         
    ##  [5] assertthat_0.2.1            latticeExtra_0.6-30        
    ##  [7] blob_1.2.3                  GenomeInfoDbData_1.2.9     
    ##  [9] Rsamtools_2.14.0            yaml_2.3.6                 
    ## [11] pillar_1.8.1                lattice_0.20-45            
    ## [13] quadprog_1.5-8              glue_1.6.2                 
    ## [15] digest_0.6.31               GenomicRanges_1.50.2       
    ## [17] RColorBrewer_1.1-3          colorspace_2.0-3           
    ## [19] htmltools_0.5.4             Matrix_1.5-3               
    ## [21] pkgconfig_2.0.3             ShortRead_1.56.1           
    ## [23] zlibbioc_1.44.0             scales_1.2.1               
    ## [25] jpeg_0.1-10                 BiocParallel_1.32.5        
    ## [27] tibble_3.1.8                generics_0.1.3             
    ## [29] cachem_1.0.6                withr_2.5.0                
    ## [31] SummarizedExperiment_1.28.0 cli_3.4.1                  
    ## [33] magrittr_2.0.3              crayon_1.5.2               
    ## [35] deldir_1.0-6                memoise_2.0.1              
    ## [37] evaluate_0.19               fansi_1.0.3                
    ## [39] nlme_3.1-161                hwriter_1.3.2.1            
    ## [41] tools_4.2.1                 lifecycle_1.0.3            
    ## [43] matrixStats_0.63.0          interp_1.1-3               
    ## [45] munsell_0.5.0               DelayedArray_0.24.0        
    ## [47] compiler_4.2.1              rlang_1.0.6                
    ## [49] grid_4.2.1                  RCurl_1.98-1.9             
    ## [51] rstudioapi_0.14             igraph_1.3.5               
    ## [53] bitops_1.0-7                rmarkdown_2.19             
    ## [55] gtable_0.3.1                codetools_0.2-18           
    ## [57] DBI_1.1.3                   reshape2_1.4.4             
    ## [59] R6_2.5.1                    GenomicAlignments_1.34.0   
    ## [61] knitr_1.41                  fastmap_1.1.0              
    ## [63] bit_4.0.5                   utf8_1.2.2                 
    ## [65] fastmatch_1.1-3             stringi_1.7.8              
    ## [67] vctrs_0.5.1                 png_0.1-8                  
    ## [69] tidyselect_1.2.0            xfun_0.36
