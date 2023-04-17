Class\_exercise
================
Corrin Bowman
3/20/2023

# Load the libraries you need

# Load functions you need “my\_class\_functions”

# load in your peak files for each replicate of each protein

# Here I am starting to analyze my data for my proteins of interest:

# proteinX, Y, Z …..

# First I will read in each replicate file

``` r
# filepath to import peaks
basepath <- "/scratch/Shares/rinnclass/CLASS_2023"
peak_path <- "2023_2_10_group3/group/results/bwa/mergedLibrary/macs/broadPeak/"
broadpeakfilepath <- file.path(basepath, peak_path)

# import peaks
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

# let's get a list of how many peaks are in each file before we create consensus peaks.
peak_num <- sapply(peak_list, length) %>% as.data.frame(row.names = T)
```

    ## Warning in as.data.frame.integer(., row.names = T): 'row.names' is not a
    ## character vector of length 22 -- omitting it. Will be an error!

``` r
# label column
names(peak_num) <- c("num_peaks")


# printing out a table of the number of peaks in each file:
peak_num <- peak_num %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")
```

# Now I am going to create consensus peaks for each protein

``` r
dbps <- unique(sapply(names(peak_list), function(x) {
   unlist(strsplit(x, "_"))[1]
}))

consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)
names(consensus_list) <- dbps

num_consensus_peaks <- sapply(consensus_list, length) %>% 
  as.data.frame() %>%
  rownames_to_column( var = "dbp") %>%
  dplyr::rename(number_consensus_peaks = ".")

peak_num <- left_join(peak_num, num_consensus_peaks)
```

    ## Joining with `by = join_by(dbp)`

``` r
# export consensus peaks to results folder

#write_csv(peak_num, "results/num_peaks_df.csv")
```

# Now I am going to make my consensus peaks compatable with UCSC genome browser

``` r
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman"
consensus_path <- "CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/consensus_peaks/"
exportpath <- file.path(basepath, consensus_path)

# now we can export each as .bed file
for(i in 1:length(consensus_list)) {
rtracklayer::export(consensus_list[[i]], paste0(exportpath, names(consensus_list)[i], "_consensus_peaks.bed") )}




consensus_file_list <- list.files("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/consensus_peaks/", full.names = T, pattern = ".bed")

peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))

names(peaks) <- dbps

# make chromosomes of interest object
canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")

# let's use lapply with filter funciton to cannonical_chr
peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))


# now that these are all nice and clean let's export:
new_filenames <- paste0("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/consensus_peaks/", names(peaks), "_consensus.bed")

for(i in 1:length(peaks)) {
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}

# we are using paste0 to print the header text
# then we are adding the names as the value the header equals with 'names' function
headers <- paste0("track type=bed name=", names(peaks))
headers
```

    ##  [1] "track type=bed name=CEBPZ" "track type=bed name=CHD2" 
    ##  [3] "track type=bed name=CTCF"  "track type=bed name=ELF1" 
    ##  [5] "track type=bed name=EP300" "track type=bed name=NA"   
    ##  [7] "track type=bed name=NA"    "track type=bed name=NA"   
    ##  [9] "track type=bed name=NA"    "track type=bed name=NA"   
    ## [11] "track type=bed name=NA"    "track type=bed name=NA"   
    ## [13] "track type=bed name=NA"    "track type=bed name=NA"   
    ## [15] "track type=bed name=NA"    "track type=bed name=NA"

``` r
# print out consensus peak files in a results/UCSC directory
# creating a path to export after we add header in for loop below
new_filenames <- paste0("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/", names(peaks), ".bed")
new_filenames
```

    ##  [1] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/CEBPZ.bed"
    ##  [2] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/CHD2.bed" 
    ##  [3] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/CTCF.bed" 
    ##  [4] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/ELF1.bed" 
    ##  [5] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/EP300.bed"
    ##  [6] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ##  [7] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ##  [8] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ##  [9] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ## [10] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ## [11] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ## [12] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ## [13] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ## [14] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ## [15] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"   
    ## [16] "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/ucsc_consensus_peaks/NA.bed"

``` r
# let's do so more indexing and only print out the first two of the three files.

for(i in 1:length(peaks)) {
  # Write the header line
  writeLines(headers[[i]], new_filenames[[i]])
  # Append the broadPeak table data
  
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}




# print out consensus peak files in a results/UCSC directory
```

# Now I want to compare a protein with a previous analysis

``` r
# goto UCSC genome browser and load in a peak file for a given protein
# load in the data for the same protein from the previous analysis
# compare how your consensus peaks are similar or different to previous analyses


#I uploaded EP300 and compared this to the EP300 from spring of 2021. I believe what I was looking at was had mostly similar peaks with minor differences. 
```

# Now I am going to determine how my peaks for each protein overlap annotations of the genome

# First I will find the overlaps between my consensus peaks with promoters of lncRNA and mRNA promoters

``` r
library(tidyverse)
library(GenomicRanges)
source("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/util/BCHM5631_my_class_functions.R")

basepath <- "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman"
peak_path <- "CLASS_2023/CLASSES/05_R_analyses//Class_exercise/Results/consensus_peaks"
consensusPeakPath <- file.path(basepath, peak_path)


 consensus_peaks_files<- list.files(consensusPeakPath, 
                                             pattern = "*.bed",
                                             full.names = TRUE)

consensus_peaks <- lapply(consensus_peaks_files, rtracklayer::import)

names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/consensus_peaks/|_consensus.bed","", consensus_peaks_files)

gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
table(gencode_gr$type)
```

    ## 
    ##           gene     transcript           exon            CDS    start_codon 
    ##          60609         227462        1372308         761508          87662 
    ##     stop_codon            UTR Selenocysteine 
    ##          79913         310193            119

``` r
rtracklayer::export(gencode_genes, "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/gene_annotations/gencode_genes.gtf")


mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"] 

rtracklayer::export(mrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/gene_annotations/mrna_genes.gtf")
table(gencode_genes$gene_type)
```

    ## 
    ##                          IG_C_gene                    IG_C_pseudogene 
    ##                                 14                                  9 
    ##                          IG_D_gene                          IG_J_gene 
    ##                                 37                                 18 
    ##                    IG_J_pseudogene                      IG_pseudogene 
    ##                                  3                                  1 
    ##                          IG_V_gene                    IG_V_pseudogene 
    ##                                144                                188 
    ##                             lncRNA                              miRNA 
    ##                              16849                               1881 
    ##                           misc_RNA                            Mt_rRNA 
    ##                               2212                                  2 
    ##                            Mt_tRNA             polymorphic_pseudogene 
    ##                                 22                                 42 
    ##               processed_pseudogene                     protein_coding 
    ##                              10171                              19965 
    ##                         pseudogene                           ribozyme 
    ##                                 18                                  8 
    ##                               rRNA                    rRNA_pseudogene 
    ##                                 52                                500 
    ##                             scaRNA                              scRNA 
    ##                                 49                                  1 
    ##                             snoRNA                              snRNA 
    ##                                942                               1901 
    ##                               sRNA                                TEC 
    ##                                  5                               1061 
    ##                          TR_C_gene                          TR_D_gene 
    ##                                  6                                  4 
    ##                          TR_J_gene                    TR_J_pseudogene 
    ##                                 79                                  4 
    ##                          TR_V_gene                    TR_V_pseudogene 
    ##                                106                                 33 
    ##   transcribed_processed_pseudogene     transcribed_unitary_pseudogene 
    ##                                495                                130 
    ## transcribed_unprocessed_pseudogene    translated_processed_pseudogene 
    ##                                923                                  2 
    ##  translated_unprocessed_pseudogene                 unitary_pseudogene 
    ##                                  2                                 98 
    ##             unprocessed_pseudogene                           vaultRNA 
    ##                               2631                                  1

``` r
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 

rtracklayer::export(lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/gene_annotations/lncrna_genes.gtf")


mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
rtracklayer::export(mrna_lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/gene_annotations/mrna_lncrna_genes.gtf")


lncrna_mrna_genes <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/gene_annotations/mrna_lncrna_genes.gtf")

# Nice that loaded so much faster -- lets see as DF
lncrna_mrna_genes_df <- lncrna_mrna_genes %>% as.data.frame()


```

``` r
rtracklayer::export(lncrna_mrna_promoters, "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/gene_annotations/lncrna_mrna_promoters.gtf")


# last handy annotation will be lncRNA and mRNA gene IDs to subset
lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
table(mrna_lncrna_genes$gene_type)
```

    ## 
    ##         lncRNA protein_coding 
    ##          16849          19965

``` r
# same for mRNAs
mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]

num_peaks_df <- data.frame("dbp" = names(consensus_peaks),
                           "num_peaks" = sapply(consensus_peaks, length))


# Now let's get the total amount of the genome covered by all the peaks for a given DBP.

num_peaks_df$total_peak_length <- sapply(consensus_peaks, function(x) sum(width(x)))

source ("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/util/my_class_functions.R")
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, type = "counts")

num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)
```

## results:

\#1) What can you determine from these overlaps?

Here we looked at the overlaps between my consensus peaks with promoters
of lncRNA and mRNA promoters.

CEBPZ- most overlaps were with mRNA promoters and very few were with
lncRNA

CHD2- 2nd highest amount of overlaps for both lncRNA and mRNA promoters-
with more being mRNA

CTCF- 3rd highest amount of overlaps for both lncRNA and mRNA promoters-
with more being mRNA

ELF1- greatest amount of overlaps with both lncRNA and mRNA promoters

EP300- least amount of overlaps with both lncRNA and mRNA promoters

# Now I want to compare the overlaps with lncRNA and mRNA promoters seperately

``` r
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids]) 

# mrna promoter overlaps
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/num_peaks_df.csv")
```

## results:

# 1) What is the difference in overlaps between mRNA and lncRNA promoters

For my proteins the majority of the peak overlaps are more mRNA compared
to lncRNA promoters- besides EP300 it has more overlaps with lncRNA.

# Now I am going to test if there is more binding over gene bodies than promoters

# I will seperate lncRNA and mRNA gene bodies to find the overlaps

``` r
source ("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/util/my_class_functions.R")
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, type = "counts")

num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                consensus_peaks, 
                                                type = "counts")

num_peaks_df$peaks_overlapping_genebody <- 
  rowSums(genebody_peak_counts)

# lncRNA gene bodies 
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA gene bodies
num_peaks_df$peaks_overlapping_mrna_genebody <- 
  rowSums(genebody_peak_counts[,mrna_gene_ids])
write_csv(num_peaks_df, "num_peaks_df.csv")
```

## results:

# 1) Do my proteins have more overlaps with promoters or genebodies?

CEBPZ- promoters

CHD2- genebody

CTCF- genebody by a lot

ELF1- genebody

EP300- genebody

# It is nice and all to find overlaps, but I am interested in how many proteins

# bind a specific promoter. I will use my handy “occurence” parameter in

# " count peaks per feature"

``` r
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, 
                                               type = "occurrence")
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

write.table(promoter_peak_occurence, "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")
```

``` r
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))

# "counts" just gave us value for the dbp of interest

write_csv(peak_occurence_df, "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/peak_occurence_dataframe.csv")
```

## results: I find the max number of proteins on a promoter to be X

10

# I am curious if my proteins are transcription factors so I will use the annotations

# in a cell paper I found and see

``` r
url <- "https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx"

destination_for_url <- "/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/TF_annotations/TF_annotations.xlsx" 

# to download we can use download.file
download.file(url, destination_for_url)


human_tfs <- readxl::read_excel("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/TF_annotations/TF_annotations.xlsx" 
,
                                sheet = 2, skip = 1)
```
  

``` r
# let's rename the 4th column to indicate if it is a TF.
names(human_tfs)[4] <- "is_tf"
num_peaks_df <- data.frame("dbp" = names(consensus_peaks),
                           "num_peaks" = sapply(consensus_peaks, length))

length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$Name)))
```

    ## [1] 0

``` r
human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(num_peaks_df$dbp), 1:4]


# adding new column names
names(human_tfs) <- c("ensembl_id",
                      "dbp",
                      "dbd",
                      "tf")
num_peaks_df <- merge(num_peaks_df, human_tfs, all.x = T) 

# Let's check how many NAs -- we should have some missing values  
dim(num_peaks_df[is.na(num_peaks_df$tf),])
```

    ## [1] 16  5

# Now I want to start plotting my results

# First I will see if there is a realtionship between peak number and total DNA covered

``` r
library(ggplot2)
library(tidyverse)


source("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/util/plotting_functions.R")
source("/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/util/_setup.R")

num_peaks_df <- read_csv('/scratch/Shares/rinnclass/CLASS_2023/corrinbowman/CLASS_2023/CLASSES/05_R_analyses/Class_exercise/Results/num_peaks_df.csv')
```

    ## Rows: 16 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): dbp
    ## dbl (5): num_peaks, total_peak_length, peaks_overlapping_promoters, peaks_ov...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ggplot(num_peaks_df, aes(x = num_peaks, 
                         y = total_peak_length)) +
  geom_point() 
```

![](Assignment_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Now I want to color my plot by wether the protein is a TF or not.

``` r
 # ggplot(num_peaks_df, aes(x = num_peaks, 
                # y = total_peak_length,
                 #color = tf == "Yes")) +
 # geom_point()
```

# I want to make a histogram of the number of peaks for each of my proteins

``` r
ggplot(num_peaks_df, aes(x = num_peaks)) +
  geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Assignment_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
hist
```

    ## function (x, ...) 
    ## UseMethod("hist")
    ## <bytecode: 0x1ddfab8>
    ## <environment: namespace:graphics>

# Now I want to facet this by the type of DNA binding domain my protein has.

``` r
#ggplot(num_peaks_df, aes(x = num_peaks, 
                # y = total_peak_length
              #  )) +
 # facet_wrap(dbd ~ .) +
 # geom_point() 
```

# Cool now I am ready to send my result to my collaborator as a

# Knitted document
