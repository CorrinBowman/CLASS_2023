# setting up the function: all what we just wrote above:
# Get's consensus peak file path, select .broadPeak and then get dbp name from file name
import_peaks <- function(consensus_file_path = broadpeakfilepath) {
  peak_files <- list.files(consensus_file_path, full.names = T, pattern = ".broadPeak")
  dbp_name <- sapply(peak_files, function(x){
    y <-  str_extract(x, "([^\\/]+$)")
    # NOTE ALTERNATIVE Gsub
    # gsub("_peaks.broadPeak", "", y)
    paste(unlist(strsplit(y, "_"))[c(1,2)], collapse = "_") 
    
  })
  
  # NOW FOR LOOP TO IMPORT FILES
  # setting an empty list "peak_list" to be filled by for loop
  peak_list <- c() 
  
  # the for loop !
  for(i in 1:length(peak_files)) {
    # Import peaks
    peaks <- rtracklayer::import(peak_files[i])
    # Append this GRanges object to the of the list.
    peak_list <- c(peak_list, peaks)
    # Name the list elements by their TF name.
    names(peak_list)[length(peak_list)] <- dbp_name[i] # need to name file so dont just have a list of them 
  }
  return(peak_list)
}












#Let's make the function for peak overlaps 

#' CREATE CONSENSUS PEAKS
#' this function will take multiple replicate .broadPeak files (also narrow)
#' find peaks that overlap in all the replicates. 
#' @description 
#' input set of chipseq replicate peak files
#' this function then creates one merged file peaks in all samples
#' @param dbp
#' This will be extracted with names(GR_list) in the lapply at end of fun
#' You will need a "dbps" or some object for the lapply that has the 
#' name of each dbp in the named GRanges list
#' 
#' @param peak_list
#' Named list of GRanges for each chipseq replicate
#' peak_list can be generated using import_peaks function above
consensus_from_reduced <- function(dbp, peak_list) {
  dbp_peaks <- peak_list[grepl(as.character(dbp), names(peak_list))]
  suppressWarnings(all_peaks <- GenomicRanges::reduce(unlist(as(dbp_peaks, "GRangesList"))))
  all_peaks <- all_peaks[grepl("chr", seqnames(all_peaks))] # taking out anything that doesnt have chr
  
  # peak_exists <- lapply(dbp_peaks, function(x) {
  #   as.numeric(countOverlaps(all_peaks, x) > 0))
  # }) %>%
  # bind_rows() OR bind_cols()
  peak_exists <- matrix(NA, nrow = length(all_peaks), ncol = length(dbp_peaks))
  for(i in 1:length(dbp_peaks)) {
    suppressWarnings(peak_exists[,i] <- as.numeric(countOverlaps(all_peaks, dbp_peaks[[i]]) > 0))
  }
  # filter to consensus requiring peaks to be in all replicates
  dbp_consensus <- all_peaks[rowSums(peak_exists) == ncol(peak_exists)]
  # Required only two replicates == dbp_consensus <- all_peaks[rowSums(peak_exists) > 1]
  return(dbp_consensus)
}




# Now running the function consensus_from_reduced


# loading peak list
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)
# getting dbp input as dbps
dbps <- unique(sapply(names(peak_list), function(x) {
  unlist(strsplit(x, "_"))[1]
}))
# Let's run it !
consensus_list <- lapply(dbps, consensus_from_reduced, peak_list) # think of lapply as another for loop to run all of the proteins 
# adding names to results
names(consensus_list) <- dbps # need to name consensus list 
sapply(consensus_list, length) # looking out length now of peaks that overlap- make sure not all the same number


