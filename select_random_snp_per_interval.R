#####################################################################
### Function to extract aleatory the SNPs on a GWASpoly "format" file
#####################################################################

select_random_snp_per_interval <- function(snp_data_chr, interval_length = 200) {
  # Make sure data is sorted by Position
  snp_data_chr <- snp_data_chr[order(snp_data_chr$Position),]
  
  # Create intervals
  max_pos <- max(snp_data_chr$Position)
  intervals <- seq(0, max_pos, by = interval_length)
  
  selected_snps <- list()
  
  for (i in seq_along(intervals)[-length(intervals)]) {
    snps_in_interval <- snp_data_chr[snp_data_chr$Position > intervals[i] & snp_data_chr$Position <= intervals[i + 1], ]
    
    if (nrow(snps_in_interval) > 0) {
      selected_snp <- snps_in_interval[sample(nrow(snps_in_interval), 1), ]
      selected_snps[[i]] <- selected_snp
    }
  }
  
  do.call(rbind, selected_snps)
}
