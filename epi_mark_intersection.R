
epi_mark_intersection <- function(mark, samples, q_value){
    
    # mark = "DNase.macs2"
    # q_value = 0.05 
    # samples <- c("E081", "E082")
    
 ### Pre-processing and filtering data ###
    peak_dir <- "Roadmap_consolidated_peaks/"
    
    E_sample_list <- lapply(samples, function(x) read.table(paste0(peak_dir, x, "-", mark, ".narrowPeak")))
    E_samples <- rbindlist(E_sample_list)
    
    colnames(E_samples) <- c("Chr", "Start", "Stop", "Name", "Score", "Strand", 
                             "SignalValue", "Nlog10pValue", "Nlog10qValue", "Peak")
    
    # Adding p q and peak width columns
    E_samples$pValue <- 10^-E_samples$Nlog10pValue
    E_samples$qValue <- 10^-E_samples$Nlog10qValue
    
    # Peak_width
    E_samples$Peak_width <- abs(E_samples$Start - E_samples$Stop)
    
    # Increases epi peak statistical stringency
    E_samples_filtered <- dplyr::filter(E_samples, qValue < q_value)
    
    # Conversion to GRanges 
    E_samples_GR <- makeGRangesFromDataFrame(E_samples_filtered)
    
    epi_peaks_filtered <- GenomicRanges::reduce(E_samples_GR)
    
 ### Adding Activity group column ###
    data.predict_w_hg19_coord$Act_group <- ifelse(data.predict_w_hg19_coord$Pvalue > 0.05, "Non_significant", NA)
    
    data.predict_w_hg19_coord$Act_group <- ifelse(data.predict_w_hg19_coord$Pvalue < 0.05, "Significant",
                                              data.predict_w_hg19_coord$Act_group)
    
    
    n_of_non_sig <-  sum(data.predict_w_hg19_coord$Act_group == "Non_significant")
    n_of_sig <- sum(data.predict_w_hg19_coord$Act_group == "Significant")
    
    # Counting overlapping amplicons
    data.predict_w_hg19_coord_GR_Non_significant <- GRanges(dplyr::filter(data.predict_w_hg19_coord, Act_group == "Non_significant"))
    data.predict_w_hg19_coord_GR_Significant <- GRanges(dplyr::filter(data.predict_w_hg19_coord, Act_group == "Significant"))
    
    
    count_overlapping_amps <- function(x){
        length(unique(as.data.frame(findOverlaps(x, epi_peaks_filtered))$queryHits))
    }
    
    n_of_non_sig_intersect <-  count_overlapping_amps(data.predict_w_hg19_coord_GR_Non_significant)
    n_of_sig_intersect <- count_overlapping_amps(data.predict_w_hg19_coord_GR_Significant)
    
    
    
    df <- as.data.frame(t(data.frame(
        "Non_significant" = c(n_of_non_sig_intersect, n_of_non_sig),
        "Significant" = c(n_of_sig_intersect, n_of_sig))))
        
    
    colnames(df) <- c("Intersecting_amplicons", "Amplicons_in_a_group")
    
    df$Fraction_intersecting <- df$Intersecting_amplicons / df$Amplicons_in_a_group  
    df$Condition <- rownames(df)
    
    df$Epi_mark <- mark
    df$Epi_peak_q_value <- q_value
    
    rownames(df) <- NULL
    df <- df[,c(4:6,1:3)]
    
    df$n_of_all_amplicons <- nrow(data.predict_w_hg19_coord)
    #df$n_of_epi_peaks <- n_of_unfiltered_peaks
    df$n_of_epi_peaks_filtered <- length(epi_peaks_filtered)
    
 ### Overlap Intersection P value - nor needed ###    
    
    # This is excellent: https://rdrr.io/bioc/regioneR/man/overlapPermTest.html
    # library(regioneR)
    # How to specify the min number of permutations: https://stats.stackexchange.com/questions/80025/required-number-of-permutations-for-a-permutation-based-p-value
    
#    ll <- list(data.predict_w_hg19_coord_GR_Low,
#               data.predict_w_hg19_coord_GR_High, 
#               data.predict_w_hg19_coord_GR_Middle)
#    
#    clust <- makeCluster(3)
#    clusterExport(clust, varlist=c("ll", "epi_peaks_filtered"), envir=environment())
#    
#    P_values_overlap <- parSapply(clust, 1:3, function(x){ 
#        
#        a <-  regioneR::overlapPermTest(ll[[x]], epi_peaks_filtered, ntimes = 200, 
#                                        force.parallel = TRUE, verbose = TRUE)
#        
#        a$numOverlaps$pval
#    })
#    stopCluster(clust)  
#    
#    df$Perm_PVal_overlap <- P_values_overlap
    
    #df
    
    ###############  Sampling permutation test #####################
    n_of_non_sig_intersect
    n_of_sig_intersect
    n_of_non_sig
    n_of_sig
    
    
    all_intersecting <- sum(n_of_non_sig_intersect, n_of_sig_intersect)
    all_amp <- nrow(data.predict_w_hg19_coord)
    
    # Random vector or 0s and 1s, all_amp values, with the number of ones =  all_intersecting
    set.seed(1234)
    random_vector <- sample(c(rep(1, all_intersecting), rep(0, all_amp - all_intersecting)))
    
    
    # Permutation test Non_sig
    set.seed(1234)
    NS_perm <- replicate(20000, sum(sample(random_vector, n_of_non_sig, FALSE)))
    
    n_of_NS_scaled <- (n_of_non_sig_intersect - mean(NS_perm)) / sd(NS_perm)
    set.seed(1234)
    NS_p <- pnorm(n_of_NS_scaled, lower.tail=FALSE)
    
    
    # Permutation test Significant
    set.seed(1234)
    S_perm <- replicate(20000, sum(sample(random_vector, n_of_sig, FALSE)))
    n_of_S_scaled <- (n_of_sig_intersect - mean(S_perm)) / sd(S_perm)
    set.seed(1234)
    S_p <- pnorm(n_of_S_scaled , lower.tail=FALSE)
    
   
    
    P_val_df <- data.frame("Act_group" = c("Non_significant", "Significant"),
                           "Z_scaled_intersections_sampling" = c(n_of_NS_scaled, n_of_S_scaled),
                           "P_values_perm" = c(NS_p, S_p))
    
    merge(df, P_val_df, by.x = "Condition", by.y = "Act_group")
    
}
