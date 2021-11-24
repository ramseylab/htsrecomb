# version 2018.01.13 - added mapping quality filter

require(Rsamtools)
require(GenomicAlignments)
require(Biostrings)
require(pbapply)
require(GenomicRanges)

nreads_use_orig <- NULL  ## for debugging, set to 100 or 1000 or whatever; for no debugging, set to NULL

pboptions(type="txt") ## display the progress bar even when R is running in batch mode

args <- commandArgs(trailing=TRUE)
# args <- c("DxL2_Hybrid.fasta",
#          "DxL2_SNPs.txt",
#          "test.bam",
#          "DNA_paired_trimmed_cut.sorted_hit_variants.bam",
#          "s001.sorted_hit_variants.bam",
#          "9",
#          "4")

if (is.na(args[1])) {
   print("Usage:  Rscript genotype_reads_with_no_indels_coords_2Snps.R DxL2_Hybrid.fasta DxL2_SNPs.txt D50L50_paired_trimmed.mka.sorted_hitsnps#.bam min_mapq_value min_alleles_per_parent")
   stop("invalid command-line")
}

genome_fasta_file_name <- args[1]
snp_file_name          <- args[2]
bam_file_name          <- args[3]
min_mapq_value         <- as.integer(args[4])
min_alleles_per_parent <- as.integer(args[5])

stopifnot(is.integer(min_mapq_value) && min_mapq_value >= 0)
stopifnot(is.integer(min_alleles_per_parent) && min_alleles_per_parent > 0)

num_max_mismatch <- 25
num_cores <- ceiling(0.75*parallel::detectCores())

## from GitHub:    https://gist.github.com/SamBuckberry/9914246
readBAM <- function(bamFile){

  bam <- Rsamtools::scanBam(bamFile)
  
  # A function for collapsing the list of lists into a single list
  # as per the Rsamtools vignette
  .unlist <- function (x){
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
  
  bam_field <- names(bam[[1]])
  
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field

  #return a list that can be called as a data frame
  return(bam_df)
}

## read the SNP data from the spreadsheet provided by Steven Carrell
snp_data <- read.table(snp_file_name,
                       header=TRUE,
                       sep="\t",
                       stringsAsFactors=FALSE,
                       quote="",
                       comment.char="")
stopifnot("integer" %in% is(snp_data$Coordinate))

bam_file_df <- readBAM(bam_file_name)

bam_file_header <- Rsamtools::scanBamHeader(bam_file_name)
assembly_length <- bam_file_header[[bam_file_name]]$targets
stopifnot(length(assembly_length) == 1)

## converts a bionconductor DNAString object to an integer vector
dnastring_to_integer <- function(DNAString_obj ) {
    ret_vec <- as.integer(DNAString_obj)
    ret_vec[ret_vec == 15] <- NA
    ret_vec[ret_vec == 16] <- NA
    log2(ret_vec) + 1
}

## takes in a character vector (like:  c("A","C","C","T")) and returns integers (like:  c(1, 2, 2, 4))
dna_char_vec_to_int <- function(char_vec) {
    ret_vec <- rep(NA, length(char_vec))
    ret_vec[char_vec == "A"] <- 1
    ret_vec[char_vec == "C"] <- 2
    ret_vec[char_vec == "G"] <- 3
    ret_vec[char_vec == "T"] <- 4
    ret_vec[char_vec == "-"] <- NA
    ret_vec[char_vec == "N"] <- NA
    ret_vec
}

## convert alleles to integers
snp_data$strain1_int <- dna_char_vec_to_int(snp_data$strain1)
snp_data$strain2_int <- dna_char_vec_to_int(snp_data$strain2)

## get "chromosome" name
scaffold_names <- unique(as.character(bam_file_df$mrnm))
stopifnot(length(scaffold_names) == 1)
stopifnot(names(assembly_length)[1] == scaffold_names)

nreads <- nrow(bam_file_df)
nreads_use <- if (is.null(nreads_use_orig)) { nreads } else { nreads_use_orig }
                  
read_pair_names <- bam_file_df$qname[seq.int(from=1,
                                             to=nreads,
                                             by=2)]

## double-check that reads are organized in pairs
stopifnot(read_pair_names == bam_file_df$qname[seq.int(from=2,
                                                    to=nreads,
                                                    by=2)])

print(sprintf("Number of reads before filtering for mapping quality: %d", nrow(bam_file_df)))
read_pair_mapq_values_R1 <- bam_file_df$mapq[seq.int(from=1, to=nreads, by=2)]
read_pair_mapq_values_R2 <- bam_file_df$mapq[seq.int(from=2, to=nreads, by=2)]
read_pairs_keep_inds <- which(read_pair_mapq_values_R1 >= min_mapq_value &
                              read_pair_mapq_values_R2 >= min_mapq_value)
reads_keep_inds <- sort(c(2*read_pairs_keep_inds, (read_pairs_keep_inds-1)*2 + 1))
bam_file_df <- bam_file_df[reads_keep_inds, ]
nreads <- nrow(bam_file_df)
nreads_use <- if (is.null(nreads_use_orig)) { nreads } else { nreads_use_orig }
print(sprintf("Number of reads after filtering for mapping quality: %d", nrow(bam_file_df)))

## make a matrix of genotypes, with nrows equal to the size of the genome
## (for bigger organisms we'd have to do something more clever here)
snp_genotype_matrix <- data.frame(strain1=rep(NA, assembly_length),
                                  strain2=rep(NA, assembly_length))
snp_genotype_matrix$strain1[snp_data$Coordinate] <- snp_data$strain1_int
snp_genotype_matrix$strain2[snp_data$Coordinate] <- snp_data$strain2_int

genotype <- function(seq_fwd_strand_vec,  ## must be a vector of integers
                     seq_start_coord,
                     seq_end_coord,
                     snp_genotype_matrix,
                     read_grange,
                     read_qname) {

    ## seq_pos_range is the coordinate range of the read, in assembly coordinates
    seq_pos_range <- seq.int(from=seq_start_coord,
                             to=seq_end_coord,
                             by=1)

    pos_strain1_snp_matches <- seq_pos_range[which(seq_fwd_strand_vec == snp_genotype_matrix$strain1[seq_pos_range])]
    pos_strain2_snp_matches <- seq_pos_range[which(seq_fwd_strand_vec == snp_genotype_matrix$strain2[seq_pos_range])]
    

    status <- 'OK'
        
    pos_strain1_genotype_matches <- sort(pos_strain1_snp_matches)
    pos_strain2_genotype_matches <- sort(pos_strain2_snp_matches)

    list(strain1=pos_strain1_genotype_matches,
         strain2=pos_strain2_genotype_matches,
         status=status)
}

## read the hybrid genome sequence
genome_sequence <- Biostrings::readDNAStringSet(filepath=genome_fasta_file_name)[[1]]

## make a parallel computing cluster
cluster <- parallel::makeForkCluster(nnodes=num_cores)

reads_list <- pbapply::pblapply(1:nreads_use,
                     function(read_ind) {
                         strand <- bam_file_df$strand[read_ind]
                         pos <- bam_file_df$pos[read_ind]
                         cigar_str <- bam_file_df$cigar[read_ind]
                         
                         ## unpack the read sequence in reference coordinates, using the CIGAR string
                         read_seq <- GenomicAlignments::sequenceLayer(bam_file_df$seq[read_ind], cigar_str)[[1]]
                         qwidth <- length(read_seq)
                         read_seq_int <- dnastring_to_integer(read_seq)

                         list(start=pos, end=pos+qwidth-1, read_seq_int=read_seq_int, qname=bam_file_df$qname[read_ind])
                     }, cl=cluster)

clusterExport(cluster, "reads_list")

reads_granges <- GRanges(seqnames="hybrid_assembly",
                         ranges=IRanges(start=sapply(reads_list, "[[", "start"),
                                        end=sapply(reads_list, "[[", "end")))
clusterExport(cluster, "reads_granges")

## genotype all of the reads in the BAM file
res_list_list <- pbapply::pblapply(1:nreads_use,
                        function(read_ind) {
                            read_list <- reads_list[[read_ind]]
                            genotype(read_list$read_seq_int,
                                     read_list$start,
                                     read_list$end,
                                     snp_genotype_matrix,
                                     reads_granges[read_ind],
                                     read_list$qname)
                        }, cl=cluster)

clusterExport(cluster, "res_list_list")

reject_overlapping_read_pair <- unlist(pbapply::pblapply(1:(as.integer(nreads_use/2)),
                                                         function(read_pair_id) {
                                                             read1_id <- 2*(read_pair_id - 1) + 1
                                                             read2_id <- 2*(read_pair_id - 1) + 2
                                                             read1_strain1_pos <- res_list_list[[read1_id]]$strain1
                                                             read1_strain2_pos <- res_list_list[[read1_id]]$strain2
                                                             read2_strain1_pos <- res_list_list[[read2_id]]$strain1
                                                             read2_strain2_pos <- res_list_list[[read2_id]]$strain2
                                                             (length(intersect(read1_strain1_pos, read2_strain2_pos)) > 0 ||
                                                              length(intersect(read1_strain2_pos, read2_strain1_pos)) > 0)
                                                         }, cl=cluster))

print(sprintf("Number of read pairs rejected because of discordant alleles at a SNP where the reads overlap: %d", length(which(reject_overlapping_read_pair))))

## shut down the parallel cluster
stopCluster(cluster)

reject_read_overlapping_read_pair <- rep(reject_overlapping_read_pair, each=2)

strain1_genotype_matches_list <- lapply(res_list_list, "[[", "strain1")
strain2_genotype_matches_list <- lapply(res_list_list, "[[", "strain2")

strain1_genotype_matches_pos_str <- unlist(lapply(strain1_genotype_matches_list, paste, collapse=";"))
strain2_genotype_matches_pos_str <- unlist(lapply(strain2_genotype_matches_list, paste, collapse=";"))

status_str <- sapply(res_list_list, "[[", "status")
status_str[status_str == "OK" & reject_read_overlapping_read_pair] <- "overlapping_read_pair"

num_strain1_genotype_matches <- unlist(lapply(strain1_genotype_matches_list, length))
num_strain2_genotype_matches <- unlist(lapply(strain2_genotype_matches_list, length))

## add a "strain1" column and a "strain2" column to the BAM file data frame
bam_file_df <- data.frame(bam_file_df[1:nreads_use, ],
                          strain1_genotype_matches=num_strain1_genotype_matches,
                          strain2_genotype_matches=num_strain2_genotype_matches,
                          strain1_genotype_matches_pos=strain1_genotype_matches_pos_str,
                          strain2_genotype_matches_pos=strain2_genotype_matches_pos_str,
                          reject_read_due_to_variant_in_both_reads=reject_read_overlapping_read_pair,
                          status_str)

## save the complete BAM file data frame along with the new columns
write.table(bam_file_df,
            file=sprintf("%s_with_genotype.txt", tools::file_path_sans_ext(bam_file_name)),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

num_read1_strain1_alleles <- bam_file_df$strain1_genotype_matches[seq.int(from=1, to=nreads_use, by=2)]
num_read2_strain1_alleles <- bam_file_df$strain1_genotype_matches[seq.int(from=2, to=nreads_use, by=2)]

num_read1_strain2_alleles <- bam_file_df$strain2_genotype_matches[seq.int(from=1, to=nreads_use, by=2)]
num_read2_strain2_alleles <- bam_file_df$strain2_genotype_matches[seq.int(from=2, to=nreads_use, by=2)]

#num_strain1_alleles <- num_read1_strain1_alleles + num_read2_strain1_alleles
#num_strain2_alleles <- num_read1_strain2_alleles + num_read2_strain2_alleles

read_ids_for_pairs <- lapply(seq.int(from=1, to=nreads_use, by=2), function(x) {c(x, x+1)})
strain1_match_pos <- lapply(lapply(read_ids_for_pairs, function(x) {unlist(unique(strain1_genotype_matches_list[x]))}), unique)
strain2_match_pos <- lapply(lapply(read_ids_for_pairs, function(x) {unlist(unique(strain2_genotype_matches_list[x]))}), unique)
num_strain1_alleles <- sapply(strain1_match_pos, length)
num_strain2_alleles <- sapply(strain2_match_pos, length)

## look for read pairs that have probable recombination events
interesting_read_pairs <- which(! reject_overlapping_read_pair &
                                num_strain1_alleles >= min_alleles_per_parent &
                                num_strain2_alleles >= min_alleles_per_parent) 

print(sprintf("We got %d read pairs that look like they might be recombination events", length(interesting_read_pairs)))

interesting_reads_rows <- sort(c(2*interesting_read_pairs, 2*(interesting_read_pairs - 1) + 1))

## save the complete information about read pairs that might contain recombination events, in a separate file
write.table(bam_file_df[interesting_reads_rows, ],
            file=sprintf("%s_with_genotype_candidate_recombinations.txt", tools::file_path_sans_ext(bam_file_name)),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

write.table(bam_file_df[-interesting_reads_rows, ],
            file=sprintf("%s_with_genotype_NOT_candidate_recombinations.txt", tools::file_path_sans_ext(bam_file_name)),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)
