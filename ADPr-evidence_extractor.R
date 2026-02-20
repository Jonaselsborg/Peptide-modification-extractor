# ------------------------------------------------------------------------------
# Quality filter for ADP-ribosylation site data from MaxQuant
#
# Input:
#   - evidence.txt (from MaxQuant, with ADPr search enabled)
#   - FASTA file used for the MaxQuant search (.fasta)
#
# Output:
#   - evidence_extracted_loc.txt (filtered and annotated ADPr sites)
#
# Author
#   Jonas D. Elsborg
#   jonas.elsborg@cpr.ku.dk
#   jonas.elsborg@path.ox.ac.uk
# ------------------------------------------------------------------------------

# clear env
rm(list = ls())

# libraries
library(tidyverse)
library(hablar)
library(seqinr)

# wd
#setwd("~/sequence extractor ADPr")

# warn
if (!file.exists("evidence.txt")) {
  stop("Could not find 'evidence.txt' in the working directory. Please ensure the file is present.")
}

fasta_files <- list.files(pattern = "\\.fasta$", full.names = TRUE)
if (length(fasta_files) == 0) {
  stop("No FASTA file found in the working directory. Please provide the FASTA used in MaxQuant.")
}

# data
evidence <- read_tsv(file = "evidence.txt") 

# define the grouping
evidence <- evidence %>% 
  mutate(experiment_group = str_remove(Experiment, pattern = "_[0-9]+$")) %>% 
  convert(fct(experiment_group))

# Define the match groups. 
# - If the same as experimental groups, then those groups define the boundaries.
# - If the match group is the same for all samples, it will match all-to-all.
# - you can define custom match groups. (match within PARPs, match within treatments etc.)

evidence <- evidence %>% 
  #mutate(match_group = experiment_group) %>%  # MBR within replicates
  mutate(match_group = "all-to-all") %>%  # if only one param group is defined, MBR between all samples (default MQ option)
  #mutate(match_group = str_extract(experiment_group, "^PARP[^-]+")) %>%  # MBR within custom group
  convert(fct(match_group))
  
# Read the FASTA file
fasta_files <- list.files(pattern = "\\.fasta$", full.names = TRUE)
fasta_file <- fasta_files[1] # take the first

# function from seqinr to read protein fasta files
sequences <- read.fasta(fasta_file, seqtype = "AA", as.string = TRUE)

# extract the uniprot ID's from the fasta, make them the names (used for functions to iter)
uniprot_ids <- sapply(names(sequences), function(x) str_split(x, "\\|")[[1]][2])
names(sequences) <- uniprot_ids
rm(uniprot_ids)

# function that returns the entire fasta sequence for a uniprot ID. Use the same fasta as search to ensure match
# Attention! this means any rows with ID's from MQ internal contamination DB is not found
get_sequence <- function(uniprot_id) {
  return(sequences[[uniprot_id]] %>% as.character(.))
}


# filter on sites. The script runs as long as the peptide is modified.
# you could also run it if you filter(Modifications != "Unmodified")
evidence_modified <- evidence %>% 
  filter(`ADP-ribosylation (CDEHKRSTY)`>0)

# for reference, without filtering.
evidence_modified_alladpr <- evidence_modified

# take just the MSMS where we could trace the MS1 peak for quantification
evidence_modified <- evidence_modified %>% 
  filter(Type == "MULTI-MSMS") 

# This function extracts the probability, position and residue per possibility within the modified sequence
extract_modifications <- function(peptide) {
  # Extract modified residues and probabilities
  matches <- str_extract_all(peptide, "[A-Z]\\(\\d+\\.?\\d*\\)")[[1]]
  
  if (length(matches) == 0) return(tibble(Position = integer(), Residue = character(), Probability = numeric()))
  
  # Extract residues and probabilities separately
  residues <- str_extract(matches, "^[A-Z]")
  probabilities <- as.numeric(str_extract(matches, "\\d+\\.?\\d*"))
  
  # Extract amino acid positions in string
  positions <- 
    str_replace_all(string = peptide,
                               pattern = "\\d+\\.?\\d*",
                               replacement = ";") %>% # puts ; for number
    str_remove_all(string = ., pattern = "\\(|\\)") %>% # removes parenthesis
    str_split_1(., pattern = ";") %>% # splits on ;
    str_length(.) %>% # calc. len
    cumsum(.) %>% # add the previous len to next (iterate the position)
    head(., -1) # the last length is the total length of peptide. But this is not a position for modification so it is removed
  
  tibble(Position = positions, Residue = residues, Probability = probabilities) # return
}

# Apply function
evidence_modified_cleaned <- evidence_modified %>%
  mutate(Extracted = map(`ADP-ribosylation (CDEHKRSTY) Probabilities`, extract_modifications)) %>%
  unnest(Extracted)


# lookup fasta for sequence
evidence_modified_cleaned_filtered <- evidence_modified_cleaned %>% 
  mutate(fasta_seq = map(`Leading razor protein`, get_sequence))

# unnest 
evidence_modified_cleaned_filtered <- evidence_modified_cleaned_filtered %>% 
  unnest_longer(fasta_seq)

# search for start and end, finding the ptm position
evidence_modified_cleaned_filtered <- evidence_modified_cleaned_filtered %>% 
  mutate(pep_start = str_locate(fasta_seq, Sequence)[, "start"],
         pep_end = str_locate(fasta_seq, Sequence)[, "end"]) %>% 
  mutate(ptm_position = pep_start+Position-1) 

# take the sites with ranked descending probability according to the number of identified modifications
evidence_modified_cleaned_filtered <- evidence_modified_cleaned_filtered %>% 
  group_by(id) %>% 
  arrange(desc(Probability), .by_group = TRUE) %>%  
  filter(row_number() <= first(`ADP-ribosylation (CDEHKRSTY)`)) %>%
  ungroup()

# MS/MS evidence with loc >= 0.9 is required for proper localisation
evidence_modified_localized <- evidence_modified_cleaned_filtered %>% 
  filter(Probability >= 0.9)

# clear memory
rm(evidence_modified_cleaned_filtered)
rm(evidence_modified_cleaned)

# map back intensities
# this finds all combinations of LOCALIZED mod_seq and which group they were found, index original df
evidence_mapped <- evidence %>% 
  semi_join(
    evidence_modified_localized %>% 
      select(`Modified sequence`, match_group) %>% distinct(),
    by = c("Modified sequence", "match_group")
  ) %>% filter(Type == "MULTI-MATCH")


# this then maps the sites within that modified sequence IF it was localized for the group
evidence_mapped <- evidence_mapped %>% 
  left_join(evidence_modified_localized %>% 
              select(`Modified sequence`, match_group, Residue, pep_start, pep_end, ptm_position) %>% distinct(), 
            by = c("Modified sequence", "match_group"),
            relationship = "many-to-many") 

# here additional stringency filters can be applied to safeguard faithful MBR transfer (e.g. trash bad matches by score)
# However, this can also easily be done with further data processing.


# combine the MS/MS with the localized sites and the MBR intensities
# note: evidences with localisation below 0.9 are not included, even for quantification.
evidence_modified_localized_mbr <- evidence_modified_localized %>% 
  select(-fasta_seq) %>% 
  bind_rows(evidence_mapped) %>% 
  arrange(id) 

# write tsv
evidence_modified_localized_mbr %>% write_tsv("evidence_extracted_loc.txt")
message("Filtered evidence written to 'evidence_extracted_loc.txt'")
