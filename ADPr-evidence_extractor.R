# coded by Jonas Elsborg
# htz513
#
# script to look at raw evidence.txt from MQ.
# 
# In an unfiltered manner, consider all evidences
# and assign intensities per experiment. 
# 
#
#
# Version 1.1
# Oxford 
# made to be generalized for any file


# clear env
rm(list = ls())

# libraries
library(tidyverse)
library(hablar)
library(seqinr)

# wd
set_wd_to_script_path()

# data
evidence <- read_tsv(file = "evidence.txt") 

# define the grouping
evidence <- evidence %>% 
  mutate(experiment_group = str_remove(Experiment, pattern = "_[0-9]+$"))

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
# nb! this means any rows with ID's from MQ internal contamination DB is not found
get_sequence <- function(uniprot_id) {
  return(sequences[[uniprot_id]] %>% as.character(.))
}


# filter on sites. The script runs as long as the peptide is modfied.
# you could also run it if you filter(Modifications != "Unmodifed")
evidence_modified <- evidence %>% 
  filter(`ADP-ribosylation (CDEHKRSTY)`>0)

# for reference, without filtering.
evidence_modified_alladpr <- evidence_modified

# take just the MSMS where we could trace the peak for quantification
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

evidence_modified_localized <- evidence_modified_cleaned_filtered %>% 
  filter(Probability >= 0.9)

# clear 
rm(evidence_modified_cleaned_filtered)
rm(evidence_modified_cleaned)

# map back intensities
# this finds all combinations of LOCALIZED mod_seq and which group they were found, index original df
evidence_mapped <- evidence %>% 
  semi_join(
    evidence_modified_localized %>% 
      select(`Modified sequence`, experiment_group) %>% distinct(),
    by = c("Modified sequence", "experiment_group")
  ) %>% filter(Type == "MULTI-MATCH")


# this then maps the sites within that mod seq if it was localized for the group
evidence_mapped <- evidence_mapped %>% 
  left_join(evidence_modified_localized %>% 
              select(`Modified sequence`, experiment_group, Residue, pep_start, pep_end, ptm_position) %>% distinct(), 
            by = c("Modified sequence", "experiment_group"),
            relationship = "many-to-many") 

# here additional stringency filters can be applied to safeguard faitful MBR transfer (e.g. trash bad matches by score)

# combine the msms with the localized sites and the mbr intensities
evidence_modified_localized_mbr <- evidence_modified_localized %>% 
  select(-fasta_seq) %>% 
  bind_rows(evidence_mapped) %>% 
  arrange(id) 

# write tsv
evidence_modified_localized_mbr %>% write_tsv("evidence_extracted_loc.txt")

 