# R-script which reads in a gff file and a file containing sRNA coordinates, and then adds the sRNAs to the gff file.
# The sRNAs are added as a new feature type, and the feature type is named "sRNA".
# Date: 2023-06-15
# Author: Jakob Jung

# Load libraries
library(tidyverse)

# Read in gff file, skip lines starting with "#". Keep spaces in strings and not \
gff <- read.csv("./data/reference_sequences/Staphylococcus_aureus_101588.gff3",
                comment.char = "#", sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
  as_tibble()

# import eggnogg annotations:
eggnogg <- read.csv("./data/reference_sequences/eggnogg_output.tsv",
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  as_tibble()

eggnogg <- eggnogg %>% mutate(genome = gsub("CDS::([^:]+).*", "\\1", query),
                              start = as.numeric(gsub("CDS::[^:]+:(\\d+)-\\d+.*", "\\1", query)),
                              end = as.numeric(gsub("CDS::[^:]+:\\d+-(\\d+).*", "\\1", query)),
                                start=start+1) %>%
  select(query, genome, start, end, Description, Preferred_name, KEGG_Pathway)

# Add gene names (Preferred_name) to gff file based on start and end coordinates:
gff <- gff %>%
  mutate(start = V4,
         end = V5,
         genome = V1) %>%
  left_join(eggnogg, by = c("genome", "start", "end")) %>%
  mutate(V9 = paste0(V9, ";prefname=", Preferred_name, ";description=", Description, ";KEGG_Pathway=", KEGG_Pathway)) %>%
  select(-query, -genome, -start, -end, -Description, -Preferred_name, -KEGG_Pathway) %>%
  # get rid of ;translation=...; in V9:
    mutate(V9 = gsub(";translation=[^;]+", "", V9)) %>%
  mutate(V9 = gsub(" CDS", "", V9))


# Read in sRNA table:
sRNAs <- read.csv("./data/reference_sequences/blast_srnas.txt",
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
  as_tibble() %>% mutate(lt = gsub("N315_", "", V1))

# read in sRNA lengths:
sRNA_lengths <- read.csv("./data/reference_sequences/lengths.tsv",
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
  mutate(lt = gsub("N315_", "", V1))


# Add sRNA lengths to sRNA table:
sRNAs <- sRNAs %>%
  # Add sRNA lengths to sRNA table:
  mutate(length = sRNA_lengths$V2[match(V1, sRNA_lengths$V1)]) %>%
  # keep only if V7 is 0 to 2 and V8 is length of sRNA +- 2
  filter(V7 <= 2 & V8 >= length - 2 & V8 <= length + 2) %>%
  # add strand from gff_old:
    mutate(strand = ifelse(V13=="minus", "-", "+"))

# find duplicates and add a number _n to the end of v1:
sRNAs$lt <- make.unique(sRNAs$lt, sep = "_")



# Add sRNAs to gff file:
sRNAs_gff <- sRNAs %>%
  mutate(sRNA = "sRNA",C3="sRNA", C6 = ".", C8 = "0", C9 = paste0("locus_tag=", lt, ";type=sRNA;Name=", lt,
                                                                  ";prefname=", lt,";description=NA;KEGG_Pathway=NA")) %>%
  select(V2, sRNA, C3, V9, V10,C6, strand, C8, C9)


colnames(sRNAs_gff) <- colnames(gff)

# bind gff and sRNAs_gff:
gff_new <- bind_rows(gff, sRNAs_gff)

# write to file:
write.table(gff_new, "./data/reference_sequences/WG_with_sRNAs_Staphylococcus_aureus_101588.gff3", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)