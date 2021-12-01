#!/usr/bin/env Rscript
# SShekarriz Dec1,2021
# comparing B vs. SHCM15


#################
library(tidyverse)
#################

path="/datastore/shekas3/CEMG_SHCM15/assembled/DonorB_DMG_coverage"
patt=".flagstat"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = F))) %>%
         unnest() %>%
  mutate(sample.id = gsub(patt, "", sample.id)) %>%
  mutate(sample.id= gsub("DonorB_D", "B", sample.id))-> df

head(df)

df %>%
  rename(Reads= X1) %>%
  filter(grepl("in total", Reads)) %>%
  mutate(Type= paste("Total")) %>%
  mutate(Reads= gsub(" .*", "", Reads))-> df_total

df %>%
  rename(Reads= X1) %>%
  filter(grepl("+ 0 mapped", Reads)) %>%
  mutate(Type= paste("Mapped")) %>%
  mutate(Reads= gsub(" .*", "", Reads)) %>%
  # binding the total reads
  bind_rows(df_total) %>%
  spread(Type, Reads) %>%
  mutate(Unmapped= as.numeric(Total) - as.numeric(Mapped)) %>%
  mutate(Mpercent= as.numeric(Mapped)/as.numeric(Total)*100) %>%
  mutate(Upercent= as.numeric(Unmapped)/as.numeric(Total)*100) %>%
  select(sample.id, Mpercent, Upercent) %>%
  gather(Type, Percentage, -sample.id) -> tbl

tbl$sample.id <- factor(tbl$sample.id, 
                        levels = rev(c("B_2013",
                                   "B_2016",
                                   "B_May17",
                                   "B_Oct17")))
tbl$Type <- factor(tbl$Type, c("Upercent", "Mpercent"))

  ggplot(tbl, aes(Percentage, sample.id, fill=Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste(round(Percentage, digits = 1), "%", sep = "")),
              position = position_stack(vjust = .5), size = 2) +
  scale_fill_manual(values = c("Mpercent" = "#31a354",
                               "Upercent" = "#bdbdbd")) +
  theme_classic() +
  theme(text = element_text(size = 7),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
ggsave("/datastore/shekas3/CEMG_SHCM15/results/DonorB_mapped_SHCM15.png", width = 5, height = 4, units = "cm")


path="/datastore/shekas3/CEMG_SHCM15/DonorB_markers/genes_from_rsCEGs"
patt="_CEGs.blastout"
COLS <- c("Genome","qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
"qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp")
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = F))) %>%
         unnest() %>%
  mutate(sample.id = gsub(patt, "", sample.id)) -> blastOuts
colnames(blastOuts) <- COLS

blastOuts %>%
    filter(pident >= 90, qcovhsp >= 90) %>%
    mutate(temp=sseqid) %>%
    separate(temp, c("Ref_genomes", "Ref_contig"), sep = ";") %>%
    filter(Ref_genomes %in% c("GCF_902497355.1_P9094",
                              "GC313_hybrid",
                              "GC568")) %>%
   rename(gene_id= qseqid) %>%
   select(gene_id, Genome, Ref_genomes, Ref_contig) %>%
   mutate(gene_id= paste("g", gene_id, sep = "_"))-> gene_in_rsCEGs


# Species markers- genomes

path="/datastore/shekas3/CEMG_SHCM15/DonorB_markers/coverage"
patt="_to_SpeciesMarker.coverage.txt"

data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(Contig= `#rname` ) %>%
         # make sure to have "__" for all contigs
         separate(Contig, c("Contig", "Cluster"), sep = "__") %>%
         # make a new variable for genome names  
         mutate(Genome = case_when(grepl("Isolate_19", Contig) ~ paste("Fus"),
                                   grepl("RefSeq_89", Contig) ~ paste("Fpr"),
                                   grepl("GC568", Contig) ~ paste("Dor"),
                                   TRUE ~ as.character("Other"))) %>%
         mutate(length= endpos) %>%
         mutate(Type=paste("Conserved"))-> cons



# Strain marker -commonly engrafted genes from phylogeny-genomes
path="/datastore/shekas3/CEMG_SHCM15/DonorB_markers/coverage"
patt="_to_CEGs.coverage.txt"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(gene_id= `#rname`) %>%
         mutate(gene_id= paste("g", gene_id, sep = "_")) %>%
         left_join(gene_in_rsCEGs, by = "gene_id") %>%
         filter(!is.na(Genome)) %>% 
         # selecting gene with minimum 100 bps
         filter(endpos >= 100) %>%
  # all the genes within a genome together
  group_by(sample.id, Genome) %>%
  summarise(length= sum(endpos),
            numreads=sum(numreads),
            covbases=sum(covbases),
            coverage=covbases/length * 100) %>%
  mutate(Type= paste("CEGs"))-> cegs


# Visualize results
strain_col <- c("Fus" = "#e41a1c",
                "Dor" = "#377eb8",
                "Fpr" = "#984ea3")
 
# merge the two tables
cons %>%
  select(colnames(cegs)) %>%
  bind_rows(cegs) %>%
  select(sample.id, Genome, coverage, Type) %>%
  spread(Type, coverage) %>%
 ggplot(aes(CEGs, Conserved, color=Genome, fill=Genome)) +
    # setting the background colors
   annotate("rect", xmin=80,xmax=Inf,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#008837") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#a6dba0") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#bababa") +
   annotate("rect", xmin=80,xmax=Inf,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#404040") +
  # main plot component
  geom_point(alpha=0.6,position="dodge", size=2) +
  geom_hline(yintercept=80, linetype="dashed", 
                color = "black") +
  geom_vline(xintercept=80, linetype="dashed", 
                color = "black") +
  scale_color_manual(values = strain_col) +
  scale_fill_manual(values = strain_col) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 100)) + 
  scale_y_continuous(limits = c(0, 100)) +
  theme(text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = "right")
ggsave("/datastore/shekas3/CEMG_SHCM15/results/marker_in_SHCM15.png", width = 7, height = 5, units = "cm")


