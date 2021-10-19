#!/usr/bin/env Rscript
#Shekarriz Jul14,2021
# Making a cumulative bar for each sample's contigs
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "quality"
  args[2] = "contig_cumulative.png"
}

#################
library(tidyverse)
#################

path=args[1]
Figure=args[2]

data.frame(sample.id = paste(dir(path, pattern = "*_headers.txt"))) %>%
mutate(file_contents = map(sample.id, ~ read_tsv(file.path(path, .), 
                                                 col_names = F))) %>%
unnest() %>%
mutate(sample.id = gsub("_contigs_edit_headers.txt", "", sample.id)) %>%
rename(Contig= "X1") %>%
mutate(Length = gsub(".*_length_", "", Contig)) %>%
mutate(Length = gsub("_cov_.*", "", Length)) %>%
  # select those that are above 2.5K
  filter(as.numeric(Length) >= 1000) %>%
  #select(sample.id, Length) %>%
  group_by(sample.id) %>% 
  mutate(Cumulative= cumsum(as.numeric(Length))) %>%
  group_by(sample.id) %>% 
  mutate(Contigid = row_number()) %>%
  mutate(Cumulative = Cumulative / 1000000) %>%
  mutate(Contigid2 = Contigid / 1000) %>%
  ggplot(aes(Contigid2, Cumulative)) +
  geom_line(size=2) +
  theme_bw() +
  facet_grid(~sample.id) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Cumulative length of assembled contigs over 1kb") +
  ylab("Cumulative length in Mbp") +
  xlab("Contigs are ordered from largest (contig #1) to smallest (x1000).") +
  theme(text = element_text(size= 10),
        axis.text.x = element_text(angle = 0))
ggsave(Figure, height = 8, width = 12, units = "cm")

