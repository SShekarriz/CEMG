#!/usr/bin/env Rscript
# SShekarriz Jun03,2021
# Making a Bin info file to be imported into anvio
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "input.txt"
  args[2] = "output.txt"
}

#################
library(tidyverse)
#################

tbl <- read.csv(args[1], sep="\t", header=FALSE)

tbl %>%
 separate(V1, c("binN", "contig"), sep = ":>") %>%
 mutate(bins= gsub("\\/.*_bins\\/", "", binN)) %>%
 mutate(bins= gsub("bin\\.", "bin_", bins)) %>%
 mutate(bins= gsub(".fa", "", bins)) %>%
 select(contig, bins) %>%
 write.table(args[2], sep="\t", quote=F, row.names=F, col.names=F)

