######################
## SShekarriz Oct23,19
##

import os
from Bio import SeqIO
import sys
import glob

##################################
# HOW TO USE THIS CODE:
#
#python extract_ContigIds_frombins.py [*.fasta] > [Trimmed.fasta]
# this script selects contigs longer than 20,000bp for a list/single fasta

#################################

long_sequence = {}
list_of_files = glob.glob(sys.argv[1])
for file_name in sys.argv[1:]:
	R = SeqIO.parse(file_name, "fasta")
	for records in R:
		if len(records.seq) >= 1000:
			with open(file_name + "_out", "a") as handle: 
				SeqIO.write(records, handle, "fasta")






	
