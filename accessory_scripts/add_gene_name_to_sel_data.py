import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:s:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

orth_file_name    = None
pos_sel_file_name = None
outprefix    = "testout"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** add_gene_name_to_positive_sel_data.py | Written by DJP, 06/05/21 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Add gene names to positive sel file")	
		print("\n**** Usage****\n")
		print("python3 accessory_scripts/add_gene_name_to_positive_sel_data.py -i data/selection/Jaron_SM_Table_4.tsv -s data/selection/timema_543_branches_with-ncat-codon-rate_sites_with_h0.tsv  \n\n")
		sys.exit(2)
	elif opt in ('-i'):
		orth_file_name  = arg
	elif opt in ('-s'):
		pos_sel_file_name = arg
	elif opt in ('-o'):
		outprefix = arg
	else:
		print("i dont know")
		sys.exit(2)

### read in file
line_N = 0
hog_to_gene_dict = {}


orth_file = open(orth_file_name)
for line in orth_file:
	if not line.startswith("#"):
		line_N = line_N + 1
		line = line.rstrip("\n").split("\t")
		if line_N == 1:
			sp_order_l = line
		else:
			HOG_name = line[0]
			for g in range(1,len(line)):
				hog_to_gene_dict[sp_order_l[g] + HOG_name] = line[g]
	
orth_file.close()

outfile = open(pos_sel_file_name.replace(".tsv", "") + "_wgenename.tsv" , "w")
pos_sel_file = open(pos_sel_file_name)
want_sp = set(["Tbi", "Tce", "Tcm", "Tpa", "Tps"])
line_N = 0
for line in pos_sel_file:
	line_N = line_N + 1
	line = line.rstrip("\n")
	if line_N == 1:
		#print(line)
		outfile.write(line + "\t" + "gene_name" + "\n")
	
	sp = line.rstrip("\n").split("\t")[2]
	if sp in want_sp:	
		gene_name = hog_to_gene_dict.get(sp + line.rstrip("\n").split("\t")[0]).split("-")[0]
		#print(gene_name)
		outfile.write(line + "\t" + gene_name + "\n")

print("\n\nDone, Bees\n\n\n")






