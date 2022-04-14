import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:c:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_dir_name = None
gene_dir = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** m8a_process.py | Written by DJP, 13/04/22 in Python 3.5 in Bangor, UK ****\n")
		print("m8a_process.py")	
		print("\n**** Usage****\n")
		print("python3 accessory_scripts/m8a_process.py -i data/selection/Selectome_timema_M8 -c data/counts/ \n\n")
		sys.exit(2)
	elif opt in ('-i'):
		in_dir_name = arg
	elif opt in ('-c'):
		gene_dir = arg
	else:
		print("i dont know")
		sys.exit(2)
		

### read selectome files

D_dict = {}
HOG_sp_dict = {}

all_sel_HOGs = set()

path = in_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		HOG = name.split(".")[0]
		if name.endswith(".godonM8.std.out"):
			curr_file = open(os.path.join(path, name))
			all_sel_HOGs.add(HOG)
			for line in curr_file:
				line = line.strip()
				if line.startswith("Final D="):
					final_D = line.split("Final D=")[1]
					D_dict[HOG ] = final_D
		if name.endswith(".fas"):
			curr_file = open(os.path.join(path, name))
			sp_in_HOG = []
			for line in curr_file:
				line = line.strip()
				if line.startswith(">"):
					sp_in_HOG.append(line.replace(">", "").upper())
			HOG_sp_dict[HOG] = sp_in_HOG


print(len(all_sel_HOGs))
all_want_HOGs = sorted(list(all_sel_HOGs))



### get chr info for each gene from each species

path = gene_dir
chr_class_dict = {}


for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith("_1000_chr_info_and_counts.csv"):
			if not name.startswith("TbiTceTcmTpaTps"):
				#print(name)
				line_N = 0
				curr_file = open(os.path.join(path, name))
				for line in curr_file:
					line_N = line_N + 1
					line = line.strip().split(",")
					if line_N == 1:
						gene_index      = line.index("gene_id")
						chr_soft_index = line.index("chr_soft")
					else:
						gene = line[gene_index]
						chr_soft = line[chr_soft_index]
						chr_class_dict[gene] = chr_soft



## gene to HOG

HOG_to_genes_dict = {}

seen_HOG = set()
gene_to_HOG_file = open(os.path.join(in_dir_name, "og2gene_orthodb_All_final.txt"))
for line in gene_to_HOG_file:
	line = line.strip().split("\t")
	#print(line)
	HOG = "HOG_" + line[0]
	if HOG not in seen_HOG:
		seen_HOG.add(HOG)
		HOG_to_genes_dict[HOG] = [line[1]]
	else:
		rec = HOG_to_genes_dict.get(HOG)
		rec.append(line[1])
		HOG_to_genes_dict[HOG] = rec
	

## output
output_file = open("data/selection/selectome_M8_wchr.csv", "w")
output_file.write("HOG,soft_chr,finalD\n")

want_chr_class = set(["XXXX", "XXXXX", "AAAA", "AAAAA"])
N_HOGs_not_classified = 0
N_HOGs_classified = 0

for h in all_want_HOGs:
	#print(h)
	genes = HOG_to_genes_dict.get(h)
	sp_in_HOG = set(HOG_sp_dict.get(h))
	final_D = D_dict.get(h)

	chr_c_all = ""
	for g in genes:
		g = g.split("-")[0]
		sp = g.split("_")[0]
		if sp in sp_in_HOG:
			chr_c = chr_class_dict.get(g)
			chr_c_all = chr_c_all + chr_c
	
	if chr_c_all in want_chr_class:
		N_HOGs_classified = N_HOGs_classified + 1
		output_file.write(HOG + "," + chr_c_all + "," + final_D + "\n")
	else:
		N_HOGs_not_classified = N_HOGs_not_classified + 1

print("\nN HOGs classified = " + str(N_HOGs_classified))
print("N HOGs not classified = " + str(N_HOGs_not_classified))

print("\n\nFinished, Snow\n\n")
