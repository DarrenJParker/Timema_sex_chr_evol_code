import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'f:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_file_name = None
outprefix    = "testout"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** filter_alignment.py | Written by DJP, 04/04/22 in Python 3.5 in Bangor, UK ****\n")
		print("Filters selectome alignments to include only sexual species. Alignments with <4 species are discarded")	
		print("\n**** Usage****\n")
		print("python3 filter_alignment.py -f [selectome fasta file] \n\n")
		sys.exit(2)
	elif opt in ('-f'):
		in_file_name  = arg
	else:
		print("i dont know")
		sys.exit(2)


#############################################
### tree topologies

tree_dict = {}
tree_dict["TbiTceTcmTpaTps"] = "(((Tps,Tcm),Tce),(Tpa,Tbi));"
tree_dict["TceTcmTpaTps"] = "(((Tps,Tcm),Tce),Tpa);"
tree_dict["TbiTcmTpaTps"] = "((Tps,Tcm),(Tpa,Tbi));"
tree_dict["TbiTceTpaTps"] = "((Tps,Tce),(Tpa,Tbi));"
tree_dict["TbiTceTcmTps"] = "(((Tps,Tcm),Tce),Tbi);"
tree_dict["TbiTceTcmTpa"] = "((Tcm,Tce),(Tpa,Tbi));"

### read in file

def fasta_to_dict(in_fasta_file_name):
	output_fasta_name = in_fasta_file_name + ".TEMP_extract_fasta_file" 
	
	output_file = open(output_fasta_name, "w")
	#print("\nUnwrapping fasta file")
	count = 0
	in_file = open(in_fasta_file_name)
	for line in in_file:
		count = count + 1
		line = line.rstrip("\n")
		if line.startswith(">") and count == 1:
			output_file.write(line + "\n")
		elif line.startswith(">") and count > 1:
			output_file.write("\n" + line + "\n")
		else: 
			output_file.write(line)	
	
	output_file.close()
	
	
	### add seqs to dictionary
	name_list = []
	seq_list = []
	seq_dict = {}
	
	done = 0
	seq_file_1 = open(output_fasta_name)
	for line in seq_file_1:
		lineA = line.rstrip("\n")
		if lineA.startswith(">"):
			lineB = lineA.lstrip(">")
			name_list.append(lineB)
		else:
			seq_list.append(lineA)
			done = done + 1
			seq_len = len(lineA)
	
	for element in range(0,len(name_list)):
		name1 = name_list[element]
		seq1 = seq_list[element].replace(" ", "") ## remove gaps if seq comes from gblocks 
		seq_dict[name1] = seq1

	## tidyup
	seq_file_1.close()
	os.remove(output_fasta_name)
	
	#print("Read " + str(done) + " sequences from " + in_fasta_file_name)
	
	return([seq_dict, name_list])


fasta_seqs_dict  = fasta_to_dict(in_file_name)[0]
fasta_seqs_names = fasta_to_dict(in_file_name)[1]



all_sexual_sp = set(["Tbi", "Tce", "Tcm", "Tpa", "Tps"])
sex_sp_in_file = []

for s in fasta_seqs_names:
	if s in all_sexual_sp:
		sex_sp_in_file.append(s)
# 
# print(fasta_seqs_names)
# print(sex_sp_in_file)

#print("N sex sp = " + str(len(sex_sp_in_file)))



if len(sex_sp_in_file) >= 4:
	sex_sp_in_file_sorted = sorted(sex_sp_in_file)
	sex_sp_in_file_str = ""
	for i in sex_sp_in_file_sorted:
		sex_sp_in_file_str = sex_sp_in_file_str + i
	tree_top = tree_dict.get(sex_sp_in_file_str)
	output_tree_file = open(in_file_name.replace(".nt_masked.fas", ".treetop.nwk"), "w")
	output_tree_file.write(tree_top)
	
	output_file = open(in_file_name.replace(".nt_masked.fas", ".nt_masked_sexonly.fas"), "w")
	for s in sex_sp_in_file:
		seq = fasta_seqs_dict.get(s)
		output_file.write(">" + s + "\n" + seq + "\n")
	




