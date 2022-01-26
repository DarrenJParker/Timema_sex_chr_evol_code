import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:Ah')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

fasta_seq_dir_name = None
outprefix    = "testout"
sp_set = set(["Tbi", "Tce", "Tcm", "Tpa", "Tps"])
min_len = 1

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** alignments_to_seq.py | Written by DJP, 16/12/21 in Python 3.5 in Lausanne, Swiss ****\n")
		print("removes masking and gaps from selectome alignments")	
		print("\n**** Usage ****\n")
		print("python3 alignments_to_seq.py -i [alignment dir] -o [outprefix] \n\n")
		print("\n**** Example ****\n")
		print("python3 accessory_scripts/alignments_to_seq.py -i data/selection/Selectome_v07_Timema-nt_masked/ -o data/output/sel_out/Ortho \n\n")
		sys.exit(2)
	elif opt in ('-i'):
		fasta_seq_dir_name  = arg
	elif opt in ('-A'):
		sp_set = set(["Tbi", "Tce", "Tcm", "Tpa", "Tps", "Tte", "Tms", "Tsi", "Tge", "Tdi"])
	elif opt in ('-i'):
		min_len = arg
	elif opt in ('-o'):
		outprefix = arg
	else:
		print("i dont know")
		sys.exit(2)



### takes dir of fasta files, unwraps them, and adds seqs to a dict with the addition of the file id

def fasta_dir_to_dict(in_seq_dict):
	seq_dict = {}
	path = in_seq_dict
	for path, subdirs, files in os.walk(path):
		for name in files:
			if name.endswith("fas"):
				#print (os.path.join(path, name))
				
				## file ID ## might need to edit this
				file_id = name.split("/")[-1].split("_to_")[0]
				#print(file_id)
				in_fasta_file_name = os.path.join(path, name)
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
					seq_dict[name1 + "__FILEID__" + file_id] = seq1
			
				## tidyup
				seq_file_1.close()
				os.remove(output_fasta_name)
	
				#print("Read " + str(done) + " sequences from " + in_fasta_file_name)
	
	return(seq_dict)


full_seq_dict = fasta_dir_to_dict(fasta_seq_dir_name)  # data/selection/Selectome_v07_Timema-nt_masked
print(len(full_seq_dict))



for sp in sp_set:
	out_file = open(outprefix + "_" + sp + ".fa", "w")
	out_file.close()
	list_file = open(outprefix + "_" + sp + "_seqnames.txt", "w")
	list_file.close()


for s in full_seq_dict:

	seq = full_seq_dict.get(s)
	#print(seq)
	cleaned_seq = seq.replace("n", "").replace("-", "")
	#print(cleaned_seq )
	
	## check we kept codons
	
	if len(cleaned_seq) >= min_len:
		if len(cleaned_seq) % 3 != 0:
			print("ERROR\n")
			sys.exit()
		
		
			
		HOG_name = s.split("__FILEID__")[1].split(".")[0]
		sp_name = s.split("__FILEID__")[0]
		
		if sp_name in sp_set:
			out_file = open(outprefix + "_" + sp_name + ".fa", "a")
			out_file.write(">" + sp_name + "_" + HOG_name + "\n" + cleaned_seq + "\n")
			out_file.close()
			list_file = open(outprefix + "_" + sp_name + "_seqnames.txt", "a")
			list_file.write(sp_name + "_" + HOG_name + "\n")
			list_file.close()
		

print("\n\nDone, Robinson \n\n")
	
	
	
	
	
	
	
	





