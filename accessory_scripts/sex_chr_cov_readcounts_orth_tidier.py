### sex_chr_cov_readcounts_tidier.py

import sys
import os
import getopt
import decimal
from decimal import *

try:
	opts, args = getopt.getopt(sys.argv[1:], 'o:r:e:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


out_prefix     = "Testttttt"
readcount_dir  = None
readcount_ext  = "1000_chr_info_and_counts.csv"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** sex_chr_cov_readcounts_orth_tidier.py | Written by DJP, 30/04/21 in Python 3.5 in Lausanne, Swiss ****\n")
		print("brings reaccounts etc for orths together from sex_chr_cov_readcounts_tidier.py")
		
		print("\n**** USAGE **** \n")
		print("python3 sex_chr_cov_readcounts_orth_tidier.py -r [read_count file] -l [gene length file] -c [M/F coverage file] -y [orth_matrix_filename] -L [linkage_filename] -g [gff file] -o [out prefix] [options] \n")
		
		print("\n**** OPTIONS ****\n")
		print("-o\toutput prefix")
		print("-r\tread_count files - output from sex_chr_cov_readcounts_tidier.py")
		print("-e\textention of read_count files\n\n")

		sys.exit(2)

	elif opt in ('-o'):
		out_prefix = arg
	elif opt in ('-r'):
		readcount_dir = arg
	elif opt in ('-e'):
		readcount_ext   = arg
	else:
		print("i dont know")
		sys.exit(2)


if readcount_dir == None:
	print("\nError: missing readcount_dir")
	sys.exit(2)
	

################################################################################################################################################
######## get read counts for all species

all_dict = {}
header_dict = {}
sp_all = set()
orth_all = set()

path = readcount_dir
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith(readcount_ext ):
			print (os.path.join(path, name))
			curr_file = open(os.path.join(path, name))
			sp = name.split("_")[0]
			sp_all.add(sp)
			print(sp)
			line_N = 0
			for line in curr_file:
				line_N = line_N + 1
				line = line.strip()
				if line_N == 1:
					header_dict[sp] = line
				else:
					orth = line.split(",")[2]
					if orth != "NA":
						orth_all.add(orth)
						a = 1
						all_dict[sp + orth] = line
						#print(line)
						
				
			
### output

out_file = open(out_prefix + "_orth_" + readcount_ext, "w")
sp_all_ls   = sorted(list(sp_all))
orth_all_ls = sorted(list(orth_all))
print(sp_all_ls)

head_out = ""
for sp in sp_all_ls:
	head = header_dict.get(sp).split(",")
	for i in head:
		if i.startswith(sp):
			head_out = head_out + "," + i
		else:
			head_out = head_out + "," + sp + "_"+ i

out_file.write("HOG" + head_out + "\n")

for orth in orth_all_ls:
	out_line = ""
	for sp in sp_all_ls:
		out_line = out_line + "," + all_dict.get(sp + orth)
	
	out_file.write(orth + out_line + "\n")
	


print("\n\nFinished, Honey.\n\n\n")




