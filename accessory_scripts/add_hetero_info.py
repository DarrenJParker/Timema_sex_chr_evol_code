import sys
import os
import getopt
import re
import decimal

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:s:c:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_dir_name = None
cov_ext = "_est.ml"
species_want = None
coverage_with_LG_filename = None


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** add_hetero_info.py| Written by DJP, 03/01/20 in Python 3.5 in Stoke-on-Trent, UK ****\n")
		print("\n**** USAGE example ****\n")
		print("python3 add_hetero_info.py -i /Users/dparker/Documents/University/Lausanne/Sex_chromosomes/v8_angsD_out/mapped_as_paired/ -s Tbi -c /Users/dparker/Documents/University/Lausanne/Sex_chromosomes/v8_cov_contig/Tbi_pairedcov_minlen=1000_contig_cov_wLGinfo.csv \n\n\n")

		sys.exit(2)
	
	elif opt in ('-i'):
		in_dir_name = arg
	elif opt in ('-s'):
		species_want = arg
	elif opt in ('-c'):
		coverage_with_LG_filename = arg
	else:
		print("i dont know")
		sys.exit(2)

#############################################################################################################################################
#### read in hetero info files into a dict 

hetero_homo_dict = {}
sample_list = []
seen_sample_scaf = set()
N_scafs_angst_D_no_calc = set()
# scaf	N_homo	N_hetero

path = in_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith(cov_ext):
			if name.startswith(species_want):
				#print(name)
				sample_name = name.split("_to_")[0]
				sample_list.append(sample_name)
				#print(sample_name)
				curr_file = open(os.path.join(in_dir_name, name))
				line_N = 0
				for line in curr_file:
					line_N = line_N + 1
					line = line.strip()
					if line_N > 1:
						line = line.replace(" ", "\t").split("\t")
						if len(line) <3: # blank lines = NA
							hetero_homo_dict[sample_name + "_" + line[0]] = "NA,NA"
							N_scafs_angst_D_no_calc.add(line[0])
						else:
							hetero_homo_dict[sample_name + "_" + line[0]] = line[1] + "," + line[2]
							if sample_name + "_" + line[0] not in seen_sample_scaf:
								seen_sample_scaf.add(sample_name + "_" + line[0])
							else:
								print("ERROR, EXITING")
								sys.exit(2)
								


sample_list_s = sorted(sample_list)
print("\nSamples:")
print(sample_list_s)

#######################################################################################################
### add to coverage info


coverage_with_LG_file = open(coverage_with_LG_filename)

#  /Users/dparker/Documents/University/Lausanne/Sex_chromosomes/v8_cov_contig/Tbi_pairedcov_minlen=1000_contig_cov_wLGinfo.csv

missing_scafs = set()
total_scafs = set()
out_file = open(coverage_with_LG_filename.replace(".csv", "_withangsD.csv"), "w")


line_N = 0
for line in coverage_with_LG_file:
	line = line.strip()
	line_N = line_N + 1
	if line_N == 1:
		#print(line)
		header = line
		for el in sample_list:
			header_bit = "homo_mlest_" + el + ",hetero_mlest_" + el
			header = header + "," + header_bit
		out_file.write(header + "\n")
	else:
		#print(line)
		scaf_name = line.split(",")[0]
		out_line = line
		for el in sample_list:
			ss = el + "_" + scaf_name
			het_homo_info = hetero_homo_dict.get(ss)
			if het_homo_info == None: # possible if giving contigs <1000
				 het_homo_info = "NA,NA"
				 missing_scafs.add(scaf_name)
			else:
				total_scafs.add(scaf_name)
			
			out_line = out_line + "," + het_homo_info
		#print(out_line)		
		out_file.write(out_line + "\n")
		
print("\nTotal number of scafs: " + str(len(total_scafs)))
print("Total number where ansgD did not calc ests for ANY sites for at least 1 sample: " + str(len(N_scafs_angst_D_no_calc)))
print("Number of scafs with missing angsD ests (i.e. in cov file but not angsD file): " + str(len(missing_scafs)))


print("\n\nFinished, Max\n\n\n")