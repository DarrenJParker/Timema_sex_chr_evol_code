### class_genes_to_lg.py

import sys
import os
import getopt
import decimal
from decimal import *

try:
	opts, args = getopt.getopt(sys.argv[1:], 'l:c:o:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

linkage_file_name = None
cov_filename = None
out_base = "TEST"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** class_genes_to_lg.py | Written by DJP, 28/11/18 in Python 3.5 in Lausanne  ****\n")
		print("Calcs rel position of genes based on Nosil's linkage map after correction by Kamil. See sex_chr_scripts/2_read_mapping_to_v8_genomes.sh and sex_chr_scripts/4_Sex_chr_expression_analyses.sh")
		
		print("\n***** USAGE *****\n")		
		print("\npython3 class_genes_to_lg.py -l [linkage alignment block file] -g [gff file] -o [output basename]\n\n")
		
		sys.exit(2)
	elif opt in ('-l'):
		linkage_file_name = arg
	elif opt in ('-c'):
		cov_filename = arg
	elif opt in ('-o'):
		out_base = arg
	else:
		print("i dont know")
		sys.exit(2)
		
		

##############################################

linkage_file = open(linkage_file_name)

linkage_info_dict = {}

seen_scaf = set()
N_Kamil_filt_lines = 0
line_N = 0
for line in linkage_file:
	line = line.strip().split("\t")
	line_N = line_N + 1
	if line_N > 1:
		UNIL_scaf = line[0]
		UNIL_coord_start = int(line[3])
		UNIL_coord_end   = int(line[4])
		UNIL_midpoint = int((UNIL_coord_start + UNIL_coord_end) / 2)
		
		
		nosil_coord_start = int(line[9])
		nosil_coord_end   = int(line[10])
		orientation = None
		
		### X chr is not ordered except within scafs - SO make multiple X linkage groups
		lg = line[8]
		if lg == "lgX":
			lg = line[5]
		
		
		
		if UNIL_coord_start > UNIL_coord_end:
			orientation = "reverse"
		else:
			orientation = "forward"
		
		### kamil extra filter
		
		block_size = int(line[2])
		block_r_start = int(line[6])
		block_r_end  =  int(line[7])
		
		if abs(block_r_end - block_r_start) / block_size > 0.5:
		
		
			### add to dict
			if UNIL_scaf not in seen_scaf:
				seen_scaf.add(UNIL_scaf)
				linkage_info_dict[UNIL_scaf] = [(lg, orientation, nosil_coord_start, nosil_coord_end, UNIL_midpoint)]
			else:
				rec = linkage_info_dict.get(UNIL_scaf)
				rec.append((lg, orientation, nosil_coord_start, nosil_coord_end,UNIL_midpoint ))
				linkage_info_dict[UNIL_scaf] = rec
		else:
			N_Kamil_filt_lines = N_Kamil_filt_lines + 1
		

print(N_Kamil_filt_lines)

### select largest alignment block. most have one large block and many little ones.
### should ask kamil about this
### If there was > 1 linkage group I make a note


scaf_to_linkage_info_dict = {}
	

N_single_linkage = 0
N_multi_linkage  = 0

for scaf in linkage_info_dict:
	rec = linkage_info_dict.get(scaf)
	#print(rec)
	lg_set = set()
	lg_set_all = ""

	for i in rec:
		lg_set.add(i[0])
		
	for el in lg_set:
		lg_set_all = lg_set_all + ";" + el
	lg_set_all = lg_set_all.strip(";")

	multi_linkage_groups = None
	if len(lg_set) == 1:
		N_single_linkage = N_single_linkage + 1
		multi_linkage_groups = "NO"
	else:
		N_multi_linkage  = N_multi_linkage + 1
		multi_linkage_groups = "YES"
		
	start_pos = []
	end_pos = []
	curr_block_len = -1
	start_biggest_block_len = None
	end_biggest_block_len   = None
	lg_biggest_block        = None
	orie_biggest_block      = None
	UNIL_mid_biggest_block  = None
	
	for i in rec:
		block_len = i[3] - i[2]
		if curr_block_len < block_len:
			curr_block_len = block_len
			start_biggest_block_len = i[2]
			end_biggest_block_len   = i[3]
			lg_biggest_block   = i[0]
			orie_biggest_block  = i[1]
			UNIL_mid_biggest_block = i[4]
	
	mid_point_biggest_block_len = int((end_biggest_block_len + start_biggest_block_len) / 2)
	
	#print(mid_point_biggest_block_len)	

	scaf_to_linkage_info_dict[scaf] = [mid_point_biggest_block_len,lg_biggest_block,orie_biggest_block,multi_linkage_groups,UNIL_mid_biggest_block,lg_set_all]	


print("\nNumber of scafs assigned to a single linkage group:  " + str(N_single_linkage))	
print("Number of scafs assigned to multiple linkage groups: " + str(N_multi_linkage))	

#print(scaf_to_linkage_info_dict)



#############################################################################################
### Add to coverage 

out_file = open(out_base + "_wLGinfopos.csv", "w")

cov_file = open(cov_filename)

line_N = 0
for line in cov_file:
	line_N = line_N + 1
	line = line.strip()
	if line_N == 1:
		out_file.write(line + ",mid_point_biggest_block_len,lg_biggest_block,orie_biggest_block,multi_linkage_groups,UNIL_mid_biggest_block,lg_set_all\n")
		#print(line)
	else:
		#print(line)
		scaf_name = line.split(",")[0]
		lg_rec = scaf_to_linkage_info_dict.get(scaf_name)
		lg_out = ""
		if lg_rec == None:
			lg_out = ",NA,NA,NA,NA,NA,NA"
		else:
			for el in lg_rec:
				lg_out = lg_out + "," + str(el)
		
		out_file.write(line + lg_out + "\n")
		
		#print(lg_out)



	
print("\nFinished, Billy Rock\n\n\n")
	
	
	

