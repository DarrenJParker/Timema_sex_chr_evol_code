# genomeCoverageBed_tidier_select_scafs.py

import sys
import os
import getopt

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:w:e:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

infile_name = None
extra_outfile_name = ""
#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** genomeCoverageBed_tidier_select_scafs.py | Written by DJP, 21/11/19 in Python 3.5 in Lausanne ****\n")
		print("This program takes the default histogram output from genomeCoverageBed (e.g. genomeCoverageBed -ibam mybam.bam -g myref.fa > mybam_coverage.out) ") 
		print("It outputs just the coverage histogram for specified contigs, which can then be plotted with plot_genome_cov.R")
		print("\n**** USAGE **** \n")
		print("genomeCoverageBed_tidier_select_scafs.py -i [name of coverage file (e.g. mybam_coverage.out)] -w [file of wanted scaf names]\n")	
		sys.exit(2)
		
	elif opt in ('-i'):
		infile_name = arg
	elif opt in ('-w'):
		wantfile_name = arg
	elif opt in ('-e'):
		extra_outfile_name = arg
	else:
		print("i dont know")
		sys.exit(2)

if infile_name == None:
	print("Please specify a input file! For more info see help with option -h")
	sys.exit(2)


want_scafs = set()

wantfile = open(wantfile_name)
for line in wantfile:
	line = line.strip().strip(">")
	want_scafs.add(line)
	
print("Number of wanted scafs:" + str(len(want_scafs)))


outfile_name = infile_name.replace(".out", "") + extra_outfile_name + "cov.txt"
outfile = open(outfile_name, "w")

found_want_scafs = set()

cov_dict = {}
all_depths = set()

seen_want_scaf = set()
total_len_wanted = 0

outfile.write("feature\tdepth\tNsites\ttotalNsites\n")
infile = open(infile_name)
for line in infile:
	line = line.strip()
	feature = line.split("\t")[0]
	depth   = int(line.split("\t")[1])
	Nsites  = int(line.split("\t")[2])
	T_sites = int(line.split("\t")[3])

	if feature in want_scafs:
		if feature not in seen_want_scaf:
			seen_want_scaf.add(feature)
			total_len_wanted = total_len_wanted + T_sites
		
		found_want_scafs.add(feature)
		rec = cov_dict.get(depth)
		all_depths.add(depth)
		if rec == None:
			cov_dict[depth] = Nsites
		else:
			cov_dict[depth] = Nsites + rec
		

cov_depth_ord = sorted(list(all_depths))
for d in cov_depth_ord:
	rec = cov_dict.get(d)
	outfile.write("sel_scaf\t" + str(d) + "\t" + str(rec) + "\t" + str(total_len_wanted) + "\n")


		
print("Number of wanted scafs found in cov file:" + str(len(found_want_scafs)))
print("Total size of wanted contigs:" + str(total_len_wanted))

print("\n\nDone, Parzival\n\n")
		
		

		

		

	
