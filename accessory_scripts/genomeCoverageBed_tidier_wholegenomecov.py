# genomeCoverageBed_tidier_wholegenomecov.py

import sys
import os
import getopt

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

infile_name = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** genomeCoverageBed_tidier_wholegenomecov.py | Written by DJP, 13/11/19 in Python 3.5 in Lausanne ****\n")
		print("This program takes the default histogram output from genomeCoverageBed (e.g. genomeCoverageBed -ibam mybam.bam -g myref.fa > mybam_coverage.out) ") 
		print("It outputs just the coverage histogram for the 'genome', which can then be plotted with plot_genome_cov.R")
		print("\n**** USAGE **** \n")
		print("genomeCoverageBed_tidier_wholegenomecov.py -i [name of coverage file (e.g. mybam_coverage.out)]\n")	
		sys.exit(2)
		
	elif opt in ('-i'):
		infile_name = arg
	else:
		print("i dont know")
		sys.exit(2)

if infile_name == None:
	print("Please specify a input file! For more info see help with option -h")
	sys.exit(2)

outfile_name = infile_name.replace(".out", "") + "_genomecov.txt"
outfile = open(outfile_name, "w")

outfile.write("feature\tdepth\tNsites\ttotalNsites\tpropsites\n")
infile = open(infile_name)
for line in infile:
	line = line.strip()
	feature = line.split("\t")[0]
	if feature == "genome":
		outfile.write(line + "\n")
		
		
print("\n\nDone, Parzival\n\n")
		
		

		

		

	
