"""
intersect_controls.py

wrapper to call bcftools isec, given input of a tab-delimited file containing a
SAMPLE and CONTROLS to isec against.

Matt Rich, 8/2021
"""

import os

def main(infile, prefix):
		for line in open(infile, "r"):
				out_l = "bcftools isec -p {} -C ".format(prefix)

				l = line.strip().split("\t")
				
				#call bcftools isec
				out_l += " ".join(l)			

				print(out_l)
				os.system(out_l)

				out_name_parts = l[0].split("/")[1].split(".")
				isec_name = ".".join([out_name_parts[0], 
							"isec", 
							".".join(out_name_parts[1:])])
				#rename the output (0000.vcf) to something meaningful
				os.system("mv {0}/0000.vcf {0}/{1}".format(prefix, isec_name.rstrip(".gz")))
		#then delete the README and sites files
		os.system("rm {}/README {}/sites.txt")

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument("INPUT_FILE", 
		help="tab-delimited file containing a SAMPLE followed by CONTROLS \
						to intersect against")
	parser.add_argument('-p', '--prefix', action = 'store', type = str, 
		dest = 'PREFIX', default="ISEC_OUTPUT", 
		help = "prefix to store output (i.e., output folder)")
	args = parser.parse_args()
	
	main(args.INPUT_FILE, args.PREFIX)	

