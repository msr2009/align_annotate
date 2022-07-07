"""
intersect_controls.py

wrapper to call bcftools isec, given input of a tab-delimited file containing a
SAMPLE and CONTROLS to isec against.

Matt Rich, 8/2021
"""

import os

def main(infile, prefix, isec_all):
		isec_commands = []
		out_l_template = "bcftools isec -p {} -C ".format(prefix) + "{}"

		if isec_all:
				#compile all the unique files to compare
				all_files = set()
				for x in open(infile, "r").readlines():
						for y in x.split("\t"):
								all_files.add(y)
				all_files = list(all_files)
				#then write all the isec commands
				for i in range(len(all_files)):
						f_list = " ".join([all_files[i]] + all_files[0:i] + all_files[i+1:])
						isec_commands.append(out_l_template.format(f_list))
				print(isec_commands)
				
		else:
			for line in open(infile, "r"):
				out_l = "bcftools isec -p {} -C ".format(prefix)

				l = line.strip().split("\t")
				
				#call bcftools isec
				isec_commands.append(out_l_template.format(" ".join(l)))		

				
		#after compiling all isec commands, run them		
		for out_l in isec_commands:
				os.system(out_l)

				out_name_parts = out_l.split(" ")[5].split(".")
				isec_name = ".".join([out_name_parts[0], 
							"isec", 
							".".join(out_name_parts[1:])])
				#remove any upstream folder names from isec_name
				try:
						isec_name = isec_name.split("/")[1]
				except IndexError:
						pass

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
	parser.add_argument('--all', action = 'store_true', 
		dest = "ISEC_ALL", default = False,
		help = "intersect each file against all files in input")
	
	args = parser.parse_args()
	
	main(args.INPUT_FILE, args.PREFIX, args.ISEC_ALL)	

