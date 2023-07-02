"""
intersect_controls.py

wrapper to call bcftools isec, given input of a tab-delimited file containing a
SAMPLE and CONTROLS to isec against.

Matt Rich, 8/2021
"""

import os
from os.path import exists

def main(infile, prefix, isec_all, skip_missing):
		isec_commands = []
		out_l_template = "bcftools isec -p {} -C ".format(prefix) + "{}"

		if isec_all:
				#compile all the unique files to compare
				all_files = set()
				for x in open(infile, "r").readlines():
						for y in x.strip().split("\t"):
								all_files.add(y)
				all_files = list(all_files)
				
				#check if the files exist:
				missing_files = []
				for f in all_files:
						if not exists(f):
								missing_files.append(f)
				#if we're skipping the missing files, remove them from list and
				#continue forward.
				if len(missing_files) != 0:
						if skip_missing:
								print("removing missing files from commands")
								for mf in missing_files:
										all_files.remove(mf)
						else:
								mf_str = "\n\t".join(missing_files)
								raise OSError("The following files are missing. To perform intersection with these files, use --skip-missing.\n{}".format(mf_str))

				#then write all the isec commands
				for i in range(len(all_files)):
						f_list = " ".join([all_files[i]] + all_files[0:i] + all_files[i+1:])
						isec_commands.append(out_l_template.format(f_list))
				#print(isec_commands)
				
		else:
			for line in open(infile, "r"):
				out_l = "bcftools isec -p {} -C ".format(prefix)

				l = line.strip().split("\t")
				
				#call bcftools isec
				isec_commands.append(out_l_template.format(" ".join(l)))		

				
		#after compiling all isec commands, run them		
		for out_l in isec_commands:
				print(out_l.split()[5])
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
		os.system("rm {}/README.txt {}/sites.txt".format(prefix, prefix))

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
	parser.add_argument('--skip-missing', action = 'store_true', 
		dest = "SKIPMISSING", default = False,
		help = "skip missing files instead of raising error.")
	
	args = parser.parse_args()
	
	main(args.INPUT_FILE, args.PREFIX, args.ISEC_ALL, args.SKIPMISSING)	

