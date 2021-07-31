"""
snpeff2tsv.py

this script takes a snpeff annotated VCF file and
expands each effect annotation to new lines (rather 
than having all potential isoform effects in a single 
line). Based on old script snpeff2human.m from Eddie
Hubjer. 

Matt Rich, 7/2021
"""

def main(vcf, append_info):
		
		ANN_HEADERS=[]
		LOF_HEADERS=[]
		NMD_HEADERS=[]
		new_header = []
		specific_headers = ["DP", "MQ", "AC", "AN"]
	
		for line in open(vcf, "r").readlines():	
			#grab snpeff headers from comments
			if line.startswith("##INFO"):
					comment_info = line.strip().split("INFO=<ID=")[1].split(",")[0] 
					if comment_info == "ANN":
							ANN_HEADERS = [ "ANN_" + "".join(x.strip().split()) for x in line.strip().split("'")[1].split("|") ]
					if comment_info == "LOF":
							LOF_HEADERS = [ "LOF_" + "".join(x.strip().split()) for x in line.strip().split("'")[1].split("|") ]
					if comment_info == "NMD":
							NMD_HEADERS = [ "NMD_" + "".join(x.strip().split()) for x in line.strip().split("'")[1].split("|") ]
					
			#grab the header line and add snpeff headers to it
			elif line.startswith("#CHROM"):
					new_header = line.strip().split("\t")[:7] + \
								specific_headers + \
								ANN_HEADERS + \
								LOF_HEADERS + \
								NMD_HEADERS + \
								line.strip().split("\t")[-2:]
					if append_info:
							new_header += ["INFO"]
					print("\t".join(new_header))
		
			#otherwise just print the comment line
			elif line[:2] == "##":
					print(line.strip())


			#after dealing with headers, we can reformat each mutation 
			#and its annotations
			else:
					l = line.strip().split("\t")
					#start building new line
					new_variant_line = l[:7]
					line_suffix = l[-2:]
					
					#parse INFO line (which contains EFF info)
					info = {}
					_info = l[7].split(";")
					if _info[0] == "INDEL":
							info = {x.split("=")[0]:x.split("=")[1] for x in _info[1:]}
					else:
							info = {x.split("=")[0]:x.split("=")[1] for x in _info}
			
					#add DP, MQ, etc... to the line
					#this line is now the same for all annotations of this
					#variant, including the FORMAT, SAMPLE, (and INFO)
					for h in specific_headers:
							new_variant_line.append(info[h])

					#how many annotations (~# of isoforms)
					for x in info["ANN"].split(","):
							ANNdat = x.split("|")
							if "LOF" in info:
									ANNdat += info["LOF"].strip("()").split("|")
							else:
									ANNdat += ["","","",""]
							if "NMD" in info:
									ANNdat += info["LOF"].strip("()").split("|")
							else:
									ANNdat += ["","","",""]
							#finally, print line
							#adding FORMAT and SAMPLE fields to end
							#and INFO, if append_info
							if append_info:
									print("\t".join(new_variant_line + ANNdat +
											line_suffix + l[7])) 
							else:
									print("\t".join(new_variant_line + ANNdat +
											line_suffix))
								
if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--vcf', action = 'store', type = str, dest = 'vcf', 
		help = "snpeff annotated vcf")
	parser.add_argument('--info', action = 'store_true', dest = 'append_info',
		help = "append full info for each variant at end of line?",
		default = False)
	args = parser.parse_args()
	
	main(args.vcf, args.append_info)	
