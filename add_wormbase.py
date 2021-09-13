"""
add_wormbase.py

Matt Rich, 08/2021
"""

def main(f, wb):
	
		#read wormbase annotations
		wb_dict = {}
		for line in open(wb, "r"):
				l = line.strip().split("\t")
				wb_dict[l[1]] = l	
				wb_dict[l[2]] = l	#add both systematic and name for each gene
				ann_length = len(l)

		#read file to be annotated
		for line in open(f, "r"):
				l = line.strip().split("\t")
				new_line = l
				if l[0] in wb_dict:
						new_line += wb_dict[l[0]][3:]
				else:
						new_line += (["n/a"] * (ann_length-3))
				print("\t".join(new_line))


if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--wb', action = 'store', type = str, dest = 'wormbase', 
		help = "tab-delimited file of annotations from wormbase")
	parser.add_argument('--in', action = 'store', type = str, dest = 'infile', 
		help = "file to be annotated; gene name must be in first column")
	args = parser.parse_args()
	
	main(args.infile, args.wormbase)	

