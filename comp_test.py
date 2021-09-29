"""
comp_test.py

Matt Rich, 08/2021
"""

def main(infiles, _filter):
		gene_mut_dict = {}

		_filter = set(_filter.split(","))

		for f in infiles:
			#remove all the excess text from sample names
			#assumes sample does not contain "."
			sample = f.split("/")[-1].split(".")[0]
		
			for line in open(f, 'r'):
				#skip header lines
				if line.startswith("#"): 
					continue
				l_strip = line.strip()
				l = l_strip.split("\t")
				
				#filter for only mutations with 
				#moderate or high predicted effects
				if l[13] in _filter or _filter == [""]:
					
					effect = l[21]
					if effect == "":
						effect = l[20]
	
					if l[14] in gene_mut_dict:
						gene_mut_dict[l[14]].append([sample, "{}:{}({})".format(sample, l[17], effect), l_strip])
					else:
						gene_mut_dict[l[14]] = [[sample, "{}:{}({})".format(sample, l[17], effect), l_strip]]
	
		#print output
		for g in gene_mut_dict:
				sample_effects = [x[1] for x in gene_mut_dict[g]]
#				sample_effects = ["{}:{}({})".format(x[0], g, x[1]) for x in gene_mut_dict[g]]
				samples = set(list(zip(*gene_mut_dict[g]))[0])
				print("\t".join([
						g, 
						str(len(samples)),
						",".join(sorted(list(set(sample_effects))))]
						))


if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('FILES', nargs='*',
		help = "files for comparing. output of snpeff2tsv.py.")
	parser.add_argument("-f", "--file_list", action='store', dest="file_list",
		help = "file containing list of files to compare, one per line")
	parser.add_argument('--filter', action='store', dest="FILTER", 
		help="comma-delimited list of filters to apply", default="")
	args = parser.parse_args()
	
	if args.file_list != None:
		#open the list of files
		list_of_files = []
		for line in open(args.file_list, "r"):
			list_of_files.append(line.strip())
#		print(len(list_of_files))
		main(list_of_files, args.FILTER)
	else:	
		main(args.FILES, args.FILTER)	

