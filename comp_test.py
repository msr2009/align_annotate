"""
comp_test.py

Matt Rich, 08/2021
"""

def main(infiles, _filter):
		gene_mut_dict = {}

		_filter = set(_filter.split(","))
		
		dat_cols = {}
	
		for f in infiles:
			#remove all the excess text from sample names
			#assumes sample does not contain "."
			sample = f.split("/")[-1].split(".")[0]
			
			for line in open(f, 'r'):
				#skip header lines
				#unless it's the header, in which case
				#id the columns we need
				if line.startswith("#CHROM"):
					header = line.strip().split("\t")
					dat_cols = {header[i]: i for i in range(len(header)) }
				elif line.startswith("#"): 
					continue

				l_strip = line.strip()
				l = l_strip.split("\t")

				#filter for only mutations with 
				#moderate or high predicted effects
				if l[dat_cols["ANN_Annotation_Impact"]] in _filter or _filter == [""]:	
					
					effect = ""
					if l[dat_cols["ANN_HGVS.p"]] != "":
						effect = l[dat_cols["ANN_HGVS.p"]]
					else:
						effect = l[dat_cols["ANN_HGVS.c"]]
	
					if l[dat_cols["ANN_Gene_Name"]] in gene_mut_dict:
						gene_mut_dict[l[dat_cols["ANN_Gene_Name"]]].append([sample, "{}:{}({})".format(sample, l[dat_cols["ANN_Feature_ID"]], effect), l_strip])
					else:
						gene_mut_dict[l[dat_cols["ANN_Gene_Name"]]] = [[sample, "{}:{}({})".format(sample, l[dat_cols["ANN_Feature_ID"]], effect), l_strip]]
	
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
		main(list_of_files, args.FILTER)
	else:	
		main(args.FILES, args.FILTER)	

