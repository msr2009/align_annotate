"""
make_background_vcf.py

script that takes list of background files as input, 
merges these into single multisample vcf, then uses pysam
to calculate the background allele frequency of every variant.

Matt Rich, 6/2023
"""

import pysam, subprocess, itertools

def calculate_bgaf(rec):
	#number of samples
	n_samples = float(len(rec.samples.values()))
	#get genotypes and calculate allele frequency
	gts = list(itertools.chain.from_iterable([s['GT'] for s in rec.samples.values()]))
	bgaf = sum(filter(None, gts))/(2*n_samples)
	return bgaf

def main(file_list, out_name):
	pwd = subprocess.run("dirname {}".format(out_name), shell=True, capture_output=True).stdout.decode().strip()
	#do merge
	subprocess.run('bcftools merge -m none -l {} | bcftools norm -m- -Oz -o {}/tmp.vcf.gz'.format(file_list, pwd), shell=True)
	subprocess.run("bcftools index {}/tmp.vcf.gz".format(pwd), shell=True)

	merged_vcf = pysam.VariantFile("{}/tmp.vcf.gz".format(pwd), 'r')
	#make new header for output
	out_header = merged_vcf.header
	#throws a ValueError if the info line already exists in header, so
	try:
		out_header.info.add('BGAF', '1', 'Float', 'Background allele frequency')
	except ValueError:
		pass
	#make output
	out_vcf = pysam.VariantFile(out_name, "w", header=out_header)

	#loop through merged vcf and calculate bgaf
	for record in merged_vcf.fetch():
		record.info["BGAF"] = calculate_bgaf(record)
		out_vcf.write(record)
	
	out_vcf.close()
	
	#index new merged/bgaf vcf
	subprocess.run("bcftools index {}".format(out_name), shell=True)
	subprocess.run("rm {}/tmp.vcf.gz".format(pwd), shell=True)
	subprocess.run("rm {}/tmp.vcf.gz.csi".format(pwd), shell=True)

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-l', '--filelist', action = 'store', type = str, dest = 'filelist', 
		help = "list of vcfs to merge. one per line with full paths.")
	parser.add_argument('-o', '--outvcf', action = 'store', type = str, dest = 'outvcf',
		help = "name for output vcf")
	args = parser.parse_args()
	
	main(args.filelist, args.outvcf)	

