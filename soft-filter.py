"""
soft-filter.py

script for performing soft-filtration of reads after base calling
uses pysam to identify heterozygous and low quality sites. can 
optionally, with a properly analyzed multisample vcf of other sequences
filter based on background allele frequency.

to add new filters, edit build_filter_sets and add new arguments (I've 
marked these places with "###ADD NEW FILTERS"

Matt Rich, 6/23
"""

import pysam, subprocess, sys


def build_filter_sets(rec, genotype, lowqual, background):
	f_dict = {"Het": False,
					"lowQual": False,
					"background": 0 }
	f_string = []
	#have to do background analysis as try-except, bc
	#unique variants in new vcf dont' get merged with 
	#bgaf field
	#also, this is optional
	if background != None:
		try:
			if rec.info["BGAF"] > background:
				f_dict["background"] = rec.info["BGAF"]
				f_string.append("background")
		except KeyError:
			pass

	#heterozygous sites
	if genotype[0] != genotype[1]: 
		f_dict["Het"] = 1
		f_string.append("Het")
    
	#low quality sites
	if [s['GQ'] for s in rec.samples.values()][-1] < 15:
		f_dict["lowQual"] = 1
		f_string.append("lowQual")
	
	###ADD NEW FILTERS

	#after looping through filters, if no flags then variant PASSes
	if len(f_string) == 0:
		f_string.append("PASS")
    
	#add record's filter string to list
	return f_string, f_dict

def main(vcf, bgvcf, bgaf, lowqual):
	
	#get working directory first
	pwd = subprocess.run("dirname {}".format(vcf), shell=True, capture_output=True).stdout.decode().strip()
	#merge vcf to-be-filtered to bgvcf (if bgvcf exists) with bcftools 
	if bgvcf != None:
		print("merging background vcf ({}) with sample vcf ({})".format(bgvcf, vcf))
		merge_call = subprocess.run("bcftools merge -m none -o {}/tmp.vcf.gz --force-samples -Oz {} {}".format(pwd, bgvcf, vcf), shell=True)
		print(merge_call)	
	else:
		subprocess.run("cp {} {}/tmp.vcf.gz".format(vcf, pwd), shell=True)
	
	#index tmp vcf
	subprocess.run("bcftools index {}/tmp.vcf.gz".format(pwd), shell=True)

	#read merged vcf
	tmpvcf = pysam.VariantFile("{}/tmp.vcf.gz".format(pwd), "r")
	
	filter_sets = []
	filter_dicts = []
	#calculate the filter string for every variant
	#that is present in the last sample (our sample of interest)
	for record in tmpvcf.fetch():
    	#if our sample (at -1) has a genotype, it isn't reference
		gt = [s['GT'] for s in record.samples.values()][-1] 
		if gt != (None, None):
			rec_filters = build_filter_sets(record, gt, lowqual, bgaf)
			filter_sets.append(rec_filters[0])
			filter_dicts.append(rec_filters[1])

	#now we build a new VCF file, in which we loop through out 
	#to-be-filtered vcf, and replace the FILTER field
	#first we make the new header, which is the same as the original header
	#but with new filters added

	#open vcf 
	in_vcf = pysam.VariantFile(vcf, "r")

	out_header = in_vcf.header

	###ADD NEW FILTERS
	new_filters = [
		["Het", None, None, "Heterozygous variant"],
		["background", None, None, "allele frequency above threshold in background samples"],
		["lowQual", None, None, "variant has low genotype quality (GQ)"],
		["BGAF", None, None, "background allele frequency"]
	]

	#if header for filter already exists, report it but don't worry about it
	for nf in new_filters:
		try:
			out_header.filters.add(nf[0], nf[1], nf[2], nf[3])
		except ValueError:
			print("{} id header already exists.".format(nf[0]), file=sys.stderr)
			
#	out_header.filters.add("lowQual", None, None, "variant has low genotype quality (GQ)")
#	out_header.filters.add("background", None, None, "allele frequency above threshold in background samples")
	out_header.info.add("BGAF", "1", "Float", "Allele frequency in background strains")

	#output name will be programmatically named by stripping ".vcf.gz" and
	#replacing with ".soft-filter.vcf.gz"
	out_name = vcf.rstrip(".vcf.gz")+".soft-filter.vcf.gz"
	out_vcf = pysam.VariantFile(out_name, "w", header=out_header)
	for record in in_vcf.fetch():	
		for f in filter_sets.pop(0):
			record.filter.add(f)
		#add bgaf info to record also
		record.info["BGAF"] = filter_dicts.pop(0)["background"]
		out_vcf.write(record)
	out_vcf.close()

	#index the output vcf
	subprocess.run("bcftools index {}".format(out_name), shell=True)

	#finally, remove the tmp vcf file and its index
	subprocess.run("rm {}/tmp.vcf.gz".format(pwd), shell=True)
	subprocess.run("rm {}/tmp.vcf.gz.csi".format(pwd), shell=True)

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-b', '--bg-vcf', action = 'store', type = str, dest = 'bgvcf', 
			default = None, help = "background multisample VCF with BGAF field")
	parser.add_argument('-v', '--in-vcf', action = 'store', type = str, 
			dest = "input_vcf", help = "single-sample vcf file needing filtration")
	parser.add_argument('--bgaf', action = 'store', type = float, dest = "bgaf", 
			default = 0.1, help = "background allele frequency threshold. must provide --bg-vcf")
	parser.add_argument('--lowqual', action = 'store', type = int, dest = "lowqual", 
			default = 15, help = "genotype quality threshold (GQ)")


	###ADD NEW FILTERS

	args = parser.parse_args()

	main(args.input_vcf, args.bgvcf, args.bgaf, args.lowqual)	

