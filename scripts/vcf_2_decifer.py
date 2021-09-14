#!/usr/bin/python -tt

"""
Please type "python3 vcf_2_decifer.py --help" for a the full list of options

This is a python script that takes as input (1) a multi-sample VCF, e.g. variants called in multiple tumor samples from a single patient, (2) a file containing copy number aberration information, e.g. the "best.seg.ucn" file that is output by the program HATCHet. The sample names must agree in both files!

This script also uses cyvcf2 to parse VCF files efficiently, and the python implementation of bedtools to find which CNA interval each SNV overlaps with. Please make a conda environment containing these python modules before running this script.
For instance, do the following to run this script:

conda create -n vcf_bedtools pybedtools cyvcf2 pandas -y
conda activate vcf_bedtools
python vcf_2_decifer.py [OPTIONS]
"""

import re
import sys
import pybedtools as pbt
from cyvcf2 import VCF
import os
import pandas as pd
import numpy as np
from collections import defaultdict
import argparse


def filterByDepth(gt_depths, gt_alt_depths, Filter):
    PASS = 1
    missing = 0
    for i in range(len(gt_depths)):
        # filter if genotype has low depth or is missing
        if (gt_depths[i] < Filter['MinDepth']): 
            missing += 1
        # filter if alt allele isn't greater than the specified threshold in at least one sample
        if not any( np.greater_equal(gt_alt_depths, Filter['MinDepthAltAllele']) ):
            missing += 1
    #(gt_alt_depths[i] < Filter['MinDepthAltAllele'])
    if missing > 0: 
        PASS = 0
    return(PASS)

def compute_ref_var_depths(vcf, FilterDP):
    ref_var_depths = defaultdict(list) # ref_var_depths[char_label] = list of (ref,alt) tuples, one for each sample, in same order as vcf.samples
    for variant in vcf:
        if len(variant.ALT) == 1 and variant.var_type == "snp":
            PASS = filterByDepth(variant.gt_depths, variant.gt_alt_depths, FilterDP)
            #print(np.greater_equal(variant.gt_alt_depths,FilterDP['MinDepthAltAllele']))
            if PASS:
                chrom = variant.CHROM
                pos = variant.end
                char_label = ".".join(map(str,[chrom, pos, variant.REF, variant.ALT[0]]))

                for i in range(len(vcf.samples)):
                    ref = variant.gt_depths[i] - variant.gt_alt_depths[i]
                    alt = variant.gt_alt_depths[i]
                    ref_var_depths[char_label].append( (ref, alt) )
    return(ref_var_depths)

def print_output(vcf, ref_var_depths, cna_overlaps, outdir):
    char_index = 0
    chars = ref_var_depths.keys() & cna_overlaps.keys() 
    header = [str(len(chars)) + " #characters"]
    header.append( str(len(vcf.samples)) + " #samples" )
    header.append( "#sample_index\tsample_label\tcharacter_index\tcharacter_label\tref\tvar" )
    print(header)
    with open(f"{outdir}/decifer.input.tsv",'w') as out:
        print("\n".join(header), file=out)
        for char_label in ref_var_depths:
            if char_label in cna_overlaps:
                for i in range(len(vcf.samples)):
                    r, v = ref_var_depths[char_label][i][0], ref_var_depths[char_label][i][1]
                    to_print = [i, vcf.samples[i], char_index, char_label, r, v]
                    cnas = cna_overlaps[char_label][i]
                    to_print.extend(cnas)
                    print("\t".join(map(str, to_print)), file=out)
                    #print(i, vcf.samples[i], char_index, char_label, r, v)
                char_index += 1

def print_purities(cna_df, sample_index, num_samples, outdir):
    purities = {}
    for i, row in cna_df.head(num_samples+1).iterrows():
        purities[row['SAMPLE']] = 1.0 - row['u_normal'] 
    with open(f"{outdir}/decifer.purity.tsv",'w') as out:
        for sample in sample_index:
            print(sample_index[sample], purities[sample], file=out, sep="\t")

def filter_high_CN_sites(cn_states_persite, max_CN):
    # returns 0 if site has CN state greater than max_CN
    # cn_states_persite is a list of 2-tuples, CN states of maternal/paternal chromosomes
    for i in cn_states_persite:
        if int(i[0]) + int(i[1]) > max_CN:
            return 0 
    return 1

def print_unique_CN_states(cn_states, max_CN, outdir):
    # print unique copy number states for sites that are below the max_CN threshold 
    cn_states = tuple(set(cn_states))
    with open(f"{outdir}/cn_states.txt", 'w') as out:
        for value in cn_states:
            PASS=1
            for i in value:
                if int(i[0]) + int(i[1]) > max_CN: 
                    PASS=0
            if PASS: 
                out.write(';'.join([','.join(i) for i in value]) + '\n')

def print_filtered_sites(filtered_sites, cna_overlaps, outdir):
    with open(f"{outdir}/filtered_sites.txt",'w') as out:
        out.write("\n".join(filtered_sites))
        print(file=out)
    with open(f"{outdir}/filtered_stats.txt",'w') as out:
        filtered = len(filtered_sites)
        total = len(cna_overlaps.keys())
        print("# sites that were filtered due to copy-number states > max_CN", file=out)
        print("filtered: ", filtered, file=out)
        print("fraction: ", float(filtered/total), file=out)

def overlap_cna_snp(vcf_samples, max_CN, out_dir):
    cna_overlaps = defaultdict(list)
    cn_states_allsites = [] # a list of tuples
    filtered_sites = set()  # sites filtered out because of high CN
    snps = pbt.BedTool(f"{out_dir}/snps.bed")
    # for each sample in VCF, intersect it's SNVs with sample-specific CNAs
    for sample in vcf_samples:
        sample_cnas = pbt.BedTool(f"{out_dir}/{sample}_cna.bed")    
        #snps.intersect(sample_cnas, wo=True).saveas(f"snps_cnas_overlap_{sample}.bed")
        bed = snps.intersect(sample_cnas, wo=True)
        for line in bed:
            line = str(line).split()
            # first 5 columns are SNP info (chr,pos_start,pos_end,REF,ALT), rest are CNA info
            # CNA info starts with CHR, START, END, which you don't want, so start at index 8
            # exclude the last index, since bedtools adds this, the number of bp of overlap
            char_label = ".".join([line[0], line[2], line[3], line[4]])
            cna_info = line[8:-1]
            cna_info_parsed = []
            cn_states_persite = []
            for i in cna_info:
                if "|" in i:
                    # this is an allele-specific CN state, int CN|int CN
                    cn = i.split("|")
                    cna_info_parsed.extend(cn)
                    # get all CN states, we'll filter out duplicates later
                    cn_states_persite.append( tuple(cn) )
                else:
                    cna_info_parsed.append(i)

            if filter_high_CN_sites(cn_states_persite, max_CN):
                cna_overlaps[char_label].append( tuple(cna_info_parsed) )
            else:
                filtered_sites.add( char_label )
            # get collection of unique CN states for this SNV site
            cn_states_allsites.append( tuple(set(cn_states_persite)) )


    return cna_overlaps, cn_states_allsites, filtered_sites

def main():

    parser = argparse.ArgumentParser(description='Generate input for Decifer using VCF file and HATCHet CNA file')
    parser.add_argument("-V","--vcf_file", required=True, type=str, help="single or multi-sample VCF file")
    parser.add_argument("-C","--cna_file", required=True, type=str, help="HATCHet CNA file: best.seg.ucn ")
    parser.add_argument("-O","--out_dir", required=True, default="./", type=str, help="directory for printing files; please make unique for each patient!")
    parser.add_argument("-M","--min_depth", required=True, type=int, help="minimum depth PER sample")
    parser.add_argument("-A","--min_alt_depth", required=True, type=int, help="minimum depth of ALT allele in at least one sample")
    parser.add_argument("-N","--max_CN", required=False, default=6, type=int, help="maximum total copy number for each observed clone")
    args = parser.parse_args()

    vcf_name = os.path.basename(args.vcf_file)
    vcf = VCF(args.vcf_file, gts012=True)
    
    # Filtering criteria
    FilterDP = {}
    FilterDP['MinDepth'] = args.min_depth
    FilterDP['MinDepthAltAllele'] = args.min_alt_depth
    
    num_samples = len(vcf.samples)
    sample_index = { vcf.samples[i] : i for i in range(len(vcf.samples)) }
    print(vcf.samples)
    print(type(vcf.samples))

    # ref_var_depths[char_label] = list of (ref,alt) tuples, one for each sample, in same order as vcf.samples
    ref_var_depths = compute_ref_var_depths(vcf, FilterDP)

    # print BED file for SNPs
    with open(f"{args.out_dir}/snps.bed", 'w') as out:
        print("chrom\tstart\tend\tREF\tALT", file=out)
        for chr_label in ref_var_depths:
            pos = chr_label.split(".")
            # subtract 1 from position to create interval in BED format
            print(pos[0], int(pos[1])-1, int(pos[1]), pos[2], pos[3],  sep="\t", file=out)

    # Load in CNA information
    cna_df = pd.read_csv(args.cna_file, sep = '\t', index_col=False)
    # print purity information
    print_purities(cna_df, sample_index, num_samples, args.out_dir)

    # prepare BED files for CNA intervals for each sample, for overlapping with SNPs
    for sample in vcf.samples:
        df = cna_df[cna_df['SAMPLE'] == sample]
        # consider subtracting 1 from start of interval to be compatible with BED format, leave end interval alone
        df.loc[:,'START'] = df['START']
        df = df.drop('SAMPLE', axis=1)
        df.to_csv(f"{args.out_dir}/{sample}_cna.bed", index=False, sep="\t")

    # overlap SNPs with CNA intervals for each sample
    # cna_overlaps[char_label] = list of tuples of CNA info (one tuple for each sample, in same order as vcf.samples)
    # this function also prints the observed CN state trees for the generatestatetrees function
    cna_overlaps, cn_states_allsites, filtered_sites = overlap_cna_snp(vcf.samples, args.max_CN, args.out_dir)
    
    # sites may have unique CN states that are duplicate; set them to find unique CN states across sites
    print_unique_CN_states(cn_states_allsites, args.max_CN, args.out_dir)
    print_filtered_sites(filtered_sites, cna_overlaps, args.out_dir)

    print_output(vcf, ref_var_depths, cna_overlaps, args.out_dir) 

    os.system(f"rm {args.out_dir}/snps.bed")
    for sample in vcf.samples:
        os.system(f"rm {args.out_dir}/{sample}_cna.bed")
    
if __name__ == '__main__':
  main()


