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


def filterByDepthAndVaf(gt_depths, gt_alt_depths, Filter):
    PASS = 1
    if any( np.less(gt_depths, Filter['MinDepth']) ):
        PASS = 0
    # filter if alt allele isn't greater than the specified threshold in at least one sample
    if not any( np.greater_equal(gt_alt_depths, Filter['MinDepthAltAllele']) ):
        PASS = 0
    # filter if VAF  isn't greater than the specified threshold in at least one sample
    if not any( np.greater_equal((gt_alt_depths/gt_depths), Filter['MinVAF']) ):
        PASS = 0
    return(PASS)

def compute_ref_var_depths(vcf, Filter):
    ref_var_depths = defaultdict(list) # ref_var_depths[char_label] = list of (ref,alt) tuples, one for each sample, in same order as vcf.samples
    for variant in vcf:
        if len(variant.ALT) == 1 and variant.var_type == "snp":
            PASS = filterByDepthAndVaf(np.asarray(variant.gt_depths), np.asarray(variant.gt_alt_depths), Filter)
            #print(np.greater_equal(variant.gt_alt_depths,Filter['MinDepthAltAllele']))
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
    #print(header)
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

def get_purities(cna_df, num_samples, min_purity):
    purities = {}
    for i, row in cna_df.head(num_samples+1).iterrows():
        purity = 1.0 - row['u_normal']
        if purity >= min_purity:
            purities[row['SAMPLE']] = purity
    return purities

def print_purities(purities, sample_index, num_samples, outdir):
    with open(f"{outdir}/decifer.purity.tsv",'w') as out:
        for sample in sample_index:
            print(sample_index[sample], purities[sample], file=out, sep="\t")

def filter_high_CN_sites(cn_states_persite, max_CN):
    # returns 0 if site has CN state greater than max_CN
    # cn_states_persite is a list of str "A|B", CN states of maternal/paternal chromosomes
    for i in cn_states_persite:
        j = i.split("|")
        if int(j[0]) + int(j[1]) > max_CN:
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

def overlap_cna_snp(vcf_samples, max_CN, snps, out_dir):
    cna_overlaps = defaultdict(list)
    cn_states_allsites = [] # a list of tuples
    filtered_sites = set()  # sites filtered out because of high CN
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
            cns = line[8:-1:2] # copy-number states, string[start:end:step]
            props = line[9:-1:2] # copy-number proportions
            if len(cns) != len(props):
               sys.exit("CNA file not formatted correctly!") 

            cn_info = defaultdict(float) # dict with cn_info["A|B"] = proportion
            for c, p in zip(cns, props):
                cn_info[c] += float(p) # this collapses nonunique CN states

            # returns 0 if CN state too high
            if filter_high_CN_sites(cn_info.keys(), max_CN):
                # store results, converting from dict to a list for later printing
                cna_info = []
                [ cna_info.extend( [c.split("|")[0], c.split("|")[1], cn_info[c]] ) for c in cn_info ]
                cna_overlaps[char_label].append( cna_info )
            else:
                filtered_sites.add( char_label )
            # get collection of unique CN states for this SNV site
            tuple_states = [(c.split("|")[0], c.split("|")[1]) for c in cn_info]
            cn_states_allsites.append( tuple(set(tuple_states)) )

    return cna_overlaps, cn_states_allsites, filtered_sites

def main():

    parser = argparse.ArgumentParser(description='Generate input for Decifer using VCF file and HATCHet CNA file')
    parser.add_argument("-V","--vcf_file", required=True, type=str, help="single or multi-sample VCF file")
    parser.add_argument("-C","--cna_file", required=True, type=str, help="HATCHet CNA file: best.seg.ucn ")
    parser.add_argument("-O","--out_dir", required=True, default="./", type=str, help="directory for printing files; please make unique for each patient!")
    parser.add_argument("-M","--min_depth", required=True, type=int, help="minimum depth PER sample")
    parser.add_argument("-A","--min_alt_depth", required=True, type=int, help="minimum depth of ALT allele in at least one sample")
    parser.add_argument("-F","--min_vaf", required=True, type=float, help="minimum VAF of ALT allele in at least one sample")
    parser.add_argument("-N","--max_CN", required=False, default=6, type=int, help="maximum total copy number for each observed clone")
    parser.add_argument("-B","--exclude_list", required=False, default=None, type=str, help="BED file of genomic regions to exclude")
    parser.add_argument("-p","--min_purity", required=False, default=0.0, type=float, help="minimum purity to consider samples")
    parser.add_argument("-S","--snp_file", required=False, default=None, type=str, help="HATCHet file containing germline SNP counts in tumor samples, baf/tumor.1bed")
    args = parser.parse_args()


    # load in vcf file
    vcf_name = os.path.basename(args.vcf_file)
    vcf = VCF(args.vcf_file, gts012=True)
    num_samples = len(vcf.samples)
    
    # Load in CNA information
    cna_df = pd.read_csv(args.cna_file, sep = '\t', index_col=False)

    # get purities and filter by min_purity
    purities = get_purities(cna_df, num_samples, args.min_purity)

    # restrict samples considered in VCF and CNA file to those that have purity > min_purity
    vcf.set_samples(list(purities.keys()))
    cna_df = cna_df.loc[cna_df['SAMPLE'].isin(list(purities.keys()))]
    # print new CNA file, filtering out samples below min_purity
    cna_df.to_csv(f"{args.out_dir}/best.seg.ucn", sep="\t", index=False)
    if args.snp_file:
        snp_df = pd.read_csv(args.snp_file, sep = '\t', index_col=False, header=None)
        snp_df = snp_df.loc[snp_df[2].isin(list(purities.keys()))]
        # rearrange columns for decifer
        snp_df = snp_df[[2, 0, 1, 3, 4]] 
        snp_df.to_csv(f"{args.out_dir}/snpfile.1bed", sep="\t", index=False, header=False) 

    num_samples = len(vcf.samples)
    # print purity information
    sample_index = { vcf.samples[i] : i for i in range(len(vcf.samples)) }
    print_purities(purities, sample_index, num_samples, args.out_dir)

    # Filtering criteria
    Filter = {}
    Filter['MinDepth'] = args.min_depth
    Filter['MinDepthAltAllele'] = args.min_alt_depth
    Filter['MinVAF'] = args.min_vaf
    

    # ref_var_depths[char_label] = list of (ref,alt) tuples, one for each sample, in same order as vcf.samples
    ref_var_depths = compute_ref_var_depths(vcf, Filter)

    # print BED file for SNPs
    with open(f"{args.out_dir}/snps.bed", 'w') as out:
        print("chrom\tstart\tend\tREF\tALT", file=out)
        # sort ref_var_depths by the first two parts of chr_label
        for chr_label in sorted(ref_var_depths, key = lambda x: (x.split('.')[0], int(x.split('.')[1]))):
            pos = chr_label.split(".")
            # subtract 1 from position to create interval in BED format
            print(pos[0], int(pos[1])-1, int(pos[1]), pos[2], pos[3],  sep="\t", file=out)

    snps = pbt.BedTool(f"{args.out_dir}/snps.bed")
    if args.exclude_list:
        blist = pbt.BedTool(f"{args.exclude_list}")
        snps = snps.subtract(blist)


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
    cna_overlaps, cn_states_allsites, filtered_sites = overlap_cna_snp(vcf.samples, args.max_CN, snps, args.out_dir)
    
    # sites may have unique CN states that are duplicate; set them to find unique CN states across sites
    print_unique_CN_states(cn_states_allsites, args.max_CN, args.out_dir)
    print_filtered_sites(filtered_sites, cna_overlaps, args.out_dir)

    print_output(vcf, ref_var_depths, cna_overlaps, args.out_dir) 

    os.system(f"rm {args.out_dir}/snps.bed")
    for sample in vcf.samples:
        os.system(f"rm {args.out_dir}/{sample}_cna.bed")
    
if __name__ == '__main__':
  main()


