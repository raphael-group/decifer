# Script for creating DeCiFer input

DeCiFer uses information about copy-number aberrations (CNAs) and single-nucleotide variants (SNVs) to compute descendant cell fractions, but it does not itself quantify CNAs or call SNVs from sequencing data. Thus, users must employ other programs to get this information, which can then be combined into correct input files for DeCiFer using the `vcf_2_decifer.py` script in this directory.

### CNA calls

We quantify CNAs using [HATCHet](https://github.com/raphael-group/hatchet), although you may use any program to identify CNAs as long as you have a file of CNAs with the same format as `best.seg.ucn` (HATCHet output file) in this directory.

### SNV calls

Our `vcf_2_decifer.py` script takes the CNA calls mentioned above along with SNV calls in the standard VCF format to generate decifer input. This VCF must be multi-sample, meaning read depth information for reference and alternate alleles is available for **all** tumor samples (from the same patient) at every polymorphic site. In other words, if an SNV was detected in one tumor sample but not in a second tumor sample from the same patient, we still need read depth for these allele in all samples at this SNV site.

We have had excellent results using [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2) and [Strelka2](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md) to call SNVs, as both programs are capable of calling SNVs for all tumor samples present from the same patient. See [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) for multi-sample (or joint) calling with Mutect2, and [here](https://github.com/Illumina/strelka/issues/59) for multi-sample calling with Strelka2. We usually consider SNVs detected by both programs for downstream analyses.

If you have your SNV data in MAF format, please use this [script](https://github.com/mskcc/vcf2maf/blob/main/maf2vcf.pl) created by the MSKCC to convert from MAF to VCF. However, as stated above, for each patient your MAF files still need depth information (for reference and alternate alleles) for all tumor samples at each SNV site, not just the samples for which the SNV was called.


## Using `vcf_2_decifer.py`

This python script uses two other python packages: pybedtools and cyvcf2. Before running, you can create a conda environment to install these:

```
conda create -n vcf_bedtools pybedtools cyvcf2 pandas -y
conda activate vcf_bedtools
python vcf_2_decifer.py [OPTIONS]
```

Use `python vcf_2_decifer.py --help` to see the options. Options `MIN_DEPTH` and `MIN_ALT_DEPTH` may be set to, for instance, 8 and 3, respectively, for WGS data. These values should be much higher for WES/WXS data, depending on your sequencing depths. We recommend setting the `MAX_CN` to 6, such that sites containing a clone with a total copy number greater than 6 will get excluded.

Lastly, the sample names in the CNA file and the VCF file (containing SNVs) must agree with one another!

### `vcf_2_decifer.py` output

1. `decifer.input.tsv`: This is the primary input file for decifer.
2. `decifer.purity.tsv`: These are the per-sample purity estimates for decifer.
3. `cn_states.txt`: The unique copy-number states observed across all subclones at a site. E.g. a line with `2,2;1,1` indicates a SNV site in which one subclone has 2 copies of both maternal and paternal alleles (WGD) whereas the other subclone is diploid with only one copy of each. We provide an example of what this file might look like in this directory.
4. `filtered_sites.txt`: A list of sites, one per line, that were filtered due to the value specified for `MAX_CN`. Each site has the format `chromosome.position.REF_allele.ALT_allele`. This list is provided to see if any important sites (e.g. for your biological story) were filtered out even before the decifer analysis.
5. `filtered_stats.txt`: Shows the total number and the fraction of SNV sites that were filtered out due to the `MAX_CN` value.


## Adressing the "Skipping mutation warning"

When you run decifer, you may see an error very early on that looks like the following:

```
decifer/src/decifer/mutation.py:158: UserWarning: Skipping mutation ###: State tree file does not contain state trees for the set of copy-number states that affect mutation ###.
 To generate state trees, see documentation for `generatestatetrees`, included in the C++ component of DeCiFer
 ```

While the state tree file we provide should accomodate many users, there could be exceptions. If you observe this message for many SNV sites, you may generate a state tree file that is tailored to the CN states and subclones observed in your data using the following command:

```
generatestatetrees cn_states.txt > my_state_trees.txt
```

Where the `cn_states.txt` is the output file mentioned above and `my_state_trees.txt` is any file name of your choosing to store the state trees. Then, re-run decifer using the `--statetrees my_state_trees.txt` option, specifying the state tree file you just created. Please note that, when generating the `cn_states.txt` with `vcf_2_decifer.py`, if the value used for `MAX_CN` is greater than 6, the `generatestatetrees` function may take a very long time.

