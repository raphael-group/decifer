# Script for creating DeCiFer input

DeCiFer uses information about copy-number aberrations (CNAs) and single-nucleotide variants (SNVs) to compute descendant cell fractions, but it does not itself quantify CNAs or call SNVs from sequencing data. Thus, users must employ other programs to get this information, which can then be combined into correct input files for DeCiFer using the `create_input.ipynb` python notebook in this directory.

Here we quantify CNAs using [HATCHet](https://github.com/raphael-group/hatchet) and call SNVs using [varscan](http://varscan.sourceforge.net/).

## Step 1

To generate SNV information in the correct format for `create_input.ipynb`, please use the following sequence of commands to use varscan's somatic function, where variables within curly brackets get replaced with the appropriate file name:
```
# In accordance with the standard varscan workflow, first generate mpileup files
samtools mpileup -f {reference_sequence.fa} -q 1 {tumor_file.bam} > {tumor.mpileup}
samtools mpileup -f {reference_sequence.fa} -q 1 {normal_file.bam} > {normal.mpileup}
# Run varscan to identify somatic SNPs in the tumor
varscan somatic {normal.mpileup} {tumor.mpileup} {base_name}
```

and varscan should create two files per tumor sample, starting with `base_name`, with extensions `.snp` and `.indel`. In this folder, we provide an example (`tumor.snp`) of what this `.snp` file should look like, as this file is used by `create_input.ipynb`.

## Step 2

Additionally, `create_input.ipynb` requires a pileup file for each tumor sample independently. This may be obtained using [bcftools](http://samtools.github.io/bcftools/bcftools.html) and using the output of varscan, to generate pileups only for the sites at which SNVs were called. To get these pileups, one for each tumor sample, use the following bcftools command:

```
bcftools mpileup {tumor_file.bam} -f {reference_sequence.fa} -T <(cut -f1-2 {tumor.snp} | grep -v position) -a INFO/AD -Ou | \
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%AD\n' > {tumor.mpileup.tsv} 
```

This should create a file with an extension `.mpileup.tsv`, and we provide an example (`tumor.mpileup.tsv`) of what this file should look like in this directory.

## Step 3

Lastly, one needs to provide information about CNAs to `create_input.ipynb`. We use the output from HATCHET, a file called `best.seg.ucn`. An example of what this file looks like is included in this directory. 

Note that if multiple tumor samples are used per patient, one needs to call SNVs and generate pileups (explained above) for each tumor sample from the patient. However, HATCHet is run on all tumor samples from the same patient simultaneously and generates a single `best.seg.ucn` file per patient. Thus, to use `create_input.ipynb` on a single patient for which *n* tumor samples are available, one needs *n* `.snp` files, *n* `.mpileup.tsv` files, and a single `best.seg.ucn` file.

## Step 4

Use `create_input.ipynb`!
