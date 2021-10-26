"""
This module reads in purity estimates from a file or calculates it from mutation data.
Purity is treated as a separate module so that it can be imported and used by both the
__main__.py top-level script and also the config.py script to rescale calculations
by purity.
"""

import sys
from decifer.parse_args import args
from decifer.fileio import read_input_file, read_purity

def calc_purity(MUTATION_DF):
    PURITY = {}
    for idx in MUTATION_DF.index:
        line = MUTATION_DF.loc[idx]
        label = line['character_label']
        index = line['character_index']
        ref = int(line['ref'])
        a = int(line['var'])
        d = a + ref
        sample = int(line['#sample_index'])

        # 6 columns for SNV info (sample index/label, character index/label, REF counts, ALT counts)
        # remainder of columns for copy number states, each state having 3 columns
        num_CNstates = (len(MUTATION_DF.columns) - 6) // 3
        CNs = []  # list of tuples of allele-specific copy numbers
        mus = []  # list of clone proportions
        for i in range(num_CNstates):
            try:
                # c = tuple of allele-specific copy numbers
                c = (int(line['c{}a'.format(i + 1)]),
                     int(line['c{}b'.format(i + 1)]))
                # mu = clone proportion
                mu = float(line['mu{}'.format(i + 1)])
                CNs.append(c)
                mus.append(mu)
            except ValueError:
                # print "ValueError"
                # TODO
                continue

        # sum all clone proportions that aren't (1,1), use to update purity[sample]
        mut_purity = sum([m for c, m in zip(CNs, mus) if c != (1, 1)])
        try:
            PURITY[sample] = max(mut_purity, PURITY[sample])
        except:
            PURITY[sample] = mut_purity

    return PURITY

def filter_MUTATION_DF(MUTATION_DF, PURITY):
    purity_new = {}

    MUTATION_DF["#sample_index"] = MUTATION_DF["#sample_index"].astype("int")
    # recode_dict is for filtering MUTATION_DF later
    recode_dict = {}
    new_index = 0
    # find which samples have inappropriate ploidy
    for old_index in range(len(PURITY.keys())):
        if PURITY[old_index] > 0.0 and PURITY[old_index] <= 1.0:
            recode_dict[old_index] = new_index
            purity_new[new_index] = PURITY[old_index]
            # new_index doesn't get incremented if bad ploidy
            new_index += 1
    # only include sample indices with ploidy > 0 and <= 1
    MUTATION_DF = MUTATION_DF.loc[ MUTATION_DF["#sample_index"].isin(list(recode_dict.keys())) ]
    # re-index samples so that they start at 0 and are consecutive integers
    MUTATION_DF = MUTATION_DF.replace({"#sample_index": recode_dict})
    return MUTATION_DF, purity_new

# Read/calculate purity information
if args['purity'] is not None:
    PURITY = read_purity(args['purity'])
else:
    # calculate purity from the mutation file as the maximum proportion of
    # non-diploid (not (1,1)) clones
    PURITY = calc_purity(MUTATION_DF)
for sample in PURITY:
    if PURITY[sample] <= 0.0 or PURITY[sample] > 1.0:
        sys.stderr.write(f"WARNING: sample with index {sample} has purity of zero, negative, or greater than one, and will be removed\n")

# read input mutation data, store as pd.DataFrame
MUTATION_DF = read_input_file(args["input"])

# filter samples with inappropriate PURITY values
# update PURITY dict to reflect re-indexed samples
MUTATION_DF, PURITY = filter_MUTATION_DF(MUTATION_DF, PURITY)
