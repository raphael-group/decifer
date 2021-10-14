"""
This module reads in purity estimates from a file or calculates it from mutation data.
Purity is treated as a separate module so that it can be imported and used by both the
__main__.py top-level script and also the config.py script to rescale calculations
by purity.
"""

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

# read input mutation data, store as pd.DataFrame
MUTATION_DF = read_input_file(args["input"])

if args['purity'] is not None:
    PURITY = read_purity(args['purity'])
else:
    # calculate purity from the mutation file as the maximum proportion of
    # non-diploid (not (1,1)) clones
    PURITY = calc_purity(MUTATION_DF)