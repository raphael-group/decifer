"""
fileio.py

author: gsatas
date: 2020-05-04
"""

from new_coordinate_ascent import *


def read_in_state_trees(filename):
    with open(filename) as f:
        lines = f.readlines()
        state_trees = {}

        i = 1
        while i < len(lines):
            states = list(set(tuple(map(int, l.split(','))) for l in lines[i].strip().split(' ')))
            count = int(lines[i+1].strip().split(' ')[0])
            sts = []

            for j in range(1, count+1):
                state_tree = [tuple(map(int, l.split(','))) for l in lines[i+1+j].strip().split(' ')]

                sts.append(state_tree)

            state_trees[tuple(states)] = sts
            i = i + 1 + count + 1

    return state_trees


import pandas as pd
def read_in_test_file(filename):

    with open(filename) as f:
        f.readline()
        f.readline()
        header = f.readline().strip().split('\t')
        lines = [line.strip().split('\t') for line in f.readlines()]
        length = max([len(line) for line in lines])
        num_cstates = (length - len(header))//3
        for i in range(num_cstates):
            header += ["c{}a".format(i+1), "c{}b".format(i+1), "mu{}".format(i+1)]
        # pad lines
        lines = [l + ['']*(length - len(l)) for l in lines]

        df = pd.DataFrame(lines)
        df.columns = header
    return df

def write_results(prefix, C, mut_cluster_assignments, mut_config_assignments, mutations, purity, bb):

    with open('{}.output'.format(prefix), 'w') as out:
        out.write('mut_index\t'+"\t".join(['VAR_{}'.format(i) for i in range(len(C))]+['TOT_{}'.format(i) for i in range(len(C))])+'\tcluster\t'+"\t".join(['clust_CF{}'.format(i) for i in range(len(C))] + ['leftb_CF{}'.format(i)  for i in range(len(C))] + ['rightb_CF{}'.format(i)  for i in range(len(C))] + ['est_VAF{}'.format(i)  for i in range(len(C))] + ['est_CF{}'.format(i)  for i in range(len(C))] + ['cmm_CF{}'.format(i)  for i in range(len(C))]) + '\tExplained\tLHs' + '\n')
        for mut, clust in zip(mutations, mut_cluster_assignments):
            label = mut.label
            CF = [c[clust] for c in C]
            var = mut.a
            VAR = [v for v in var]
            tot = mut.d
            TOT = [v for v in tot]
            vaf = [float(v)/t if t > 0 else 0.5 for v,t in zip(var, tot)]
            config = mut.assigned_config
            leftbC = [config.cf_bounds(i)[0] for i in range(len(vaf))]
            rightbC = [config.cf_bounds(i)[1] for i in range(len(vaf))]
            estC = [config.v_to_cf(vaf[i], i, truncate = False) for i in range(len(vaf))]
            cmmC = [mut.compute_cmm_ccf(vaf[i], purity, i) for i in range(len(vaf))]

            explained = []
            lhs = []
            for cl in xrange(len(C[0])):
                if bb is None:
                    form = (lambda cf, sam : (cf.cf_to_v(C[sam][cl], sam), mut.a[sam]+1, (mut.d[sam] - mut.a[sam]) + 1))
                else:
                    form = (lambda cf, sam : (cf.cf_to_v(C[sam][cl], sam), mut.a[sam], mut.d[sam] - mut.a[sam], bb[sam]))
                grow = (lambda sample : compute_pdfs(*zip(*[form(cb, sample) for cb in mut.configs])))
                best = np.min(-np.sum(np.array([grow(sample) for sample in xrange(len(C))]), axis=0))
                if best < np.inf:
                    explained.append(cl)
                    lhs.append(best)
            explained = ';'.join(map(str, explained))
            lhs = ';'.join(map(str, lhs))
            
            out.write('\t'.join(map(str, [label] + VAR + TOT + [clust] + CF + leftbC + rightbC + vaf + estC + cmmC + [explained, lhs]))+'\n')

    with open('{}.C'.format(prefix), 'w') as out:
        for c in C:
            out.write('\t'.join(map(str, c)) + '\n')
    with open('{}.cluster_assignments'.format(prefix), 'w') as out:
        out.write('\n'.join(map(str,mut_cluster_assignments)))
    with open('{}.config_assignments'.format(prefix), 'w') as out:
        out.write('\n'.join(map(str,mut_config_assignments)))

def write_results_decifer_format(mutations, mut_cluster_assignments, prefix, num_clusters, num_samples, C):
    #SNV	cluster0	cluster1	cluster2	cluster	dcf0	dcf1	dcf2
    #0	1	0	0.00019	0	0.53	0.14	0.21
    #1	0	0	1	2	0.73	0.19	0.29

    # SNV index, cluster posteriors, assigned cluster, cluster dcfs
    with open('{}.tsv'.format(prefix), 'w') as out:
        out.write("\t".join(['SNV'] + ['cluster{}'.format(i) for i in range(num_clusters)] + ['cluster'] + ['dcf{}'.format(i) for i in range(num_samples)]) + '\n')

        for mut, clust in zip(mutations, mut_cluster_assignments):
            label = mut.index
            CF = [c[clust] for c in C]
            out.write('\t'.join(map(str, [label] + [0,]*num_clusters + [clust] + CF))+'\n')

def write_results_BIC(bic, ll, num_clusters, prefix):
    with open('{}.bic'.format(prefix), 'a') as out:
        out.write('\t'.join(map(str,[num_clusters, ll, bic])) + '\n')
