"""
fileio.py

author: gsatas
date: 2020-05-04
"""

from decifer.new_coordinate_ascent import compute_pdfs
import numpy as np
import pandas as pd

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


def read_input_file(filename):

    with open(filename) as f:
        f.readline()
        f.readline()
        header = f.readline().strip().split('\t')
        lines = [line.strip().split('\t') for line in f.readlines()]
        # copy number states vary across SNVs, find max number of observed states
        length = max([len(line) for line in lines])
        num_cstates = (length - len(header))//3

        # copy number states are in unlabelled columns, give them column names here
        # allows easy access to these columns later
        for i in range(num_cstates):
            header += ["c{}a".format(i+1), "c{}b".format(i+1), "mu{}".format(i+1)]
        # pad lines
        lines = [l + ['']*(length - len(l)) for l in lines]

        df = pd.DataFrame(lines)
        df.columns = header
    return df

def read_purity(purity_file):
    PURITY = {}
    with open(purity_file) as f:
        for line in f:
            line = line.strip().split('\t')
            PURITY[int(line[0])] = float(line[1])
    return PURITY

def write_results(prefix, C, CIs, mut_cluster_assignments, mut_config_assignments, mutations, purity, bb, kind, printallk, k):
    if printallk:
        name = f"{prefix}_output_K{k}.tsv"
    else:
        name = f"{prefix}_output.tsv"
    with open(name, 'w') as out:
        out.write('mut_index\t'+"\t".join(['VAR_{}'.format(i) for i in range(len(C))]+['TOT_{}'.format(i) for i in range(len(C))])+'\tcluster\tstate_tree\t'+"\t".join(['true_cluster_{}{}'.format(kind, i) for i in range(len(C))] + ['point_estimate_{}{}'.format(kind, i)  for i in range(len(C))] + ['cmm_CCF{}'.format(i) for i in range(len(C))]) + '\tExplained\tLHs' + '\n')
        #for mut, clust in zip(mutations, mut_cluster_assignments):
        for mut, clust in zip(mutations, mut_cluster_assignments):
            label = mut.label
            #CF = [c[clust] for c in C]
            CF = [ ";".join([str(C[i][clust]), str(CIs[i][clust]).replace(" ", "")]) for i in range(len(C))]
            var = mut.a
            VAR = [v for v in var]
            tot = mut.d
            TOT = [v for v in tot]
            vaf = [float(v)/t if t > 0 else 0.5 for v,t in zip(var, tot)]
            config = mut.assigned_config
            assigned_tree = mut.assigned_tree()
            tree = ';'.join(map(lambda p : '({},{},{})->({},{},{})'.format(*(p[0]+p[1])), zip(assigned_tree[:-1:2], assigned_tree[1::2])))
            leftbC = [config.cf_bounds(i)[0] for i in range(len(vaf))]
            rightbC = [config.cf_bounds(i)[1] for i in range(len(vaf))]
            estC = [config.v_to_cf(vaf[i], i, truncate = False) for i in range(len(vaf))]
            cmmC = [mut.compute_cmm_ccf(vaf[i], purity, i) for i in range(len(vaf))]

            explained = []
            lhs = []
            for cl in range(len(C[0])):
                if bb is None:
                    form = (lambda cf, sam : (cf.cf_to_v(C[sam][cl], sam), mut.a[sam]+1, (mut.d[sam] - mut.a[sam]) + 1))
                else:
                    form = (lambda cf, sam : (cf.cf_to_v(C[sam][cl], sam), mut.a[sam], mut.d[sam] - mut.a[sam], bb[sam]))
                grow = (lambda sample : compute_pdfs(*zip(*[form(cb, sample) for cb in mut.configs])))
                best = np.min(-np.sum(np.array([grow(sample) for sample in range(len(C))]), axis=0))
                if best < np.inf:
                    explained.append(cl)
                    lhs.append(best)
            explained = ';'.join(map(str, explained))
            lhs = ';'.join(map(str, lhs))
            
            #out.write('\t'.join(map(str, [label] + VAR + TOT + [clust] + [tree] + CF + leftbC + rightbC + vaf + estC + cmmC + [explained, lhs]))+'\n')
            out.write('\t'.join(map(str, [label] + VAR + TOT + [clust] + [tree] + CF + estC + cmmC + [explained, lhs]))+'\n')

    # with open('{}.C'.format(prefix), 'w') as out:
    #     for c in C:
    #         out.write('\t'.join(map(str, c)) + '\n')
    # with open('{}.cluster_assignments'.format(prefix), 'w') as out:
    #     out.write('\n'.join(map(str,mut_cluster_assignments)))
#    with open('{}.config_assignments'.format(prefix), 'w') as out:
#        out.write('\n'.join(map(str,mut_config_assignments)))

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

def write_results_CIs(prefix, num_samples, clus, sample_ids, CIs, printallk, k):
    if printallk:
        name = f"{prefix}_clusterCIs_K{k}.tsv"
    else:
        name = f"{prefix}_clusterCIs.tsv"
    with open(name, 'w') as f:
        f.write(" ".join( [str(num_samples), "#anatomical sites", "\n"] ))
        f.write(" ".join( [str(num_samples), "#samples", "\n"] ))
        f.write(" ".join( [str(len(set(clus))), "#mutation clusters", "\n"] ))
        header = ["#sample_index", "sample_label", "anatomical_site_index", "anatomical_site_label", "cluster_index", "cluster_label", "f_lb", "f_ub", "\n"]
        f.write("\t".join(header))
        for sample_index, sample_label in sample_ids.items():
            for cluster_index, cluster_name in enumerate(set(clus)):
                info = [sample_index, sample_label]
                info.extend([sample_index, sample_label]) # using these as anatomical index and anatomical name, for now
                info.extend( [cluster_index, cluster_name] )
                info.extend( [CIs[sample_index][cluster_name][0], CIs[sample_index][cluster_name][1], "\n"])
                f.write("\t".join(list(map(str, info))))

def write_model_selection_results( k, mink, maxk, objs, elbow, selected, prefix ):
    with open(f"{prefix}_model_selection.tsv", 'w') as f:
        f.write('\t'.join(['#NUM_CLUSTERS', 'BEST_OBJ', 'ELBOW_SCORE', 'SELECTED', "\n"]))
        for k in range(mink, maxk + 1):
            f.write('\t'.join(map(str, [k, objs[k], elbow[k] if k < maxk else 'NaN', selected == k])))
            f.write("\n")