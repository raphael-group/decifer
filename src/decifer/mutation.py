"""
mutation.py

author: gsatas
date: 2020-05-04
"""
from config import *
import warnings

class mutation:
    ''' A mutation has a set of configurations, and read counts'''
    def __init__(self, configs, a, d, label="NA", index = 0):
        self.configs = configs
        self.a = a
        self.d = d
        self.assigned_config = self.configs[0]
        self.assigned_cluster = 0
        self.label = label
        self.index = index

    def __str__(self):
        return "\t".join(map(str, [self.label, self.index, self.a, self.d, self.assigned_cluster, self.assigned_config]))

    def add_sample(self, a, d, cn_props):
        self.a.append(a)
        self.d.append(d)
        for config in self.configs:
            config.add_sample(cn_props, self.label)

    def compute_cmm_ccf(self, vaf, purity, sample):
        # Maybe CN proportions and F should live here in the mutation class. As they're
        # the same for all configs. not sure
        config = self.configs[0]
        cn_props = config.cn_props
        F = config.F(sample)
        # m = round(VAF * F / purity)
        m = max(1, round(vaf * F/purity[sample]))
        try:
            c = vaf * F / m
        except ZeroDivisionError:
            print(F, vaf, purity[sample], m, vaf * F/purity[sample])
            raise
        return c

def get_edge_list(state_tree):
    edge_list = []
    for i in range(0, len(state_tree),2):
        edge_list.append((state_tree[i], state_tree[i+1]))
    return edge_list




def get_desc_set(state_tree, vertex):

    edge_list = get_edge_list(state_tree)

    desc_set = []

    vertex = list(vertex)
    if len(vertex) == 2:
        vertex.append(1)

    for edge in [edge for edge in edge_list if list(edge[0]) == vertex]:
        desc_set.append(edge[1])
        desc_set += get_desc_set(state_tree, edge[1] )

    return desc_set

def get_configs(state_trees, Cs, mus, dcf_mode): #c1,mu1, c2 = None, mu2 = None):


    if len(Cs) == 1:
        # Config is simple
        configuration = config(mut_state=Cs[0], other_states=[], cn_props={Cs[0]:[mus[0],]}, desc_set = [], dcf_mode=dcf_mode)
        return [configuration,]

    else:


	state_trees = state_trees[tuple(set(Cs))]


        configurations = []
        for state_tree_full in state_trees:
            state_tree = set([tuple(s) for s in state_tree_full])

            # Identify mutation state
            for c in Cs:

                try:
                    num_states_c = sum([1 for s in state_tree if s[0] == c[0] and s[1] == c[1]])
                except:
                    print(state_tree)
                    print(Cs)
                    print(mus)
                    raise

                if num_states_c == 0: continue
                elif num_states_c == 2:
                    mut_state = c

                    # Don't break out of this here because we still want to skip this
                    # state tree if not all states are in it

            # Identify the descendant set of mut_state

            desc_set = get_desc_set(state_tree_full, mut_state)

            other_states = [s for s in state_tree if s[0] != mut_state[0] or s[1] != mut_state[1]]
            cn_props = {c:[mu,] for c,mu in zip(Cs, mus)}
            configuration = config(mut_state=mut_state, other_states=other_states, cn_props = cn_props, desc_set=desc_set, dcf_mode=dcf_mode)

            configurations.append(configuration)
        return configurations



import pandas as pd
def create_mutations(mutation_data, state_trees, dcf_mode = True):
    # From a mutation data file, create a set of mutations
    mutations = {}
    purity = {}

    for idx in mutation_data.index:
        line = mutation_data.loc[idx]
        label = line['character_label']
        index = line['character_index']
        ref = int(line['ref'])
        a = int(line['var'])
        d = a + ref


        sample = int(line['#sample_index'])

        c1 = (int(line['c1a']), int(line['c1b']))

        mu1 = float(line['mu1'])

        num_cstates = (len(mutation_data.columns) - 6)//3
        Cs = []
        mus = []
        for i in range(num_cstates):
            try:
                c = (int(line['c{}a'.format(i+1)]), int(line['c{}b'.format(i+1)]))
                mu = float(line['mu{}'.format(i+1)])
                Cs.append(c)
                mus.append(mu)
            except ValueError:
                #print "ValueError"
                #TODO
                continue

        mut_purity = sum([m for c,m in zip (Cs, mus) if c != (1,1)])
        try:
            purity[sample] = max(mut_purity, purity[sample])
        except:
            purity[sample] = mut_purity

        if index in mutations:
            cn_props = {c:mu for c, mu in zip(Cs, mus)}
            mutations[index].add_sample(a,d,cn_props)
        else:
            try:
                configs = get_configs(state_trees, Cs, mus, dcf_mode)
            except KeyError:
                
		msg="Skipping mutation {}: State tree file does not contain state trees for the set of copy-number states that affect mutation {}.\n To generate state trees, see documentation for `generatestatetrees`, included in the C++ component of DeCiFer.".format(idx, idx)
                warnings.warn(msg)
                continue

            mut = mutation(configs, [a,], [d,], label, index)
            mutations[index] = mut

    return [mutations[m] for m in mutations], purity
