"""
mutation.py

author: gsatas
date: 2020-05-04
"""
from decifer.config import Config
import warnings


class Mutation:
    ''' A Mutation has a set of configurations, and read counts'''

    def __init__(self, configs, alltrees, a, d, label="NA", index=0):
        self.configs = configs
        self.trees = alltrees
        self.a = a
        self.d = d
        self.assigned_config = self.configs[0]
        self.assigned_cluster = 0
        self.label = label
        self.index = index

    def __str__(self):
        return "\t".join(map(str, [self.label, self.index, self.a, self.d, self.assigned_cluster, self.assigned_config]))

    def add_sample(self, a, d, cn_props):
        # cn_props dict: keys=cn_state, vals=proportions, e.g. cn_props[(1,1)] = 1.0 for diploid segment
        self.a.append(a)
        self.d.append(d)
        for config in self.configs:
            config.add_sample(cn_props, self.label)

    def compute_cmm_ccf(self, vaf, purity, sample):
        # Maybe CN proportions and F should live here in the Mutation class. As they're
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

    def assigned_tree(self):
        for c, t in zip(self.configs, self.trees):
            if c == self.assigned_config:
                return t


def get_edge_list(state_tree):
    # pairs off genotypes (x,y,m) in state tree into tuples
    # 2nd element of tuple is descendent of the first
    edge_list = []
    for i in range(0, len(state_tree), 2):
        edge_list.append((state_tree[i], state_tree[i+1]))
    return edge_list


def get_desc_set(state_tree, vertex):
    # vertex is mut_state, SSCN state
    edge_list = get_edge_list(state_tree)
    desc_set = []
    # convert tuple to list
    vertex = list(vertex)
    if len(vertex) == 2:
        vertex.append(1)
    for edge in [edge for edge in edge_list if list(edge[0]) == vertex]:
        desc_set.append(edge[1])
        desc_set += get_desc_set(state_tree, edge[1])
    return desc_set


# c1,mu1, c2 = None, mu2 = None):
def get_configs(_state_trees, CNs, mus, dcf_mode):
    if len(CNs) == 1:
        # Config is simple
        configuration = Config(mut_state=CNs[0], other_states=[], cn_props={
                               CNs[0]: [mus[0], ]}, desc_set=[], dcf_mode=dcf_mode)
        return [configuration, ], [[(1, 1, 0), (1, 1, 1)], ]
    else:
        # get state trees for observed CN states covering this mutation
        # state_trees (below) is list of lists, each list a particular state tree
        state_trees = _state_trees[tuple(set(CNs))]
        configurations = []
        alltrees = []

        # identify mut_state, the 2-tuple that indicates the CN state in which the mutation occurs
        # this CN state is in exactly two genotype states (x*,y*,0), (x*,y*,m)
        for state_tree_full in state_trees:
            saved_tree = [tuple(s) for s in state_tree_full]
            state_tree = set([tuple(s) for s in state_tree_full])
            # Identify mutation state, the CN state that the mutation occurs in
            for c in CNs:
                # for each observed CN state c, count how many times it occurs in state_tree set
                try:
                    num_states_c = sum(
                        [1 for s in state_tree if s[0] == c[0] and s[1] == c[1]])
                except:
                    print(state_tree)
                    print(CNs)
                    print(mus)
                    raise
                if num_states_c == 0:
                    continue
                elif num_states_c == 2:
                    mut_state = c
                    # Don't break out of this here because we still want to skip this
                    # state tree if not all states are in it

            # Identify the descendant genotypes (x,y,m) of mut_state
            desc_set = get_desc_set(state_tree_full, mut_state)
            # other_states include any genotype (x,y,m) with CN state != mut_state
            other_states = [s for s in state_tree if s[0]
                            != mut_state[0] or s[1] != mut_state[1]]
            # cn_props dict of lists, one proportion (mu) per sample
            cn_props = {c: [mu, ] for c, mu in zip(CNs, mus)}
            configuration = Config(mut_state=mut_state, other_states=other_states,
                                   cn_props=cn_props, desc_set=desc_set, dcf_mode=dcf_mode)
            configurations.append(configuration)
            alltrees.append(saved_tree)
        return configurations, alltrees


def create_mutations(mutation_data, state_trees, dcf_mode=True):
    # From a mutation data file, create mutation dict of Mutation objects
    mutations = {}

    for idx in mutation_data.index:
        line = mutation_data.loc[idx]
        label = line['character_label']
        index = line['character_index']
        ref = int(line['ref'])
        a = int(line['var'])
        d = a + ref
        sample = int(line['#sample_index'])

        # 6 columns for SNV info (sample index/label, character index/label, REF counts, ALT counts)
        # remainder of columns for copy number states, each state having 3 columns
        num_CNstates = (len(mutation_data.columns) - 6)//3
        CNs = [] # list of tuples of allele-specific copy numbers
        mus = [] # list of clone proportions
        for i in range(num_CNstates):
            try:
                # c = tuple of allele-specific copy numbers
                c = (int(line['c{}a'.format(i+1)]),
                     int(line['c{}b'.format(i+1)]))
                # mu = clone proportion
                mu = float(line['mu{}'.format(i+1)])
                CNs.append(c)
                mus.append(mu)
            except ValueError:
                # print "ValueError"
                # TODO
                continue

        # if multiple samples, encounter mutation index again
        # for newly encountered sample, append new a, d, and for each CN state, append the clone proportion
        if index in mutations:
            # mutation already encountered, add info for an additional sample, cn_props gets added to self.config
            # cn_props dict: keys=cn_state, vals=proportions, e.g. cn_props[(1,1)] = 1.0 for diploid segment
            cn_props = {c: mu for c, mu in zip(CNs, mus)}
            mutations[index].add_sample(a, d, cn_props)
        else:
            try:
                configs, alltrees = get_configs(state_trees, CNs, mus, dcf_mode)
            except KeyError:

                msg = "Skipping mutation {}: State tree file does not contain state trees for the set of copy-number states that affect mutation {}.\n To generate state trees, see documentation for `generatestatetrees`, included in the C++ component of DeCiFer.".format(
                    idx, idx)
                warnings.warn(msg)
                continue

            mut = Mutation(configs, alltrees, [a, ], [d, ], label, index)
            mutations[index] = mut

    return [mutations[m] for m in mutations]
