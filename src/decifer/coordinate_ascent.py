"""
coordinate_ascent.py

author: gsatas
date: 2020-05-04
"""



EPSILON=10e-10

from fileio import *
from mutation import *
from config import *
from mutation import *
#from visualize import *
from random import choice, random, seed
seed(6)

import math
#from matplotlib import pyplot

#NUM_CLUSTERS = 3

def optimize_cluster_centers(mutations, num_samples, C_old):

    from scipy.optimize import minimize_scalar

    C_new = []
    obj_total = 0
    for sample in range(num_samples):
        C_sample = []
        for i in range(NUM_CLUSTERS):
            mutations_i = [m for m in mutations if m.assigned_cluster == i]
            if len(mutations_i) > 0:
                v = minimize_scalar(objective, args = (mutations_i,sample), method = 'bounded', bounds = [0,1]).x

                # Minimize scalar doesn't always find the global minima. This prevents
                # the objective from getting worse than we already had
                obj= objective(v, mutations_i, sample)
                if C_old != None:
                    obj_old = objective(C_old[sample][i], mutations_i, sample)
                    if obj_old < obj:
                         v = C_old[sample][i]
                         obj = obj_old
                obj_total += obj
            else:
                v = random()

            C_sample.append(v)
        C_new.append(C_sample)



    return C_new, obj_total

def optimize_assignments(mutations, C, num_samples):
    # For each mutation iterate through possible cluster assignments and configs

    total_obj = 0
    for m in mutations:
        min_post = float('inf')
        min_config = None
        min_clust = None
        for clust in range(NUM_CLUSTERS):
            for j, config in enumerate(m.configs):
                v = -sum([posterior_c(C[sample][clust], m.a, m.d, config, sample) for sample in range(num_samples)])
                if v < min_post:
                    min_post = v
                    min_config = config
                    min_clust = clust
        total_obj += min_post
        m.assigned_config = min_config
        m.assigned_cluster = min_clust

    print("--------")
    return mutations, total_obj

def coordinate_ascent(mutations, num_samples, max_iters):
    ### Initialize
    for m in mutations:
        # Randomly initialize config
        chosen_config = choice(range(len(m.configs)))
        m.assigned_config = m.configs[chosen_config]
        m.assigned_cluster = choice(range(NUM_CLUSTERS))

    tolerance = 0.001
    C = None
    C_old = None

    def converged(C, C_old):
        if C_old == None: return False
        return all([all([abs(d-d2) < tolerance for d, d2 in zip (c, c2)])for c,c2 in zip(C, C_old)])

    objs = []
    iterN = 0
    while not converged(C,C_old) and iterN < max_iters:
        iterN += 1

        C_old = C
        C, obj = optimize_cluster_centers(mutations, num_samples, C_old)
        objs.append(obj)
        print ("After optimize cluster centers", [["{:.2f}".format(v) for v in c] for c in C],obj)
        mutations, obj = optimize_assignments(mutations, C, num_samples)
        objs.append(obj)
        #for mutation in mutations:
        #    print(mutation.assigned_cluster)
        print("After optimize assignments    ", [["{:.2f}".format(v) for v in c] for c in C], obj)

    mut_cluster_assignments = [m.assigned_cluster for m in mutations]
    mut_config_assignments = [m.assigned_config for m in mutations]
    #plot_objective(objs)
    return C, mutations, mut_cluster_assignments, mut_config_assignments, objs[-1]

#def plot_objective(objs):
#    pyplot.plot(range(len(objs)), objs)

from scipy.stats import beta
def posterior_c(c, a, d, config, sample, mut=None):
    try:
        v = config.c_to_v(c, sample)
    except IndexError:
        print("Mutation {}".format(mut.label))
        raise
    if v:
        post = max(beta.pdf(v, a[sample]+1, (d[sample]-a[sample]) + 1), EPSILON)
    else:
        post = EPSILON

    return math.log(post)

def objective(c, mutations, sample):
    obj = sum([posterior_c(c, mut.a, mut.d, mut.assigned_config, sample, mut) for mut in mutations])
    return -obj

from numpy import prod

def compute_BIC(LL, num_mutations, num_clusters, num_samples):
    return (-2 * LL) + (math.log(num_mutations) * num_clusters * num_samples)


if __name__ == "__main__":

    from fileio import *

    import sys

    if len(sys.argv) < 4:
        print("Usage: coordinate_ascent.py input_filename output_prefix num_clusters")
        exit(1)
    filename = sys.argv[1]
    global NUM_CLUSTERS
    prefix = sys.argv[2]
    NUM_CLUSTERS = int(sys.argv[3])

    max_iters = 100

    if len(sys.argv) > 4:
        print("Purity file provided")
        purity_file = sys.argv[4]
        

    mutation_data = read_in_test_file(filename)

    num_samples = len(mutation_data['sample_label'].unique())

    script_dir=sys.path[0]
    state_trees = read_in_state_trees(script_dir+'/state_trees.txt')
    mutations, purity = create_mutations(mutation_data, state_trees)

    if purity_file:
        purity = {}
        with open(purity_file) as f:
            for line in f:
                line = line.strip().split('\t')
                purity[int(line[0])] = float(line[1])
       

    min_obj = float('inf')
    min_mutations = None
    from copy import deepcopy
    for i in range(5):
        C, best_mutations, mut_cluster_assignments, mut_config_assignments, obj = coordinate_ascent(mutations, num_samples, max_iters)
        if obj < min_obj:
            C_best = C
            mut_cluster_assignments_best = mut_cluster_assignments
            mut_config_assignments_best = mut_config_assignments
            min_obj = obj
            min_mutations = deepcopy(best_mutations)

        print("------------------------------------------"*3)



    #prefix = 'temp.1000'
    write_results(prefix,C_best,mut_cluster_assignments_best,mut_config_assignments_best, min_mutations,purity)
    write_results_decifer_format(min_mutations, mut_cluster_assignments_best, prefix, NUM_CLUSTERS, num_samples, C_best)
    write_results_BIC(compute_BIC(-min_obj, len(mutations), NUM_CLUSTERS, num_samples), -min_obj, NUM_CLUSTERS, prefix)
    #pyplot.xlabel('iteration')
    #pyplot.ylabel('Log posterior probability')
    #pyplot.show()

    #plot_ccfs(C_best, min_mutations, purity)
