"""
decifer.py

author: Simone Zaccaria
date: 2020-05-21
"""

import sys, os
import warnings
import datetime
import traceback
import multiprocessing as mp
import random as rand
from collections import defaultdict

from copy import deepcopy
from multiprocessing import Lock, Value, Pool, Manager
import numpy as np
import math
from bisect import bisect_left

# decifer
from decifer.parse_args import args
from decifer.fileio import write_results, write_results_CIs, read_in_state_trees, write_model_selection_results
from decifer.new_coordinate_ascent import coordinate_descent, objective
from decifer.mutation import create_mutations
from decifer.process_input import PURITY, MUTATION_DF
from decifer.progress_bar import ProgressBar
from decifer.generator import fit_betabinom


def main():
    sys.stderr.write('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n')

    # create dictionary of sample indices and labels for printing later
    sample_ids = { int(i[0]) : i[1] for i in zip(MUTATION_DF['#sample_index'].unique(), MUTATION_DF['sample_label'].unique()) }
    for i in sample_ids:
        print(i, sample_ids[i])
    num_samples =  len(sample_ids)

    if args['mink'] < 2 + num_samples:
        args['mink'] = 2 + num_samples
        sys.stderr.write('## The minimum number of clusters has been increased to {} to account for fixed clusters!\n'.format(args['mink']))
    if args['maxk'] < args['mink']:
        args['maxk'] = args['mink']
        sys.stderr.write('## The maximum number of clusters has been increased to {} to be higher than the minimum!\n'.format(args['maxk']))
    # state_trees dict: keys are potential CN observed CN states (x,y), e.g. ((1, 0), (1, 1), (0, 0))
    # values are lists of lists, where each list contains all possible genotypes (x,y,m), e.g.
    # [(1, 0, 0), (0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 0, 0)]
    # these lists are always even, and each pair (i,i+1) is an edge in the state tree
    state_trees = read_in_state_trees(args['statetrees'])

    # store info from MUTATION_DF pd.DataFrame in mutations list, containing Mutation objects as elements
    # infer purity from SNV data, but overwrite purity dict if purity file provided (next line)
    mutations = create_mutations(MUTATION_DF, state_trees, not args['ccf'])

    if args['betabinomial']:
        betabinom_out = fit_betabinom(args)
        # invert sample_ids to get index for each sample name, according to input file
        get_index = {v : k for k,v in sample_ids.items()}
        betabinomial = {get_index[sam] : betabinom_out[sam][0] for sam in betabinom_out.keys()}

    if args['record']:
        manager = mp.Manager()
        record = manager.list()
    if not args['iterative']:
        print("Using binary-search model selection")
        run_coordinator_binary(mutations, num_samples, PURITY, args,
                               record if args['record'] else None,
                               betabinomial if args['betabinomial'] else None)
    else:
        print("Using iterative model selection")
        run_coordinator_iterative(mutations, sample_ids, num_samples, PURITY, args,
                                  record if args['record'] else None,
                                  betabinomial if args['betabinomial'] else None)
    if args['record']:
        with open('record.log.tsv', 'w') as o:
            o.write('#NUM_CLUSTERS\tRESTART\tSEED\tITERATION\tOBJECTIVE\n')
            for r in record:
                o.write('{}\t{}\t{}\t{}\t{}\n'.format(r[0], r[1], r[2], r[3], r[4]))

def run_coordinator_iterative(mutations, sample_ids, num_samples, PURITY, args, record, betabinomial):
    mink, maxk, maxit, prefix, restarts, ubleft, J = unpck(args)
    jobs = [(x, k, np.random.randint(low=0, high=2**10)) for x in range(restarts) for k in range(mink, maxk+1)]
    # run in single-thread mode for development/debugging
    if args['debug']:
        shared = defaultdict(dict)
        # make objects global
        init_descent(mutations, num_samples, maxit, shared, record, betabinomial, PURITY)
        for job in jobs:
            run_descent(job)
    else:
        manager, shared = setup_shared()
        initargs = (mutations, num_samples, maxit, shared, record, betabinomial, PURITY)
        pool = Pool(processes=min(J, len(jobs)), initializer=init_descent, initargs=initargs)
        bar = ProgressBar(total=len(jobs), length=30, verbose=False, lock=Lock(), counter=Value('i', 0))
        bar.progress(advance=False, msg="Started")
        report = (lambda r : bar.progress(advance=True, msg="Completed {} for k={} [Iterations: {}]".format(r[0], r[1], r[3])))
        list(map(report, pool.imap_unordered(run_descent, jobs)))
    # best[cluster number k] = min objective across runs/restarts
    best = {k : min(list(filter(lambda j : j[1] == k, jobs)), key=(lambda j : shared['objs'][j])) for k in range(mink, maxk+1)}

    #ubleft = .25 * len(mutations) * num_samples * 10
    objs = {k : shared['objs'][best[k]] for k in best}
    for k in range(mink+1, maxk+1):
        if objs[k - 1] < objs[k]:
            best[k] = best[k - 1]
            objs[k] = objs[k - 1]
    chk = (lambda v : v if v != 0.0 else 0.01)
    left = (lambda k : min((objs[k - 1] - objs[k]) / abs(chk(objs[k - 1])), ubleft) if k > mink else ubleft)
    right = (lambda k : (objs[k] - objs[k+1]) / abs(chk(objs[k])))

    elbow = {k : left(k) - right(k) for k in range(mink, maxk)}
    if mink < maxk:
        selected = max(range(mink, maxk), key=(lambda k : elbow[k]))
    else:
        selected = mink

    write_model_selection_results( mink, maxk, objs, elbow, selected, prefix )

    if args['printallk']:
        k_to_print = [ k for k in range(mink, maxk+1) ]
    else:
        k_to_print = [ selected ]
    for k in k_to_print:
        C, bmut, clus, conf, objs = map(lambda D : shared[D][best[k]], ['C', 'bmut', 'clus', 'conf', 'objs'])
        # C is list of lists; rows are samples, columns are cluster IDs, values are CCFs
        #CIs = [[()]*len(C[i]) for i in range(len(C))] # list of lists to store CIs, same structure as C
        CIs, PDFs = compute_CIs_mp(set(clus), bmut, num_samples, betabinomial, J, C)

        write_results_CIs(prefix, num_samples, clus, sample_ids, CIs, args['printallk'], k)
        write_results(prefix, C, CIs, clus, conf, bmut, PURITY, betabinomial, 'CCF' if args['ccf'] else 'DCF', args['printallk'], k)
        #write_results_decifer_format(bmut, clus, prefix, selected, num_samples, C)

    """
    # FOR TESTING
    #print_PDF(set(clus), bmut, num_samples, args['betabinomial'], C)
    #print_feasibleVAFs(set(clus), bmut, num_samples, args['betabinomial'], C)
    with open("pdfs.txt", 'w') as f:
        for c in set(clus):
            for s in range(num_samples):
                i = str(c) + "_" + str(s) + " "
                f.write(i)
                f.write(" ".join( list(map(str, PDFs[s][c]))))
                f.write("\n")
    with open("max_dcfs.txt", 'w') as f:
        for c in set(clus):
            for s in range(num_samples):
                print c, s, C[s][c], CIs[s][c][0], CIs[s][c][1]
                f.write(" ".join( list(map(str, [c, s, C[s][c], CIs[s][c][0], CIs[s][c][1]] ))))
                f.write("\n")
    """



def print_feasibleVAFs(cluster_ids, muts, num_samples, bb, C):
    with open("feasibleVAFs.txt", 'w') as f:
        for i in cluster_ids:
            mut = list(filter(lambda m : m.assigned_cluster == i, muts))
            for s in range(0,num_samples):
                lowers = [m.assigned_config.cf_bounds(s)[0] for m in mut]
                uppers = [m.assigned_config.cf_bounds(s)[1] for m in mut]
                f.write(" ".join(list(map(str, [i, s, max(lowers), min(uppers)]))))
                f.write("\n")

def print_PDF(cluster_ids, muts, num_samples, bb, C):
    with open("pdfs.txt", 'w') as f:
        for i in cluster_ids:
            mut = list(filter(lambda m : m.assigned_cluster == i, muts))
            for s in range(0,num_samples):
                max_dcf = C[s][i] # dcf value that maximizes posterior for this sample and cluster ID
                delta = (-1*objective(max_dcf, mut, s, bb))-2
                #prob = (lambda x: math.exp(-1*(x+delta))) # convert neg log to probability
                prob = (lambda x: math.exp(-1*(x))) # convert neg log to probability
                l = []
                for j in np.linspace(0, 1, 1000):
                   l.append(prob(objective(j, mut, s, bb)))
                total = np.sum(l)
                l = [x/total for x in l]
                f.write(" ".join(list(map(str, [i,s] + l))))
                f.write("\n")


def compute_CIs_mp(cluster_ids, muts, num_samples, bb, J, C):
    CIs = [[()]*len(C[i]) for i in range(len(C))] # list of lists to store CIs, same structure as C
    PDFs = [[()]*len(C[i]) for i in range(len(C))] # list of lists to store CIs, same structure as C
    num_tests = float(len(cluster_ids)*num_samples) # bonferroni correction for multiple hypothesis testing
    #C[s][i] is the putative mode of the pdf
    jobs = [(c, s, muts, num_tests, bb) for c in cluster_ids for s in range(num_samples)]
    pool = Pool(processes=min(J, len(jobs)))
    results = pool.imap_unordered(CI, jobs)
    pool.close()
    pool.join()
    for i in results:
        clust, samp, lower, upper, pdf = i[0], i[1], i[2], i[3], i[4]
        CIs[samp][clust] = (lower,upper)
        PDFs[samp][clust] = pdf

    return CIs, PDFs 

def CI(job):
    """
    Computes CIs for a sample-cluster combination

    There have been two issues in dealing with the "objective" function to characterize the PDF of CCF/DCF values, needed for obtaining CIs.
    1.) converting -log(unnormalized probability) from objective to an unnormalized probability produced prohibitively large numbers (-log numbers 
    are very negative such that e^-x huge), so we rescaled all values based on the most negative -log value observed from many samples from objective
    (previously we tried using cluster centers, but these do not succeed in finding mode especially when PDF is extremely narrow/disjoint from
    infeasible VAF values truncating distribution).
    2.) across the support of the CCF/DCF PDF distribution, some sample-cluster combinations have 0 values everywhere except an extremely narrow range,
    so trying to integrate approximately with functions like scipy.integrate.quad to get a normalization constant produces 0, because all sampled
    points yield 0.

    Thus, we have used a brute force method where we sample num_pts from the objective function, and use these to create and cahracterize the PDF of
    the CCF/DCF distribution.
    """
    c, s, muts, num_tests, bb = job # c is cluster, s is sample
    mut = list(filter(lambda m : m.assigned_cluster == c, muts))
    
    num_pts = 10000
    grid = [objective(j, mut, s, bb) for j in np.linspace(0, PURITY[s], num_pts)]
    min_log = min(grid)
    delta = (-1*min_log)-2      # constant to make -log(pdf) values less negative
    prob = (lambda x: math.exp(-1*(x+delta)))           # convert -log(pdf) to unnormalized probability
    total = sum([prob(x) for x in grid])                # unnormalized probabilities across support
    pdf = [prob(x)/total for x in grid]                 # normalized probabilities across support
    cdf = np.cumsum(pdf)

    low_ci = 0.025/num_tests                            # divide the desired CI quantile by the number of tests, bonferonni correction
    high_ci = 1 - low_ci

    low_index = take_closest(cdf, low_ci)
    high_index = take_closest(cdf, high_ci)

    l = float(low_index)/num_pts
    u = float(high_index)/num_pts

    return (c, s, l, u, pdf)

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.
    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos <= 0:
        return 0
    if pos >= len(myList):
        return pos-1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return pos 
    else:
        return pos-1

def run_coordinator_binary(mutations, num_samples, PURITY, args, record, betabinomial):
    L, R, maxit, prefix, restarts, ubleft, J = unpck(args)
    results = {}
    def evaluate(V):
        if V not in results:
            sys.stderr.write('[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + 'Computing for {} clusters...\n'.format(V))
            results[V] = run(mutations, num_samples, V, maxit, prefix, PURITY, restarts, ubleft, J, record, betabinomial)
        sys.stderr.write('[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + 'Objective for {} clusters: {}\n'.format(V, results[V][-1]))
        return

    MAXR = R
    M = L
    evaluate(MAXR)

    while(R - L > 1):
        M = int(math.floor(float(R + L) / 2.0))
        assert M not in {L, R}
        evaluate(M)
        if float(results[M][-1] - results[MAXR][-1]) / abs(results[MAXR][-1]) > ubleft:
            L = M
        else:
            R = M

    evaluate(L)
    evaluate(R)
    selected = L if float(results[L][-1] - results[MAXR][-1]) / abs(results[MAXR][-1]) <= ubleft else R
    C, bmut, clus, conf, objs = results[selected]
    write_results(prefix, C, clus, conf, bmut, PURITY, betabinomial, 'CCF' if args['ccf'] else 'DCF')
    #write_results_decifer_format(bmut, clus, prefix, selected, num_samples, C)


def run(mutations, num_samples, K, maxit, prefix, PURITY, restarts, ubleft, J, record, betabinomial):
    jobs = [(x, K, np.random.randint(low=0, high=2**10)) for x in range(restarts)]
    manager, shared = setup_shared()
    initargs = (mutations, num_samples, maxit, shared, record, betabinomial, PURITY)
    pool = Pool(processes=min(J, len(jobs)), initializer=init_descent, initargs=initargs)
    bar = ProgressBar(total=len(jobs), length=30, verbose=False, lock=Lock(), counter=Value('i', 0))
    bar.progress(advance=False, msg="Started")
    report = (lambda r : bar.progress(advance=True, msg="Completed {} for k={} [Iterations: {}]".format(r[0], r[1], r[3])))
    map(report, pool.imap_unordered(run_descent, jobs))
    best = min(jobs, key=(lambda j : shared['objs'][j]))
    C, bmut, clus, conf, objs = map(lambda D : shared[D][best], ['C', 'bmut', 'clus', 'conf', 'objs'])
    return C, bmut, clus, conf, objs


def setup_shared():
    manager = Manager()
    shared = {}
    shared['C'] = manager.dict()
    shared['bmut'] = manager.dict()
    shared['clus'] = manager.dict()
    shared['conf'] = manager.dict()
    shared['objs'] = manager.dict()
    return manager, shared


def init_descent(_mutations, _num_samples, _maxit, _shared, _record, _betabinomial, _purity):
    global gmutations, num_samples, maxit, shared, record, betabinomial, purity
    gmutations = _mutations
    num_samples = _num_samples
    maxit = _maxit
    shared = _shared
    record = _record
    betabinomial = _betabinomial
    purity = _purity


def run_descent(job):
    try:
        x, k, seed = job
        rand.seed(seed)
        np.random.seed(seed)
        mutations = deepcopy(gmutations)
        #with warnings.catch_warnings() as w:
        #   warnings.simplefilter("ignore")
        np.seterr(all='warn')
        C, best_mutations, mut_cluster_assignments, mut_config_assignments, obj, it = coordinate_descent(x, seed, mutations, num_samples, k, maxit, record, betabinomial, PURITY)
        shared['C'][job] = C
        shared['bmut'][job] = best_mutations
        shared['clus'][job] = mut_cluster_assignments
        shared['conf'][job] = mut_config_assignments
        shared['objs'][job] = obj
        return job + (it, obj)
    except:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def unpck(args):
    return args['mink'], args['maxk'], args['maxit'], args['output'], args['restarts'], args['elbow'], args['J']





if __name__ == "__main__":
    main()
