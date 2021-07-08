"""
decifer.py

author: Simone Zaccaria
date: 2020-05-21
"""


import sys, os
import argparse
import warnings
import datetime
import traceback
import multiprocessing as mp
import random as rand

from pkg_resources import resource_filename
from copy import deepcopy
from multiprocessing import Lock, Value, Pool, Manager
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import bisect

from fileio import *
from new_coordinate_ascent import *



def main():
    sys.stderr.write('Parsing and checking arguments\n')
    args = parse_args()
    sys.stderr.write('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n')
    mutation_data = read_in_test_file(args["input"])
    
    #create dictionary of sample indices and labels for printing later
    sample_ids = { int(i[0]) : i[1] for i in zip(mutation_data['#sample_index'].unique(), mutation_data['sample_label'].unique()) }
    for i in sample_ids:
        print i, sample_ids[i]
    num_samples =  len(sample_ids)

    if args['mink'] < 2 + num_samples:
        args['mink'] = 2 + num_samples
        sys.stderr.write('## The minimum number of clusters has been increased to {} to account for fixed clusters!\n'.format(args['mink']))
    if args['maxk'] < args['mink']:
        args['maxk'] = args['mink']
        sys.stderr.write('## The maximum number of clusters has been increased to {} to be higher than the minimum!\n'.format(args['maxk']))        
    script_dir = sys.path[0]
    state_trees = read_in_state_trees(args['statetrees'])
    mutations, purity = create_mutations(mutation_data, state_trees, not args['ccf'])
    if args['purity'] is not None:
        purity = read_purity(args['purity'], purity)
    if args['record']:
        manager = mp.Manager()
        record = manager.list()
    if not args['iterative']:
        print "Using binary-search model selection"
        run_coordinator_binary(mutations, num_samples, purity, args, record if args['record'] else None)
    else:
        print "Using iterative model selection"
        run_coordinator_iterative(mutations, sample_ids, num_samples, purity, args, record if args['record'] else None)
    if args['record']:
        with open('record.log.tsv', 'w') as o:
            o.write('#NUM_CLUSTERS\tRESTART\tSEED\tITERATION\tOBJECTIVE\n')
            for r in record:
                o.write('{}\t{}\t{}\t{}\t{}\n'.format(r[0], r[1], r[2], r[3], r[4]))


def parse_args():
    description = "DeCiFer."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", type=str, help="Input file in DeCiFer format.")
    parser.add_argument("-p","--purityfile", type=str, required=True, help="File with purity of each sample (TSV file in two columns`SAMPLE PURITY`)")
    parser.add_argument("-b","--betabinomial", type=str, required=False, default=None, help="File with precisions for betabinomial fit (default: binomial likelihood)")
    parser.add_argument("--ccf", required=False, default=False, action='store_true', help="Run with CCF instead of DCF (default: False)")
    parser.add_argument("-k","--mink", type=int, required=False, default=2, help="Minimum number of clusters, which must be at least 2 (default: 2)")
    parser.add_argument("-K","--maxk", type=int, required=False, default=12, help="Maximum number of clusters (default: 12)")
    parser.add_argument("-r","--restarts", type=int, required=False, default=100, help="Number of restarts (default: 100)")
    parser.add_argument("-t","--maxit", type=int, required=False, default=200, help="Maximum number of iterations per restart (default: 200)")
    parser.add_argument("-e","--elbow", type=float, required=False, default=0.06, help="Elbow sensitivity, lower values increase sensitivity (default: 0.06)")
    parser.add_argument("--binarysearch", required=False, default=False, action='store_true', help='Use binary-search model selection (default: False, iterative is used; use binary search when considering large numbers of clusters')
    parser.add_argument("--record", required=False, default=False, action='store_true', help='Record objectives (default: False')
    parser.add_argument("-j","--jobs", type=int, required=False, default=mp.cpu_count(), help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-o","--output", type=str, required=False, default="./decifer", help="Output prefix (default: ./decifer)")
    parser.add_argument("--statetrees", type=str, required=False, default=None, help="Filename of state-trees file (default: use state_trees.txt in the package)")
    parser.add_argument("--seed", type=int, required=False, default=None, help="Random-generator seed (default: None)")
    args = parser.parse_args()

    if not os.path.isfile(args.INPUT):
        raise ValueError("INPUT file does not exist!")
    if args.mink < 2:
        raise ValueError("The minimum number of clusters must be at least 2 since the first two clusters are always used to represent mutations that are either present or absent in all samples")

    if args.statetrees is not None:
        statetrees = args.statetrees
    else:
        statetrees = resource_filename(__name__, 'state_trees.txt')
    if not os.path.isfile(statetrees):
        raise ValueError("State tree file does not exist:\n{}".format(statetrees))
    
    if args.betabinomial is not None:
        with open(args.betabinomial, 'r') as i:
            betabinomial = {int(l.split()[0]) : float(l.split()[1]) for l in i if '#' not in l}
    else:
        betabinomial = None

    if args.seed:
        rand.seed(args.seed)
        np.random.seed(seed=args.seed)

    return {
        "input" : args.INPUT,
        "mink" : args.mink,
        "maxk" : args.maxk,
        "maxit" : args.maxit,
        "purity" : args.purityfile,
        "restarts" : args.restarts,
        "elbow" : args.elbow,
        "iterative" : not args.binarysearch,
        "record" : args.record,
        "J" : args.jobs,
        "output" : args.output,
        "ccf" : args.ccf,
        "betabinomial" : betabinomial,
        "statetrees" : statetrees
    }


def read_purity(purity_file, purity):
    purity = {}
    with open(purity_file) as f:
        for line in f:
            line = line.strip().split('\t')
            purity[int(line[0])] = float(line[1])
    return purity


def run_coordinator_iterative(mutations, sample_ids, num_samples, purity, args, record):
    mink, maxk, maxit, prefix, restarts, ubleft, J = unpck(args)
    jobs = [(x, k, np.random.randint(low=0, high=2**10)) for x in xrange(restarts) for k in xrange(mink, maxk+1)]
    manager, shared = setup_shared()
    initargs = (mutations, num_samples, maxit, shared, record, args['betabinomial'], purity)
    pool = Pool(processes=min(J, len(jobs)), initializer=init_descent, initargs=initargs)
    bar = ProgressBar(total=len(jobs), length=30, verbose=False, lock=Lock(), counter=Value('i', 0))
    bar.progress(advance=False, msg="Started")
    report = (lambda r : bar.progress(advance=True, msg="Completed {} for k={} [Iterations: {}]".format(r[0], r[1], r[3])))
    map(report, pool.imap_unordered(run_descent, jobs))
    best = {k : min(filter(lambda j : j[1] == k, jobs), key=(lambda j : shared['objs'][j])) for k in xrange(mink, maxk+1)}

    #ubleft = .25 * len(mutations) * num_samples * 10
    objs = {k : shared['objs'][best[k]] for k in best}
    for k in xrange(mink+1, maxk+1):
        if objs[k - 1] < objs[k]:
            best[k] = best[k - 1]
            objs[k] = objs[k - 1]
    chk = (lambda v : v if v != 0.0 else 0.01)
    left = (lambda k : min((objs[k - 1] - objs[k]) / abs(chk(objs[k - 1])), ubleft) if k > mink else ubleft)
    right = (lambda k : (objs[k] - objs[k+1]) / abs(chk(objs[k])))

    elbow = {k : left(k) - right(k) for k in xrange(mink, maxk)}
    if mink < maxk:
        selected = max(xrange(mink, maxk), key=(lambda k : elbow[k]))
    else:
        selected = mink
    print '\t'.join(['#NUM_CLUSTERS', 'BEST_OBJ', 'ELBOW_SCORE', 'SELECTED'])
    for k in xrange(mink, maxk+1):
        print '\t'.join(map(str, [k, objs[k], elbow[k] if k < maxk else 'NaN', selected==k]))

    C, bmut, clus, conf, objs = map(lambda D : shared[D][best[selected]], ['C', 'bmut', 'clus', 'conf', 'objs'])
    
    # C is list of lists; rows are samples, columns are cluster IDs, values are CCFs
    #CIs = compute_CIs_OLD(set(clus), bmut, num_samples, args['betabinomial'], C)
    #CIs = [[()]*len(C[i]) for i in range(len(C))] # list of lists to store CIs, same structure as C
    CIs, PDFs = compute_CIs_mp(set(clus), bmut, num_samples, args['betabinomial'], J, C)
    #print_PDF(set(clus), bmut, num_samples, args['betabinomial'], C)
    print_feasibleVAFs(set(clus), bmut, num_samples, args['betabinomial'], C)

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

    write_results_machina(num_samples, clus, sample_ids, CIs) 
    """
    with open("for_machina.txt", 'w') as f:
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
                #print c, s, C[s][c], CIs[s][c][0], CIs[s][c][1]
                #f.write(" ".join( list(map(str, [c, s, C[s][c], CIs[s][c][0], CIs[s][c][1]] ))))
                #f.write("\n")
    """
    

    write_results(prefix, C, CIs, clus, conf, bmut, purity, args['betabinomial'], 'CCF' if args['ccf'] else 'DCF')
    #write_results_decifer_format(bmut, clus, prefix, selected, num_samples, C)

def print_feasibleVAFs(cluster_ids, muts, num_samples, bb, C):
    with open("feasibleVAFs.txt", 'w') as f:
        for i in cluster_ids:
            mut = filter(lambda m : m.assigned_cluster == i, muts) 
            for s in range(0,num_samples):
                lowers = [m.assigned_config.cf_bounds(s)[0] for m in mut]
                uppers = [m.assigned_config.cf_bounds(s)[1] for m in mut]
                f.write(" ".join(list(map(str, [i, s, max(lowers), min(uppers)]))))
                f.write("\n")

def print_PDF(cluster_ids, muts, num_samples, bb, C):
    with open("pdfs.txt", 'w') as f:
        for i in cluster_ids:
            mut = filter(lambda m : m.assigned_cluster == i, muts) 
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
    mut = filter(lambda m : m.assigned_cluster == c, muts) 
    
    # narrow bounds for operaations; helped prevent errors induced by pdf=0 for majority of [0,1]
    lowerb = max(max([m.assigned_config.cf_bounds(s)[0] for m in mut]) - 0.05, 0.0) # max lower feasible VAF, 0 as limit
    upperb = min(min([m.assigned_config.cf_bounds(s)[1] for m in mut]) + 0.05, 1.0) # min upper feasible VAF, 1 as limit

    # more refined search for efficiency, but lose x-value info for pdf values
    #steps = int( (upperb-lowerb)/0.001 )                # grid search for min -log(pdf) by steps of 0.001
    #grid = [objective(j, mut, s, bb) for j in np.linspace(lowerb, upperb, steps)]
    num_pts = 10000

    grid = [objective(j, mut, s, bb) for j in np.linspace(0, 1, num_pts)]
    min_log = min(grid)
    delta = (-1*min_log)-2      # constant to make -log(pdf) values less negative
    prob = (lambda x: math.exp(-1*(x+delta)))           # convert -log(pdf) to unnormalized probability
    total = sum([prob(x) for x in grid])                # unnormalized probabilities across support
    pdf = [prob(x)/total for x in grid]                 # unnormalized probabilities across support
    quant = (lambda x,q: sum(pdf[0:int(x)+1]) - q)
    low_ci = 0.025/num_tests                            # divide the desired CI quantile by the number of tests, bonferonni correction

    # get indices of pdf that correspond to quantiles of interest
    # if first pdf point already greater than low ci, assign 0 (otherwise bisect fails; f(a) and f(b) have same sign)
    l = 0 if quant(1, low_ci) > 0 else int((bisect(quant, a=1, b=num_pts, args=(low_ci), xtol=0.001)))
    high_ci = 1 - low_ci
    # if last point already lower than high ci, assign num_pts
    u = num_pts if quant(num_pts, high_ci) < 0 else int((bisect(quant, a=1, b=num_pts, args=(high_ci), xtol=0.001)))
    l = float(l)/num_pts
    u = float(u)/num_pts

    """
    # integrate unormalized dcf/ccf PDF from lowerb -> x; note quad func takes it's own lambda func and return tuple of (area,error)
    area0x = (lambda x: integrate.quad(lambda dcf: prob(objective(dcf, mut, s, bb)), lowerb, x)) 
    Z = area0x(upperb)[0] # compute normalization constant to get probabilities
    if Z == 0:
        print "Z is 0"
        print "min_log ", min_log
        print "delta ", delta
        print "lower/upper ", lowerb, upperb 
    # Bisection method
    quant = (lambda x,q: (area0x(x)[0])/Z - q) # q is the desired CI quantile to find via bisection
    low_ci = 0.025/num_tests # divide the desired CI quantile by the number of tests, bonferonni correction
    l = bisect(quant, a=lowerb, b=upperb, args=(low_ci), xtol=0.001)
    high_ci = 1 - low_ci
    u = bisect(quant, a=lowerb, b=upperb, args=(high_ci), xtol=0.001)
    """

    return (c, s, l, u, pdf)

def compute_CIs_OLD(cluster_ids, muts, num_samples, bb, C):
    CIs = [[()]*len(C[i]) for i in range(len(C))] # list of lists to store CIs, same structure as C
    num_tests = float(len(cluster_ids)*num_samples) # bonferroni correction for multiple hypothesis testing
    print num_tests
    with open("max_dcfs.txt", 'w') as f:
        for i in cluster_ids:
            mut = filter(lambda m : m.assigned_cluster == i, muts) 
            for s in range(0,num_samples):
                max_dcf = C[s][i] # dcf value that maximizes posterior for this sample and cluster ID
                #delta = (-1*objective(max_dcf, mut, s, bb))-2
                #prob = (lambda x: math.exp(-1*(x+delta))) # convert neg log to probability
                prob = (lambda x: math.exp(-1*x)) # convert neg log to probability
                lowerb = max(max([m.assigned_config.cf_bounds(s)[0] for m in mut]) - 0.05, 0.0) # max lower bound, 0 as limit
                upperb = min(min([m.assigned_config.cf_bounds(s)[1] for m in mut]) + 0.05, 1.0) # min upper bound, 1 as limit
                # integrate unormalized dcf/ccf PDF from 0 -> x; note quad func takes it's own lambda func and return tuple of (area,error)
                area0x = (lambda x: integrate.quad(lambda dcf: prob(objective(dcf, mut, s, bb)), lowerb, x)) 
                Z = area0x(upperb)[0] # compute normalization constant to get probabilities

                quant = (lambda x,q: (area0x(x)[0])/Z - q) # divide the desired CI quantile by the number of tests

                f.write(" ".join( list(map(str, [i, s, max_dcf, objective(max_dcf, mut, s, bb), delta, Z] ))))
                f.write("\n")

                # Bisection method
                """
                low_ci = 0.025/num_tests
                l = bisect(quant, a=lowerb, b=upperb, args=(low_ci), xtol=0.001)
                high_ci = 1 - low_ci
                u = bisect(quant, a=lowerb, b=upperb, args=(high_ci), xtol=0.001)
                print i, s, l, u
                """




                """
                # OLD CODE
                area0x = (lambda x: integrate.quad(lambda dcf: prob(objective(dcf, mut, s, bb)), 0.0, x)) 
                Z = area0x(1)[0] # compute normalization constant to get probabilities
                #mid = (lambda x, mut, s, bb: (integrate.quad(lambda dcf: prob(objective(dcf, mut, s, bb)), lowerb, x) - 0.5)/Z )
                minim = minimize_scalar(mid, method='bounded', bounds=[lowerb,upperb], tol=TOLERANCE).x
                print "bounded"
                print type(minim)
                minim = minimize(mid, method='SLSQP', x0=((lowerb+upperb)/2.0), bounds=((lowerb,upperb),), tol=TOLERANCE).x[0]
                print "slsqp"
                print type(minim)

                print "Delta: ", delta
                print "obj at max: ", objective(max_dcf, mut, s, bb)
                print "Normalization const: ", Z
                print
                if Z > 0:
                    a,b = max_dcf, max_dcf # initialize variables for CI search
                    cdf_low, cdf_high = (area0x(a)[0]/Z), (area0x(b)[0]/Z)
                    # CI search, brute force as some dcf/ccf distributions truncated
                    while cdf_low > 0.025: 
                        a = a-0.005 if cdf_low > 0.1 else a-0.001 # smaller steps as approach desired probability
                        cdf_low = (area0x(a)[0]/Z) 
                    while cdf_high < 0.975: 
                        b = b+0.005 if cdf_high < 0.9 else b+0.001
                        cdf_high = (area0x(b)[0]/Z) 
                    a,b = max(0.0, a),min(1.0, b)
                else:
                    a,b = "NA", "NA"
                    check = []
                    for j in np.linspace(0, 1, 500):
                        check.append(str(objective(j, mut, s, bb)))

                    print i, s, " ".join(check)
                CIs[s][i] = (a,b)
                """
    return CIs

def run_coordinator_binary(mutations, num_samples, purity, args, record):
    L, R, maxit, prefix, restarts, ubleft, J = unpck(args)
    results = {}
    def evaluate(V):
        if V not in results:
            sys.stderr.write('[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + 'Computing for {} clusters...\n'.format(V))
            results[V] = run(mutations, num_samples, V, maxit, prefix, purity, restarts, ubleft, J, record, args['betabinomial'])
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
    write_results(prefix, C, clus, conf, bmut, purity, args['betabinomial'], 'CCF' if args['ccf'] else 'DCF')
    #write_results_decifer_format(bmut, clus, prefix, selected, num_samples, C)


def run(mutations, num_samples, K, maxit, prefix, purity, restarts, ubleft, J, record, betabinomial):
    jobs = [(x, K, np.random.randint(low=0, high=2**10)) for x in xrange(restarts)]
    manager, shared = setup_shared()
    initargs = (mutations, num_samples, maxit, shared, record, betabinomial, purity)
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
        C, best_mutations, mut_cluster_assignments, mut_config_assignments, obj, it = coordinate_descent(x, seed, mutations, num_samples, k, maxit, record, betabinomial, purity)
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


class ProgressBar:

    def __init__(self, total, length, lock, counter, verbose=False, decimals=1, fill=unichr(9608), prefix = 'Progress:', suffix = 'Complete'):
        self.total = total
        self.length = length
        self.decimals = decimals
        self.fill = fill
        self.prefix = prefix
        self.suffix = suffix
        self.lock = lock
        self.counter = counter
        assert lock is not None or counter == 0
        self.verbose = verbose

    def progress(self, advance=True, msg=""):
        flush = sys.stderr.flush
        write = sys.stderr.write
        if advance:
            with self.counter.get_lock():
                self.counter.value += 1
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (self.counter.value / float(self.total)))
        filledLength = int(self.length * self.counter.value // self.total)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        rewind = '\x1b[2K\r'
        result = '%s |%s| %s%% %s' % (self.prefix, bar, percent, self.suffix)
        msg = '[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + msg
        if not self.verbose:
            toprint = rewind + result + " [%s]" % (msg)
        else:
            toprint = rewind + msg + "\n" + result
        with self.lock:
            write(toprint.encode('utf-8'))
            flush()
            if self.counter.value == self.total:
                write("\n")
                flush()
        return True


if __name__ == "__main__":
    main()
