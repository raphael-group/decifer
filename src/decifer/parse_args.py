import sys
import os
import argparse
import multiprocessing as mp
from pkg_resources import resource_filename
import random as rand
import numpy as np

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
    parser.add_argument("--record", required=False, default=False, action='store_true', help='Record objectives (default: False)')
    parser.add_argument("-j","--jobs", type=int, required=False, default=mp.cpu_count(), help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-o","--output", type=str, required=False, default="./decifer", help="Output prefix (default: ./decifer)")
    parser.add_argument("--statetrees", type=str, required=False, default=None, help="Filename of state-trees file (default: use state_trees.txt in the package)")
    parser.add_argument("--seed", type=int, required=False, default=None, help="Random-generator seed (default: None)")
    parser.add_argument("--debug", required=False, default=False, action='store_true', help='single-threaded mode for development/debugging')
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
        "statetrees" : statetrees,
        "debug" : args.debug
    }

args = parse_args()




