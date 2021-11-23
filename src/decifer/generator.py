"""
decifer.py

author: Simone Zaccaria
date: 2020-06-04
"""

import sys, os
import argparse
import warnings
import datetime
import multiprocessing as mp
import random as rand
import bisect
from copy import deepcopy
from multiprocessing import Lock, Value, Pool, Manager, Queue

import numpy as np
from scipy.special import comb
from scipy.special import betaln
from scipy.special import gammaln
from scipy.optimize import minimize_scalar

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns

from collections import defaultdict

from decifer.fileio import *
from decifer.progress_bar import ProgressBar


SEQERROR = 1e-03


def main():
    args = parse_args()
    sys.stderr.write('> Reading true BAF from SEG file of copy numbers and proportions\n')
    baf, bsamples = read_baf(args['segfile'], args['threshold'])
    sys.stderr.write('> Reading SNPs\n')
    snps, osamples = read_snps(args['snpfile'])
    assert bsamples == osamples, 'Found different samples in the SEG and SNP files'
    samples = bsamples.union(osamples)
    sys.stderr.write('> Overlapping SNPs and copy-number segments\n')
    overlap, tot, sel = get_overlap(baf, snps, args['J'])
    sys.stderr.write('The number of retained SNPs is {} over a total of {} SNPs ({:.2%})\n'.format(sel, tot, sel / float(tot)))
    sys.stderr.write('> Fitting Beta-Binomial\n')
    betabinom = fit(baf, overlap, samples, args['restarts'], args['skip'], args['J'])
    with open(f"{args['output']}_betabinom.tsv", 'w') as out:
        out.write('#SAMPLE\tPRECISION\tNLH\n')
        for sam in sorted(betabinom):
            out.write('\t'.join(map(str, [sam, betabinom[sam][0], betabinom[sam][1], "\n"])))
    sys.stderr.write('> Plotting binomial fit\n')
    set_style()
    plot_binomial(overlap, baf, snps, samples)
    sys.stderr.write('> Plotting beta-binomial fit\n')
    plot_betabinomial(overlap, baf, snps, samples, betabinom)
    

def parse_args():
    description = "DeCiFer."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i","--snpfile", type=str, required=True, help="File with read counts for germline SNP from all samples")
    parser.add_argument("-s","--segfile", type=str, required=True, help="SEG file with allele-specific copy numbers")
    parser.add_argument("-e","--sensitivity", type=float, required=False, default=0.1, help="Sensitivity E to exclude SNPs with 0.5 - E <= BAF < 0.5 (default: 0.1)")
    parser.add_argument("-j","--jobs", type=int, required=False, default=mp.cpu_count(), help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-r","--restarts", type=int, required=False, default=10000, help="Maximum size of brute-force search (default: 1e4)")
    parser.add_argument("-k","--skip", type=int, required=False, default=10, help="Numbers to skip in the brute-force search (default: 10)")
    parser.add_argument("-o","--output", type=str, required=False, default="./", help="Output prefix (default: ./)")
    parser.add_argument("--seed", type=int, required=False, default=None, help="Random-generator seed (default: None)")
    args = parser.parse_args()

    if not os.path.isfile(args.segfile):
        raise ValueError("SEG file does not exist!")
    if not os.path.isfile(args.snpfile):
        raise ValueError("SNP file does not exist!")

    if args.seed:
        rand.seed(args.seed)
        np.random.seed(seed=args.seed)

    return {
        "segfile" : args.segfile,
        "snpfile" : args.snpfile,
        "threshold" : args.sensitivity,
        "J" : args.jobs,
        "restarts" : args.restarts,
        "skip" : args.skip,
        "output" : args.output
    }


def read_baf(segfile, threshold):
    cns = defaultdict(lambda : defaultdict(lambda : dict()))
    samples = set()
    mirror = (lambda v : min(v, 1.0 - v))
    def check(L):
        assert 0.98 <= sum(i[1] for i in L) <= 1.02, '{}'.format(L)
        return True
    getbaf = (lambda L : mirror(sum(e[0][1] * e[1] for e in L) / sum(sum(e[0]) * e[1] for e in L)) if check(L) else None)
    getcns = (lambda e : tuple(map(int, e.split('|'))))
    with open(segfile, 'r') as i:
        form = (lambda p : (p[0], (int(p[1]), int(p[2])), p[3], getbaf([(getcns(e), float(p[4:][x+1])) for x, e in enumerate(p[4:]) if x%2==0])))
        for l in (g for g in i if '#' not in g):
            chro, seg, sam, baf = form(l.strip().split())
            assert 0.0 <= baf <= 0.5, 'Error in the input:\n\tRead line: {}\tSegment: {}:{}-{}\n\tBAF: {}'.format(l, chro, seg[0], seg[1], baf)
            if not (0.5 - threshold <= baf < 0.5):
                cns[chro][seg][sam] = baf
                samples.add(sam)
    return {chro : {seg : cns[chro][seg] for seg in cns[chro]} for chro in cns}, samples


def read_snps(snpfile):
    samples = set()
    snps = defaultdict(lambda : defaultdict(lambda : dict()))
    with open(snpfile, 'r') as i:
        form = (lambda p : (p[0], p[1], int(p[2]), (int(p[3]), int(p[4]))))
        for l in (g for g in i if '#' not in g):
            sam, chro, pos, cnt = form(l.strip().split())
            snps[chro][pos][sam] = cnt
            samples.add(sam)
    return snps, samples


def get_overlap(baf, snps, J):
    extract = (lambda chro, D, L, R : {o : {sam : snps[chro][o][sam] for sam in snps[chro][o]} for o in D[L:R]} if L < R else dict())
    getoverlap = (lambda chro, D : {seg : extract(chro, D, bisect.bisect_left(D, seg[0]), bisect.bisect_right(D, seg[1])) for seg in baf[chro]})
    overlap = {chro : getoverlap(chro, sorted(snps[chro].keys())) for chro in baf if chro in snps}
    tot = sum(len(snps[c]) for c in snps)
    sel = sum(len(overlap[c][seg]) for c in overlap for seg in overlap[c])
    return overlap, tot, sel


def fit(baf, overlap, samples, R, skip, J):
    res = {}
    for sam in sorted(samples):
        sys.stderr.write('>> Fitting sample {}...\n'.format(sam))
        getphased = (lambda cnt, baf : min(cnt) if baf < 0.5 else cnt[rand.randint(0, 1)])
        form = (lambda c, s, o : (getphased(overlap[c][s][o][sam], baf[c][s][sam]), sum(overlap[c][s][o][sam]), baf[c][s][sam]))
        KS, NS, BAFS = map(np.array, zip(*[form(c, seg, o) for c in overlap for seg in overlap[c] for o in overlap[c][seg] if sam in overlap[c][seg][o] and sam in baf[c][seg]]))
        BAFS[BAFS < SEQERROR] = SEQERROR
        BAFS[BAFS > (1 - SEQERROR)] = 1 - SEQERROR
        CS = gammaln(NS + 1) - gammaln(NS - KS + 1) - gammaln(KS + 1) #comb(NS, KS)
        NSminusKS = NS - KS
        minusBAFS = 1 - BAFS
        LH = (lambda S : -np.sum(CS + betaln(KS + np.abs(S) * BAFS, NSminusKS + np.abs(S) * minusBAFS) - betaln(np.abs(S) * BAFS, np.abs(S) * minusBAFS)) if S != 0 else np.inf)
        pool = Pool(processes=J, initializer=init_brute, initargs=(CS, NS, KS, NSminusKS, BAFS, minusBAFS))
        bar = ProgressBar(total=R/skip, length=30, verbose=False, lock=Lock(), counter=Value('i', 0))
        bar.progress(advance=False, msg="Started")
        report = (lambda r : bar.progress(advance=True, msg="Completed {} [OBJ: {}]".format(r[0], r[1])))
        bbopt = min((x for x in pool.imap_unordered(run_brute, range(0, R, skip)) if report(x)), key=(lambda x : x[1]))
        pool.close()
        pool.join()
        opt = minimize_scalar(LH, bracket=(1, R))
        if not opt.success:
           sys.stderr.write('\tStatus: {}\n\tMessage: {}\n'.format(opt.status, opt.message))
        sys.stderr.write('Best brute-force fit at {} with obj {}\nBest brent fit at {} with obj {}\n'.format(bbopt[0], bbopt[1], np.abs(opt.x), opt.fun))
        res[sam] = (np.abs(opt.x), opt.fun) if opt.fun <= bbopt[1] else bbopt
    return res


def init_brute(CS, NS, KS, NSminusKS, BAFS, minusBAFS):
    global LH
    LH = (lambda S : -np.sum(CS + betaln(KS + np.abs(S) * BAFS, NSminusKS + np.abs(S) * minusBAFS) - betaln(np.abs(S) * BAFS, np.abs(S) * minusBAFS)) if S != 0 else np.inf)


def run_brute(S):
    return (S, LH(S))


def set_style(args=None):
    plt.style.use('ggplot')
    sns.set_style("whitegrid")
    #plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["axes.grid"] = True
    plt.rcParams["axes.edgecolor"] = "k"
    plt.rcParams["axes.linewidth"]  = 1.5


def plot_binomial(overlap, baf, snps, samples):
    fig, axes = plt.subplots(len(samples))
    for x, sam in enumerate(sorted(samples)):
        obsbaf = (lambda c, o : float(snps[c][o][sam][rand.randint(0, 1)]) / float(sum(snps[c][o][sam])))
        form = (lambda c, s, o : {'CN' : baf[c][s][sam], 'Observed BAF' : obsbaf(c, o)})
        df = [form(c, s, o) for c in overlap for s in overlap[c] for o in overlap[c][s] if sam in overlap[c][s][o] and sam in baf[c][s] and sum(snps[c][o][sam]) > 0]
        #g = distplot_fig(data=pd.DataFrame(df), x='Observed BAF', hue='CN', ax=axes[x])
        g = distplot_fig(data=pd.DataFrame(df), x='Observed BAF', hue='CN')
    plt.savefig('test.png', dpi=600, bbox_inches='tight')
        


def plot_betabinomial(overlap, baf, snps, samples, betabinom):
    return


def distplot_fig(data, x, hue=None, row=None, col=None, legend=True, hist=False, **kwargs):
    """A figure-level distribution plot with support for hue, col, row arguments."""
    bins = kwargs.pop('bins', None)
    if (bins is None) and hist: 
        # Make sure that the groups have equal-sized bins
        bins = np.histogram_bin_edges(data[x].dropna())
    g = sns.FacetGrid(data, hue=hue, row=row, col=col)
    g.map(sns.distplot, x, bins=bins, hist=hist, **kwargs)
    if legend and (hue is not None) and (hue not in [x, row, col]):
        g.add_legend(title=hue) 
    return g


if __name__ == "__main__":
    main()
