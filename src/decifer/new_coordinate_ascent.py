"""
new_coordinate_ascent.py

author: Simone Zaccaria
date: 2020-05-25
"""

import sys
import warnings

import numpy as np
from scipy.stats import beta
from scipy.special import betaln
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
import time

TOLERANCE = 1e-03
EPSILON = -1e40
SEQERROR = 1e-40


def coordinate_descent(restart, seed, mutations, num_samples, num_clusters, MAX_IT=5, record=None, bb=None, purity=None):
    # rest = number of variably present clusters subject to model selection, fixed clusters include truncal, absent, and sample specific SNVs
    rest = num_clusters - (2 + num_samples)
    assert rest >= 0, 'Number of clusters is too low!'
    # initialize cluster centers C; a list of lists, where each list is the cluster centers for a sample:
    # [absent, truncal, sample specific clusters (truncal if sample index, 0 otherwise),
    # variably present clusters with random initializations]
    form = ( lambda L, sam : [0.0, 1.0] + [purity[sam] if s == sam else 0.0 for s in range(num_samples)]
                            + list(map(lambda v : v / 10.0, list(L))) )
    C = [ form(np.random.randint(low=0, high=11, size=rest) if rest > 0 else [], sam) for sam in range(num_samples) ]
    V, C_old, V_old = None, None, None
    IT = 0
    objs = []

    # convergence if the difference between the last and 3rd-to-last objective less than TOLERANCE
    isconverged = (lambda O : False if len(O) < 3 else abs(O[-1] - O[-3]) < TOLERANCE)
    while not isconverged(objs) and IT < MAX_IT:
        mutations, obj = optimize_assignments(mutations, C, num_samples, num_clusters, bb)
        assert len(objs) == 0 or obj - objs[-1] <= TOLERANCE, 'Objective is not decreasing after assignments: {}'.format(objs + [obj])
        objs.append(obj)
        C_old, V_old = C, update_objs(C, mutations, num_samples, num_clusters, bb)
        C, V, obj = optimize_cluster_centers(mutations, num_samples, C_old, V_old, num_clusters, bb, purity)
        assert len(objs) == 0 or obj - objs[-1] <= TOLERANCE, 'Objective is not decreasing after centers: {}'.format(objs + [obj])
        objs.append(obj)
        IT += 1
    #C[0][5] = 0.0
    mutations, obj = optimize_assignments(mutations, C, num_samples, num_clusters, bb, last=True)
    assert len(objs) == 0 or obj - objs[-1] <= TOLERANCE, 'Objective is not decreasing after assignments: {}'.format(objs + [obj])    
    objs.append(obj)
    
    mut_cluster_assignments = [m.assigned_cluster for m in mutations]
    mut_config_assignments = [m.assigned_config for m in mutations]
    if record is not None:
        for x, obj in enumerate(objs):
            record.append((num_clusters, restart, seed, x, obj))
    return C, mutations, mut_cluster_assignments, mut_config_assignments, objs[-1], IT


def optimize_assignments(mutations, C, num_samples, num_clusters, bb, last=False):

    def select(m):
        # m in a single element from mutations

        # i *think* each element of m.configs list is Config object for a particular cn_state, cn_prop
        combs = [(config, clust) for config in m.configs for clust in range(num_clusters)]
        # form returns a 3-tuple (VAF, a, d-a), where
        # VAF is computed by converting from DCF/CCF cluster center for a sample, a is ALT depth, d is total depth
        if bb is None:
            form = (lambda cf, cl, sam : (cf.cf_to_v(C[sam][cl], sam), m.a[sam]+1, (m.d[sam] - m.a[sam]) + 1))
        else:
            form = (lambda cf, cl, sam : (cf.cf_to_v(C[sam][cl], sam), m.a[sam], m.d[sam] - m.a[sam], bb[sam]))
        # given a sample, returns log prob of mutation for each cluster???
        # zip takes list of 3-tuples (asterisk w/in zip makes each tuple separate arg)
        # zip returns 3 tuples/lists, one for each return value of form lambda: VAF, a, d-a
        # asterisk preceding zip makes each tuple/list a separate arg for compute_pdfs
        grow = (lambda sample : compute_pdfs( *zip(*[form(cb[0], cb[1], sample) for cb in combs]) ))
        # grow(sample) returns array of log prob of mutation (a, d-a) for each cluster
        # obj is array of summed log prob(mutation) across samples, 1 element for each cluster
        objs = -np.sum(np.array([grow(sample) for sample in range(num_samples)]), axis=0)
        # taking minimum of negative log probability maximizes probability, finding best cluster
        best = objs.argmin()
        assert objs[best] < np.inf, 'A mutation cannot be assigned to any cluster: {},{}\n\t{}\n'.format(m.index, m.label, map(lambda x : '{},{}'.format(x.mut_state, x.other_states), m.configs)) + str([m.configs[0].cf_bounds(sam) for sam in range(num_samples)]) + '\n' + str([C[sam][1] for sam in range(num_samples)]) + '\n'
        # compare best
        found = combs.index((m.assigned_config, m.assigned_cluster))
        before = -sum(compute_pdfs(*zip(*[form(combs[found][0], combs[found][1], sample) for sample in range(num_samples)])))
        assert objs[best] <= before, 'Non scende: {}'.format(C)
        return (float(objs[best]), combs[best][0], combs[best][1])

    def update(m, best):
        m.assigned_config = best[1]
        m.assigned_cluster = best[2]
        for sam in range(num_samples):
            assert m.assigned_config.cf_bounds(sam)[0] - 0.05 <= C[sam][m.assigned_cluster] <= m.assigned_config.cf_bounds(sam)[1] + 0.05, (C[sam][m.assigned_cluster], m.assigned_config.cf_bounds(sam), m.assigned_config.d_to_lam(C[sam][m.assigned_cluster], sam))
        return best[0]

    return mutations, sum(map(lambda m : update(m, select(m)), mutations))


def optimize_cluster_centers(mutations, num_samples, C_old, V_old, num_clusters, bb, purity):
    # vmin and vmax find the minimum and maximum feasible VAFs, respectively
    vmin = (lambda muti, sam: max( [max( [m.assigned_config.cf_bounds(sam)[0] for m in muti] ) - 0.05, 0.0]))
    vmax = (lambda muti, sam: min( [min( [m.assigned_config.cf_bounds(sam)[1] for m in muti] ) + 0.05, 1.0]))
    # given a sample and mutations, minim finds minimum -log(prob) of each DCF/CCF value (ci in objective)
    # i.e. it finds cluster center
    minim = (lambda muti, sam : minimize_scalar(objective, args=(muti, sam, bb), method='bounded', bounds=[vmin(muti,sam),vmax(muti,sam)], options={'xatol' : TOLERANCE}).x)
    # getmi reports the DCF/CCF value and the objective value (-log(prob)) for a particular DCF/CCF value (x)
    # returns tuple (DCF/CCF, -log(prob))
    getmi = (lambda muti, sam, x : (x, objective(x, muti, sam, bb) if len(muti) > 0 else 0.0))
    # caseg implements getmi, substituting random number for minim cluster center if no mutations
    caseg = (lambda muti, sam : getmi(muti, sam, minim(muti, sam) if len(muti) > 0 else np.random.rand()))
    # (i - 2) == sam is true if i is the sample-specific cluster for that sample
    cases = (lambda muti, sam, i : getmi(muti, sam, minim(muti, sam) if (i - 2) == sam and len(muti) > 0 else 0.0))

    # cas01 is same as getmi
    cas01 = (lambda muti, sam, x : (x, objective(x, muti, sam, bb) if len(muti) > 0 else 0.0))
    # caseb runs cas01 if i in {0,1}, absent or truncal cluster, o.w. runs cases
    caseb = (lambda muti, sam, i : cas01(muti, sam, 0.0 if i == 0 else 1.0) if i in {0, 1} else cases(muti, sam, i))
    # obj_i runs caseb if cluster is absent, truncal, or sample specific, o.w. runs caseg
    obj_i = (lambda muti, sam, i : caseb(muti, sam, i) if i < 2 + num_samples else caseg(muti, sam))
    selec = (lambda sam, i, R : (C_old[sam][i], V_old[sam][i]) if V_old[sam][i] < R[1] else R)
    R_sam = (lambda sam : [selec(sam, i, obj_i(list(filter(lambda m : m.assigned_cluster == i, mutations)), sam, i)) for i in range(num_clusters)])
    R_new = [R_sam(sam) for sam in range(num_samples)]
    C_new, V_new = map(list, zip(*[map(list, zip(*r)) for r in R_new]))
    return C_new, V_new, sum(el for r in V_new for el in r)


def objective(ci, muti, sample, bb):
    if bb is None:
        form = (lambda mut : (mut.assigned_config.cf_to_v(ci, sample), mut.a[sample] + 1, (mut.d[sample] - mut.a[sample]) + 1))
    else:
        form = (lambda mut : (mut.assigned_config.cf_to_v(ci, sample), mut.a[sample], mut.d[sample] - mut.a[sample], bb[sample]))
    pdfs = compute_pdfs(*zip(*[form(m) for m in muti]), check=True)
    return -EPSILON if np.any(np.isneginf(pdfs)) else -np.sum(pdfs)


def update_objs(C, mutations, num_samples, num_clusters, bb):
    iobj = (lambda sam, muti, c : objective(c, muti, sam, bb) if len(muti) > 0 else 0.0)
    comp = (lambda sam, i : iobj(sam, list(filter(lambda m : m.assigned_cluster == i, mutations)), C[sam][i]))
    return [[comp(sam, i) for i in range(num_clusters)] for sam in range(num_samples)]


def compute_pdfs(_VS, _AS, _BS, _BB=None, check=False):
     M = np.array(list(map(lambda v : False if v is False else True, _VS)))
     if _BB is None:
         VS, AS, BS = map(np.array, [_VS, _AS, _BS])
     else:
         VS, AS, BS, BB = map(np.array, [_VS, _AS, _BS, _BB])
     VS[VS < SEQERROR] = SEQERROR
     VS[VS > (1.0 - SEQERROR)] = 1.0 - SEQERROR
     RS = np.full(M.shape[0], np.NINF)
     if _BB is None:
        """
        x = beta.logpdf(VS[M], AS[M], BS[M])
        y = VS[M]
        print np.ndim(x)
        print str(x)
        print str(y)
        if np.ndim(x) > 1:
            tmp = str(y) + "\n" + str(x)
            sys.exit(tmp)
        """
        RS[M] = beta.logpdf(VS[M], AS[M], BS[M])
     else:
         NS = AS[M] + BS[M]
         KS = AS[M]
         CS = gammaln(NS + 1) - gammaln(NS - KS + 1) - gammaln(KS + 1) #comb(NS, KS)
         NSminusKS = NS - KS
         BAFS = VS[M]
         minusBAFS = 1 - BAFS
         SS = BB[M]
         RS[M] = CS + betaln(KS + SS * BAFS, NSminusKS + SS * minusBAFS) - betaln(SS * BAFS, SS * minusBAFS)
     RS[M & (RS < EPSILON)] = EPSILON
     if check:
         assert np.NINF not in RS[M]
     return RS

