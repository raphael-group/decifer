"""
config.py

author: gsatas
date: 2020-05-04
"""

from decifer.process_input import PURITY

THRESHOLD=0.05


class Config:
    def __init__(self, mut_state, other_states, cn_props, desc_set, dcf_mode = True):
        '''
        mut_state: 2-tuple that indicates the CN state that the
            mutation occurs in
        other_states: a list of 3-tuples indicating all other mutation
            states
        cn_props: a dict of lists, mapping 2-tuple CN states (keys) to proportions for each sample (list value)
        '''
        self.mut_state = mut_state
        self.other_states = other_states
        self.cn_props = cn_props
        self.desc_set = desc_set
        self.dcf_mode = dcf_mode

    def F(self, sample):
        '''returns fractional copy number F'''
        return self.cn_props[self.mut_state][sample] * sum(self.mut_state) \
            + sum([self.cn_props[s[:2]][sample] * sum(s[:2]) for s in self.other_states])

    def cf_bounds(self, sample):
        if self.dcf_mode:
            return self.d_bounds(sample)
        else:
            return self.c_bounds(sample)

    def c_bounds(self, sample):
        M = sum([self.cn_props[s[:2]][sample] for s in self.other_states if s[2] > 0])
        return M, M + self.cn_props[self.mut_state[:2]][sample]

    def d_bounds(self, sample):
        M = sum([self.cn_props[s[:2]][sample] for s in self.other_states if s in self.desc_set])
        return M, M + self.cn_props[self.mut_state[:2]][sample]

    #def M(self):
    #    return sum([self.cn_props[s[:2]] * (s[2]-1) for s in self.other_states if s[2] > 0])

    def c(self, lam, sample):
        return self.cn_props[self.mut_state][sample] * lam \
            + sum([self.cn_props[s[:2]][sample] for s in self.other_states if s[2] >= 1])

    def v(self, lam, sample):
            F = self.F(sample)
            val = self.cn_props[self.mut_state][sample] * lam \
            + sum([self.cn_props[s[:2]][sample] * s[2] for s in self.other_states])
            return 1./F * val

    def d(self, lam, sample):
        # multiplies lam by SSCN CN proportion in which mutation arose, canceling out earlier division in v_to_lam
        return self.cn_props[self.mut_state][sample] * lam \
                + sum([self.cn_props[s[:2]][sample] for s in self.other_states if s in self.desc_set])

    def c_to_lam(self, c, sample):
            if self.cn_props[self.mut_state][sample] == 0: return 0
            lam = (c - sum([self.cn_props[s[:2]][sample] for s in self.other_states if s[2] >= 1]))/self.cn_props[self.mut_state][sample]
            return lam

    def v_to_lam(self, v, sample):
            if self.cn_props[self.mut_state][sample] == 0: return 0
            # sum term iterates across other_states in which m (of (x,y,m)) >=1
            # multiplies CN proportions for these states by m
            # and divides by the SSCN CN proportion in which mutation arose
            # note: other_states include any genotype (x,y,m) with CN state != mut_state
            lam = (v*self.F(sample) - sum([self.cn_props[s[:2]][sample]*s[2] for s in self.other_states if s[2] >= 1]))/self.cn_props[self.mut_state][sample]
            return lam

    def d_to_lam(self, d, sample):
            if self.cn_props[self.mut_state][sample] == 0: return 0
            lam = (d - sum([self.cn_props[s[:2]][sample] for s in self.other_states if s in self.desc_set]))/self.cn_props[self.mut_state][sample]
            return lam

    def v_to_cf(self, v, sample, truncate = True):
        # calls d_to_v or c_to_v depending on dcf_mode
        if self.dcf_mode:
            cf = self.v_to_d(v, sample, truncate)
        else:
            cf = self.v_to_c(v, sample, truncate)
        return min(max(cf, 0.0), PURITY[sample])

    def cf_to_v(self, c, sample, truncate = True):
        # calls d_to_v or c_to_v depending on dcf_mode
        if not (self.cf_bounds(sample)[0] - THRESHOLD <= c <= self.cf_bounds(sample)[1] + THRESHOLD):
            return False
        if self.dcf_mode:
            v = self.d_to_v(c, sample)
        else:
            v = self.c_to_v(c, sample)
        assert v is False or 0.0 <= v <= 1.0
        return v

    def c_to_v(self, c, sample):
        lam = self.c_to_lam(c, sample)
        v = self.v(lam, sample)
        if (lam > -THRESHOLD and lam < 1 + THRESHOLD): 
            if v < 0: return 0.0
            if v > 1: return 1.0
            return v
        else: return False

    def v_to_c(self, v, sample, truncate = True):
        # If truncrate is true, then it returns False when there is no feasible
        # ccf that would result in v. If truncate is False, then it returns the
        # nearest feasible ccf. The latter is used for visualization purposes mainly
        lam = self.v_to_lam(v, sample)
        c = self.c(lam, sample)
        if truncate:
            if (lam > -THRESHOLD and lam < 1 + THRESHOLD): 
                if c < 0: return 0.0
                if c > PURITY[sample]: return PURITY[sample]
                return c
            else: return False
        else:
            return c

    def d_to_v(self, d, sample):
        lam = self.d_to_lam(d, sample)
        v = self.v(lam, sample)
        if (lam > -THRESHOLD and lam < 1 + THRESHOLD): 
            if v < 0: return 0.0
            if v > 1: return 1.0
            return v
        else: return False

    def v_to_d(self, v, sample, truncate = True):
        # If truncrate is true, then it returns False when there is no feasible
        # ccf that would result in v. If truncate is False, then it returns the
        # nearest feasible ccf. The latter is used for visualization purposes mainly
        lam = self.v_to_lam(v, sample)
        d = self.d(lam, sample)
        if truncate:
            if (lam > -THRESHOLD and lam < 1 + THRESHOLD): 
                if d < 0: return 0.0
                if d > PURITY[sample]: return PURITY[sample]
                return d
            else: return False
        else:
            return d

    def add_sample(self, cn_props, mut_label = None):
        # NOTE, this assumes samples are always in the same order in the input file
        # and that all copy-number states are included for all samples even if the prop. is 0
        try:
            # ensure same number of CN states for this sample compared to previous samples
            assert(len(cn_props) == len(self.cn_props))
            # ensure the CN states are also identical
            assert(set(cn_props.keys()) == set(self.cn_props.keys()) )
        except AssertionError:
            print(cn_props)
            print(self.cn_props)
            raise Exception("The same copy-number states and proportions must be provided across samples for a mutation, even when the proportion is 0. Mutation {}".format(mut_label))

        # append CN proportions for this sample
        for c in cn_props:
            self.cn_props[c].append(cn_props[c])
