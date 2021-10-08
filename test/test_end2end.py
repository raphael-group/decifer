import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
import subprocess


def test_end2end():
    mutations = "test/data/mutations.tsv"
    purities = "test/data/purity.tsv"
    k = "5"
    K = "8"
    r = "20"
    seed = "17"
    j = "4"
    out = subprocess.run(["decifer", mutations, 
                    "-p", purities, 
                    "-k", k,
                    "-K", K,
                    "-r", r,
                    "--seed", seed,
                    "-j", j])
    test_df = pd.read_csv("decifer.output.tsv", "\t")
    truth_df = pd.read_csv("data/decifer.output.tsv", "\t")
    assert_frame_equal(test_df, truth_df)

    test_df = pd.read_csv("decifer.cluster.CIs.tsv", "\t", skiprows=3)
    truth_df = pd.read_csv("data/decifer.cluster.CIs.tsv", "\t", skiprows=3)
    assert_frame_equal(test_df, truth_df)
