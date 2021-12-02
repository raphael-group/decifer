import os
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
import subprocess


this_dir = os.path.dirname(__file__)

def test_end2end():
    mutations = f"{this_dir}/data/decifer_input.tsv"
    purities = f"{this_dir}/data/purity.tsv"
    segfile = f"{this_dir}/data/segfile.tsv"
    snpfile = f"{this_dir}/data/snpfile.tsv"
    k = "6"
    K = "6"
    seed = "17"
    j = "2"
    test_out = this_dir + "/decifer"
    out = subprocess.run(["decifer", mutations, 
                    "-p", purities, 
                    "--betabinomial",
                    "--segfile", segfile,
                    "--snpfile", snpfile,
                    "-k", k,
                    "-K", K,
                    "--seed", seed,
                    "-j", j,
                    "-o", test_out,
                    "--debug"])
    """
    test_df = pd.read_csv(f"{test_out}.output.tsv", "\t").sort_values(by=['mut_index'])
    truth_df = pd.read_csv(f"{this_dir}/data/decifer.output.tsv", "\t").sort_values(by=['mut_index'])
    assert_frame_equal(test_df.reset_index(drop=True), 
                        truth_df.reset_index(drop=True),
                        check_less_precise=3)
    """

    test_df = pd.read_csv(f"{test_out}_clusterCIs.tsv", "\t", skiprows=3)
    truth_df = pd.read_csv(f"{this_dir}/data/decifer_clusterCIs.tsv", "\t", skiprows=3)
    assert_frame_equal(test_df, truth_df)
