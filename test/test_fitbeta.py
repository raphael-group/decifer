import os
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
import subprocess


this_dir = os.path.dirname(__file__)

def test_fitbeta():
    seg_file = f"{this_dir}/data/fitbeta.best.seg.ucn"
    snp_file = f"{this_dir}/data/fitbeta_tumor.1bed"
    seed = "17"
    file_prefix = "decifer"
    test_out = f"{this_dir}/{file_prefix}"

    test_file = f"{this_dir}/{file_prefix}_betabinom.tsv"
    truth_file = f"{this_dir}/data/{file_prefix}_betabinom.tsv"

    out = subprocess.run(["fitbeta", 
                    "-i", snp_file, 
                    "-s", seg_file,
                    "--seed", seed,
                    "-o", test_out])

    test_df = pd.read_csv(test_file, "\t")
    truth_df = pd.read_csv(truth_file, "\t")
    assert_frame_equal(test_df, truth_df)
