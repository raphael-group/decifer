# Complete demo for DeCiFer on prostate cancer patient A12
: ex: set ft=markdown ;:<<'```shell' #

The following DeCiFer demo represents a guided example of running DeCiFer on somatic mutations inferred from prostate cancer patient A12 from [Gundem et al., Nature, 2015](https://doi.org/10.1038/nature14347).
The exemplary input files are already available from [DeCiFer data repository](https://github.com/raphael-group/decifer-data), where copy numbers as well as tumour purity in every sample have been inferred using [HATCHet](https://doi.org/10.1038/s41467-020-17967-y).
From this directory, simply run this file through BASH as a standard script to run the complete demo.
The demo can also be considered as a guided example of a complete execution and is correspondingly commented.

## Requirements and set up

The demo requires that DeCiFer has been succesfully installed with conda. If the custom installation was used, please make sure that you can succesfully run the command `decifer`. The demo includes the downloading of all the required files and will terminate in <20 minutes on machine with minimum requirements satisfied.

We gurantee that the running directory in the same directory of the demo and we remove previous results.

```shell
cd $( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
:<<'```shell' # Ignore this line
```

We also ask the demo to terminate in case of errors and to print a trace of the execution by the following commands
```shell
set -e
set -o xtrace
PS4='[\t]'
:<<'```shell' # Ignore this line
```

## Downloading of data

The demo auomatically downloads the required barcoded single-cell and matched-normal BAM files in `data` folder.

```shell
# Creating data folder
mkdir -p data

# Downloading mutation input file
curl -L 'https://raw.githubusercontent.com/raphael-group/decifer-data/main/input/prostate/mutations/A12.decifer.input.tsv' > data/mutations.tsv

# Downloading purity input file
curl -L 'https://raw.githubusercontent.com/raphael-group/decifer-data/main/input/prostate/purity/A12.purity.txt' > data/purity.tsv
:<<'```shell' # Ignore this line
```

## Run DeCiFer

We now run the command `decifer`, by providing the required input files and by fixing the number of clusters between 5 and 8 as well as the number or restarts to 20 for the sake of testing in short time. Note that we also fix the random seed to 17 for deterministic reproducibility.

```shell
decifer data/mutations.tsv -p data/purity.tsv -k 5 -K 8 -r 20 --seed 17
exit $?
```
