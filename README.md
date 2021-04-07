# DeCiFer

DeCiFer is an algorithm that simultaneously selects mutation multiplicities and clusters SNVs by their corresponding descendant cell fractions (DCF), a statistic that quantifies the proportion of cells which acquired the SNV or whose ancestors acquired the SNV. DCF is related to the commonly used cancer cell fraction (CCF) but further accounts for SNVs which are lost due to deleterious somatic copy-number aberrations (CNAs), identifying clusters of SNVs which occur in the same phylogenetic branch of tumour evolution.

The full description of the algorithm and its application on published cancer datasets are described in

[Gryte Satas†, Simone Zaccaria†, Mohammed El-Kebir†,\* and Ben Raphael\*, 2021](https://doi.org/10.1101/2021.02.27.429196)\
† Joint First Authors\
\* Corresponding Authors

The results of the related paper are available at:

[DeCiFer data](https://github.com/raphael-group/decifer-data)

This repository includes detailed instructions for installation and requirements, demos and tutorials of DeCiFer, a list of current issues, and contacts.
This repository is currently in a preliminary release and improved versions are released frequently.
During this stage, please keep checking for updates.

## Contents ##

1. [Algorithm](#algorithm)
2. [Installation](#installation)
3. [Usage](#usage)
    - [Required input data](#requireddata)
    - [Optional input data](#optionaldata)
    - [Output](#output)
    - [System requirements](#requirements)
    - [Demos](#demos)
    - [Reccomendations and quality control](#reccomendations)
4. [Development](#development)
5. [Contacts](#contacts)

<a name="algorithm"></a>
## Algorithm

<img src="doc/decifer.png" width="500">

DeCiFer uses the Single Split Copy Number (SSCN) assumption and evolutionary constraints to enumerate potential genotype sets.
This allows DeCiFer to exclude genotype sets with constant-mutation multiplicity (CMM) that are not biologically likely (red crosses) and include additional genotype sets (green star) that are.
DeCiFer simultaneously selects a genotype set for each SNV and clusters all SNVs based on a probabilistic model of DCFs, which summarize both the prevalence of the SNV and its evolutionary history.

<a name="installation"></a>
## Installation

DeCiFer is mostly written in Python 2.7 and has an optional component in C++. The recommended installation is through conda but we also provide custom instructions to install DeCiFer in any Python evironment.

### Automatic installation 

The recommended installation is through [bioconda](https://bioconda.github.io/recipes/decifer/README.html) and requires `conda`, which can be easily and locally obtained by installing one of the two most common freely available distributions: [anaconda](https://www.anaconda.com/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). Please make sure to have executed the [required channel setup](https://bioconda.github.io/user/install.html#set-up-channels) for bioconda. Thus, the following one-time commands are sufficient to fully install DeCiFer within a virtual conda envirnoment called `decifer`:

```shell
conda install -n decifer decifer -y
```

After such one-time installation, DeCiFer can be executed in every new session after activating the `decifer` environment as follows:

```shell
conda activate decifer
```

### Manual installation

DeCiFer can also be installed in a conda envirnoment directly from this repo. Thus, the following one-time commands are sufficient to fully install DeCiFer within a virtual conda envirnoment called `decifer` from this Git repo:

```shell
git clone https://github.com/raphael-group/decifer.git && cd decifer/
conda create -c anaconda -c conda-forge -n decifer python=2.7 numpy scipy matplotlib-base pandas seaborn -y
pip install .
```

### Custom installation

DeCiFer can be installed with `pip` by the command `pip install .` in any Python2.7 envirnoment with the following packages or compatible versions:

| Package | Tested version | Comments |
|---------|----------------|----------|
| [numpy](https://numpy.org/) | 1.16.1 | Efficient scientific computations |
| [scipy](https://www.scipy.org/) | 1.2.1 | Efficient mathematical functions and methods |
| [pandas](https://pandas.pydata.org/) | 0.20.1 | Dataframe management |
| [matplotlib](https://matplotlib.org/) | 2.0.2 | Basic plotting utilities |
| [seaborn](https://seaborn.pydata.org/) | 0.7.1 | Advanced plotting utilities |


### Installation of C++ component

DeCiFer includes C++ code to enumerate state/genotype trees. The
dependencies for this code are as follows.

| Package | Tested version | Comments |
|---------|----------------|----------|
| [cmake](https://cmake.org/) | >= 2.8 | Build environment |
| [lemon](https://lemon.cs.elte.hu/trac/lemon) | 1.3.1 | C++ graph library |
| [boost](https://www.boost.org/) | >= 1.69.0 | C++ library for scientific computing |

To build this code, enter the following commands from the root of the
repository:

```shell
mkdir build
cd build
# OPTIONAL: specify lemon and/or Boost paths if not detected automatically.
cmake ../src/decifer/cpp/ -LIBLEMON_ROOT=/usr/local/ -DBOOST_ROOT=/scratch/software/boost_1_69_0/
make
```

<a name="usage"></a>
## Usage

DeCiFer can be executed using the command `decifer`, whose [manual](man/man-decifer.md) describes the available parameters and argument options. See more details below.

1. [Required input data](#requireddata)
2. [Optional input data](#optionaldata)
3. [Output](#output)
4. [System requirements](#requirements)
5. [Demos](#demos)
6. [Reccomendations and quality control](#reccomendations)

<a name="requireddata"></a>
### Required input data

DeCiFer requires two input data:

1. Input mutations with nucleotide counts and related copy numbers in a tab-separated file (TSV) with three header lines ((1) The first specifies the number of mutations; (2) The second specifies the number of samples; and (3) The third is equal to: `#sample_index	sample_label	character_index	character_label	ref	var`) and where every other row has the following values for every mutation in every sample:

| Name | Description | Mandatory |
|------|-------------|-----------|
| Sample index | a unique number identifying the sample | Yes |
| Sample label | a unique name for the sample | Yes |
| Mutation index | a unique number identifying the mutation | Yes |
| Mutation label | a unique name identifying the mutation | Yes |
| REF | Number of reads with reference allele for the mutation | Yes |
| ALT | Number of reads with alternate allele for the mutation | Yes |
| Copy numbers and proportions | Tab-separated `A  B  U` where `A,B` are the inferred allele-specific copy numbers for the segment harboring the mutationa and `U` is the corresponding proportion of cells (normal and tumour) with those copy numbers | Yes |
| Additional copy numbers | An arbitrary number of fields with the same format as of `Copy numbers and proportions` describing the proportions of cells with different copy numbers. Note that all proporions should always sum up to 1. | No |

2. Input tumour purity in a two-column tab-separated file where every row `SAMPLE-INDEX   TUMOUR-PURITY` defines the tumour purity `TUMOUR-PURITY` of a sample with index `SAMPLE-INDEX`.

<a name="optionaldata"></a>
### Optional input data

DeCiFer can use the following additional and optional input data:

1. Precision parameters for fitting beta-binomial distributions when clustering mutations. Specifically, this is a two-column tab-separated file where every row `SAMPLE-INDEX   PRECISION` defines the precision parameter `PRECISION` of the beta binomial distribution of a sample with index `SAMPLE-INDEX`. When this optional file is not provided, DeCiFer will fit a binomial distribution instead. This file can be estimated by using the command `decifer_fit`, whose usage is described in the corresponding [manual](man/man-decifer-fit.md).

2. File containing the set of all possible state trees evaluated by DeCiFer. State trees have been generated for the set of most common copy numbers, however a dataset might have a combination of copy numbers which has not been included. In this case, the user can use the command `decifer_statetrees` to generate all the state trees needed for its dataset, following the instructions in the corresponding [manual](man/man-decifer-manual.md).

<a name="output"></a>
### Output

DeCiFer's output corresponds to a single TSV file encoding a dataframe where every row corresponds to an input mutation and with the following fields:

| Name | Description |
|------|-------------|
| `mut_index` | Unique identified for a mutation |
| `VAR_{SAMPLE}` | Variant sequencing read count of the mutation for every sample with index `{SAMPLE}` |
| `TOT_{SAMPLE}` | Total sequencing read count of the mutation for every sample with index `{SAMPLE}` |
| `VAR_{SAMPLE}` | Variant sequencing read count of the mutation for every sample with index `{SAMPLE}` |
| `cluster` | Unique identifier of the inferred mutation cluster |
| `state_tree` | Inferred state tree defined as a `->`-separated edge list of genotypes |
| `cluster` | Unique identifier of the inferred mutation cluster |
| `true_cluster_DCF{SAMPLE}` | Inferred true cluster DCF of the mutation in every sample with index `{SAMPLE}`; when execute in CCF-mode, DCF will be CCF instead |
| `point_estimate_DCF{SAMPLE}` | Point estimate of the mutation DCF in every sample with index `{SAMPLE}`; when execute in CCF-mode, DCF will be CCF instead |
| `cmm_CCF{SAMPLE}` | Inferred CCF of the mutation under the previous CMM assumption in every sample with index `{SAMPLE}` |
| `Explained` | `;`-separated list of all the clusters to which the mutation could be assigned |
| `LHs` | `;`-separated list of the negative-log likelihoods of assigned the mutation to all clusters in `Explained` |

<a name="requirements"></a>
### System requirements

DeCiFer is highly parallelized in order to make efficient the extensive computations needed for clustering under a probabilistic model thousands of mutations across multiple tumour samples from the same patient. We recommend executing DeCiFer on multi-processing computing machines as the running time will scale down nearly proportionally with the number of parallel jobs, which can be specified with the argument `-j`. If the parameter is not specified, then DeCiFer will attempt to use all available CPUs; however, when using a computing cluster, we strongly recommend the user to always specifies `-j` in order to match the number of requested CPUs and avoid computing competition. Finally, note that also required memory also scales with the number of parallel processes; however in all previous tests on thousands of mutations with high number of parallel processses, DeCiFer never required more than 80GB of RAM. Please lower `-j` in case of exceeding memory.

<a name="demos"></a>
### Demos

Each demo is an exemplary and guided execution of a DeCiFer.Each demo is simultaneously a guided description of the entire example and a BASH script which can be directly executed to run the complete demo from this repository. As such, the user can both read the guided description as a web page and run the same script to execute the demo. At this time the following demos are available (more demos will be available soon):

| Demo | Description |
|------|-------------|
| [A12](demos/demo-A12.sh) | Demo of DeCiFer basic command on prostate cancer patient A12 |

<a name="reccomendations"></a>
### Reccomendations and quality control

- We reccomend to initially select a reasonably high maximum number of clusters with option `-K` and then further increase it if the selected best number of clusters is close to the maximum limit.
- DeCiFer outputs the decreasing objective function which is used to select the number of clusters based on the Elbow criterion; if the function is still substantially decreasing near the selecte maximum number of clusters please try to further increase this value.
- You can adapt the sensitivity of the Elbow criterion by adjustiv the corresponding parameter.

<a name="development"></a>
## Development

DeCiFer is in active development, please report any issue or question as this could help the development and improvement of DeCiFer. Known issues with current version are reported here below.

<a name="contacts"></a>
## Contacts

DeCiFer has been developped and actively mantained by three previous Postdoctoral Research Associates at Princeton University in the research group of prof. Ben Raphael:

- [Simone Zaccaria](https://simozacca.github.io/) (most of algorithmic core, likelihood computation, and Python implementation), group leader of the [Computational Cancer Genomics research group](https://www.ucl.ac.uk/cancer/zaccaria-lab) at UCL Cancer Institute, London, UK

- [Gryte Satas](linkedin.com/in/gryte-satas-23a74844) (most of input/output interfact and CF computations for single mutations), Postdoctoral Research Fellow at Memorial Sloan Kettering Cancer Center, NY, USA

- [Mohammed El-Kebir](http://www.el-kebir.net/) (generation of state trees), Assistant Professor in the Computer Science department at the University of Illinois at Urbana-Champaign (UIUC), IL, USA.
