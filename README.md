# [SequenceLab](https://github.com/CMU-SAFARI/SequenceLab)

SequenceLab enables comprehensive benchmarks of computational methods for comparing genomic sequences.
SequenceLab is written in C, but compiled in using `g++` to increase compatibility with filters written in C++.  It was designed to be incredibly lightweight and only add minimal overhead to running the filters themselves.

SequenceLab includes implementations and comparisons of the following pre-alignment filtering and alignment algorithms, and new algorithms can be added easily if a C/C++ implementation is available:

0. Adjacency
1. Base Counting
2. Shouji
3. Q-Gram
4. MAGNET
5. Hamming Distance
6. GateKeeper
7. SneakySnake
8. GRIM
9. KSW2
10. Edlib

## Installation
Download SequenceLab and compile it using `g++`:
```bash
git clone https://github.com/CMU-SAFARI/SequenceLab && cd SequenceLab
make
```

## Obtaining Test Data
We provide `.tss` files, containing sequence pairs derived from real sequencing data.
They can be downloaded from [Zenodo](https://zenodo.org/records/10028978) using the following commands:
```bash
testdata/download_hifi.sh
testdata/download_ont.sh
testdata/download_illumina.sh
```

## Usage
At a high level, SequenceLab accepts a `.tss` file containing genomic sequence pairs as input and runs each filter on each pair. It then outputs a binary value for each filter and pair, indicating whether the filter accepted or rejected the pair. Acceptance or rejection depends on the edit distance threshold, specified as a parameter sweep.

For example, the following code runs all filters on the Illumina dataset, using the default parameters:
```bash
./main --file testdata/illumina/Illumina_chained_10M.tss
```
This will generate a `results` directory with three types of output files:
- `results/grid_X.csv`: A CSV file containing a line for each sequence pair. Each evaluated filter has a column, and entries are boolean values indicating if a given filter accepted or rejected a given sequence pair for the edit distance threshold `X` (in percent).
- `results/histogram.csv`: A histogram of the evaluated read error rates (in percent).
- `results/runtime.csv`: A CSV file containing the runtime of each filter for each edit distance threshold (average).
These files can be used to e.g. calculate the accuracy of each individual filter, combinations of filters, or their sensitivity to the edit distance threshold using simply pandas-based python scripts.

By default, SequenceLab runs each filter on only the first 100 sequence pairs. This value can (and for full experiment runs should) be overridden. E.g., to use all sequence pairs in `Illumina_chained_10M.tss`, use the following commands:
```bash
NUM_LINES=$(wc -l < testdata/illumina/Illumina_chained_10M.tss)
./main --file testdata/illumina/Illumina_chained_10M.tss --read-count $NUM_LINES
```

You can find a full list of parameters, their descriptions, and default values by running `./main --help`.


## Output Formats
The suite outputs three types of `.csv` files.

### Histogram
The histogram describes the edit distance distribution over the supplied sequence pairs. Suppose a pair has an edit distance $e$, with sequence lengths $s_1$ and $s_2$. Then the bucket with index $\lfloor\frac{e}{min(s_1, \; s_2)}\rfloor$ is incremented by one. In lay-terms, we count how often pairs of a given edit distance *percentage* occur.

If `full-edlib` is not set, the histogram only measures edit distances up to the `max-edit-distance` specified. Pairs with a higher edit distance are reported as having an edit distance of 100%.

The prefix – without the `.csv` ending – can be specified using the `--histogram-file-prefix` argument. The default file name is `histogram.csv`.

### Grid
For *each* edit distance threshold, we output a grid `.csv` that records the accept/reject decision for each filter and sequence pair. Additionally, we record the length of the read, the *absolute* threshold & Edlib's edit distance estimate.

An example grid file for a handful of sequences and filters is shown below. An entry is $1$, if the filter accepts the pair and $0$ if it rejects the pair. The edlib_estimate is either a positive number if the actual edit distance is below the error_threshold or $-1$ if it is above it.

| Shouji | Q-Gram | SneakySnake | GRIM | edlib | edlib_estimate | error_threshold | read_length |
|--------|--------|-------------|------|-------|----------------|-----------------|-------------|
| 0      | 0      | 0           | 1    | 0     | -1             | 10              | 151         |
| 1      | 0      | 0           | 0    | 0     | -1             | 10              | 151         |
| 0      | 0      | 0           | 0    | 0     | -1             | 10              | 151         |
| 0      | 0      | 0           | 0    | 0     | -1             | 10              | 151         |
| 1      | 0      | 0           | 0    | 0     | -1             | 10              | 151         |
| 1      | 0      | 0           | 0    | 0     | -1             | 10              | 151         |
| 1      | 0      | 0           | 1    | 0     | -1             | 10              | 151         |
| 0      | 0      | 0           | 0    | 0     | -1             | 10              | 151         |

This grid can then be used to determine false accepts, sensitivity, or any other binary metric. Furthermore, it can be used to understand which filters do similar work – i.e. discard the same sequence pairs.

The prefix – without the `.csv` ending – can be specified using the `--grid-prefix` argument. The default grid names are `grid_X.csv`, where `X` corresponds to the edit distance threshold in percent.

### Execution Time
This `.csv` records the execution time for each of the filters. The first column shows the edit distance in percent. The next set of columns show the total execution time for each filter in seconds. This can then be used to calculate the execution time per sequence pair for each filter.

The filename – **with** the `.csv` ending – can be specified using the `--runtime-filename` argument. The default file name is `runtime.csv`.

## TSS Format

The TSS (Tab Separated Sequences) format contains a pair of genomic sequences on each line, separated by tabs.
This simplified format enables evaluating genomic tools little overhead on real datasets.

### Specification

- Each line consists of a pair of genomic basepair sequences, separated by a tab character.
- Each line is terminated by a single newline character, i.e. in UNIX style. Windows style linebreaks (carriage return + newline) are *not* permitted.
- Basepair sequences may consist of uppercase and lowercase nucleic or amino acid codes, as allowed in the [FASTA format](https://zhanggroup.org/FASTA/).
- If the dataset is for a readmapping usecase, the first sequence is the read or query, the second is the reference or target.

Below is an example file in TSS format:
```
ACGTACGT        ACGTACGT
ACGTACGT        TGCATGCA
GTACGT          TGCATGCA
ACGTACGT        ATG
```

## ONT Parameters
There are a few *fast* filters on ONT reads and a few very slow ones (mostly because of their quadratic space complexity):

*Fast filters*:
- Base counting
- Q-gram
- Hamming Distance
- GRIM

*Slow filters:*
- Adjacency
- Shouji
- MAGNET
- GateKeeper
- SneakySnake
- KSW2
- Edlib

To run only the fast filters, use the following flag:
```bash
./main --activate-individually 0101010010001 [...]
```
The slow filters can be run individually by activating their corresponding index.

## Citing SequenceLab

If you use SequenceLab in your work, please cite:

> Maximilian-David Rumpf, Mohammed Alser, Arvid E. Gollwitzer, Joel Lindegger, Nour Almadhoun, Can Firtina, Serghei Mangul, Onur Mutlu. 
> "SequenceLab: A Comprehensive Benchmark of Computational Methods for Comparing Genomic Sequences" 
> (2023). https://doi.org/10.48550/arXiv.2310.16908)

You may use the following BibTeX:

```bibtex
@misc{rumpf2023sequencelab,
      title={SequenceLab: A Comprehensive Benchmark of Computational Methods for Comparing Genomic Sequences}, 
      author={Maximilian-David Rumpf and Mohammed Alser and Arvid E. Gollwitzer and Joel Lindegger and Nour Almadhoun and Can Firtina and Serghei Mangul and Onur Mutlu},
      year={2023},
      eprint={2310.16908},
      archivePrefix={arXiv},
      primaryClass={q-bio.GN}
}
```
