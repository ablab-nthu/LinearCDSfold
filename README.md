# LinearCDSfold (version 1.1.0)

LinearCDSfold is a tool for designing coding sequences by jointly optimizing their secondary structure stability and codon usage preferences. 
It supports exact, beam, and Pareto-optimal search. 

## Installation Instructions

### Prerequisites

* Compiler: `gcc/g++` (version 7.0 or higher recommended)
* Build Tool: `make`

### Building from Source

#### 1. Clone the repository:

```
git clone https://github.com/ablab-nthu/LinearCDSfold.git
cd LinearCDSfold
```

#### 2. Compile the project:

```
make
```

#### 3. Verify Installation:

```
./LinearCDSfold -h
```

## Usage Instructions

```
./LinearCDSfold -o <output_file> [other options] <input_file>
```

### Input File

The `input_file` is an amino acid sequence file in FASTA format, and an example is shown below.

```
> P15421
MYGKIIFVLLLSGIVSISASSTTGVAMHTSTSSSVTKSYISSQTNGITLINWWAMARVIFEVMLVVVGMIILISYCIR
```

### Output File

The `-o` option (lowercase o) specifies the output file name (`output_file`) containing LinearCDSfold results in plain text. If `-o` is not provided, the default output file is `result.txt`.

### Other Options

#### 1. Search Mode:

```
-m <search_mode> or --mode <search_mode>
```

`search_mode` selects the search strategy in LinearCDSfold (default: `exact`):

* `exact`: performs an exact search for the optimal CDS.
* `beam`: uses beam search for fast approximate CDS generation.
* `pareto`: generates Pareto-optimal CDSs via `exact` search with the DERNA objective.

If no mode is specified, the default is `exact`.

#### 2. Objective Function:

```
-O <objective_function>
```

The `-O` option (uppercase letter O) specifies the objective function for joint optimization of MFE (Minimum Free Energy) and CAI (Codon Adaptation Index). The parameter `objective_function` must be set to either `LD` (default) or `DN`.

* `LD`: LinearDesign objective function `MFE − lambda × l × log(CAI)`, where `l` is the amino acid sequence length and `lambda` is a real value ranging from 0 to ∞.

* `DN`: DERNA objective function `lambda × MFE − (1 − lambda) × l × log(CAI)`, where `lambda` is a real value ranging from 0 to 1.

The parameter `lambda` is used to balance the relative contribution of MFE and CAI.

#### 3. Codon Usage File:

```
-c <codon_usage_file> 
```

`codon_usage_file` is a CSV file containing three fields: codon triplet, amino acid, and relative frequency. The first line may begin with `#` to indicate a comment. An example is provided below.

```
# codon triplet, amino acid, relative frequency
GCU,A,0.26
GCC,A,0.4
GCA,A,0.23
GCG,A,0.11
...
```

If `-c` is not specified, `codon_usage_file` defaults to `codon_usage_freq_table_human.csv`.

#### 4. Lambda Value:

```
-l <lambda>
```

`lambda` is a non-negative real-valued scaling parameter that balances MFE and CAI in the optimization objective.

* With LinearDesign objective (`-O LD`): 
    * Default: `lambda = 0`. 
    * Weight: MFE = 1 and CAI = `lambda`. 
    * Effect: Larger `lambda` values increase CAI contribution.

* With DERNA objective (`-O DN`): 
    * Default: `lambda = 1`. 
    * Weight: MFE = `lambda` and CAI = 1 - `lambda`. 
    * Effect: Larger `lambda` values increase MFE contribution.

#### 5. Beam Size:

```
-b <beam_size>
```

`beam_size` is a positive integer specifying the beam size used in beam search (`-m beam`) and its default value is 500.

#### 6. Termination Thresholds (Pareto-optimal search):

```
-t <tau1>   or   --tau1 <tau1>
-u <tau2>   or   --tau2 <tau2>
```

Both `tau1` and `tau2` are termination thresholds for Pareto-optimal search (`-m pareto`).

* `tau1` (default: 0.0025): smaller values generate more Pareto-optimal CDSs but increase runtime.
* `tau2` (default: 0.00075): must be smaller than `tau1`; smaller values explore more CDSs with higher runtime.

#### 7. Additional Output File:

```
-f <FILE_NAME>
```

`-f` specifies an additional CSV file containing only MFE and CAI results for each `lambda`. The file should use the .csv extension. If `-f` is not specified, the default output file is `result.csv`.

## Examples

### 1. Exact search using LinearDesign objective function:

```
> ./LinearCDSfold -m exact -O LD -l 3 -o P15421_exact_LD.txt example/P15421.fasta
```

Output: `cat P15421_exact_LD.txt`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: LinearDesign
Search mode: Exact search
Lambda: 3.000
Processing: [==================================================]  100%
Coding sequence and its secondary structure:
AUGUAUGGCAAGAUCAUCUUUGUGCUGCUGCUGAGCGGGAUUGUGUCAAUCUCCGCCAGCAGCACCACAGGGGUGGCCAUGCACACCUCCACCAGCAGCAGCGUCACCAAGAGCUACAUCUCUAGCCAGACCAAUGGCAUCACCCUGAUCAACUGGUGGGCCAUGGCUAGAGUGAUCUUCGAGGUGAUGCUGGUGGUGGUGGGGAUGAUCAUCCUGAUCAGCUACUGCAUCAGG
.((((((((.....((((((((((.((((((((.(((((((((...))))).)))))))))))).))))))))))))))))))((.((((((((.((.((((((((((..(((...((.(((((((((......((((.(((((..........))))))))))))))))))))...)))..)))))))))).)).)))))))).)).....(((((((((...))).))))))
Folding free energy: -126.600 kcal/mol
CAI: 0.947
Total runtime: 2.907 s
```

### 2. Exact search using DERNA objective function:

```
> ./LinearCDSfold -m exact -O DN -l 0.002 -o P15421_exact_DN.txt example/P15421.fasta
```

Output: `cat P15421_exact_DN.txt`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: DERNA
Search mode: Exact search
Lambda: 0.002
Processing: [==================================================]  100%
Coding sequence and its secondary structure:
AUGUAUGGCAAGAUCAUCUUUGUGCUGCUGCUGAGCGGGAUCGUGAGCAUCUCCGCCAGCAGCACCACAGGGGUGGCCAUGCACACCUCCACCAGCAGCAGCGUCACCAAGAGCUACAUCUCUAGCCAGACCAAUGGCAUCACCCUGAUCAACUGGUGGGCCAUGGCUAGAGUGAUCUUCGAGGUGAUGCUGGUGGUGGUGGGGAUGAUCAUCCUGAUCAGCUACUGCAUCAGG
.((((((((.....((((((((((.((((((((.((((((((....).))).)))))))))))).))))))))))))))))))((.((((((((.((.((((((((((..(((...((.(((((((((......((((.(((((..........))))))))))))))))))))...)))..)))))))))).)).)))))))).)).....(((((((((...))).))))))
Folding free energy: -124.300 kcal/mol
CAI: 0.956
Total runtime: 2.934 s

```

### 3. Beam search using LinearDesign objective function:

```
> ./LinearCDSfold -m beam -O LD -l 3 -b 500 -o P15421_beam_LD.txt example/P15421.fasta
```

Output: `cat P15421_beam_LD.txt`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: LinearDesign
Search mode: Beam search
Beam size: 500
Lambda: 3.000
Processing: [==================================================]  100%
Coding sequence and its secondary structure:
AUGUAUGGCAAGAUCAUCUUUGUGCUGCUGCUGAGCGGAAUUGUGAGCAUUUCCGCCAGCAGCACCACAGGGGUGGCCAUGCACACCUCCACCAGCAGCAGCGUCACCAAGAGCUACAUCUCUAGCCAGACCAAUGGCAUCACCCUGAUCAACUGGUGGGCCAUGGCUAGAGUGAUCUUCGAGGUGAUGCUGGUGGUGGUGGGGAUGAUCAUCCUGAUCAGCUACUGCAUCAGG
.((((((((.....((((((((((.((((((((.((((((.((....)).)))))))))))))).))))))))))))))))))((.((((((((.((.((((((((((..(((...((.(((((((((......((((.(((((..........))))))))))))))))))))...)))..)))))))))).)).)))))))).)).....(((((((((...))).))))))
Folding free energy: -125.700 kcal/mol
CAI: 0.949
Total runtime: 0.702 s
```

### 4. Beam search using DERNA objective function:

```
> ./LinearCDSfold -m beam -O DN -l 0.002 -b 500 -o P15421_beam_DN.txt example/P15421.fasta
```

Output: `cat P15421_beam_DN.txt`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: DERNA
Search mode: Beam search
Beam size: 500
Lambda: 0.002
Processing: [==================================================]  100%
Coding sequence and its secondary structure:
AUGUAUGGCAAGAUCAUCUUUGUGCUGCUGCUGAGCGGAAUUGUGAGCAUCUCCGCCAGCAGCACCACAGGGGUGGCCAUGCACACCUCCACCAGCAGCAGCGUCACCAAGAGCUACAUCUCUAGCCAGACCAAUGGCAUCACCCUGAUCAACUGGUGGGCCAUGGCUAGAGUGAUCUUCGAGGUGAUGCUGGUGGUGGUGGGGAUGAUCAUCCUGAUCAGCUACUGCAUCAGG
.((((((((.....((((((((((.((((((((.(((((..((....))..))))))))))))).))))))))))))))))))((.((((((((.((.((((((((((..(((...((.(((((((((......((((.(((((..........))))))))))))))))))))...)))..)))))))))).)).)))))))).)).....(((((((((...))).))))))
Folding free energy: -124.500 kcal/mol
CAI: 0.953
Total runtime: 0.705 s
```

### 5. Pareto-optimal search using default termination thresholds:

```
./LinearCDSfold -m pareto -o P15421_pareto.txt -f P15421_pareto.csv example/P15421.fasta
```

Output: `cat P15421_pareto.csv`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: DERNA
Search mode: Pareto-optimal search
Lambda: 0.000
Coding sequence and its secondary structure:
AUGUACGGCAAGAUCAUCUUCGUGCUGCUGCUGAGCGGCAUCGUGAGCAUCAGCGCCAGCAGCACCACCGGCGUGGCCAUGCACACCAGCACCAGCAGCAGCGUGACCAA
GAGCUACAUCAGCAGCCAGACCAACGGCAUCACCCUGAUCAACUGGUGGGCCAUGGCCCGGGUGAUCUUCGAGGUGAUGCUGGUGGUGGUGGGCAUGAUCAUCCUGAUCA
GCUACUGCAUCCGG
(((((.(((..((((((((((((((((((((((.(..((...(((.((((..(((((.(........).)))))....)))).)))..)).)))))))))))))))...)
)).....((((((.(((..((((.(((((((((((.((((.(((..(((((....))))))))))))...).)))))))))).)))))))..)).)))).....))))).
)))..)))))....
Folding free energy: -88.000 kcal/mol
CAI: 1.000
Total runtime: 2.911 s
Lambda: 1.000
Coding sequence and its secondary structure:
AUGUAUGGCAAGAUCAUCUUCGUCCUGCUGCUCUCCGGGAUCGUGUCGAUCUCGGCGAGCAGCACGACGGGGGUGGCCAUGCAUACGAGUACCAGCAGUAGCGUGACUAA
GAGUUAUAUAUCCUCACAGACCAACGGCAUCACCUUGAUAAAUUGGUGGGCGAUGGCCCGCGUAAUUUUCGAGGUGAUGCUGGUGGUCGUGGGGAUGAUAAUUCUUAUCA
GCUACUGCAUUCGU
(((((((((.....((((((((((.((((((((.(((((((((...))))))))).)))))))).)))))))))))))))))))((((((....((((((((.(((.(((
(((((((.(((((((((.(((((.((((((((((((((.((((((((((((....)))))).)))))))))))))))))))).)))))))))))))))))))))))))))
))))))))))))))
Folding free energy: -148.700 kcal/mol
CAI: 0.734
Total runtime: 2.912 s
...
```

Output: `cat P15421_pareto.csv`

```
lambda,MFE,CAI
0.000010,-88.000000,1.000000
0.999990,-148.700000,0.733673
0.500000,-148.700000,0.733673
...
```

## Contact Information

Corresponding author: Prof. Chin Lung Lu (Email: cllu@cs.nthu.edu.tw)
