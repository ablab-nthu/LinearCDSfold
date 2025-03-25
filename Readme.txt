LinearCDSfold Software

Description:
LinearCDSfold is a tool for designing a coding sequence by jointly optimizing its secondary structure stability and codon usage.

Download:
Please click this link (http://genome.cs.nthu.edu.tw/LinearCDSfold) to download the zipped file of LinearCDSfold.

To Compile:
Run the following command to compile LinearCDSfold.

make

To Run:
Use the following command to run LinearCDSfold.

./LinearCDSfold [OPTIONS] SEQFILE

SEQFILE:
SEQFILE is an amino acid sequence file in FASTA format. The following is an example of SEQFILE.

> 72_1175.fa
MGQSKSKEEKGISGTSRAEILPDTTYLGPLNCKSCWQKFDSFSKCHDHYLCRHCLNLLLTSSDRCPLCKYPL

OPTIONS:

-b <BEAMSIZE>
BEAMSIZE is a parameter of beam size used by LinearCDSfold when running with beam pruning (the default of BEAMSIZE is 500). To run LinearCDSfold without beam pruning, use the -b 0 option.

-cai <CAIFILE>
CAIFILE is a file of codon usage frequencies, which contains contains the relative adaptiveness values of all the 64 codons arranged in the following order: AUU, AUA, AUC, CUA, CUC, CUG, CUU, UUA, UUG, GUU, GUA, GUC, GUG, UUU, UUC, AUG, UGU, UGC, GCA, GCC, GCG, GCU, GGU, GGC, GGA, GGG, CCU, CCC, CCA, CCG, ACU, ACC, ACA, ACG, UCU, UCC, UCA, UCG, AGU, AGC, UAU, UAC, UGG, CAA, CAG, AAU, AAC, CAU, CAC, GAA, GAG, GAU, GAC, AAA, AAG, CGU, CGC, CGA, CGG, AGA, AGG, UAA, UAG, UGA. The following is an example of CAIFILE.

1.00, 0.59, 0.57, 0.48, 0.21, 0.38, 0.45, 0.97, 1.00, 1.00, 0.54, 0.54, 0.49, 1.00, 0.69, 1.00, 1.00, 0.59, 0.76, 0.58, 0.29, 1.00, 1.00, 0.40, 0.47, 0.26, 0.74, 0.36, 1.00, 0.29, 1.00, 0.63, 0.86, 0.40, 1.00, 0.62, 0.81, 0.38, 0.62, 0.42, 1.00, 0.79, 1.00, 1.00, 0.45, 1.00, 0.69, 1.00, 0.56, 1.00, 0.43, 1.00, 0.54, 1.00, 0.72, 0.29, 0.13, 0.15, 0.08, 1.00, 0.44, 1.00, 0.49, 0.64

Note that the default of CAIFILE is yeast_relative_adaptiveness.txt.

-l <LAMBDA>
LAMBDA is a scaling parameter of a non-negative real number to balance the contributions of MFE (minimum free energy) of a coding sequence and its CAI (Codon Adaptation Index) value to the joint optimization objective. If LAMBDA is 0, then only MFE is considered in the joint optimization objective; otherwise, both MFE and CAI are considered and the larger the value of LAMBDA, the greater the contribution of CAI. The default of LAMBDA is 3.

Example:
> ./LinearCDSfold -b 200 -cai human_relative_adaptiveness.txt -l 2 example/72_1175.fa
Amino acid file: example/72_1175.fa
Beam size: 200
Codon usage table: human_relative_adaptiveness.txt
Lambda: 2
Coding sequence and its secondary structure:
AUGGGUCAGAGCAAGUCCAAGGAGGAGAAGGGCAUCAGCGGCACCAGUCGGGCAGAGAUUCUGCCCGACACCACCUACCUGGGGCCGCUGAAUUGCAAGAGCUGUUGGCAGAAGUUCGACAGCUUUAGCAAGUGCCACGACCACUACCUGUGUCGGCACUGCCUGAAUCUCCUCCUGACUAGCUCUGACCGGUGCCCGCUGUGCAAGUACCCCCUG
...(((((((((.((((..(((((((((.(((((((((((((.(((((((((((((...))))))))))..........))).)))))))).((((.((((((((((((....)).)))))))))).))))(((((..((.(((.....)))))))))))))))...))))))))))))).)))))))))(((((..((...))..))))).....
Folding free energy: -114.7 kcal/mol
CAI: 0.903

Contact Information:
Corresponding author: Prof. Chin Lung Lu (Email: cllu@cs.nthu.edu.tw)
