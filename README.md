# GenHap

GenHap is a novel computational method based on Genetic Algorithms for haplotype assembly, developed to manage long reads like those produced by SMRT sequencing technologies (PacBio RS II).

GenHap tackles the computational complexity of the haplotyping problem by exploiting Genetic Algorithms to solve the weighted Minimum Error Correction problem (wMEC), a variant of the well-known MEC problem.
GenHap can efficiently solve large instances of the wMEC problem, yielding optimal solutions by means of a global search process, without any a priori hypothesis about the sequencing error distribution in reads. Moreover, GenHap is capable of processing datasets composed of both long reads and coverages up to 60x on personal computers.

  1. [Reference](#ref) 
  2. [Compilation](#comp) 
  3. [Input parameters](#inp)
  4. [Data](#data)
  5. [License](#lic)
  6. [Contacts](#cont)


## <a name="ref"></a>Reference ##

A detailed description of GenHap, as well as a complete experimental comparison with an efficient state-of-the-art algorithm for haplotype phasing by using the dataset described below ([Data](#data)), can be found in:

Tangherloni A., Spolaor S., Rundo L., Nobile M.S., Cazzaniga P., Mauri G., Liò P., Merelli I., and Besozzi D.: "GenHap: A Novel Computational Method Based on Genetic Algorithms for Haplotype Assembly", _submitted_.

## <a name="comp"></a>Compilation ##

GenHap has been developed and tested on Ubuntu Linux but should work also on MacOS X.

GeneHap depends on:
- Message Passing Interface (MPI),

which is used to distribute the optimizations, exploiting a Master-Slave distributed programming paradigm.
For unix like operating systems, we suggest to download `mpich` as high performance and widely portable implementation of MPI.

In order to compile GenHap, we provide a bash script (`compile.sh`).

The resulting executable file (named `GenHap`) is a standalone executable program.

## <a name="inp"></a>Input parameters ##

In oder to run GenHap, the following parameter must be provided:

- `-i` to specify the input file containing the reads in input (in WIF format).
  
Optional parameters could be provided:

- `-o` to specify the output folder;
- `-n` to indicate the file name containing the estimated haplotypes;
- `-g` to specify the number of reads composing each sub-matrix of each sub-problem;
- `-v` to enable or disable the verbose modality;
- `-s` to save or not save the files containing the partitions and fragment matrix;
- `-x` to enable or disable the masking of the ambiguous positions in the output haplotypes with a X;
- `-p` to specify the file containing the GA settings.

For example, GenHap can be executed with the following command:

    mpiexec -np 4 ./GenHap -i Models/Freq_100/10000SNPs_60x/PacBio/Model_0.wif -o Output/Freq_100/10000SNPs_60x/PacBio

which will save the estimated haplotypes solution in the file `haplotypes` into the folder Output/Freq_100/10000SNPs_60x/PacBio.
In this case, `GenHap` exploits 4 cores to perform the required optimizations.

By running `GenHap` without specifying any parameter, all the above parameters will be listed.

## <a name="data"></a>Data ##

GenHap takes as input a dataset in WIF format, introduced in Patterson, M. et al. "WhatsHap: weighted haplotype assembly for future-generation sequencing reads", J Comput Biol 22(6), 2015, doi:10.1089/cmb.2014.0157.
This format can be summarized as follows:
- each line represents a read that must be ordered depending on the starting position. The line must terminate with the comment symbol `#`.
- Each read is composed of SNPs, which must be separated by using `:` as separator.
- Each SNP is codified in the format `POS ALL BIN WEI`: i) `POS` is the position of the SNP corresponding to its genomic coordinate; ii) `ALL` is the allele of the corresponding SNP, equal to a value in {A, C, G, T}; iii) `BIN` is the binary representation of the allele (i.e., 0 is the wild-type and 1 is the SNP); iv)`WEI` is the phred score (or weight) of the SNP reported in the numeric format.

For instance, the following example represents a correct instance:

    4918 G 1 12 : 4959 A 1 32 : # 60 : NA 
    4959 A 1 61 : 5063 A 0 42 : # 60 : NA 
    6127 A 0 31 : 6181 G 1 16 : 6182 C 0 25 : # 60 : NA
    19205 T 0 44 : 19391 T 0 15 : 19437 T 1 27 : 19440 A 0 17 : # 60 : NA 

### Dataset ###

In order to test the performances of GenHap, we generated two synthetic (yet realistic) datasets, each one consisting of instances obtained from a specific sequencing technology.
In particular, we considered the Roche/454 genome sequencer (Roche AG, Basel, Switzerland), representing one of the next-generation sequencing (NGS) systems able to produce long and precise reads, and the PacBio RS II sequencer, which is an emerging third-generation sequencing technology.

All the tested models can be downloaded at <https://github.com/andrea-tango/GenHap/blob/master/Models.zip>.

For each model, we provide a txt file containing the correct haplotypes (i.e., the ground truth).

## <a name="lic"></a>License ##

GenHap is licensed under the terms of the GNU GPL v3.0

## <a name="cont"></a>Contacts ##

For questions or support, please contact <andrea.tangherloni@disco.unimib.it>
or <simone.spolaor@disco.unimib.it>
