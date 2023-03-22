# LCAT
long-read error correction algorithm for transcriptome sequencing data
## 1.Introduction
LCAT (An isoform-sensitive error correction for transcriptome sequencing long reads) is a wrapper algorithm of MECAT, to reduce the loss of isoform diversity while keeping MECAT's error correction performance. The experimental results show that LCAT not only can improve the quality of transcriptome sequencing long reads, but also keeps the diversity of isoforms.
## 2.Installation
### Install LCAT
```
git clone https://github.com/luckylyw/LCAT.git
cd LCAT
make
cd ..
export PATH=/home/luoluo/tool/LCAT/Linux-amd64/bin:$PATH
After installation, all the executables are found in LCAT/ Linux-amd64/bin.
```
### Install HDF5
```
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.15-patch1/src/hdf5-1.8.15-patch1.tar.gz
tar xzvf hdf5-1.8.15-patch1.tar.gz
mkdir hdf5
cd hdf5-1.8.15-patch1
./configure --enable-cxx --prefix=/home/luoluo/tool/hdf5
make
make install
cd ..
export HDF5_INCLUDE=/home/wjzhang/lyw/tool/hdf5/include
export HDF5_LIB=/home/wjzhang/lyw/tool/hdf5/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/luoluo/tool/hdf5/lib
```
The header files of HDF5 are in hdf5/include. The library files of HDF5 are in hdf5/lib
### Install dextract
```
git clone https://github.com/PacificBiosciences/DEXTRACTOR.git
cp LCAT/dextract_makefile DEXTRACTOR
cd DEXTRACTOR
export PATH=/home/luoluo/tool/DEXTRACTOR:$PATH
edit the dextractor_makefile (line 7) :
${CC} $(CFLAGS) -I$(HDF5_INCLUDE) -L$(HDF5_LIB) -o dextract dextract.c sam.c bax.c expr.c DB.c QV.c -lhdf5 -lz
make -f dextract_makefile
cd ..
```
## 3.Quick Start
We can use LCAT to correct RNA long reads produced by PacBio or Nanopore. Options of each command will be given in next section.
### Correcting Pacbio Data
* step 1, using lcat2pw to detect overlapping candidates
```
lcat2pw x 0 -d SRR6238555.fastq -o SRR6238555.fastq.pm.can  -w wrk_dir -t 40 -n 100 -a 100 -k 4 -g 0
```
* step 2, using lcat2cns to correct the noisy RNA reads based on their pairwise overlapping candidates.
```
lcat2cns -x 0 -t 40 -p 100000 -a 100 -l 100 -r 0.6  -c 4  -k 10 SRR6238555.fastq.pm.can SRR6238555.fastq corrected_reads.fastq
```
### Correcting Nanopore Data
* step 1, using lcat2pw to detect overlapping candidates
```
lcat2pw -x 1 -d ERR2401483_proccessed_normalid.fasta  -o candidatex.txt -w wrk_dir -t 40 -n 100 -a 100 -k 4 -g 0
```
* step 2, using lcat2cns to correct the noisy RNA reads based on their pairwise overlapping candidates.
```
lcat2cns -x 0 -t 40 -p 100000 -a 100 -l 100 -r 0.6  -c 4  -k 10 candidatex.txt ERR2401483_proccessed_normalid.fasta corrected_reads.fastq
```
## 4.Program Descriptions
We describe in detail each module of LCAT, including their options and output formats.
### lcat2pw
#### Input Format
LCAT is capable of processing FASTA, FASTQ, format files.
#### Options
The command for running lcat2pw is
```
lcat2pw [-j task] [-d dataset] [-o output] [-w working dir] [-t threads] [-n candidates] [-g 0/1]
```
The options are:
```
-j <integer>    job: 0 = seeding, 1 = align
       default: 0
-d <string>    reads file name
-o <string>    output file name
-w <string>    working folder name, will be created if not exist
-t <integer>    number of cput threads
       default: 1
-n <integer>    number of candidates for gapped extension
       Default: 100
-a <integer>    minimum size of overlaps
       Default: 2000 if x = 0, 500 if x = 1
-k <integer>    minimum number of kmer match a matched block has
       Default: 4 if x = 0, 2 if x = 1
-g <0/1>    whether print gapped extension start point, 0 = no, 1 = yes
       Default: 0
-x <0/x>    sequencing technology: 0 = pacbio, 1 = nanopore
       Default: 0
```
#### Output Format
the results are output in can format, each result of which occupies one line and 9 fields:
```
[A ID] [B ID] [A strand] [B strand] [A gapped start] [B gapped start] [voting score] [A length] [B length]
```
If the -g option is set to 1, two more fields indicating the extension starting points are given:
```
[A ID] [B ID] [% identity] [voting score] [A strand] [A start] [A end] [A length] [B strand] [B start] [B end] [B length] [A ext start] [B ext start]
```
In the strand field, 0 stands for the forward strand and 1 stands for the reverse strand. All the positions are zero-based and are based on the forward strand, whatever which strand the sequence is mapped.
### lcat2cns
lcat2cns is RNA long reads self error correction tool.
#### Input Format
inputs to lcat2cns can be can format files.
#### Options
The command for running lcat2cns is
```
lcat2cns [options] input reads output
```
The options are:
```
-x <0/1>    sequencing platform: 0 = PACBIO, 1 = NANOPORE
       default: 0
-t <Integer>    number of threads (CPUs)
-p <Integer>    batch size that the reads will be partitioned
-r <Real>    minimum mapping ratio
-a <Integer>    minimum overlap size
-c <Integer>    minimum coverage under consideration
-l <Integer>    minimum length of corrected sequence
-k <Integer>    number of partition files when partitioning overlap results (if < 0, then it will be set to system limit value)
-d <Real>    identity threshold
-w <Integer>    slide window length
-m <Real>    minimum coverage rate of modify region
-h        print usage info.
```
If 'x' is set to be '0' (pacbio), then the other options have the following default values:
-t 1 -p 100000 -r 0.9 -a 2000 -c 6 -l 5000 -k 10 -d 0.65 -w 75 -m 0.05
If 'x' is set to be '1' (nanopore), then the other options have the following default values:
-t 1 -p 100000 -r 0.4 -a 400 -c 6 -l 2000 -k 10 -d 0.65 -w 75 -m 0.05
#### Output Format
The corrected sequences are given in FASTA format. The header of each corrected sequence consists of three components seperated by underlines:
```
>A_B_C_D
```
where
A is the original read id
B is the left-most effective position
C is the right-most effective position
D is the length of the corrected sequence
by effective position we mean the position in the original sequence that is covered by at least c (the argument to the option -c) reads.
## 5.Citation
LCAT: long-read error correction algorithm for transcriptome sequencing data

