# Version
This software was last updated 05/19/2014. To report bugs, please send emails
to nick.zhangyuan@gmail.com.

# Installation
1. Untar the source file called "SAT-Assembler.tar.gz".
2. g++ compiler is required in your Unix system. To install SAT-Assembler, run the Makeme file using the following command:
make 
3. Dependencies:
  1) HMMER3 (http://hmmer.janelia.org/). The bin file hmmsearch should be in the path.  
  2) Python version 2.5 or later.  
  3) Python libraries of NetworkX (http://networkx.lanl.gov/) and Biopython (http://biopython.org/wiki/Main_Page).  

# Run SAT-Assembler
To run SAT-Assembler, use the following command:
```
./SAT-Assembler.sh -m <HMM file> -f <fasta file> [options]
  Options:
    -h:  show this message
    -t:  alignment overlap threshold, default: 20;
    -d:  relative overlap difference threshold: 0.15;
    -o:  output file name, default: stdandard error
```

1. The hmm file can contain multiple hmm models and should be in HMMER3.0's hmm file format. All the hmm files of Pfam database can be downloaded from Pfam's website.
2. If you build the hmm file yourself using hmm-build in HMMER3, please make sure you have a accession number (the line that begins with ACC) as its unique identifier. Otherwise, please manual add it. 
3. The nucleotide sequence file should be in fasta format. All the reads should be in a single fasta file. 
4. The format of paired-end reads is should be in ".1" and ".2" notation. An example of a paired-end read will be gnl|SRA|SRR360147.1.1 and gnl|SRA|SRR360147.1.2.
 
# Output
The output includes a contig file in fasta format and a scaffold file. The name of the fasta file and the scaffold file include the name of the family. For example, the contig file and scaffold file for PF00005 are PF00005_contigs.fa and PF00005_scaffolds.txt. Each line of the scaffold file indicates the contigs that are from the same scaffold.

# License
Copyright (C) 2014 Yuan Zhang, Yanni Sun, and James Cole.

