#!/bin/bash
# Copyright (c) 2013 Yuan Zhang, Yanni Sun.
# You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
# pipeline of MetaDomain. 

# input: 
# -m: hmm file;
# -f: fasta file;
# -t: alignment overlap threshold, default: 20;
# -d: relative overlap difference threshold: 0.15;

# output:
# -o: output_contig_file

# get the installation path.
usage() {
  echo "SAT-Assembler.sh -m <HMM file> -f <fasta file> -o <output folder> [options] 
  Options:
    -h:  show this message
    -t:  alignment overlap threshold, default: 20;
    -d:  relative overlap difference threshold: 0.15;
    -o:  output file name, default: stdandard error"
}

hmm=
fasta=
t=20
d=0.15
out=
# installation dir.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while getopts "hm:f:t:d:o:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    m)
      hmm=$OPTARG
      ;;
    f)
      fasta=$OPTARG
      ;;
    t)
      t=$OPTARG
      ;;
    d)
      d=$OPTARG
      ;;
    o)
      out=$OPTARG     
      ;;
    esac
done

if [ "$hmm" == "" ];then
  echo "Please specify the hmm file."
  usage
  exit
fi

if [ "$fasta" == "" ];then
  echo "Please specify the input fasta file."
  usage
  exit
fi

if [ "$out" == "" ];then
  echo "Please specify the output folder."
  usage
  exit
fi


if [ `which hmmsearch 2> /dev/null | wc -l` -eq 0 ]; then
  echo "hmmsearch is not found.";
  usage
  exit
fi

if [ `$DIR/check_python_packages.py` -eq 1 ];then
  echo "Biopython or NetworkX is not found."
  usage
  exit
fi 

# generate a temporary folder.
if [ ! -d $out/ ];then
  mkdir $out/
else
  echo 'Output folder exists. Please specify another output folder.'
  exit
fi
tmp="$(cd $out && pwd)"
base_fasta=`echo $fasta | awk '{split($1,a,"/"); print a[length(a)]}'`

$DIR/DNA2Protein 1-6 $fasta $tmp/${base_fasta} 
# generate a list of domains in the input hmm file.
python $DIR/parse_hmm_files.py $hmm $tmp/HMMs
ls $tmp/HMMs | while read line
do
  hmm_acc=`echo $line | awk '{print substr($1,1,7)}'`
  cat /dev/null >$tmp/${base_fasta}_${hmm_acc}.hmmer
  for i in {1..6}
  do
    bash $DIR/hmmer3_pipeline_strand.sh $tmp/HMMs/$line $tmp/${base_fasta}.frame${i} $i >$tmp/${base_fasta}_${hmm_acc}.hmmer
  done
  python $DIR/assembler.py $tmp/${base_fasta}_${hmm_acc}.hmmer $fasta ${hmm_acc} $t $d $out 
done
rm -r $tmp/HMMs
rm $tmp/${base_fasta}.frame?
