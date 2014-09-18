#!/bin/bash
# this file is to install SAT-Assembler in the current folder.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
chmod 755 hmmer3_pipeline_strand.sh
chmod 755 SAT-Assembler.sh
chmod 755 parse_hmm_files.py
chmod 755 assembler.py
chmod 755 check_python_packages.py
