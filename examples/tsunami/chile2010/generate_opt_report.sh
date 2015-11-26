#!/bin/bash

echo -n "file="
read file
echo -n "-qopt-report="
read n
echo -n "-qopt-report-file="
read reportfile
echo -n "-qopt-report-phase="
read phase
#filter="\"$filter\""
echo
#-align array64byte 
read -r -d '' cmd << EOF
ifort -cpp -qopenmp-simd -xAVX -fopenmp -c $file
-S -fsource-asm -masm=intel
-I$CLAW/amrclaw/src/2d
-I$CLAW/geoclaw/src/2d/shallow
-I$TACC_PAPI_INC
-qopt-report=$n
-qopt-report-file=$reportfile
-qopt-report-phase=$phase
EOF

set -o xtrace
exec $cmd
