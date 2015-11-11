#!/bin/bash

echo -n "file="
read file
echo -n "-qopt-report="
read n
echo -n "-qopt-report-file="
read reportfile
echo -n "-qopt-report-phase="
read phase
echo -n "-qopt-report-filter="
read filter
#filter="\"$filter\""
echo

read -r -d '' cmd << EOF
ifort -cpp -c $file
-I$CLAW/amrclaw/src/2d
-I$CLAW/geoclaw/src/2d/shallow
-I$TACC_PAPI_INC
-qopt-report=$n
-qopt-report-file=$reportfile
-qopt-report-phase=$phase
-qopt-report-filter=$filter
EOF

set -o xtrace
exec $cmd
