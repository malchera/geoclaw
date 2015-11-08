#!/bin/bash

COMPILERS=("ifort" "gfortran")
#FLAGS=("" "-O3" "-O3 -xhost")
FLAGS=("-O2" "-O3")
RESOLUTIONS=("60" "120" "360")
# TODO: TRESOLUTIONS=(...)
ORIGINAL_FC=$FC


function usage
{
    echo "Usage: ./script [all|compile|runs]"
    exit 0
}

function checkparams
{
    if [[ $# -gt 0 ]]
    then
        if [ "$1" == "compile" ]
        then
            echo "compile"
        elif [ "$1" == "runs" ]
        then
            echo "runs"
        elif [ "$1" == "all" ]
        then
            echo "all"
        else
            usage
        fi
    else
        usage
    fi
}

function main
{
    checkparams $1
# TODO: Check if files exist: Makefile, setrun.py, job.sh and maketopo.py
    for compiler in "${COMPILERS[@]}";
    do
        for flags in "${FLAGS[@]}";
        do
            # Replace dashes by underscores and trim whitespaces for file names
            flagstring=`echo $flags | sed -r "s/ *-/_/g"`
            # Binary name containing compiler name and flags
            binname="xgeoclaw_${compiler}${flagstring}"
            # Create temporary directory, copy Makefile to it, sed executable name and flags and build binary
            tmpdir="_tmp_${compiler}_${flagstring}"
            mkdir -p $tmpdir
            sed -r -e "s/^EXE.*/EXE = ${binname}/" \
                   -e "s/^FFLAGS.*/FFLAGS = ${flags}/" \
                   -e "s/(^CLAW_PKG.*$)/FC=${compiler}\t# Set by script, do not change\n\1/" \
                Makefile > $tmpdir/Makefile
            cd $tmpdir
            if [[ ! -f $binname ]]
            then
                echo "=== BUILDING NEW EXECUTABLE! THIS MIGHT TAKE SOME TIME ==="
                echo -n "Will build $binname with flags \"$flags\" in "
                for i in {3..1}
                do
                    echo -n "$i... "
                    sleep 1
                done
                echo
                echo "Building $binname..."
                make new &> /dev/null
                echo "Done!"
            else
                echo "$binname already exists. skipping..."
            fi
            cd ..

            for res in "${RESOLUTIONS[@]}";
            do
                dirname="run_${compiler}${flagstring}_${res}"
                mkdir -p $dirname
                cp $tmpdir/* $dirname # Copy modified Makefile and corresponding binary
                cp maketopo.py $dirname # This is needed to create topography files
                # TODO: replace setrun.py, job.sh copies by sed commands
                # Change grid resolution in setrun.py and set AMR levels to 1
                sed -r -e "s/( *clawdata\.num_cells\[[01]\] *= *).*$/\1${res}/g" \
                       -e "s/( *amrdata.amr_levels_max *= *).*$/\11/g" \
                    setrun.py > $dirname/setrun.py
                # Change job file's job name and error and output file names, respectively.
                jobname="chile_${compiler}${flagstring}_${res}"
                sed -r -e "s/(^#SBATCH +-J +)(\w*)(#*.*$)/\1${jobname}\3/" \
                       -e "s/(^#SBATCH +-o +)(\w*)(#*.*$)/\1${jobname^^}_OUT\3/" \
                       -e "s/(^#SBATCH +-e +)(\w*)(#*.*$)/\1${jobname^^}_ERR\3/" \
                    job.sh > $dirname/job.sh

                # SLURM Job submit
                echo -n "Submit job ${jobname}? 'y' for yes, other key to skip: "
                read yesno
                if [[ $yesno == "y" ]]; then
                    echo "SUBMITTING JOB \"${jobname}\"..."
                    # Submit job, hide stdout
                    cd $dirname
                    sbatch job.sh 1> SBATCH_${jobname} & 
                    cd ..
                else
                    echo "Skipping..."
                fi 
            done
        done
    done
}

main $1
