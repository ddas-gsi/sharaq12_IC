#!/bin/bash

### Usage: ./run_parallel.sh 1052 1054 1055 1056

# Check for at least one RUN number
if [ "$#" -lt 1 ]; then
    echo "Usage: ./run_parallel.sh RUN_Nbr1 [RUN_Nbr2 ...]"
    exit 1
fi

MAX_PARALLEL=3   # Maximum number of parallel jobs

# Function to run all scripts for a single RUN
run_single() {
    local RUN_Nbr=$1

    echo "--------------------------------------------"
    echo "Starting full analysis for RUN Number: $RUN_Nbr"
    echo "--------------------------------------------"

    echo "Running icChargeCorrection.C..."
    root -l -q "icChargeCorrection.C($RUN_Nbr)"
    sleep 10

    echo "Running slowBeamCharge2Energy.C..."
    root -l -q "slowBeamCharge2Energy.C($RUN_Nbr)"
    sleep 10

    echo "Running fitSpline.py..."
    python3 fitSpline.py $RUN_Nbr
    sleep 10

    echo "Running splineAnaCorrection24.py..."
    python3 splineAnaCorrection24.py $RUN_Nbr
    sleep 10

    echo "Running analyzer_spline.C..."
    root -l -q "analyzer_spline.C($RUN_Nbr)"
    sleep 10

    echo "Running zetACalSpline.py..."
    python3 zetACalSpline.py $RUN_Nbr
    sleep 10

    echo "✅ Finished processing RUN Number: $RUN_Nbr"
    echo ""
}

# Counter for running jobs
job_count=0

for RUN_Nbr in "$@"; do
    run_single "$RUN_Nbr" &   # Run in background
    ((job_count++))

    # Wait if max parallel jobs reached
    if [[ $job_count -ge $MAX_PARALLEL ]]; then
        wait -n              # Wait for any job to finish
        ((job_count--))      # Decrease job count
    fi
done

# Wait for all remaining background jobs
wait

echo "✅✅ All RUN numbers processed."
