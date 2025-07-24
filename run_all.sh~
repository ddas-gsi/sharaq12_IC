#!/bin/bash

### Usage: ./run_all.sh 1052 1054 1055 1056

# Check for at least one RUN number
if [ "$#" -lt 1 ]; then
    echo "Usage: ./run_all.sh RUN_Nbr1 [RUN_Nbr2 ...]"
    exit 1
fi

for RUN_Nbr in "$@"
do
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
done

echo "✅✅ All RUN numbers processed."



















# ========================================================

# #!/bin/bash

# # Check for RUN_Nbr argument
# if [ -z "$1" ]; then
#     echo "Usage: ./run_all.sh RUN_Nbr"
#     exit 1
# fi

# RUN_Nbr=$1

# echo "Starting full analysis for RUN Number: $RUN_Nbr"

# echo "Running icChargeCorrection.C..."
# root -l -q "icChargeCorrection.C($RUN_Nbr)"
# sleep 10

# echo "Running slowBeamCharge2Energy.C..."
# root -l -q "slowBeamCharge2Energy.C($RUN_Nbr)"
# sleep 10

# echo "Running fitSpline.py..."
# python3 fitSpline.py $RUN_Nbr
# sleep 10

# echo "Running splineAnaCorrection24.py..."
# python3 splineAnaCorrection24.py $RUN_Nbr
# sleep 10

# echo "Running analyzer_spline.C..."
# root -l -q "analyzer_spline.C($RUN_Nbr)"
# sleep 10

# echo "Running zetACalSpline.py..."
# python3 zetACalSpline.py $RUN_Nbr
# sleep 10

# echo "All scripts executed for RUN $RUN_Nbr."
# ========================================================


# ========================================================

# For LOG FILE logging, uncomment the following lines:

# RUN_Nbr=$1
# LOGFILE="/u/ddas/software/work/artemis-oedo/output/logs/run_all.log"

# # Start log
# echo "--------------------------------------------" | tee -a $LOGFILE
# echo "[`date`] Starting analysis for RUN $RUN_Nbr" | tee -a $LOGFILE
# echo "--------------------------------------------" | tee -a $LOGFILE

# echo "[`date`] Running icChargeCorrection.C..." | tee -a $LOGFILE
# root -l -q "icChargeCorrection.C($RUN_Nbr)" 2>&1 | tee -a $LOGFILE
# sleep 10

# echo "[`date`] Running fitSpline.py..." | tee -a $LOGFILE
# python3 fitSpline.py $RUN_Nbr 2>&1 | tee -a $LOGFILE
# sleep 10

# echo "[`date`] Running splineAnaCorrection24.py..." | tee -a $LOGFILE
# python3 splineAnaCorrection24.py $RUN_Nbr 2>&1 | tee -a $LOGFILE
# sleep 10

# echo "[`date`] Running analyzer_spline.C..." | tee -a $LOGFILE
# root -l -q "analyzer_spline.C($RUN_Nbr)" 2>&1 | tee -a $LOGFILE
# sleep 10

# echo "[`date`] Running zetACalSpline.py..." | tee -a $LOGFILE
# python3 zetACalSpline.py $RUN_Nbr 2>&1 | tee -a $LOGFILE
# sleep 10

# echo "[`date`] All scripts completed for RUN $RUN_Nbr." | tee -a $LOGFILE

# ========================================================