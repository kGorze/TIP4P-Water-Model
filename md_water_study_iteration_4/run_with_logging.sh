#!/bin/bash
# Wrapper script to run the workflow with logging

# Define log file
LOG_DIR="/home/konrad_guest/Documents/research/cursor/md_water_study_iteration_4/logs"
LOG_FILE="${LOG_DIR}/workflow_$(date +%Y%m%d_%H%M%S).log"

# Create log directory if it doesn't exist
mkdir -p "${LOG_DIR}"

# Function to check log for errors
check_log_for_errors() {
    echo "====== LOG STATUS CHECK at $(date) ======"
    
    # Count errors in log
    ERROR_COUNT=$(grep -i "error" "${LOG_FILE}" | wc -l)
    WARNING_COUNT=$(grep -i "warning" "${LOG_FILE}" | wc -l)
    
    echo "Found ${ERROR_COUNT} errors and ${WARNING_COUNT} warnings so far."
    
    # Show the last error if any
    if [ ${ERROR_COUNT} -gt 0 ]; then
        echo "Latest error:"
        grep -i "error" "${LOG_FILE}" | tail -1
    fi
    
    # Show latest log entries
    echo "Latest log entries:"
    tail -5 "${LOG_FILE}"
    echo "====================================="
}

# Announce start
echo "Starting workflow with logging to ${LOG_FILE}"
echo "Will check log status every 60 seconds."
echo "You can manually check the log with: tail -f ${LOG_FILE}"

# Run the workflow and redirect all output to log
{
    echo "===== WORKFLOW STARTED at $(date) ====="
    echo "Running on host: $(hostname)"
    echo ""
    
    # Run the actual workflow script
    /home/konrad_guest/Documents/research/cursor/md_water_study_iteration_4/run_workflow.sh
    
    WORKFLOW_EXIT_CODE=$?
    
    echo ""
    echo "===== WORKFLOW FINISHED at $(date) with exit code ${WORKFLOW_EXIT_CODE} ====="
} > "${LOG_FILE}" 2>&1 &

# Get the PID of the workflow process
WORKFLOW_PID=$!

# Periodically check the log while the workflow is running
while kill -0 ${WORKFLOW_PID} 2>/dev/null; do
    # Wait for 60 seconds
    sleep 60
    
    # Check log status
    check_log_for_errors
    
    # Check if there are specific errors we're looking for
    if grep -q "C-rescale does not support pressure coupling type Anisotropic" "${LOG_FILE}"; then
        echo "DETECTED KNOWN ISSUE: C-rescale does not support anisotropic pressure coupling."
        echo "This is an expected error. The workflow will continue but NPT equilibration may fail."
    fi
    
    if grep -q "Simulation" "${LOG_FILE}" | grep -q "not defined"; then
        echo "DETECTED KNOWN ISSUE: Simulation not defined in Packmol. Check if Packmol API is correctly used."
    fi
    
    echo "Workflow still running (PID: ${WORKFLOW_PID}). Next check in 60 seconds..."
done

# Final check after workflow completes
check_log_for_errors

# Summarize the workflow
echo "Workflow completed. Full log available at: ${LOG_FILE}"
echo "Summary of important events:"
echo "--------------------------"
echo "Errors:"
grep -i "error" "${LOG_FILE}" | head -10  # Show first 10 errors

echo "Steps completed:"
grep -A 1 "Step" "${LOG_FILE}" | grep -v -- "--" | tail -10  # Show last completed steps

# Check final outcome
if grep -q "Workflow completed!" "${LOG_FILE}"; then
    echo "OVERALL STATUS: WORKFLOW COMPLETED SUCCESSFULLY"
else
    echo "OVERALL STATUS: WORKFLOW ENCOUNTERED ISSUES"
    # Try to identify where it failed
    grep -A 2 "failed" "${LOG_FILE}" | tail -3
fi 