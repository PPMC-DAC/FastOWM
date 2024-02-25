#!/bin/bash

# run through all the scripts in the directory with "run_" in the name
for script in run_*.py; do
    # Skip specific scripts
    if [[ "$script" == "run_par_baseline.py" || "$script" == "run_par_rev1.py" ]]; then
        echo "Skipping $script"
    else
        echo "Running $script"
        python "$script"
        echo "$script execution completed"
    fi
done