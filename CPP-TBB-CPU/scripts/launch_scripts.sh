#!/bin/bash

# generate a list of scripts to run
list_of_scripts="run_par_baseline.py run_par_optim1_qtree.py run_par_optim2.py run_par_rev1.py"
# iterate the list of scripts and run them
for script in $list_of_scripts; do
    echo "Running $script"
    python "$script"
    echo "$script execution completed"
done