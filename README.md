# FastOWM

This repo includes different implementations of the [OWM algorithm](https://www.mdpi.com/2072-4292/12/7/1051). There are four main directories:

* ``RandC-Baseline`` includes the original R implementation and the porting to language C that we consider as the baseline. This directory should be taken just as a reference of the seminal work implemented in R and C, but it is old code, significantly superseded in the other directories.
* ``CPP-CPU-Hetero`` includes a C++ version with several optimizations and a heterogeneous implementation for CPU+GPU that uses a dynamic heterogeneous scheduler and SYCL for the GPU implementation.
* ``SYCL-CUDA`` includes the SYCL implementation (that can run on CPU or GPU) and a pure CUDA version.
* ``Results`` includes the Jupyter Notebooks that process the execution times and other results obtained when executing the different implementations with different optimizations.

With more details, the **CPP-CPU-Hetero directory**, includes several implementations and optimizations that are validated (using the python scripts in the ``scripts`` directory) and reported in the ``Results`` directory. The summary of the implementations and optimizations are described in this [paper](https://www.overleaf.com/project/615584a8f2c4278161fc2b94):

* Baseline: see [process_baseline_results.ipynb](Results/process_baseline_results.ipynb)
* Optimization 1: 2D-space bipartition instead of 3D-space bipartition: see [process_opt1_results.ipynb](Results/process_opt1_results.ipynb), [o2and3_minradius.ipynb](Results/o2and3_minradius.ipynb) and [o2and3_maxnumber.ipynb](Results/o2and3_maxnumber.ipynb).
* Optimization 2: tune the granularity and load balance: see [process_o2and3_results.ipynb](Results/process_o2and3_results.ipynb).
* Optimization 3: reduce the number of accesses by memoization: see [sycl_minrad_maxnum.ipynb](Results/sycl_minrad_maxnum.ipynb)
* Parallelization: see the results with more than one core in [sycl_minrad_maxnum.ipynb](Results/sycl_minrad_maxnum.ipynb)
* Comparison of all optimizations: see [process_all_optimizations.ipynb](Results/process_all_optimizations.ipynb)