# FastOWM

This repo includes different implementations of the [OWM algorithm](https://www.mdpi.com/2072-4292/12/7/1051). 

Please, cite [this paper](https://authors.elsevier.com/sd/article/S1877-7503(24)00110-8) if you use these implementations for your own developments (BibTeX format at the end):

```
Felipe Muñoz, Rafael Asenjo, Angeles Navarro, J. Carlos Cabaleiro,
CPU and GPU oriented optimizations for LiDAR data processing,
Journal of Computational Science, 2024. 
```
[https://doi.org/10.1016/j.jocs.2024.102317](https://doi.org/10.1016/j.jocs.2024.102317)

There are four main directories:

* ``RandC-Baseline`` includes the original R implementation and the porting to language C that we consider as the baseline. This directory should be taken just as a reference of the seminal work implemented in R and C, but it is old code, significantly superseded in the other directories.
* ``CPP-TBB-CPU`` includes a C++ version with several optimizations that uses TBB to exploit parallelism on multicore CPUs.
* ``SYCL-CUDA-GPU`` includes the SYCL implementation (that can run on CPU or GPU) and a pure CUDA version.
* ``Results`` includes the Jupyter Notebooks that process the execution times and other results obtained when executing the different implementations with different optimizations.

With more details, the **CPP-TBB-CPU directory**, includes several implementations and optimizations that are validated (using the python scripts in the ``scripts`` directory) and reported in the ``Results`` directory. The summary of the implementations and optimizations are described in this [paper](https://www.overleaf.com/project/615584a8f2c4278161fc2b94):

* Baseline: see [process_baseline_results.ipynb](Results/process_baseline_results.ipynb)
* Optimization 1: 2D-space bipartition instead of 3D-space bipartition: see [process_o1quadtree_results.ipynb](Results/process_o1quadtree_results.ipynb).
* Optimization 2: CPU parallel optimizations: see [process_o2partree_results.ipynb](Results/process_o2partree_results.ipynb).
* Optimization 3: reduce the number of accesses by memoization: see [process_o3memo_results.ipynb](Results/process_o3memo_results.ipynb).
* Optimization 4: tune de granularity: see [process_o4minrad_maxnum.ipynb](Results/process_o4minrad_maxnum.ipynb)
* Comparison of all optimizations: see [process_all_optim.ipynb](Results/process_all_optim.ipynb)

On the other hand, the **SYCL-CUDA-GPU directory** includes the SYCL and CUDA implementations. Have a look at the ``Makefile`` and ``scripts``directory to see how to compile the different implementations (with different optimizations) and how to collect the execution data. The results have been analyzed in these Jupyter Notebooks:

* Study of all optimizations in [process_sycl_cuda_results.ipynb](Results/process_sycl_cuda_results.ipynb)
* Comparison of all optimizations and with the baseline and best CPU version in [process_all_optimSYCL-CUDA.ipynb](Results/process_all_optimSYCL-CUDA.ipynb)

Finally, to compare the results across the different hosts, we have used [process_all_servers.ipynb](Results/process_all_servers.ipynb).

The LiDAR clouds used in our experiments can be downloaded from [here](https://www.dropbox.com/s/0vpr8gow624ngqz/nubesLidar.tgz?dl=0) (2.49GB).

BibTeX format:

```
@article{Felipe-JCS-2024-FastOWM,
title = {CPU and GPU oriented optimizations for LiDAR data processing},
author = {Felipe Muñoz and Rafael Asenjo and Angeles Navarro and J. Carlos Cabaleiro},
journal = {Journal of Computational Science},
pages = {102317},
year = {2024},
issn = {1877-7503},
doi = {https://doi.org/10.1016/j.jocs.2024.102317},
url = {https://www.sciencedirect.com/science/article/pii/S1877750324001108},
keywords = {LiDAR data processing, Digital Terrain Model, Tree data structures, Parallel optimization, GPU, SYCL, CUDA, oneAPI}}
```
