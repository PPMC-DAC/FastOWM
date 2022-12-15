# FastOWM

This repo includes different implementations of the [OWM algorithm](https://www.mdpi.com/2072-4292/12/7/1051).

* RandC-Baseline includes the original R implementation and the porting to language C that we consider as the baseline.
* CPP-CPU-Hetero includes a C++ version with several optimizations and a heterogeneous implementation for CPU+GPU that uses a dynamic heterogeneous scheduler and SYCL for the GPU implementation
* SYCL-CUDA includes the SYCL implementation (that can run on CPU or GPU) and a pure CUDA version