//============================================================================
// Name			: Dynamic.h
// Author		: Denisa C.
// Version		: 1.0
// Date			: 02 / 01 / 2014
// Copyright	: Department. Computer's Architecture (c)
// Description	: Dynamic scheduler implementation
//============================================================================

//#define DEEP_GPU_REPORT

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "tbb/pipeline.h"
#include "tbb/tick_count.h"
// #include "tbb/atomic.h"
#include <atomic>
#include <math.h>

#ifdef ENERGYCOUNTERS
#ifdef Win32
#include "PCM_Win/windriver.h"
#else
#include "cpucounters.h"
#endif
#endif

#ifdef PJTRACER
#include "pfortrace.h"
#endif

#include "SchedulerOneApi.h"

#include "CL/sycl.hpp"

using namespace cl::sycl;

using namespace std;
using namespace tbb;

/*****************************************************************************
 * Defines
 * **************************************************************************/
#define CPU 0
#define GPU 1
#define GPU_OFF -100 //Arbitrary value
#define SEP "\t"

//#define  GET_KERNEL_TIME

/*****************************************************************************
 * Global variables
 * **************************************************************************/
std::atomic<int> gpuStatus;
int chunkCPU;
int chunkGPU;
float fG;
float gpuThroughput, cpuThroughput;
// To calculate scheduling partition overhead
// tick_count end_tc;
// cl::sycl::event event_kernel, event_send, event_receive; // unable to create

#ifdef DEEP_CPU_REPORT
ofstream deep_cpu_report;
#endif

#ifdef DEEP_GPU_REPORT
ofstream deep_gpu_report;
#endif

#ifdef OVERHEAD_ANALYSIS
// overhead accumulators
double overhead_sp = 0.0;
double overhead_h2d = 0.0;
double overhead_kl = 0.0;
double kernel_execution = 0.0;
double overhead_d2h = 0.0;
double overhead_td = 0.0;
#endif

long totalIterationsGPU = 0;
long totalIterations = 0;

/*****************************************************************************
 * Heterogeneous Scheduler
 * **************************************************************************/
/*Bundle class: This class is used to store the information that items need while walking throught pipeline's stages.*/
class Bundle
{
public:
	ulong begin;
	ulong end;
	int type; //GPU = 0, CPU=1

	// Bundle(){};
};

/*My serial filter class represents the partitioner of the engine. This class selects a device and a subrange of iterations*/
class MySerialFilter : public filter
{
private:
	int begin;
	int end;
	int nCPUs;

public:
	/*Class constructor, it only needs the first and last iteration indexes.*/
	MySerialFilter(int b, int e, int ncpus) : filter(true)
	{
#ifdef DEBUG
		printf("\nDEBUG MySerialFilter: begin = %d\n", b);
		printf("\nDEBUG MySerialFilter: end = %d\n", e);
#endif
		begin = b;
		end = e;
		nCPUs = ncpus;
	}

	/*Mandatory operator method, TBB rules*/
	void *operator()(void *)
	{
		Bundle *bundle = new Bundle();

		/*If there are remaining iterations*/
		if (begin < end /*&& begin > -1 && end > -1*/)
		{
#ifdef DEBUG
			/*if (begin<0)*/ printf("\nDEBUG operator: begin = %d\n", begin);
			/*if (end<0)*/ printf("\nDEBUG operator: end = %d\n", end);
#endif
			//Checking which resources are available
			if (--gpuStatus >= 0)
			{
				//GPU WORK
				int auxEnd = begin + chunkGPU;
				if (nCPUs < 1)
				{
					auxEnd = (auxEnd > end) ? end : auxEnd;
				}

				if (auxEnd <= end)
				{
					//GPU CHUNK
					bundle->begin = begin;
					bundle->end = auxEnd;
					begin = auxEnd;
					bundle->type = GPU;
				}
				else
				{
					//CPU CHUNK
					gpuStatus = GPU_OFF;
					if (fG > 0)
					{
						chunkCPU = std::min(chunkCPU, std::max((end - begin) / (nCPUs + 1), 1));
					}
					bundle->begin = begin;
					bundle->end = begin + chunkCPU;
					begin = begin + chunkCPU;
					bundle->type = CPU;
				}
#ifdef DEBUG
				if (bundle->begin < 0)
					printf("\nDEBUG operator: if ( --gpuStatus >= 0 ) begin = %d\n", bundle->begin);
				if (bundle->end < 0)
					printf("\nDEBUG operator: if ( --gpuStatus >= 0 ) end = %d\n", bundle->begin);
#endif

				return bundle;
				/*! HERE IS THE PROBLEM */
			}
			else
			{
				//CPU WORK
				gpuStatus++;

				/*Calculating next chunkCPU*/
				if (fG > 0.01 /*0*/)
				{ /*! HERE IS THE PROBLEM */
#ifdef DEBUG
					printf("\nDEBUG operator: else  chunkCPU init = %d \n", chunkCPU);
#endif
					chunkCPU = chunkGPU / fG; /*fG 0 0.000..1 => a chunk size that is not accheptable*/
#ifdef DEBUG
					printf("\nDEBUG operator: else  chunkCPU/fG = %d, fg = %f \n", chunkCPU, fG);
#endif
					chunkCPU = std::min(chunkCPU, std::max((end - begin) / (nCPUs), 1));
				}

				/*Taking a iteration chunk for CPU*/
				int auxEnd = begin + chunkCPU;
#ifdef DEBUG
				printf("\nDEBUG operator: else auxEnd = %d, chunkCPU = %d \n", auxEnd, chunkCPU);
#endif
				auxEnd = (auxEnd > end) ? end : auxEnd;
				bundle->begin = begin;
				bundle->end = auxEnd;
				begin = auxEnd;
				bundle->type = CPU;
#ifdef DEBUG
				if (bundle->begin < 0)
					printf("\nDEBUG operator: else begin = %d\n", bundle->begin);
				if (bundle->end < 0)
					printf("\nDEBUG operator: else end = %d\n", bundle->begin);
#endif
				return bundle;
			}
		}
		return NULL;
	} // end operator
};

/*MyParallelFilter class is the executor component of the engine, it executes the subrange onto the device selected by SerialFilter*/
template <class B>
class MyParallelFilter : public filter
{
private:
	B *body;

public:
	/*Class' constructor*/
	//template <class B>
	MyParallelFilter(B *b) : filter(false)
	{
		body = b;
	}

	/*Operator function*/
	void *operator()(void *item)
	{
		//variables
		Bundle *bundle = (Bundle *)item;

		totalIterations += bundle->end - bundle->begin;

		if (bundle->type == GPU)
		{
			// GPU WORK
#ifdef DEBUG
			std::cerr << "\nDEBUG::launchGPU(): begin: " << bundle->begin << " end: " << bundle->end << std::endl;
#endif

			tick_count start_tc = tick_count::now();

			// body->sendObjectToGPU(bundle->begin, bundle->end);

			cl::sycl::event event_kernel;

			body->OperatorGPU(bundle->begin, bundle->end, &event_kernel);

			// body->getBackObjectFromGPU(bundle->begin, bundle->end);

			//todo oneAPI
			event_kernel.wait();

#ifdef GET_KERNEL_TIME
			ulong startk = event_kernel.template get_profiling_info<cl::sycl::info::event_profiling::command_start>();
			ulong endk = event_kernel.template get_profiling_info<cl::sycl::info::event_profiling::command_end>();
			float kernel_time = (float)(endk - startk) * 1e-6f; //to milliseconds

			std::cout << "Kernel time: " << kernel_time << " ms\n";
#endif

			// end_tc = tick_count::now();

			float time = (tick_count::now() - start_tc).seconds() * 1000;
			gpuThroughput = (bundle->end - bundle->begin) / time;

			totalIterationsGPU += bundle->end - bundle->begin;

#ifdef DEEP_GPU_REPORT
			deep_gpu_report << bundle->end - bundle->begin << "\t" << gpuThroughput << std::endl;
#endif
#ifdef DEBUG
			std::cerr << "\nDEBUG::launchGPU(): chunk: " << bundle->end - bundle->begin << " TH: " << gpuThroughput << std::endl;
#endif

			/*If CPU has already computed some chunk, then we update fG (factor GPU)*/
			if (cpuThroughput > 0)
			{
				fG = gpuThroughput / cpuThroughput;
			}

			/*To release GPU token*/
			gpuStatus++;
		}
		else
		{ // CPU WORK

#ifdef DEBUG
			std::cerr << "\nDEBUG::launchCPU(): begin: " << bundle->begin << " end: " << bundle->end << std::endl;
#endif
			tick_count start = tick_count::now();
			body->OperatorCPU(bundle->begin, bundle->end);
			tick_count end = tick_count::now();
			float time = (end - start).seconds() / 1000;
			cpuThroughput = (bundle->end - bundle->begin) / time;

			/*If GPU has already computed some chunk, then we update fG (factor GPU)*/
			if (gpuThroughput > 0)
			{
				fG = gpuThroughput / cpuThroughput;
			}
		}
		/*Deleting bundle to avoid memory leaking*/
		delete bundle;
		return NULL;
	}
};
//end class

/*Oracle Class: This scheduler version let us to split the workload in two subranges, one for GPU and one for CPUs*/
class Dynamic : public SchedulerOneApi<Dynamic>
{
	Params *pars;

public:
	/*This constructor just call his parent's contructor*/
	Dynamic(void *params) : SchedulerOneApi(params)
	{
		Params *p = (Params *)params;
		pars = p;

		chunkCPU = 0;
		chunkGPU = p->gpuChunk;
		fG = 0.0;
		gpuThroughput = 0.0;
		cpuThroughput = 0.0;

		//Initializing library PJTRACER
		initializePJTRACER();
	}

	/*Initializes PJTRACER library*/
	void initializePJTRACER()
	{
#ifdef PJTRACER
		char traceFname[1024];
		sprintf(traceFname, "DYNAMIC_C_%d_G_%d.trace", nCPUs, nGPUs);
		tracer = new PFORTRACER(traceFname);
		tracer->beginThreadTrace();
#endif
	}

	/*The main function to be implemented*/
	template <class T>
	void heterogeneous_parallel_for(int begin, int end, T *body)
	{
		/*El tipo T, según aparece en el mainOneApiSchedulers.cpp, es del tipo ValueIterationSchedulerOneApi*/
#ifdef DEBUG
		std::cerr << "Heterogeneous Parallel For Dynamic " << nCPUs << " , " << chunkCPU << " , "
														<< nGPUs << " , " << chunkGPU << std::endl;
#endif
		/*Preparing pipeline*/
		pipeline pipe;
		MySerialFilter serial_filter(begin, end, nCPUs); //create bundle: begin, end, type of device to execute on (CPU or GPU)
		MyParallelFilter<T> parallel_filter(body);		 // <T> = ValueIterationSchedulerOneApi> if bundle type is GPU, execute body->OperadorGPU, else, body->OperatorCPU
		pipe.add_filter(serial_filter);
		pipe.add_filter(parallel_filter);
		chunkCPU = chunkGPU * 0.2;
		fG = 0.0;
		gpuStatus = nGPUs;
		//body->firsttime = true; //todo: firsttime is true after every iteration, not sure I want this, so I initialize it in the ValueIterationSchedulerOneApi constructor

#ifdef DEEP_CPU_REPORT
		char version[50];
		sprintf(version, "_Dynamic_deep_CPU_report_grain_%d.txt", chunkCPU);
		strcat(pars->benchName, version);
		deep_cpu_report.open(pars->benchName, ios::out | ios::app);
#endif
#ifdef DEEP_GPU_REPORT
		char version[100];
		sprintf(version, "%s_Dynamic_deep_GPU_report_grain_%d.txt", pars->benchName, pars->gpuChunk);
		deep_gpu_report.open(version, ios::out | ios::app);
#endif

		/*Seeting a mark to recognize a timestep*/
#ifdef PJTRACER
		tracer->newEvent();
#endif

#ifdef OVERHEAD_ANALYSIS
		end_tc = tick_count::now();
#endif
		/*Run the pipeline*/
		pipe.run(nCPUs + nGPUs);
		pipe.clear();
#ifdef DEEP_CPU_REPORT
		deep_cpu_report.close();
#endif
#ifdef DEEP_GPU_REPORT
		deep_gpu_report.close();
#endif
	}

	/*this function print info to a Log file*/
// 	void saveResultsForBench()
// 	{

// 		char *execution_name = (char *)malloc(sizeof(char) * 100);
// 		sprintf(execution_name, "_Dynamic_OneAPI_%d_%d.txt", nCPUs, nGPUs);
// 		strcat(pars->benchName, execution_name);

// 		/*Checking if the file already exists*/
// 		bool fileExists = isFile(pars->benchName);
// 		ofstream file(pars->benchName, ios::out | ios::app);
// 		if (!fileExists)
// 		{
// 			printHeaderToFile(file);
// 		}
// 		file << nCPUs << "\t" << nGPUs << "\t" << runtime << "\t"
// 			 << getPP0ConsumedJoules(sstate1, sstate2) << "\t" << getPP1ConsumedJoules(sstate1, sstate2) << "\t"
// 			 << getConsumedJoules(sstate1, sstate2) - getPP0ConsumedJoules(sstate1, sstate2) - getPP1ConsumedJoules(sstate1, sstate2) << "\t" << getConsumedJoules(sstate1, sstate2) << "\t"
// 			 << getL2CacheHits(sktstate1[0], sktstate2[0]) << "\t" << getL2CacheMisses(sktstate1[0], sktstate2[0]) << "\t" << getL2CacheHitRatio(sktstate1[0], sktstate2[0]) << "\t"
// 			 << getL3CacheHits(sktstate1[0], sktstate2[0]) << "\t" << getL3CacheMisses(sktstate1[0], sktstate2[0]) << "\t" << getL3CacheHitRatio(sktstate1[0], sktstate2[0]) << "\t"
// 			 << getCyclesLostDueL3CacheMisses(sstate1, sstate2)
// 			 << "\t" << chunkGPU << std::endl;
// 		file.close();
// #ifdef DEBUG
// 		std::cerr << nCPUs << "\t" << nGPUs << "\t" << runtime << "\t"
// 			 << getPP0ConsumedJoules(sstate1, sstate2) << "\t" << getPP1ConsumedJoules(sstate1, sstate2) << "\t"
// 			 << getConsumedJoules(sstate1, sstate2) - getPP0ConsumedJoules(sstate1, sstate2) - getPP1ConsumedJoules(sstate1, sstate2) << "\t" << getConsumedJoules(sstate1, sstate2) << "\t"
// 			 << getL2CacheHits(sktstate1[0], sktstate2[0]) << "\t" << getL2CacheMisses(sktstate1[0], sktstate2[0]) << "\t" << getL2CacheHitRatio(sktstate1[0], sktstate2[0]) << "\t"
// 			 << getL3CacheHits(sktstate1[0], sktstate2[0]) << "\t" << getL3CacheMisses(sktstate1[0], sktstate2[0]) << "\t" << getL3CacheHitRatio(sktstate1[0], sktstate2[0]) << "\t"
// 			 << getCyclesLostDueL3CacheMisses(sstate1, sstate2) << std::endl;
// #endif
// 	}

	//#define PRINTCPUSGPUS
	/*this function print info to a Log file*/
// 	void saveResultsForBenchDen(int inputId, int nrIter, int platformId, float runtime_k1, float energy_k1)
// 	{

// 		char *execution_name = (char *)malloc(sizeof(char) * 100);
// 		sprintf(execution_name, "_GD_OneApi.txt");
// 		strcat(pars->benchName, execution_name);

// 		/*Checking if the file already exists*/
// 		bool fileExists = isFile(pars->benchName);
// 		ofstream file(pars->benchName, ios::out | ios::app);
// 		if (!fileExists)
// 		{
// 			printHeaderToFileDen(file);
// 		}
// 		double eG, eC, eM, eT;
// 		eG = getPP1ConsumedJoules(sstate1, sstate2);
// 		eC = getPP0ConsumedJoules(sstate1, sstate2);
// 		eT = getConsumedJoules(sstate1, sstate2);
// 		eM = eT - eG - eC;
// 		//print implementation id and size

// 		//file << (int)DYNAMIC << "\t" << inputId << "\t" << platformId << "\t"
// 		file << (int)DYNAMIC_ONE_API << "\t" << inputId << "\t"
// #ifdef PRINTCPUSGPUS
// 			 << pars->numcpus << "\t"
// 			 << pars->numgpus << "\t"
// #endif
// 			 << platformId << "\t"
// 			 << chunkGPU << "\t"
// 			 << eG << "\t"
// 			 << eC << "\t"
// 			 << eM << "\t"
// 			 << eT << "\t"
// 			 << runtime << "\t"
// 			 << nrIter << "\t"
// 			 << runtime_k1 << "\t"
// 			 << energy_k1
// 			 << std::endl;
// 		file.close();

// 		cout << "\nDYNAMIC-oneAPI"
// 			 << "\t MDP_SIZE: " << inputId << "\t PLATFORM_ID: " << platformId
// 			 << "\t ratioGPU: " << chunkGPU
// 			 << "\t totalEnergy(Joule): " << eT
// 			 << "\t totalRuntimeValueIteration(s): " << runtime
// 			 << "\tnrIter: " << nrIter
// 			 << "\tgpuRatio: " << (float)totalIterationsGPU / totalIterations
// 			 << "\t runtimeHeterFor(s):  " << runtime_k1
// 			 << "\t energyHeterFor(Joule):  " << energy_k1
// 			 << std::endl;

// #ifdef DEBUG
// 		std::cerr << nCPUs << "\t" << nGPUs << "\t" << runtime << "\t"
// 			 << getPP0ConsumedJoules(sstate1, sstate2) << "\t" << getPP1ConsumedJoules(sstate1, sstate2) << "\t"
// 			 << getConsumedJoules(sstate1, sstate2) - getPP0ConsumedJoules(sstate1, sstate2) - getPP1ConsumedJoules(sstate1, sstate2) << "\t" << getConsumedJoules(sstate1, sstate2) << "\t"
// 			 << getL2CacheHits(sktstate1[0], sktstate2[0]) << "\t" << getL2CacheMisses(sktstate1[0], sktstate2[0]) << "\t" << getL2CacheHitRatio(sktstate1[0], sktstate2[0]) << "\t"
// 			 << getL3CacheHits(sktstate1[0], sktstate2[0]) << "\t" << getL3CacheMisses(sktstate1[0], sktstate2[0]) << "\t" << getL3CacheHitRatio(sktstate1[0], sktstate2[0]) << "\t"
// 			 << getCyclesLostDueL3CacheMisses(sstate1, sstate2) << std::endl;
// #endif
// #ifdef OVERHEAD_ANALYSIS //blows up if this is not defined
// 		cout << "\n KernelEx(ms): " << kernel_execution
// 			 << "\t ThDispatch(ms): " << overhead_td
// 			 << "\t SchedPartitioning(ms): " << overhead_sp
// 			 << "\t H2D(ms): " << overhead_h2d
// 			 << "\t KernelLauch(ms): " << overhead_kl
// 			 << "\t D2H(ms): " << overhead_d2h
// 			 << "\t Total overhead(s): " << (kernel_execution + overhead_td + overhead_sp + overhead_h2d + overhead_kl + overhead_d2h) / 1000
// 			 << std::endl;
// #endif
// 	}

	void printHeaderToFile(ofstream &file)
	{
		file << "N. CPUs" << SEP << "N. GPUs" << SEP << "Time (ms)" << SEP
			 << "CPU Energy(J)" << SEP << "GPU Enegy(J)" << SEP << "Uncore Energy(J)" << SEP << "Total Energy (J)" << SEP
			 << "L2 Cache Hits" << SEP << "L2 Cache Misses" << SEP << "L2 Cache Hit Ratio" << SEP
			 << "L3 Cache Hits" << SEP << "L3 Cache Misses" << SEP << "L3 Cache Hit Ratio" << SEP << "Cycles lost Due to L3 Misses" << SEP
			 << "chunkGPU" << std::endl;
	}
	void printHeaderToFileDen(ofstream &file)
	{
		file << "ALG_ID" << SEP << "INPUT_ID" << SEP << "PLATFORM_ID" << SEP
			 << "GPU_RATIO" << SEP << "E_GPU(J)" << SEP << "E_CPU(J)" << SEP << "E_MEM(J)" << SEP << "E_TOT(J)" << SEP
			 << "TIME(s)" << SEP
			 << "NR_ITER" << SEP
			 << "HET_FOR_TIME(s)" << SEP
			 << "HET_FOR_ENERGY(Joule)"
			 << std::endl;
	}
};
