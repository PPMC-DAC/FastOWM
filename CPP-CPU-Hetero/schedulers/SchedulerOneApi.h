//============================================================================
// Name			: Scheduler.h
// Author		: Antonio Vilches
// Version		: 1.0
// Date			: 26 / 12 / 2014
// Copyright	: Department. Computer's Architecture (c)
// Description	: Main scheduler interface class
//============================================================================

// #define ENERGYCOUNTERS

#include <cstdlib>
#include <iostream>
#include <fstream>

#include <CL/cl.h>

#ifdef ENERGYCOUNTERS
#ifdef Win32
#include "PCM_Win/windriver.h"
#else
#include "cpucounters.h"
#endif
#endif

// #include "../utils/cl_helper.h"

// #include "tbb/task_scheduler_init.h"
#include "tbb/global_control.h"
#include "tbb/tick_count.h"


#ifdef PJTRACER
#include "pfortrace.h"
#endif


#include "CL/sycl.hpp"

using namespace cl::sycl;

using namespace std;
using namespace tbb;




/*****************************************************************************
 * types
 * **************************************************************************/
typedef struct{
	int numcpus;
	int numgpus;
	int gpuChunk;
	// char benchName[256];
//	char kernelName[50];
	// float ratioG;
}Params;

/*****************************************************************************
 * Global Variables For OpenCL
 * **************************************************************************/

// int error;
// cl_uint num_max_platforms;
// cl_uint num_max_devices;
// cl_uint num_platforms;
// cl_uint num_devices;
// cl_platform_id* platforms_id;
// cl_device_id device_id;
// cl_context context;
// cl_command_queue command_queue;
// cl_program program;
// cl_kernel kernel;
// int computeUnits = 48;
// int minChunkGPU;
// size_t vectorization;


//profiler
#ifdef PJTRACER
PFORTRACER * tracer;
#endif

/****************************************************************************
 * Base Scheduler class
 * **************************************************************************/

/*This Scheduler Base class implementation follows a Singleton pattern*/
template <class T>
class SchedulerOneApi{
protected:
// Class members
	//Scheduler Itself
	static T *instance;
	// task_scheduler_init *init;
	int nCPUs;
	int nGPUs;

	//Energy Counters and timing
#ifdef ENERGYCOUNTERS
	PCM * pcm;
	vector<CoreCounterState> cstates1, cstates2;
	vector<SocketCounterState> sktstate1, sktstate2;
	SystemCounterState sstate1, sstate2;
	//timing
	tick_count start, end;
	float runtime;
#endif

	//OneApi

    //auto context = command_queue.get_context();

//End class members

	/*Scheduler Constructor, forbidden access to this constructor from outside*/
    SchedulerOneApi(void * params) {
		Params * p = (Params *) params;
		nCPUs = p->numcpus;
		nGPUs = p->numgpus;

		// global_control c(global_control::max_allowed_parallelism, nCPUs + nGPUs);
		// init = new task_scheduler_init(nCPUs + nGPUs);

#ifdef DEBUG
		cerr << "TBB scheduler is active " << "(" << nCPUs << ", " << nGPUs << ")" << endl;
#endif

// #ifdef DEBUG
// 		cerr << "iNITIALIZING OPENCL" << endl;
// #endif

#ifdef ENERGYCOUNTERS
	#ifdef DEBUG
		cerr << "INITIALIZING pcm" << endl;
	#endif
		initializePCM();
		runtime = 0.0;
#endif

#ifdef HOST_PRIORITY  // rise GPU host-thread priority
		initializeHOSTPRI();
#endif		

	}


//     SchedulerOneApi() {

// #ifdef ENERGYCOUNTERS
// 	#ifdef DEBUG
// 		cerr << "INITIALIZING pcm" << endl;
// 	#endif
// 		initializePCM();
// #endif

// #ifdef HOST_PRIORITY  // rise GPU host-thread priority
// 		initializeHOSTPRI();
// #endif		
// 		runtime = 0.0;
// 	}



#ifdef HOST_PRIORITY  // rise GPU host-thread priority

	void initializeHOSTPRI(){

	#ifdef DEBUG
		cerr << "INITIALIZING HOSTPRIORITY" << endl;
	#endif

		unsigned long dwError, dwThreadPri;
		//dwThreadPri = GetThreadPriority(GetCurrentThread());
		//printf("Current thread priority is 0x%x\n", dwThreadPri);
		if(!SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL)){
			dwError = GetLastError();
			if(dwError){
				cerr << "Failed to set hight priority to host thread (" << dwError << ")" << std::endl;
			}
		}
	//dwThreadPri = GetThreadPriority(GetCurrentThread());
	//printf("Current thread priority is 0x%x\n", dwThreadPri);
	}

#endif


public:
	/*Class destructor*/
	~SchedulerOneApi(){
		// init->~task_scheduler_init();
		delete instance;
		instance = NULL;
	}

	/*This function creates only one instance per process, if you want a thread safe behavior protect the if clausule with a Lock*/
	static T * getInstance(void * params){
		if(! instance){
			instance = new T(params);
		}
		return instance;
	}

	/*Initializing PCM library*/
#ifdef ENERGYCOUNTERS
	void initializePCM(){
		/*This function prints lot of information*/
		pcm = PCM::getInstance();
		pcm->resetPMU();
		if (pcm->program() != PCM::Success){
			cerr << "Error in PCM library initialization" << std::endl;
			exit(-1);
		}
	}

	/*Sets the start mark of energy and time*/
	void startTimeAndEnergy(){

		pcm->getAllCounterStates(sstate1, sktstate1, cstates1);

		start = tick_count::now();
	}

	/*Sets the end mark of energy and time*/
	void endTimeAndEnergy(){
		end = tick_count::now();

		pcm->getAllCounterStates(sstate2, sktstate2, cstates2);

		runtime = (end-start).seconds();
	}
#endif

	/*Checks if a File already exists*/
	bool isFile(char *filename){
		//open file
		ifstream ifile(filename);
		return !ifile.fail();
	}
};

template <class T>
T* SchedulerOneApi<T>::instance = NULL;
