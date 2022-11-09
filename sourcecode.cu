#define Original_isGPU 1//1-Gpu, 0-cpu
#define Original_isFull 1 //1-full, 0-red

#define Original_isRev 0  //if(!full): then, 1-reverse(Dropoff first, then pickup),0-NotReverse.
#define SecondInsertionInCPU 1 //1-Yes, 0-No

// Run out of memory is possible due to two reasons: Expectation & Edge[2*n][2*n] array. Interactive gpu.
// NOTE: change TW in generate_samples() for each problem.

#define TotalRequests 120
#define TotalVehicles 11

#define Start_Window 200//330
#define End_Window 512//590

#define Str1 "./problems/general/pr05.txt"//a5-50.txt"
#define DeclaredProblem Str1
#define Defined_MaxTemperature 20

#define ExpectedScenarios 20
#define DegreeofDynamism 0.5 //(0,1)

#define Expectation 8500

#define TotalThreads 32 //max:1024
#define Expected_Blocks ((Expectation+T-1)/T)

#define Relocation_Strategy 0

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <ctime>
#include <algorithm>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/count.h>
#include <thrust/remove.h>
#include <iterator>
#include <cuda.h>
#include "cuda_runtime.h"
#include "cuda_occupancy.h"
#include "device_launch_parameters.h"
#include "cuda_fp16.h"
#include <cuda_runtime.h>
#include <memory>
#include <random>
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
//#define RandomNumber(Min,Max) ((float(rand()) / float(RAND_MAX)) * (Max - Min)) + Min;
#define RandomNumber(min,max) rand() % (max - min + 1) + min;
#define RandomPick(b,a) ((b - a) * ((float)rand() / RAND_MAX)) + a //b=MAX, a=MIN // returns in-between numbers as well

//#include <cooperative_groups.h>
//using namespace cooperative_groups;

using namespace std;

#define n TotalRequests
#define m TotalVehicles
#define Var	2*(TotalRequests+TotalVehicles)
#define TotalNodes 2*n

//Class generation and Kernel call that maps __device__ function inside class to be called for every element of an array (threads)
//Initial commit
#define CHECK(r) {_check((r), __LINE__);}
#define Blocks 1
#define Threads 1
#define ExpectedPath 2*TotalRequests

#define w1 1//8
#define w2 0//3
#define w3 0//1
#define w4 0//1
#define w5 0//n

#define WarpSize 32
#define AutoBlock(a,t) (round)(a/t)+1

unsigned int time_now = 0;

#ifdef __CUDACC__
#define LaunchBound(x,y) __launch_bounds__(x,y)
#else
#define LaunchBound(x,y)
#endif

#if defined(__CUDACC__) // NVCC
#define MY_ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
#define MY_ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
#define MY_ALIGN(n) __declspec(align(n))
#else
#error "Please provide a definition for MY_ALIGN macro for your host compiler!"
#endif

#ifdef __CUDACC__
#define KERNEL_ARGS2(grid, block) <<< grid, block >>>
#define KERNEL_ARGS3(grid, block, sh_mem) <<< grid, block, sh_mem >>>
#define KERNEL_ARGS4(grid, block, sh_mem, stream) <<< grid, block, sh_mem, stream >>>
#else
#define KERNEL_ARGS2(grid, block)
#define KERNEL_ARGS3(grid, block, sh_mem)
#define KERNEL_ARGS4(grid, block, sh_mem, stream)
#endif

#pragma region CudaErrorCheck Function
void _check(cudaError_t r, int line) {
	if (r != cudaSuccess) {
		printf("CUDA error on line %d: %s\n", line, cudaGetErrorString(r), line);
		exit(0);
	}
}
#pragma endregion

struct Managed
{
	void *operator new(size_t len) {
		void *ptr;
		CHECK(cudaMallocManaged(&ptr, len));
		return ptr;
	}
	void operator delete(void *ptr) {
		CHECK(cudaFree(ptr));
	}
};

int CurrentTime = 0;
void PAUSER(int HardStop = 0)
{
	if (0)//HardStop == 4000)//(CurrentTime > 170 && CurrentTime < 172) || HardStop > 8000)//HardStop > (n-2))// || HardStop == 3000)// || (CurrentTime > 170 && CurrentTime < 172) || (CurrentTime > 220))
	{
		printf("\nSystem is PAUSED: Enter 1 to continue:\n");
		int a;
		cin >> a;

		if (a == 0) // exit-code
			exit(0);
	}
}

//////////////////////////////////////////////////////

#include "utility.cuh"

#include "objectiveFunction.cuh"

#include "problem.cuh"

#include "cudaObjects.cuh"

#include "solution.cuh"

//#include "MED_solution.cuh"

#include "bootProblem.cuh"

#include "cudakernels.cuh"

//#include "operators.cuh"
#include "operators_simultaneous.cuh"

//#include "explore_neighborhood.cuh"
#include "explore_neighborhood_simultaneous.cuh"

//#include "lns.cuh"
#include "lns_simultaneous.cuh"

//////////////////////////////////////////////////////

int main()
{
	//File creation for data tracking
	//ofstream out_data("filename.txt");

	printf(".......THE PROGRAM HAS STARTED......\n");
	printf("Original_isGPU=%d, Original_isFull=%d, Original_isRev=%d\n",
			Original_isGPU, Original_isFull, Original_isRev);
	printf("Expectation: %d, ExpectedScenarios: %d\n", Expectation, ExpectedScenarios);

	for (int i = 0; i < 1; i++)
	{

		GenerateScenarios GS;
		GS.Creation();
		GS.Optimize(10000, 6000);
		GS.GPUInformation();
		GS.Destruction();
		CHECK(cudaDeviceReset());

		printf("\nABOVE IS THE RUN: %d\n\n", i + 1);
	}
	printf("The Program Ends here..\n");


	return 0;
}
