
struct GenerateScenarios
{
	clock_t TotalTime;

	bool DO_intra = false;
	bool Trigger_restart_mechanism = false; // CAUTION: THIS GIVES RUN-TIME ERROR FOR SOME INSTANCES!!!

	// parameters
	float *temperature;
	float *maxTemperature;
	float *minTemperature;
	float *coolingRate;
	int *noImprovement; //multiple
	int *noImprovementLimit; //multiple
	FSM *fsm;

	int *existence[ExpectedScenarios]; //multiple
	float *best_cost; //multiple
	solution *Best_Solution; //multiple
	problem *Problem; //multiple
	solution *Solution; //multiple
	coefficient *Coefficient; //multiple
	explore_neighborhood *Explore_Neighborhood; //multiple

	reinsertion *Reinsertion;
	removal *Removal;
	worstcostRemoval *WorstCostRemoval;
	noise *Noise;
	insertion **Insertion;

	myList *AvailableRequests;
	myList *UnAvailableRequests;
	myList *ServedRequests;
	myList *ServedRequests_BoundaryMark;
	myList *RejectedRequests;
	myList *SampledRequests;

	solution *Master_Solution;
	problem *Master_Problem;
	problem *UnModified_Problem;

	int **req_list; //utility
	int **temp_req_list; //utility

	int *RequestType[ExpectedScenarios]; //0-sampled, 1-real, 2-unmaterialized
	int Total_Static_Requests;
	int Total_Accepted_Dynamic_Requests = 0;
	int Total_Rejected_Dynamic_Requests = 0;
	int Idle_Vehicle_List[m];

    int *Do_Not_Modify_Scenario;

	//int CurrentTime = 0;
	int Next_Injection_Time = 0;
	int Next_Dynamic_Request;

    int All_Dynamic_Requests_Inserted = 0;
    int Stop_Checking = 0;

	void Creation()
	{

		// import Problem
		Problem = new problem[ExpectedScenarios];
		Master_Problem = new problem;
		UnModified_Problem = new problem;
		Import_Problem();

		AvailableRequests = new myList[ExpectedScenarios];
		UnAvailableRequests = new myList[ExpectedScenarios];
		ServedRequests = new myList;
		ServedRequests_BoundaryMark = new myList;
		RejectedRequests = new myList;
		SampledRequests = new myList[ExpectedScenarios];

		Master_Solution = new solution;


		req_list = new int*[ExpectedScenarios];
		temp_req_list = new int*[ExpectedScenarios];
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			req_list[s] = new int[n];
			temp_req_list[s] = new int[n];
		}


		Reinsertion = new reinsertion;
		Removal = new removal[ExpectedScenarios];
		WorstCostRemoval = new worstcostRemoval;
		WorstCostRemoval[0].create();
		Noise = new noise;

		Insertion = new insertion*[ExpectedScenarios];
		for (int s = 0; s < ExpectedScenarios; s++)
			Insertion[s] = new insertion[n];

		Explore_Neighborhood = new explore_neighborhood;
		Explore_Neighborhood[0].creation(); //cuda-mallocManaged is used here!! BEWARE!!

		// parameters
		temperature = new float[ExpectedScenarios];
		maxTemperature = new float[ExpectedScenarios];
		minTemperature = new float[ExpectedScenarios];
		coolingRate = new float[ExpectedScenarios];
		fsm = new FSM;
		noImprovement = new int[ExpectedScenarios];
		noImprovementLimit = new int[ExpectedScenarios];
        Do_Not_Modify_Scenario = new int[ExpectedScenarios];
		Clear_Do_Not_Modify_Scenario();

		for (int s = 0; s < ExpectedScenarios; s++)
			existence[s] = new int[n];
		for (int s = 0; s < ExpectedScenarios; s++)
			for (int i = 0; i < n; i++)
				existence[s][i] = 0;
		best_cost = new float[ExpectedScenarios];
		Best_Solution = new solution[ExpectedScenarios];
		Solution = new solution[ExpectedScenarios];
		Coefficient = new coefficient[ExpectedScenarios];


		for (int s = 0; s < ExpectedScenarios; s++)
			RequestType[s] = new int[n];


		Parameterization();
		Construct_Initial_Solution();

		///////////////////////////////////////////////////////
		/////// Just printing the total GPU memory usage //////
		///////////////////////////////////////////////////////

		size_t free, total;
		CHECK(cudaMemGetInfo(&free,&total));
		printf("GPU Memory usage info:\n Total = %0.02f GB, Used = %0.02f MB\n", ((float)total/1000000000), 1000* (((float)total/1000000000) - ((float)free)/(1000000000)));
		printf("Disclaimer: Make sure Expectation is greater than .allinsertions!\n");
	}

	void Destruction()
	{
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			delete[] req_list[s];
			delete[] temp_req_list[s];
		}

		delete[] AvailableRequests;
		delete[] UnAvailableRequests;
		delete ServedRequests;
		delete ServedRequests_BoundaryMark;
		delete RejectedRequests;
		delete[] SampledRequests;

		delete Master_Solution;
		delete Master_Problem;
		delete UnModified_Problem;

		delete[] req_list;
		delete[] temp_req_list;

		delete Reinsertion;
		delete[] Removal;
		WorstCostRemoval[0].destroy();
		delete WorstCostRemoval;
		delete Noise;

		for (int s = 0; s < ExpectedScenarios; s++)
			delete[] Insertion[s];
		delete[] Insertion;

		//below-three must be in this order.
		Explore_Neighborhood[0].destroy(); //cuda-mallocManaged is used here!! BEWARE!!
		delete Explore_Neighborhood;

		delete[] temperature;
		delete[] maxTemperature;
		delete[] minTemperature;
		delete[] coolingRate;
		delete fsm;
		delete[] noImprovement;
		delete[] noImprovementLimit;
        delete[] Do_Not_Modify_Scenario;

		for(int s = 0; s < ExpectedScenarios; s++)
			delete[] existence[s];

		delete[] best_cost;
		delete[] Best_Solution;
		delete[] Problem;
		delete[] Solution;
		delete[] Coefficient;

		for (int s = 0; s < ExpectedScenarios; s++)
			delete[] RequestType[s];

	}

	void Clear_Do_Not_Modify_Scenario()
	{
	    for (int s = 0; s < ExpectedScenarios; s++)
	        Do_Not_Modify_Scenario[s] = 0;
	}

	void GPUInformation()
	{
		printf("%s Starting...\n\n");
		printf(" CUDA Device Query (Runtime API) version (CUDART static linking)\n\n");

		int deviceCount = 0;
		cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

		if (error_id != cudaSuccess)
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
			printf("Result = FAIL\n");
			exit(EXIT_FAILURE);
		}

		// This function call returns 0 if there are no CUDA capable devices.
		if (deviceCount == 0)
		{
			printf("There are no available device(s) that support CUDA\n");
		}
		else
		{
			printf("Detected %d CUDA Capable device(s)\n", deviceCount);
		}

		int dev, driverVersion = 0, runtimeVersion = 0;

		for (dev = 0; dev < deviceCount; ++dev)
		{
			cudaSetDevice(dev);
			cudaDeviceProp deviceProp;
			cudaGetDeviceProperties(&deviceProp, dev);

			printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

			// Console log
			cudaDriverGetVersion(&driverVersion);
			cudaRuntimeGetVersion(&runtimeVersion);
			printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n", driverVersion / 1000, (driverVersion % 100) / 10, runtimeVersion / 1000, (runtimeVersion % 100) / 10);
			printf("  CUDA Capability Major/Minor version number:    %d.%d\n", deviceProp.major, deviceProp.minor);

			char msg[256];
			printf("  Total amount of global memory:                 %.0f GB (%llu bytes)\n",
				(float)deviceProp.totalGlobalMem / 1048576.0f / 1024, (unsigned long long) deviceProp.totalGlobalMem);
			//printf("%s", msg);

			printf("  Streamings Multiprocessors                     (%2d)\n", deviceProp.multiProcessorCount);

			//printf("  GPU Max Clock rate:                            %.0f MHz (%0.2f GHz)\n", deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);


#if CUDART_VERSION >= 5000
			// This is supported in CUDA 5.0 (runtime API device properties)
			/*
			printf("  Memory Clock rate:                             %.0f Mhz\n", deviceProp.memoryClockRate * 1e-3f);
			printf("  Memory Bus Width:                              %d-bit\n", deviceProp.memoryBusWidth);

			if (deviceProp.l2CacheSize)
			{
			printf("  L2 Cache Size:                                 %d bytes\n", deviceProp.l2CacheSize);
			}
			*/
#else
			// This only available in CUDA 4.0-4.2 (but these were only exposed in the CUDA Driver API)
			int memoryClock;
			getCudaAttribute<int>(&memoryClock, CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE, dev);
			printf("  Memory Clock rate:                             %.0f Mhz\n", memoryClock * 1e-3f);
			int memBusWidth;
			getCudaAttribute<int>(&memBusWidth, CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH, dev);
			printf("  Memory Bus Width:                              %d-bit\n", memBusWidth);
			int L2CacheSize;
			getCudaAttribute<int>(&L2CacheSize, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, dev);

			if (L2CacheSize)
			{
				printf("  L2 Cache Size:                                 %d bytes\n", L2CacheSize);
			}

#endif

			/*printf("  Maximum Texture Dimension Size (x,y,z)         1D=(%d), 2D=(%d, %d), 3D=(%d, %d, %d)\n",
			deviceProp.maxTexture1D, deviceProp.maxTexture2D[0], deviceProp.maxTexture2D[1],
			deviceProp.maxTexture3D[0], deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
			printf("  Maximum Layered 1D Texture Size, (num) layers  1D=(%d), %d layers\n",
			deviceProp.maxTexture1DLayered[0], deviceProp.maxTexture1DLayered[1]);
			printf("  Maximum Layered 2D Texture Size, (num) layers  2D=(%d, %d), %d layers\n",
			deviceProp.maxTexture2DLayered[0], deviceProp.maxTexture2DLayered[1], deviceProp.maxTexture2DLayered[2]);


			printf("  Total amount of constant memory:               %lu bytes\n", deviceProp.totalConstMem);
			printf("  Total amount of shared memory per block:       %lu bytes\n", deviceProp.sharedMemPerBlock);
			printf("  Total number of registers available per block: %d\n", deviceProp.regsPerBlock);
			printf("  Warp size:                                     %d\n", deviceProp.warpSize);
			*/
			printf("  Maximum number of threads per multiprocessor:  %d\n", deviceProp.maxThreadsPerMultiProcessor);
			/*
			printf("  Maximum number of threads per block:           %d\n", deviceProp.maxThreadsPerBlock);
			printf("  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n",
			deviceProp.maxThreadsDim[0],
			deviceProp.maxThreadsDim[1],
			deviceProp.maxThreadsDim[2]);
			printf("  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n",
			deviceProp.maxGridSize[0],
			deviceProp.maxGridSize[1],
			deviceProp.maxGridSize[2]);
			printf("  Maximum memory pitch:                          %lu bytes\n", deviceProp.memPitch);
			printf("  Texture alignment:                             %lu bytes\n", deviceProp.textureAlignment);
			printf("  Concurrent copy and kernel execution:          %s with %d copy engine(s)\n", (deviceProp.deviceOverlap ? "Yes" : "No"), deviceProp.asyncEngineCount);
			printf("  Run time limit on kernels:                     %s\n", deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
			printf("  Integrated GPU sharing Host Memory:            %s\n", deviceProp.integrated ? "Yes" : "No");
			printf("  Support host page-locked memory mapping:       %s\n", deviceProp.canMapHostMemory ? "Yes" : "No");
			printf("  Alignment requirement for Surfaces:            %s\n", deviceProp.surfaceAlignment ? "Yes" : "No");
			*/
			printf("  Device has ECC support:                        %s\n", deviceProp.ECCEnabled ? "Enabled" : "Disabled");
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
			printf("  CUDA Device Driver Mode (TCC or WDDM):         %s\n", deviceProp.tccDriver ? "TCC (Tesla Compute Cluster Driver)" : "WDDM (Windows Display Driver Model)");
#endif
			printf("  Device supports Unified Addressing (UVA):      %s\n", deviceProp.unifiedAddressing ? "Yes" : "No");
			/*
			printf("  Device PCI Domain ID / Bus ID / location ID:   %d / %d / %d\n", deviceProp.pciDomainID, deviceProp.pciBusID, deviceProp.pciDeviceID);

			const char *sComputeMode[] =
			{
			"Default (multiple host threads can use ::cudaSetDevice() with device simultaneously)",
			"Exclusive (only one host thread in one process is able to use ::cudaSetDevice() with this device)",
			"Prohibited (no host thread can use ::cudaSetDevice() with this device)",
			"Exclusive Process (many threads in one process is able to use ::cudaSetDevice() with this device)",
			"Unknown",
			NULL
			};
			printf("  Compute Mode:\n");
			printf("     < %s >\n", sComputeMode[deviceProp.computeMode]);
			*/
		}
	}

	void Parameterization()
	{
		// 1) Parameterization

				for (int s = 0; s < ExpectedScenarios; s++)
				{
					Coefficient[s].initialize();

					// SA-parameters
					coolingRate[s] = 0.99975;
					minTemperature[s] = 0.001;
					maxTemperature[s] = Defined_MaxTemperature; //20
					temperature[s] = maxTemperature[s];

					// restart parameter
					noImprovement[s] = 0;
					noImprovementLimit[s] = 200;
				}

				// fsm parameter
				fsm[0].initialize();
				fsm[0].FSM_Enable = false;
	}

	void Import_Problem()
	{
		// 2) Import the problem instance

				for (int s = 0 ; s < ExpectedScenarios; s++)
					BootProblem(Problem[s]);
				BootProblem(Master_Problem[0]);
				BootProblem(UnModified_Problem[0]);

				/*if (1)
				{
					//printf("%d %d %d %d %d\n", m, 2*n, temp[0], temp[1], temp[2]);

					//printf("0	%0.02f	%0.02f	0	0	%d	%d\n", Depot[0].x, Depot[0].y, temp[3], temp[4]);

					for (int i = 0; i < TotalRequests; i++)
					{
						Problem[0].Request[i].pickup.Location.ID = i+1;
						printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Problem[0].Request[i].pickup.Location.ID, Problem[0].Request[i].pickup.Location.x, Problem[0].Request[i].pickup.Location.y, Problem[0].Request[i].pickup.ServiceTime,
								Problem[0].Request[i].load, Problem[0].Request[i].pickup.EarliestTime, Problem[0].Request[i].pickup.LatestTime);
					}

					for (int i = 0; i < TotalRequests; i++)
					{
						Problem[0].Request[i].dropoff.Location.ID = n+i+1;
						printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Problem[0].Request[i].dropoff.Location.ID, Problem[0].Request[i].dropoff.Location.x, Problem[0].Request[i].dropoff.Location.y, Problem[0].Request[i].dropoff.ServiceTime,
							0, Problem[0].Request[i].dropoff.EarliestTime, Problem[0].Request[i].dropoff.LatestTime);
					}
				}*/

				if (1)
				{
					printf("--------MASTER PROBLEM---------\n");
					for (int i = 0; i < TotalRequests; i++)
					{
						Master_Problem[0].Request[i].pickup.Location.ID = i+1;
						printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Master_Problem[0].Request[i].pickup.Location.ID, Master_Problem[0].Request[i].pickup.Location.x, Master_Problem[0].Request[i].pickup.Location.y, Master_Problem[0].Request[i].pickup.ServiceTime,
								Master_Problem[0].Request[i].load, Master_Problem[0].Request[i].pickup.EarliestTime, Master_Problem[0].Request[i].pickup.LatestTime);
					}

					for (int i = 0; i < TotalRequests; i++)
					{
						Master_Problem[0].Request[i].dropoff.Location.ID = n+i+1;
						printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Master_Problem[0].Request[i].dropoff.Location.ID, Master_Problem[0].Request[i].dropoff.Location.x, Master_Problem[0].Request[i].dropoff.Location.y, Master_Problem[0].Request[i].dropoff.ServiceTime,
							0, Master_Problem[0].Request[i].dropoff.EarliestTime, Master_Problem[0].Request[i].dropoff.LatestTime);
					}
					printf("-------------------------------\n");
				}


				// WHAT IF I SIMPLY MODIFY THE PROBLEM?

				//Store Problem in GPU's Constant Memory
				//Update_ConstantMemory(Problem);
				//CHECK(cudaDeviceSynchronize());

	}

	void Construct_Initial_Solution()
	{
		// 2) Construct an initial solution

				//printf("Randomized_construction_heuristic is used!!\n");

				for (int s = 0 ; s < ExpectedScenarios; s++)
				{
					Solution[s].Problem.Copy(Problem[s]);

					if (0)
						IH_Construction_heuristic(s);
					else
						Solution[s].RandomShuffle();

					Solution[s].IdentifyCurrvehicle();
					Solution[s].OfflineEvaluateCost();

					// Update existence
					Update_Existence(Solution, s);
				}

				for (int s = 0 ; s < ExpectedScenarios; s++)
				{
					// Update
					int dummy;
					best_cost[s] = 0;
					//Test_UpdateBestSolution(0, Solution, dummy, -1, -1);
				}


		// noise mechanism (can be set only after booting problem)
				Noise[0].update_noise_parameters(Solution);

		// Fill candidate edges
				//Explore_Neighborhood[0].fill_Candidate_Edges(Problem);

	}

	void Maintain_Master_Solution_for_FirstTime(int Master_k, int firstTime = 0)
	{
		//update for the first time
		// 1) find best solution for static requests
		int best_s = 0; float best_value = 0;
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			float current_value = Coefficient[s].cost_function(Best_Solution[s].Cost);
			if (current_value != 0)
				if (best_value == 0 || current_value < best_value)
				{
					best_value = current_value;
					best_s = s;
				}
		}

		printf("BEFORE:\n");
		for (int s = 0; s < ExpectedScenarios; s++)
			printf("best_cost[%d] = %0.02f (n=%d)\n", s, Coefficient[s].cost_function(Best_Solution[s].Cost), Coefficient[s].cost_function(Best_Solution[s].Cost) == 0 ? 0 : Best_Solution[s].total_req_served);


		// 2) copy best static solution to all scenarios
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			best_cost[s] = best_value;

			Best_Solution[s].Copy(Best_Solution[best_s]);
			Best_Solution[s].OfflineEvaluateCost();

			Solution[s].Copy(Best_Solution[best_s]);
			Solution[s].OfflineEvaluateCost();

			//printf("s = %d cost = %0.02f\n", s, Coefficient[s].cost_function(Best_Solution[s].Cost));
		}

		printf("AFTER:\n");
		for (int s = 0; s < ExpectedScenarios; s++)
			printf("best_cost[%d] = %0.02f (n=%d)\n", s, Coefficient[s].cost_function(Solution[s].Cost), Coefficient[s].cost_function(Solution[s].Cost) == 0 ? 0 : Solution[s].total_req_served);


		// 3) Initialization of ServedRequests_BoundaryMark
		for (int k = 0; k < m; k++)
		{
			ServedRequests_BoundaryMark[0].push(0);
			Idle_Vehicle_List[k] = 0;
		}

		// 4) Update Served Requests: dispatch vehicles to first customers
		for (int k = 0; k < m; k++)
		{
			Master_Solution[0].Vehicle[k].size = 0;
			if (Best_Solution[best_s].Vehicle[k].size > 0)
			{
				int req = Best_Solution[best_s].Vehicle[k].path[0] - 1;

				// 1) Update schedules in Master solution
				Master_Solution[0].Vehicle[k].path[0] = req + 1;
				Master_Solution[0].Vehicle[k].size++;

				Master_Solution[0].Request[req].B = Best_Solution[best_s].Request[req].B;
				Master_Solution[0].Request[req].C = Best_Solution[best_s].Request[req].C;

				ServedRequests_BoundaryMark[0].list[k] += 1;
				ServedRequests[0].push(req);

				for (int s = 0; s < ExpectedScenarios; s++)
					AvailableRequests[s].popout(req);


				// 2) Adjust Problems in effect with the scheduled service locations
				if (req < n)
				{
					Master_Problem[0].Request[req].pickup.EarliestTime = Master_Solution[0].Request[req].B ;
					Master_Problem[0].Request[req].pickup.LatestTime = Master_Solution[0].Request[req].C;
				}
				else
				{
					Master_Problem[0].Request[req - n].dropoff.EarliestTime = Master_Solution[0].Request[req].B ;
					Master_Problem[0].Request[req - n].dropoff.LatestTime = Master_Solution[0].Request[req].C;
				}
			}
		}




		// print utility
		for (int s = 0; s < ExpectedScenarios; s++)
			Best_Solution[s].printinput();

		for (int j = 0; j < ServedRequests[0].size; j++)
			printf("%d - ", ServedRequests[0].list[j]);
		printf("\n");

		for (int k = 0; k < m; k++)
			printf("%d | ", ServedRequests_BoundaryMark[0].list[k]);
		printf("\n");

		printf("Master Solution:");
		Master_Solution[0].printinput();
		Master_Solution[0].printinput_with_details(Master_Problem);



	}

	void Maintain_Master_Solution(int Master_k)
	{
		// 1) retrieve scenario solutions
		for (int s = 0; s < ExpectedScenarios; s++)
		{

			if (best_cost[s] > 0)
			{
				Solution[s].Copy(Best_Solution[s]);
				Solution[s].OfflineEvaluateCost();
			}
			else
				printf("SCENARIO %d is In-Feasible for Maintain_Master_Solution at Line 576!!\n", s);

			//printf("s = %d cost = %0.02f\n", s, Coefficient[s].cost_function(Best_Solution[s].Cost));
		}

		// 2) Invoke Consensus Algorithm: to select the next service location for Master_k
		Consensus_Function(Master_k);
		Reschedule_Master_Solution();

		printf("CHECKING UPDATE: MASTER SOLUTION\n");
		Master_Solution[0].printinput_with_details(Master_Problem);
		PAUSER();

		printf("Master Problem:\n");
		Master_Problem[0].print();

		// 3) Update the Scenario solutions from the Master Solution
		for (int s = 0; s < ExpectedScenarios; s++)
		    Modify_Scenarios_From_Master_Solution(s);
		Clear_Do_Not_Modify_Scenario();

		//print utility
		for (int s = 0; s < ExpectedScenarios; s++)
			Solution[s].printinput();

		printf("SERVED REQUESTS: "); ServedRequests[0].simple_print();
		PAUSER(3000);


	}

    void Update_Do_Not_Modify_Scenario(int best_node, int k)
    {
        int Best_UnModified_Scenario = 0;
        float Min_Cost = 0;
        if (1)//isRealRequest)
            for (int s = 0; s < ExpectedScenarios; s++)
                if (best_cost[s] > 0) // DON'T consider infeasible solutions for score update
                {
                    int index = ServedRequests_BoundaryMark[0].list[k];
                    if (index > 0 && index < Best_Solution[s].Vehicle[k].size)
                    {
                        int node = Best_Solution[s].Vehicle[k].path[index] - 1;

                        if (node == best_node)
                        {
                            Do_Not_Modify_Scenario[s] = 1;
                            if (Min_Cost == 0 || Coefficient[s].cost_function(Best_Solution[s].Cost) < Min_Cost)
                            {
                                Min_Cost = Coefficient[s].cost_function(Best_Solution[s].Cost);
                                Best_UnModified_Scenario = s;
                            }

                            printf("UnModified Scenario %d:\n", s); Best_Solution[s].printinput();
                            printf("Cost: %0.02f:\n", Coefficient[s].cost_function(Best_Solution[s].Cost));
                        }
                        else
                        {
                            printf("OTHER Scenario %d:\n", s); Best_Solution[s].printinput();
                            printf("Cost: %0.02f:\n", Coefficient[s].cost_function(Best_Solution[s].Cost));

                        }
                    }
                }

        for (int s = 0; s < ExpectedScenarios; s++)
            if (Do_Not_Modify_Scenario[s] == 1)
                printf("Do_Not_Modify_Scenario[%d] = %d\n", s, Do_Not_Modify_Scenario[s]);

        printf("Best UnModified Scenario %d:\n", Best_UnModified_Scenario);

        PAUSER(3000);
    }

	void Consensus_Function(int Master_k)
	{
		// Select next service location for $k$ and put into Master_Solution

		int *Score = new int[2*n];
		int k = Master_k;
        int best_score = 0;
        int best_node = 0;
		//for (int k = 0; k < TotalVehicles; k++)
		{
			if (Idle_Vehicle_List[k] == 1)
			{
				// 1) reset score for each service location
				for (int i = 0; i < 2*n; i++)
					Score[i] = 0;

				// 2) update score
				for (int s = 0; s < ExpectedScenarios; s++)
				{
					if (best_cost[s] > 0) // DON'T consider infeasible solutions for score update
					{
						int index = ServedRequests_BoundaryMark[0].list[k];
						if (index > 0 && index < Best_Solution[s].Vehicle[k].size)
						{
							int node = Best_Solution[s].Vehicle[k].path[index] - 1;

							// 2.2) ensure nodes are followed by served node of that vehicle
							// unnecessary if statement
							if (Master_Solution[0].Vehicle[k].path[index-1] == Best_Solution[s].Vehicle[k].path[index-1])
							{
								if (node < n)
								{
									if (Relocation_Strategy)
										Score[node] += 1;
									else if (RequestType[s][node] == 1) //0-sampled,1-real,2-unmaterialized
										Score[node] += 1;
								}
								else
								{
									if (Relocation_Strategy)
										Score[node] += 1;
									else if (RequestType[s][node-n] == 1) //0-sampled,1-real,2-unmaterialized
										Score[node] += 1;
								}
							}
						}
					}
				}


				// 3) find best scored service location for this vehicle k
				//printf("FOR VEHICLE k=%d:\n", k);
				for (int i = 0; i < 2*n; i++)
				{
					printf("Score[%d]=%d\n", i, Score[i]);
					if (best_score == 0 || Score[i] > best_score)
					{
						best_score = Score[i];
						best_node = i;
					}
				}

				// 4) check whether best_node is sampled or real request
				/*int isRealRequest = 0;
				if (best_node < n)
				{
					if (RequestType[s][node] == 1)//0-sampled,1-real,2-unmaterialized
						isRealRequest = 1;
				}
				else
				{
					if (RequestType[s][node - n] == 1)//0-sampled,1-real,2-unmaterialized
						isRealRequest = 1;
				}*/

				// 5) update the best_node in the Master Solution (only if it is real request)
				if (1)//isRealRequest)
					if (best_score > 0)
					{
						int index = Master_Solution[0].Vehicle[k].size;
						Master_Solution[0].Vehicle[k].path[index] = best_node + 1;
						Master_Solution[0].Vehicle[k].size++;
					}
					else
					{
						printf("No feasible next service location!\n");
						PAUSER(3000);
					}

				// print utility
				printf("best_node: %d, best_score: %d\n", best_node, best_score);
				/*for (int s = 0; s < ExpectedScenarios; s++)
					Best_Solution[s].printinput();*/
				printf("***\n");
			}
		}


        // 6) Update Do_Not_Modify_Scenario
        if (best_score > 0)
            Update_Do_Not_Modify_Scenario(best_node, k);

		// 7) update Served Requests
		Update_ServedRequests(); // Make sure to update only real requests

		printf("CHECKING SERVED REQUESTS:\n"); ServedRequests[0].simple_print();
		for (int k = 0; k < m; k++)
			printf("%d | ", ServedRequests_BoundaryMark[0].list[k]);
		printf("\n");


		// delete memory
		delete[] Score;
	}

	void Find_Next_DynamicRequest_Details()
	{
		int found = 0;
		UnAvailableRequests[0].sort();
		for (int i = 0; i < UnAvailableRequests[0].size; i++)
		{
			int req = UnAvailableRequests[0].list[i];
			if (RejectedRequests[0].isFound(req))
				continue;
			else
			{
				Next_Dynamic_Request = req;
				Next_Injection_Time = (int)(Master_Problem[0].Request[req].pickup.EarliestTime) - 90;
				found = 1;
				break;
			}
		}
		if (!found)
		{
			printf("THERE IS NO MORE DYNAMIC REQUEST !!!!!!!!!!!\n");
			//Master_Solution[0].printinput_with_details(Master_Problem);
			//Destruction();
			//printf("Exited the program!\n");
			//exit(0);
			All_Dynamic_Requests_Inserted = 1;
		}

		printf("UPDATED: Next Dynamic Req: %d, Next_InjectionTime: %d\n", Next_Dynamic_Request, Next_Injection_Time);
	}

	int Feasible_Insert(int Dynamic_req, solution *Temp_Solution)
	{

	    Clear_Do_Not_Modify_Scenario();
		//solution *Temp_Solution = new solution[ExpectedScenarios];
		int feasible_s, feasible_start, feasible_gap;
		int Accept = 0;

		for (int s = 0; s < ExpectedScenarios; s++)
		{
            //copy scenario from its best scenario
            //Temp_Solution[s].Copy(Best_Solution[s]);
            Temp_Solution[s].Problem.Copy(Master_Problem[0]);
            Temp_Solution[s].OfflineEvaluateCost();

            // 1) remove all sampled requests
            printf("enter 1\n");
            for (int i = 0; i < SampledRequests[s].size; i++)
            {
                int req = SampledRequests[s].list[i];
                Temp_Solution[s].RemoveRequest(req, Temp_Solution[s].Request[req].currvehicle);
                Temp_Solution[s].OfflineEvaluateCost();
            }

			if (Temp_Solution[s].Cost.getFeasibility()) //check if the sampled request equivalent has any feasible solution
			{
				printf("enter 2: scenario %d is feasible upon removing sampled request\n", s);

				// 2) insert dynamic request in all possible positions

				int found = 0;
				int indices = 0;
				int tempstart;
				for (int k = 0; k < TotalVehicles; k++)
				{
					int size = Temp_Solution[s].Vehicle[k].size;
					//printf("veh %d size: %d\n", k, size);
					for (int limiter = size; limiter >= 0; limiter--)
					{
						for (tempstart = 0; tempstart <= limiter; tempstart++, indices++)
						{
							int start = tempstart;
							int gap = size - limiter;

							//printf("	start= %d, gap= %d\n", start, gap);

							// ensure insertion after served requests
							if (start > (ServedRequests_BoundaryMark[0].list[k] - 1) )
							{
								// reinsertion
								Temp_Solution[s].ReInsertRequest(Dynamic_req, k, start, start+gap+1);
								Temp_Solution[s].OfflineEvaluateCost();

								// feasibility check
								if (Temp_Solution[s].Cost.getFeasibility())
								{
									Accept = 1;
									found = 1;
									printf("For SCENARIO: %d - Feasible solution can be found at inserting various pos!\n", s);
									Temp_Solution[s].printinput();
									Temp_Solution[s].Cost.print();
									Solution[s].printinput();
									Solution[s].Cost.print();
									printf("-------------\n");
									PAUSER(3000);
								}
								else
								{
									Temp_Solution[s].RemoveRequest(Dynamic_req, k);
									Temp_Solution[s].OfflineEvaluateCost();
								}
							}
							if (found) break;
						}
						if (found) break;
					}
					if (found) break;
				}
				if (found) Do_Not_Modify_Scenario[s] = 1;

				// 3) accept or reject

			}
		}


        for (int s = 0; s < ExpectedScenarios; s++)
	    {
            printf("Do_Not_Modify_Scenario[%d] = %d\n", s, Do_Not_Modify_Scenario[s]);
            Temp_Solution[s].printinput();
            printf("Cost: %0.02f\n", Coefficient[s].cost_function(Temp_Solution[s].Cost));
	    }
		//Clear_Do_Not_Modify_Scenario();
		//delete[] Temp_Solution;

		return Accept;
	}

	void Inject_DynamicRequests()
	{

		// 1) Accept/Reject the new dynamic request
		int Dynamic_req = Next_Dynamic_Request;
		int Accept = 0;
		for (int s = 0; s < ExpectedScenarios; s++)
			printf("best_cost[%d]=%0.0f\n", s, best_cost[s]);
		PAUSER();

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			if (best_cost[s] > 0) //check if the sampled request equivalent has any feasible solution
			{
				//copy scenario from its best scenario
				Solution[s].Copy(Best_Solution[s]);
				Solution[s].OfflineEvaluateCost();

				for (int i = 0; i < AvailableRequests[s].size; i++)
					if (AvailableRequests[s].list[i] == Dynamic_req)
					{
						printf("FOR Scenario %d:\n", s);
						printf("BEFORE CORRECTION:\n");
						//Solution[s].Problem.print();
						//Solution[s].print_with_battery_level();
						Solution[s].Cost.print();
						//PAUSER();

						// check feasibility with the original constraints
						Solution[s].Problem.Copy(Master_Problem[0]);
						Solution[s].OfflineEvaluateCost();

						printf("AFTER CORRECTION:\n");
						//Solution[s].Problem.print();
						//Solution[s].print_with_battery_level();
						Solution[s].Cost.print();
						//PAUSER();

						if (Solution[s].Cost.getFeasibility())
							Accept = 1;
						else
						{
							printf("The Dynamic req is there, but not feasible with the original TW constraint!\n");
							PAUSER();
						}
						if (Accept)
							break;
					}
			}

			if (Accept)
				break;
		}

		printf("TRYING TO INJECT DYNAMIC REQ %d with Pickup TW [%0.0f,%0.0f]\n", Dynamic_req, Master_Problem[0].Request[Dynamic_req].pickup.EarliestTime, Master_Problem[0].Request[Dynamic_req].pickup.LatestTime);

		PAUSER();

		printf("Accept = %d; Now, trying FEASIBLE_INSERT()!\n", Accept);
		PAUSER();
        // NEW MODIFICATION:
        // 1) Take a copy of Scenarios in a Temporary_Solution
		solution *Temporary_Solution = new solution[ExpectedScenarios];
        for (int s = 0; s < ExpectedScenarios; s++)
        {
            if (best_cost[s] > 0)
                Temporary_Solution[s].Copy(Best_Solution[s]);
            else
                Temporary_Solution[s].Copy(Solution[s]);
            Temporary_Solution[s].OfflineEvaluateCost();
		}
		printf("Enter Feasible_Insert()\n");
		Accept = Feasible_Insert(Dynamic_req, Solution);
        printf("Exited Feasible_Insert()\n");
        if (Accept == 0)
			printf("Dynamic request %d is rejected despite trying the FEASIBLE_INSERT():!\n", Dynamic_req);
	    PAUSER(3000);

        // 2) Accept/Reject
        if (Accept)
        {
            Total_Accepted_Dynamic_Requests++;
            if (ServedRequests[0].isFound(Dynamic_req))
			{
				printf("ERR: Dynamic Request is already served! Line 733\n");
				Destruction(); exit(0);
			}

            // 1) Adjust Available and Sampled lists
            printf("Available: "); AvailableRequests[0].simple_print();
			printf("UnAvailable: ");UnAvailableRequests[0].simple_print();
			printf("Sampled: ");SampledRequests[0].simple_print();
			printf("Rejected: ");RejectedRequests[0].simple_print();
			printf("\n-----\n");
            for (int s = 0; s < ExpectedScenarios; s++)
			{
			    // 1.1) add/remove sampled request in Available/Unavailable lists
				SampledRequests[s].simple_print();
				for (int index = 0; index < SampledRequests[s].size; index++)
				{
					int i = SampledRequests[s].list[index];
					if (AvailableRequests[s].isFound(i))
					{
					    AvailableRequests[s].popout(i);
					    UnAvailableRequests[s].push(i);
					    RequestType[s][i] = 2; //0-sampled,1-real,2-unavailable
					}

				}
				SampledRequests[s].clear();

		        // 1.1) add/remove dynamic request in Available/Unavailable lists
		        AvailableRequests[s].push(Dynamic_req);
		        UnAvailableRequests[s].popout(Dynamic_req);
		        RequestType[s][Dynamic_req] = 1; //0-sampled,1-real,2-unavailable
			}
			printf("Available: "); AvailableRequests[0].simple_print();
			printf("UnAvailable: ");UnAvailableRequests[0].simple_print();
			printf("Sampled: ");SampledRequests[0].simple_print();
			printf("Rejected: ");RejectedRequests[0].simple_print();
            PAUSER(9000);

            // 2) Update the Scenario solutions from the Master Solution
            for (int s = 0; s < ExpectedScenarios; s++)
            {
                if (Do_Not_Modify_Scenario[s]) //copy the scenario
                {
                    Best_Solution[s].Copy(Solution[s]);
                    Best_Solution[s].OfflineEvaluateCost();

                    best_cost[s] = Coefficient[s].cost_function(Best_Solution[s].Cost);
                    Update_Existence(Solution, s);

                    Solution[s].printinput();
                }
                else //rebuild the scenario
                    Modify_Scenarios_From_Master_Solution(s);
            }
            Clear_Do_Not_Modify_Scenario();
		}
		else //REJECTION
		{
		    Total_Rejected_Dynamic_Requests++;
			RejectedRequests[0].push(Dynamic_req);
		}
		printf("Attempting to delete Temporary_Solution\n");
		delete[] Temporary_Solution;
        printf("Deleted Successfully!\n");


		// 2) If accepted: update scenarios and lists
		/*if (Accept)
		{
			Total_Accepted_Dynamic_Requests++;
			if (ServedRequests[0].isFound(Dynamic_req))
			{
				printf("ERR: Dynamic Request is already served! Line 733\n");
				exit(0);
			}

			// Part-I: Re-update AvailableRequests, UnAvailableRequests (exclude ServedRequests)
			for (int s = 0; s < ExpectedScenarios; s++)
			{
				AvailableRequests[s].clear();
				UnAvailableRequests[s].clear();
				SampledRequests[s].clear();

				// update requestTypes
				for (int i = 0; i < n; i++)
					if (ServedRequests[0].isFound(i) == 0)
					{
						if (!RejectedRequests[0].isFound(i) && i <= Dynamic_req) // real requests (exclude served reqs)
						{
							AvailableRequests[s].push(i);
							RequestType[s][i] = 1; //0-sampled,1-real,2-unavailable
						}
						else // unavailable requests
						{
							UnAvailableRequests[s].push(i);
							RequestType[s][i] = 2; //0-sampled,1-real,2-unavailable
						}
					}
			}
			AvailableRequests[0].print();
			UnAvailableRequests[0].print();
			RejectedRequests[0].print();

			// Part-II: Rebuild scenario solutions from Master Solution that includes the dynamic req
			Rebuild_Scenarios_From_Master_Solution();

		}
		else
		{
			Total_Rejected_Dynamic_Requests++;
            RejectedRequests[0].push(Dynamic_req);

			printf("Dynamic request is rejected!\n");
			//PAUSER(1);
		}*/

		// print utility
		if (Accept)
		{
			printf("DYNAMIC REQUEST %d is ACCEPTED! :\n)", Dynamic_req);
			PAUSER(1);
			PAUSER(Dynamic_req);

			printf("Accept: Removing existing samples from Scenarios!\n");
            Remove_All_Samples_From_Scenarios();
		}
		else
		{
			printf("DYNAMIC REQUEST %d is REJECTED!\n", Dynamic_req);
			PAUSER(1);
			PAUSER(Dynamic_req);

			printf("Reject: Removing existing samples from Scenarios!\n");
            Remove_All_Samples_From_Scenarios();

		}

		PAUSER(4000);

	}

	void Generate_Samples()
	{
		//int ADD_REQUEST = 35;
		for (int s = 0 ; s < ExpectedScenarios; s++)
		{
			SampledRequests[s].clear();

			////////////////////////////////////////
			//Code: to generate sampled requests////
			////////////////////////////////////////
			std::random_device rd; // uniformly-distributed integer random number generator
			std::mt19937 rng (rd ()); // mt19937: Pseudo-random number generation
			// Samples generate between [sumArrivalTimes, firstPartof_averageArrivalTimes]
			//double PeriodStart = 300;
			double sumArrivalTimes = max(Start_Window, CurrentTime);//max(229, CurrentTime);
			double averageArrival = End_Window - sumArrivalTimes;//500 - sumArrivalTimes;
			int totalRequests = UnAvailableRequests[s].size;//10
			double lamda = (double)totalRequests/averageArrival;//((24 * 1) / 100);
			std::exponential_distribution<double> exp (lamda);
			double newArrivalTime;
			////////////////////////////////////////

			for (int j = 0; j < 1; j++)
			{
				if (UnAvailableRequests[s].size > 0)
				{

					//ensure generated samples are newly dynamic ones

					/*int index = rand() % UnAvailableRequests[s].size;
					int i = UnAvailableRequests[s].list[index];*/

					//REJECTED REQUESTS SHOULD NOT BE SAMPLED!
					int found = 0; int i;
					for (int idx = 0; idx < 500; idx++)
					{
						int index = rand() % UnAvailableRequests[s].size;
						i = UnAvailableRequests[s].list[index];

						if (RejectedRequests[0].isFound(i) == 0)
							break;

						if (idx > 400)
						{
							printf("Exiting program: error? making samples in Generate_Samples();\n");
							Destruction(); exit(0);
						}
					}
					SampledRequests[s].push(i);

					//fill the constraints
					newArrivalTime = exp.operator() (rng); // generates the next random number in the distribution
					sumArrivalTimes = sumArrivalTimes + newArrivalTime;
					Solution[s].Problem.Request[i].pickup.EarliestTime = sumArrivalTimes;
					Solution[s].Problem.Request[i].pickup.LatestTime = sumArrivalTimes + 35;
					Solution[s].Problem.Request[i].dropoff.EarliestTime = 0;
					Solution[s].Problem.Request[i].dropoff.LatestTime = 1440;

					printf("Scenario %d: ", s);
					printf("%d:[%0.0f,%0.0f]\n---\n", i, Solution[s].Problem.Request[i].pickup.EarliestTime, Solution[s].Problem.Request[i].pickup.LatestTime);

					//print utility
					/*printf("BEFORE:\n");
					printf("SCENARIO %d\n size %d\n", s, AvailableRequests[s].size);
					for (int i = 0; i < AvailableRequests[s].size; i++)
						printf("AvailableRequests[%d].list[%d]=%d\n", s, i, AvailableRequests[s].list[i]);
					printf("Unavailable size %d\n", s, UnAvailableRequests[s].size);
					for (int i = 0; i < UnAvailableRequests[s].size; i++)
						printf("UnAvailableRequests[%d].list[%d]=%d\n", s, i, UnAvailableRequests[s].list[i]);
					printf("BEFORE\n");
					Solution[s].printinput();
					*/

					best_cost[s] = 0;
					int k = 0;
					Solution[s].ReInsertRequest(i, k, Solution[s].Vehicle[k].size, Solution[s].Vehicle[k].size+1);
					Solution[s].OfflineEvaluateCost();

					//remove it from unavailable request
					RequestType[s][i] = 0; //sampled-requests
					AvailableRequests[s].push(i);
					UnAvailableRequests[s].popout(i);
					Update_Existence(Solution, s);

					//print utility
					/*printf("AFTER\n");
					Solution[s].printinput();
					printf("AFTER:\n");
					printf("SCENARIO %d\n size %d\n", s, AvailableRequests[s].size);
					for (int i = 0; i < AvailableRequests[s].size; i++)
						printf("AvailableRequests[%d].list[%d]=%d\n", s, i, AvailableRequests[s].list[i]);
					printf("Unavailable size %d\n", s, UnAvailableRequests[s].size);
					for (int i = 0; i < UnAvailableRequests[s].size; i++)
						printf("UnAvailableRequests[%d].list[%d]=%d\n", s, i, UnAvailableRequests[s].list[i]);
					 */
				}
			}


			////////////////////////////////////////

		}
	}

	void Update_StaticRequests()
	{
		// Declare here
		int Total_Dynamical_Request = round(DegreeofDynamism * n);
		Total_Static_Requests = n - Total_Dynamical_Request;

		for (int s = 0 ; s < ExpectedScenarios; s++)
		{
			for (int i = 0; i < n; i++)
				if (i < Total_Static_Requests)// || i > 20)//ADD_REQUEST)
				{
					AvailableRequests[s].push(i);
					RequestType[s][i] = 1; //real-requests
				}
				else
				{
					UnAvailableRequests[s].push(i);
					RequestType[s][i] = 2; //hasn't materialized
					Solution[s].RemoveRequest(i, Solution[s].Request[i].currvehicle);
				}

			/*printf("SCENARIO %d\n size %d\n", s, AvailableRequests[s].size);
			for (int i = 0; i < AvailableRequests[s].size; i++)
				printf("AvailableRequests[%d].list[%d]=%d\n", s, i, AvailableRequests[s].list[i]);
			printf("Unavailable size %d\n", s, UnAvailableRequests[s].size);
			for (int i = 0; i < UnAvailableRequests[s].size; i++)
				printf("UnAvailableRequests[%d].list[%d]=%d\n", s, i, UnAvailableRequests[s].list[i]);
			 */

			////////////////////////////////////////

			SampledRequests[s].clear();

		}
		RejectedRequests[0].clear();

		AvailableRequests[0].print();
		UnAvailableRequests[0].print();
		RejectedRequests[0].print();

	}

	void Update_ServedRequests(int firstTime = 0)
	{
		if (1)//firstTime)
		{
			// 1) Update first requests of each vehicle as "served".
			for (int k = 0; k < TotalVehicles; k++)
			{
				int found = 0; //found_next_real_request
				for (int j = 0; j < Master_Solution[0].Vehicle[k].size; j++)
				{
					int node = Master_Solution[0].Vehicle[k].path[j] - 1;
					if (node < n)
					{
						if (ServedRequests[0].isFound(node) == 0)
						{

							// 1.1) update boundary mark
							//if (RequestType[s][node] == 1) //0-sampled, 1-real, 2-not materialized yet
							{
								//ServedRequests_BoundaryMark[0].list[k] += 1;
							}

							// 1.2) update served/available request list
							ServedRequests[0].push(node);
							for (int s = 0; s < ExpectedScenarios; s++)
								AvailableRequests[s].popout(node);
						}
					}
					else
					{
						// 1.3) update boundary mark

							//if (RequestType[s][node - n] == 1) //0-sampled, 1-real, 2-unmaterialized
							{
								//ServedRequests_BoundaryMark[0].list[k] += 1;
							}

					}

					ServedRequests_BoundaryMark[0].list[k] = Master_Solution[0].Vehicle[k].size;
				}
			}

			//print utility
			/*for (int s = 0; s < ExpectedScenarios; s++)
				Best_Solution[s].printinput();

			for (int j = 0; j < ServedRequests[0].size; j++)
				printf("%d - ", ServedRequests[0].list[j]);
			printf("\n");

			for (int k = 0; k < m; k++)
				printf("%d | ", ServedRequests_BoundaryMark[0].list[k]);
			printf("\n");*/

		}


	}

	void Check_for_Next_Idle_Vehicle(float CurrentTime)
	{
		for (int k = 0; k < TotalVehicles; k++)
			if (Master_Solution[0].Vehicle[k].size > 0)
			{
				int index = Master_Solution[0].Vehicle[k].size - 1;
				int node = Master_Solution[0].Vehicle[k].path[index] - 1;

				// making sure it is not served already
				if (node < n)
				{
					if (ServedRequests[0].isFound(node))
					{
						float B = Master_Solution[0].Request[node].B;
						float C = Master_Solution[0].Request[node].C;
						if (CurrentTime > B)// && CurrentTime < C) // whether k is idle?
							Idle_Vehicle_List[k] = 1;
					}
				}
				else
				{
					if (ServedRequests[0].isFound(node-n))
					{
						float B = Master_Solution[0].Request[node].B;
						float C = Master_Solution[0].Request[node].C;
						if (CurrentTime > B)// && CurrentTime < C) // whether k is idle?
							Idle_Vehicle_List[k] = 1;
					}
				}
			}
	}

	void Reschedule_Master_Solution()
	{
		for (int k = 0; k < m; k++)
		{
			if (Master_Solution[0].Vehicle[k].size > 1)
				for (int j = 1; j < Master_Solution[0].Vehicle[k].size; j++)
				{
					int prev_node = Master_Solution[0].Vehicle[k].path[j-1] - 1;
					int node = Master_Solution[0].Vehicle[k].path[j] - 1;
					float ArrivalTime = Master_Solution[0].Request[prev_node].C + Master_Problem[0].TravelTime[prev_node][node];

					if (node < n)
					{
						Master_Solution[0].Request[node].B = max(ArrivalTime, Master_Problem[0].Request[node].pickup.EarliestTime);
						Master_Solution[0].Request[node].C = Master_Solution[0].Request[node].B + Master_Problem[0].Request[node].pickup.ServiceTime;

						// Adjust Problems in effect with the scheduled service locations
						Master_Problem[0].Request[node].pickup.EarliestTime = Master_Solution[0].Request[node].B ;
						Master_Problem[0].Request[node].pickup.LatestTime = Master_Solution[0].Request[node].C;

					}
					else
					{
						//Master_Solution[0].Request[node].B = min(ArrivalTime, Master_Problem[0].Request[node-n].dropoff.EarliestTime);
						Master_Solution[0].Request[node].B = ArrivalTime;
						Master_Solution[0].Request[node].C = Master_Solution[0].Request[node].B + Master_Problem[0].Request[node-n].dropoff.ServiceTime;

						// Adjust Problems in effect with the scheduled service locations
						Master_Problem[0].Request[node - n].dropoff.EarliestTime = Master_Solution[0].Request[node].B ;
						Master_Problem[0].Request[node - n].dropoff.LatestTime = Master_Solution[0].Request[node].C;

					}
				}
		}
	}

    void Modify_Scenarios_From_Master_Solution(int s)
    {

        if (Do_Not_Modify_Scenario[s])
        {
            // 3) copy the problems from Master_Problem
			//Solution[s].Problem.Copy_Except_SampledRequests(Master_Problem[0], s, SampledRequests); //HERE-is-the-problem
			for (int k = 0; k < m; k++)
				for (int j = 0; j < Master_Solution[0].Vehicle[k].size; j++)
				{
					int i = Master_Solution[0].Vehicle[k].path[j] - 1;
					if ( i < n && RequestType[s][i] != 0 )
						Solution[s].Problem.Request[i].Copy(Master_Problem[0].Request[i]);
				}
        }
        else
		{
			// 1) Copy plan for served requests from Master Solution
			Fill_ServedRequests_from_Master_Solution(s);

			// 2) put respective dropoff points for served requests if not there
			Fill_Dropoffs_into_Solution(Solution, s);

			// 3) copy the problems from Master_Problem
			//Solution[s].Problem.Copy_Except_SampledRequests(Master_Problem[0], s, SampledRequests); //HERE-is-the-probklem
			for (int k = 0; k < m; k++)
				for (int j = 0; j < Master_Solution[0].Vehicle[k].size; j++)
				{
					int i = Master_Solution[0].Vehicle[k].path[j] - 1;
					if ( i < n && RequestType[s][i] != 0 )
						Solution[s].Problem.Request[i].Copy(Master_Problem[0].Request[i]);

				}

			// 4) put all available (real+sampled) requests into scenario solutions
			Fill_Available_Requests_into_Solution(Solution, s);
			Update_Existence(Solution, s);
			best_cost[s] = 0;

			float current_cost = Coefficient[s].cost_function(Solution[s].Cost);
			if (Solution[s].Cost.getFeasibility())
			{
				best_cost[s] = current_cost;
				Best_Solution[s].Copy(Solution[s]);
				Best_Solution[s].OfflineEvaluateCost();
			}
		}


    }

	void Rebuild_Scenarios_From_Master_Solution()
	{
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			// 1) Copy plan for served requests from Master Solution
			Fill_ServedRequests_from_Master_Solution(s);

			// 2) put respective dropoff points for served requests if not there
			Fill_Dropoffs_into_Solution(Solution, s);

			// 3) copy the problems from Master_Problem
			//Solution[s].Problem.Copy_Except_SampledRequests(Master_Problem[0], s, SampledRequests); //HERE-is-the-probklem
			for (int k = 0; k < m; k++)
				for (int j = 0; j < Master_Solution[0].Vehicle[k].size; j++)
				{
					int i = Master_Solution[0].Vehicle[k].path[j] - 1;
					if ( i < n && RequestType[s][i] != 0 )
						Solution[s].Problem.Request[i].Copy(Master_Problem[0].Request[i]);

				}

			// 4) put all available requests into scenario solutions
			Fill_Available_Requests_into_Solution(Solution, s);
			Update_Existence(Solution, s);
			best_cost[s] = 0;

			float current_cost = Coefficient[s].cost_function(Solution[s].Cost);
			if (Solution[s].Cost.isFeasible)
			{
				best_cost[s] = current_cost;
				Best_Solution[s].Copy(Solution[s]);
				Best_Solution[s].OfflineEvaluateCost();
			}
		}
	}

	void Remove_All_Samples_From_Scenarios()
	{
            for (int s = 0; s < ExpectedScenarios; s++)
			{
				SampledRequests[s].simple_print();

				printf("Control enters the scenario %d\n", s);

				for (int index = 0; index < SampledRequests[s].size; index++) {
					int i = SampledRequests[s].list[index];

				    printf("-- Before removing request %d from vehicle %d\n", i, Solution[s].Request[i].currvehicle);
                    Solution[s].printinput();
                    printf("Available: "); AvailableRequests[0].simple_print();
                    printf("UnAvailable: ");UnAvailableRequests[0].simple_print();
                    printf("Sampled: ");SampledRequests[0].simple_print();
                    printf("Rejected: ");RejectedRequests[0].simple_print(); cout<<endl;


					if (AvailableRequests[s].isFound(i)) {
						AvailableRequests[s].popout(i); printf("- pop success\n");
						UnAvailableRequests[s].push(i); printf("- push success\n");
						Solution[s].RemoveRequest(i, Solution[s].Request[i].currvehicle); printf("- remove success\n");
						Solution[s].OfflineEvaluateCost(); printf("- eval success\n");
						Update_Existence(Solution, s); printf("- updateExistence success\n");
					} else {
					    printf("- particular request is in the unavailable!\n");
					}

					printf("-- After removing request %d from vehicle %d\n", i, Solution[s].Request[i].currvehicle);
					Solution[s].printinput();
                    printf("Available: "); AvailableRequests[0].simple_print();
                    printf("UnAvailable: ");UnAvailableRequests[0].simple_print();
                    printf("Sampled: ");SampledRequests[0].simple_print();
                    printf("Rejected: ");RejectedRequests[0].simple_print();
                    printf("\n---------------\n\n");
				}
				SampledRequests[s].clear(); printf("- sampleClear success\n");

				printf("--- Control leaves the scenario %d\n\n", s);
			}
	}

	void Fill_Dropoffs_into_Solution(solution *Solution, int s)
	{
		// 1) fill dropoffs
		for (int i = 0; i < ServedRequests[0].size; i++)
		{
			int pickup = ServedRequests[0].list[i];

			// 3.1.1) find if there is a drop off node already!
			int locatedPickUp = 0;
			int locatedDropOff = 0;
			int k;
			for (k = 0; k < TotalVehicles; k++)
			{
				for (int j = 0; j < Solution[s].Vehicle[k].size; j++)
				{
					int P_node = Solution[s].Vehicle[k].path[j] - 1;
					if (pickup == P_node) // located-pickup node
					{
						locatedPickUp = 1;
						//now, locate the dropoff-node
						for (int l = j; l < Solution[s].Vehicle[k].size; l++)
							if (pickup + n == Solution[s].Vehicle[k].path[l] - 1)
								locatedDropOff = 1;
					}
					if (locatedPickUp)
						break;
				}
				if (locatedPickUp)
					break;
			}

			// 2.1.2) add the drop off node since it's not there already
			if (locatedDropOff == 0)
			{
				int index = Solution[s].Vehicle[k].size;
				Solution[s].Vehicle[k].path[index] = pickup + n + 1;
				Solution[s].Vehicle[k].size++;
			}

		}
	}

	void Fill_ServedRequests_from_Master_Solution(int s)
	{
		for (int k = 0; k < TotalVehicles; k++)
		{
			Solution[s].Vehicle[k].size = 0;
			for (int j = 0; j < Master_Solution[0].Vehicle[k].size; j++)
			{
				int index = Solution[s].Vehicle[k].size;
				Solution[s].Vehicle[k].path[index] = Master_Solution[0].Vehicle[k].path[j];
				Solution[s].Vehicle[k].size++;
			}
		}
	}

	void Fill_Available_Requests_into_Solution(solution *Solution, int s)
	{
		for (int i = 0; i < AvailableRequests[s].size; i++)
		{
			int k = rand() % m;

			int index = Solution[s].Vehicle[k].size;
			Solution[s].Vehicle[k].path[index] = AvailableRequests[s].list[i] + 1;
			Solution[s].Vehicle[k].size++;

			index = Solution[s].Vehicle[k].size;
			Solution[s].Vehicle[k].path[index] = AvailableRequests[s].list[i] + 1 + n;
			Solution[s].Vehicle[k].size++;
		}
		Solution[s].OfflineEvaluateCost();
		Solution[s].IdentifyCurrvehicle();
	}

	void Interrupt_Check(int iteration, float *currentbestcost, int Master_k)
	{
		if (AvailableRequests[0].size < 1)
		{
			Fill_Dropoffs_into_Solution(Master_Solution, 0);
			Reschedule_Master_Solution();
			Master_Solution[0].printinput();
			Master_Solution[0].printinput_with_details(Master_Problem);
			printf("Exiting the Program since no more dynamic requests!\n");
			exit(0);
		}

		if (iteration == 40) // Static over
		{

			printf("-----Maintaining Master Solution-----\n");
			//update for the first time
			Maintain_Master_Solution_for_FirstTime(1,1);

			PAUSER();

			ServedRequests[0].simple_print();
			AvailableRequests[0].simple_print();
			UnAvailableRequests[0].simple_print();
			RejectedRequests[0].simple_print();
			SampledRequests[0].simple_print();

			PAUSER();

			printf("-----Generating Samples at Itr %d-----\n", iteration);
			Generate_Samples();
			for (int s = 0; s < ExpectedScenarios; s++)
				currentbestcost[s] = best_cost[s];

			for (int s = 0; s < ExpectedScenarios; s++)
			{
				int vehSize = 0;
				for (int k = 0; k < m; k++)
					vehSize += Solution[s].Vehicle[k].size;
				printf("Solution size[%d]=%d\n", s, vehSize);
				//Solution[s].printinput();
			}

			ServedRequests[0].simple_print();
			AvailableRequests[0].simple_print();
			UnAvailableRequests[0].simple_print();
			RejectedRequests[0].simple_print();
			SampledRequests[0].simple_print();

			PAUSER();
		}
		else if (iteration > 40 && iteration % 6 == 0 && All_Dynamic_Requests_Inserted == 0) // dynamic in-progress
		{
			if (CurrentTime >= Next_Injection_Time)
			{

				//////////////////////////////////////////
				/////// 1) Inject Dynamic Requests ///////
				//////////////////////////////////////////

				printf("-------------------------------------\n");
				printf("BEFORE injecting dynamic request:\n");
				ServedRequests[0].simple_print();
				AvailableRequests[0].simple_print();
				UnAvailableRequests[0].simple_print();
				RejectedRequests[0].simple_print();
				SampledRequests[0].simple_print();

				printf("-------------------------------------\n");
				printf("-----Injecting a Dynamic Request at Itr %d-----\n", iteration);
				Inject_DynamicRequests();

				printf("-------------------------------------\n");
				printf("AFTER injecting dynamic request:\n");
				ServedRequests[0].simple_print();
				AvailableRequests[0].simple_print();
				UnAvailableRequests[0].simple_print();
				RejectedRequests[0].simple_print();
				SampledRequests[0].simple_print();

				for (int s = 0; s < ExpectedScenarios; s++)
				{
					int vehSize = 0;
					for (int k = 0; k < m; k++)
						vehSize += Solution[s].Vehicle[k].size;
					printf("Solution size[%d]=%d\n", s, vehSize);
					//Solution[s].printinput();
				}

				//PAUSER(1);

                //////////////////////////////////////////
                /////// 2) Update Next Injection Details ///////
                //////////////////////////////////////////

                Find_Next_DynamicRequest_Details();

				//////////////////////////////////////////
				/////// 3) Generate Future Samples ///////
				//////////////////////////////////////////

                if (All_Dynamic_Requests_Inserted == 0)
                {
                    printf("-----Generating Samples at Itr %d-----\n", iteration);
				    Generate_Samples();
                }
                else
                {
                    printf("-----Stopped generating samples since all dynamic reqs are tried!-----\n", iteration);
                }
				for (int s = 0; s < ExpectedScenarios; s++)
					currentbestcost[s] = best_cost[s];
				printf("-------------------------------------\n");

				for (int s = 0; s < ExpectedScenarios; s++)
				{
					int vehSize = 0;
					for (int k = 0; k < m; k++)
						vehSize += Solution[s].Vehicle[k].size;
					printf("Solution size[%d]=%d\n", s, vehSize);
					//Solution[s].printinput();
				}
				printf("Request info:\n Static: %d\n, Accept: %d\n, Reject: %d\n, Unavailable: %d\n, Served: %d\n, Available: %d\n, Sampled: %d\n\n", Total_Static_Requests, Total_Accepted_Dynamic_Requests, Total_Rejected_Dynamic_Requests, UnAvailableRequests[0].size, ServedRequests[0].size, AvailableRequests[0].size, SampledRequests[0].size);

				if ((ServedRequests[0].size + AvailableRequests[0].size - SampledRequests[0].size) == (Total_Static_Requests + Total_Accepted_Dynamic_Requests))
				    printf("Summation counts fine!!\n");
			    else
			        printf("SUM DOES NOT TALLY !!!!!!!!!\n");

                if ((ServedRequests[0].size + AvailableRequests[0].size + UnAvailableRequests[0].size) == n)
                    printf("Summation counts fine!!\n");
                else
                    printf("SUM DOES NOT TALLY !!!!!!!!!\n");
				//PAUSER(1);
			}
		}
		else if (iteration > 40 && iteration % 20 == 0) // Dynamic & Static over.
		{
		        if (Stop_Checking == 0)
		            if (All_Dynamic_Requests_Inserted)
                    {
                        printf("Since all dynamic requests are tried, let's remove all samples from Scenarios!\n");
                        Remove_All_Samples_From_Scenarios();
                        Stop_Checking = 1;
                    }

				printf("-------------------------------------\n");
				printf("-----Maintaining Master Solution at Itr %d-----\n", iteration);

				CurrentTime += 3;
				Check_for_Next_Idle_Vehicle(CurrentTime);

				myList NewList;
				for (int k = 0; k < m; k++)
				    NewList.list[k] = k;
				NewList.randomize();

                int Triggered = 0;
				for (int l = 0; l < m; l++)
				{
				    int k = NewList.list[l];
				    if (Idle_Vehicle_List[k]) //trigger if any idle vehicle
					{
						printf("CurrentTime: %d\n", CurrentTime);
						printf("Idle_Vehicle_List: ");
						for (int index = 0; index < m; index++)
							printf("%d ", Idle_Vehicle_List[index]);
						printf("\n");

						Master_k = k;
						Maintain_Master_Solution(Master_k);

						Master_Solution[0].printinput();
						Master_Solution[0].printinput_with_details(Master_Problem);
						printf("-------------------------------------\n");

						//AvailableRequests[0].print();
						printf("CurrentTime: %d\n", CurrentTime);

						Triggered = 1;
						//break; // THIS CAN BE PROBLEMATIC !!!!!!!!!!!!!
					}
                }

				if (Triggered == 0)
					printf("No Idle Vehicles at Current Time %d\n", CurrentTime);
				else
					printf("Triggered at Current Time %d\n", CurrentTime);

				//clear idle vehicle list
				for (int k = 0; k < m; k++)
					Idle_Vehicle_List[k] = 0;

				for (int s = 0; s < ExpectedScenarios; s++)
					printf("best_cost[%d]=%0.0f\n", s, best_cost[s]);
				/*for (int s = 0; s < ExpectedScenarios; s++)
				{
					int vehSize = 0;
					for (int k = 0; k < m; k++)
						vehSize += Solution[s].Vehicle[k].size;
					printf("Solution size[%d]=%d\n", s, vehSize);
					//Solution[s].printinput();
				}*/

				ServedRequests[0].simple_print();
				AvailableRequests[0].simple_print();
				UnAvailableRequests[0].simple_print();
				RejectedRequests[0].simple_print();
				SampledRequests[0].simple_print();

				printf("Request info:\n Static: %d\n, Accept: %d\n, Reject: %d\n, Unavailable: %d\n, Served: %d\n, Available: %d\n, Sampled: %d\n\n", Total_Static_Requests, Total_Accepted_Dynamic_Requests, Total_Rejected_Dynamic_Requests, UnAvailableRequests[0].size, ServedRequests[0].size, AvailableRequests[0].size, SampledRequests[0].size);
				if ((ServedRequests[0].size + AvailableRequests[0].size - SampledRequests[0].size) == (Total_Static_Requests + Total_Accepted_Dynamic_Requests))
                    printf("Summation counts fine!!\n");
			    else
			        printf("SUM DOES NOT TALLY !!!!!!!!!\n");

                if ((ServedRequests[0].size + AvailableRequests[0].size + UnAvailableRequests[0].size) == n)
                    printf("Summation counts fine!!\n");
                else
                    printf("SUM DOES NOT TALLY !!!!!!!!!\n");

				PAUSER(3000);
		}
	}

	void Optimize(int num_of_iteration, int TimeLimit)
	{
		// 1) Update StaticRequests
		Update_StaticRequests();
		Find_Next_DynamicRequest_Details();

		TotalTime = clock();

		// Optimization begins

		solution *Temp_Solution = new solution[ExpectedScenarios];
		float *currentbestcost = new float[ExpectedScenarios];
		for (int s = 0; s < ExpectedScenarios; s++)
			currentbestcost[s] = best_cost[s];

		int Cycles = 0;
		clock_t InjectionTime = clock();
		clock_t BatchingTime = clock();

		int Master_k = 0;
		CurrentTime = 90;

		int intensity = 1, reduction_choice = 0, iteration = 1;
		for (iteration = 1; iteration <= num_of_iteration; iteration++)
		{

			// OSCO Interrupt_Check
				Looper:
					Interrupt_Check(iteration, currentbestcost,Master_k);
					if (AvailableRequests[0].size == 0)
					{
						printf("No more available requests here!\n");
						iteration++;
						goto Looper;
					}

			// Exit if all requests are served
					if (ServedRequests[0].size == n)
					{
						printf("All requests are served! So, exiting program!\n");
						Master_Solution[0].printinput();
						exit(0);
					}


			// Backup a copy of the scenarios into Temp_Solution
					for (int s = 0; s < ExpectedScenarios; s++)
					{
						Temp_Solution[s].Copy(Solution[s]);
						Temp_Solution[s].OfflineEvaluateCost();
					}

			// utility: trigger exit, print itr time
					if ((float)((clock() - TotalTime) / CLOCKS_PER_SEC) > TimeLimit)
						goto displayResult;


					if (iteration % 500 == 0)
						printf("At Iter %d, TimeElapsed: %0.02f\n", iteration, (float)((clock() - TotalTime) / CLOCKS_PER_SEC));

if (CurrentTime > 115) {
    printf("\nReached 1");
}

			// 1) update parameters
					clock_t CoefficientUpdateTime = clock();
					for (int s = 0; s < ExpectedScenarios; s++)
					{
						Coefficient[s].update(Solution[s].Cost);
						if (iteration % 5 == 0)
							Coefficient[s].reset();
					}
					//printf("CoefficientUpdateTime: %0.04f ms\n-----\n", (float)(clock() - CoefficientUpdateTime) / CLOCKS_PER_SEC * 1000);

					reduction_choice = 0;//rand() % 3; // 0 - reduce normal cost_function, 1 - reduce cost_function_with_batteryViols, 2 - reduce cost_function_with_batteryViols_2
					intensity = 1 + rand() % 5;//1 + rand() % 10;//1 + rand() % 5;
					intensity = min(intensity, AvailableRequests[0].size);

					int min = 0;
					for (int s = 0; s < ExpectedScenarios; s++)
						if (AvailableRequests[s].size < min || min == 0)
							min = AvailableRequests[s].size;
					if (intensity == 0 || intensity > n || intensity > min)
						printf("ERRORRR!!!! intensity is greater than n or AvailableReqSize\n");

					int select = rand() % 3; //0-shaw, 1-random, 2-worst
					int pick = 0;//rand() % 3; //0-greedy, 1-regret, 2-sameOrder
					int grade = rand() % (m-1);


					if (CurrentTime > 115) {
    printf("Reached 2");
}

			// 2) explore the neighborhood

					select = 1 + rand() % 2;

					if (0)//iteration % 5 == 0)
						IntraRoute_Insertion(Solution, intensity, iteration, ServedRequests_BoundaryMark);
					else if (0)//pick == 0 || pick == 1)//iteration % 10 == 0)
						Remove_and_Reinsert_Thoroughly(select, pick, grade, Solution, intensity, iteration, reduction_choice, ServedRequests_BoundaryMark);
					else
					{
						for (int i=1; i<=2; i++)
							Remove_and_Reinsert_Loosely(i, Solution, intensity, iteration, reduction_choice, ServedRequests_BoundaryMark);
					}


		if (CurrentTime > 115) {
    printf("Reached 3");
}

			// 3) restart mechanism
					clock_t TriggerRestartTime = clock();
					for (int s = 0; s < ExpectedScenarios; s++)
						Trigger_Restart_Mechanism(s, currentbestcost, iteration);
					//printf("TriggerRestartTime: %0.04f ms\n-----\n", (float)(clock() - TriggerRestartTime) / CLOCKS_PER_SEC * 1000);

			// 4) SA: acceptance criterion
					clock_t SAacceptTIme = clock();
					for (int s = 0; s < ExpectedScenarios; s++)
					{
						if (Solution[s].Cost.getFeasibility())
						{
							if (Coefficient[s].cost_function(Solution[s].Cost) < currentbestcost[s])
								currentbestcost[s] = Coefficient[s].cost_function(Solution[s].Cost);
							// 1) New best solution

							//int intensity_dummy = 0;
							//Test_UpdateBestSolution(0, Solution, intensity_dummy, iteration, 1);

						}
						else if ((Coefficient[s].cost_function(Solution[s].Cost) < Coefficient[s].cost_function(Temp_Solution[s].Cost)))
						{
							// 2) New improving solution
						}
						else
							{
								float p = exp((-1)*(Coefficient[s].cost_function(Solution[s].Cost) - Coefficient[s].cost_function(Temp_Solution[s].Cost)) / temperature[s]);

								float random_num = rand() % 10;
								random_num = random_num / 10;

								if (random_num < p)
								{
									// 3) New SA-accepted solution
								}
								else
								{
									Solution[s].Copy(Temp_Solution[s]);
									Solution[s].OfflineEvaluateCost();
								}
							}

					}
					//printf("AcceptTime: %0.04f ms\n-----\n", (float)(clock() - SAacceptTIme) / CLOCKS_PER_SEC * 1000);


			// SA: cooling down & temperature-check

					for (int s = 0; s < ExpectedScenarios; s++)
					{
						temperature[s] *= coolingRate[s];
						if (temperature[s] < minTemperature[s])
							temperature[s] = maxTemperature[s];
					}


					if (CurrentTime > 115) {
    printf("Reached 4\n");
}
		}

	displayResult:
		{

			// deleting Temp_Solution created for SA-acceptance
			delete[] Temp_Solution;
			delete[] currentbestcost;


			/*if (best_cost[0] > 0)
			{
				Best_Solution[0].OfflineEvaluateCost();
				Best_Solution[0].print_with_battery_level();

				printf("Solution is %s\n", Best_Solution[0].Cost.isFeasible == true ? "Feasible" : " NOT Feasible");
				//if (Best_Solution[0].Cost.getFeasibility())
				//printf("-->> Cost = %f\n", Best_Solution[0].Cost.travel_cost);
			}
			else
				printf("CANNOT FIND FEASIBLE SOLUTION WITHIN GIVEN TIME!!\n-------\n");*/


/*
			for (int s = 0; s < ExpectedScenarios; s++)
				printf("best_cost[%d] = %0.02f (n=%d)\n", s, Coefficient[s].cost_function(Best_Solution[s].Cost), Coefficient[s].cost_function(Best_Solution[s].Cost) == 0 ? 0 : Best_Solution[s].total_req_served);
*/


			for (int s = 0; s < ExpectedScenarios; s++)
				Best_Solution[s].printinput();

			printf("Master Solution:");
			Master_Solution[0].printinput();

			for (int s = 0; s < ExpectedScenarios; s++)
				printf("best_cost[%d] = %0.02f (n=%d)\n", s, best_cost[s], Coefficient[s].cost_function(Best_Solution[s].Cost) == 0 ? 0 : Best_Solution[s].total_req_served);

			//print utility
			for (int s = 0; s < ExpectedScenarios; s++)
				Best_Solution[s].printinput();

			printf("Total Iterations: %d\n", iteration);
			printf("Execution time: %0.0f\n\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000);
		}
	}

	void Remove_and_Reinsert_Thoroughly(int select, int pick, int grade, solution *TempSolution, int &intensity, int iteration, int choice, myList *ServedRequests_BoundaryMark)
	{


		int isGPU = Original_isGPU;
		int isFull = Original_isFull;
		int isRev = Original_isRev;

		// 1.1) choose operator
		for (int s = 0; s < ExpectedScenarios; s++)
			Removal[s].operator_selection(select, intensity, &TempSolution[s], &AvailableRequests[s]);

		// 1.2) Remove many requests together
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			for (int p = 0; p < intensity; p++)
			{
				req_list[s][p] = Removal[s].D_list[p];

				if (existence[s][req_list[s][p]])
				{
					TempSolution[s].RemoveRequest(req_list[s][p], TempSolution[s].Request[req_list[s][p]].currvehicle);
					existence[s][req_list[s][p]] = 0;
				}
			}
		}

		// 2) Sequentially reinsert those requests

		///////////////////////////////////////////////////////////////////////////
		for (int total_req_on_bank = intensity; total_req_on_bank > 0; total_req_on_bank--)
		{
			//printf("CheckPt 2.1\n");


			// Activate current bank si
			for (int s = 0; s < ExpectedScenarios; s++)
			{
				int index = 0;
				for (int j = 0; j < intensity; j++)
					if (existence[s][req_list[s][j]] == 0)
						temp_req_list[s][index++] = req_list[s][j];
			}


			// 2.1) Evaluate insertions of each request

			for (int p = 0; p < total_req_on_bank; p++)
			{
				// 1) explore neighborhood
				solution *NewSolution = new solution[ExpectedScenarios];

				for (int s = 0; s < ExpectedScenarios; s++)
				{
					NewSolution[s].Copy(TempSolution[s]);
					NewSolution[s].OfflineEvaluateCost();
				}

				int *currvehicleList = new int[ExpectedScenarios];

				Explore_Neighborhood[0].Exploration_Decision_Maker(isGPU, isFull, isRev, p, NewSolution, Problem, temp_req_list, Insertion, Coefficient, fsm, choice, currvehicleList, ServedRequests_BoundaryMark);

				delete[] currvehicleList;
				delete[] NewSolution;
			}


			int noise_application = 0;// grade = rand() % (m-1);//0;
			for (int s = 0; s < ExpectedScenarios; s++)
			{
				// 2.2) prepare the request_bank
				Reinsertion[0].operator_selection(pick, Removal, Insertion[s], grade, noise_application, Noise, total_req_on_bank);


				// 2.3) re-insert only the best request from request bank

				int best_p = Reinsertion[0].sequence[0];
				TempSolution[s].ReInsertRequest(Insertion[s][best_p].request, Insertion[s][best_p].vehid, Insertion[s][best_p].start, Insertion[s][best_p].start + Insertion[s][best_p].gap + 1);
				TempSolution[s].OfflineEvaluateCost();

				// Update existence
				Update_Existence(TempSolution, s);

			}

		}
		///////////////////////////////////////////////////////////////////////////

		// Update the best solution
		for (int s = 0; s < ExpectedScenarios; s++)
			Test_UpdateBestSolution(s, TempSolution, intensity, iteration, 1, select, pick, grade);

		// Deallocation of memory



	}

	void Remove_and_Reinsert_Loosely(int removalOp, solution *TempSolution, int &intensity, int iteration, int choice, myList *ServedRequests_BoundaryMark)
	{
		int isGPU = Original_isGPU;
		int isFull = Original_isFull;
		int isRev = Original_isRev;

		//select ==> 0 shaw, 1 random, 2 worst
		int select = removalOp;

		clock_t Op_selection_time = clock();
		if (select == 2)
		{

					//clock_t New_op_selection_time = clock();
					//WorstCostRemoval[0].execute_in_host(TempSolution, AvailableRequests, Coefficient);
					//printf("New_op_selection_time: %0.02f\n\n", (float)(clock() - New_op_selection_time) / CLOCKS_PER_SEC * 1000);

					//clock_t New_op_selection_time = clock();
					WorstCostRemoval[0].execute(TempSolution, AvailableRequests, Coefficient);
					//printf("New_op_selection_time 2: %0.02f\n\n", (float)(clock() - New_op_selection_time) / CLOCKS_PER_SEC * 1000);

		}
		else
		{
			// 1.1) choose operator
			for (int s = 0; s < ExpectedScenarios; s++)
				Removal[s].operator_selection(select, intensity, &TempSolution[s], &AvailableRequests[s]);
		}
		//printf("Op_selection_time: %0.02f\n\n", (float)(clock() - Op_selection_time) / CLOCKS_PER_SEC * 1000);


    if (CurrentTime > 115) {
        printf("Reached 2.1\n");
}

		clock_t Removal_Time = clock();
		// 1.2) Remove many requests together
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			for (int p = 0; p < intensity; p++)
			{

				if (select == 2)
					req_list[s][p] = WorstCostRemoval[0].h_reqList[p + (n*s)];
				else
					req_list[s][p] = Removal[s].D_list[p];


				if (UnAvailableRequests[s].isFound(req_list[s][p]) == 1)
					printf("SOMETHING WRONG AT LINE 651??: disabled req is being removed!\n");

				//if (ServedRequests[0].isFound(req_list[s][p]) == 1)
					//printf("SOMETHING WRONG AT LINE 652??: served req is being removed!\n");


				if (existence[s][req_list[s][p]])
				{
					TempSolution[s].RemoveRequest(req_list[s][p], TempSolution[s].Request[req_list[s][p]].currvehicle);
					existence[s][req_list[s][p]] = 0;
				}
			}
		}
		//printf("Removal_Time: %0.02f\n\n", (float)(clock() - Removal_Time) / CLOCKS_PER_SEC * 1000);


    if (CurrentTime > 115) {
        printf("Reached 2.2\n");
}

		// 2) Sequentially reinsert those requests
		//////////////////////////////////////////////////////////////////////////
		int *currvehicleList = new int[ExpectedScenarios];
		clock_t ExplorationTime = clock();
		for (int p = 0; p < intensity; p++)
			Explore_Neighborhood[0].Exploration_Decision_Maker(isGPU, isFull, isRev, p, TempSolution, Problem, req_list, Insertion, Coefficient, fsm, choice, currvehicleList, ServedRequests_BoundaryMark);
		//printf("ExplorationTime: %0.02f\n\n", (float)(clock() - ExplorationTime) / CLOCKS_PER_SEC * 1000);
		delete[] currvehicleList;

    if (CurrentTime > 115) {
        printf("Reached 2.3\n");
}

		clock_t UpdateExistence_and_BestSolutionTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			TempSolution[s].OfflineEvaluateCost();
			Update_Existence(TempSolution, s);

			// Update the best solution
			Test_UpdateBestSolution(s, TempSolution, intensity, iteration, select+1);
		}
		//printf("UpdateExistence_and_BestSolutionTime: %0.02f\n\n", (float)(clock() - UpdateExistence_and_BestSolutionTime) / CLOCKS_PER_SEC * 1000);
		///////////////////////////////////////////////////////////////////////////


		//printf("Total_Per_Move_Time: %0.02f\n-----\n\n", (float)(clock() - Op_selection_time) / CLOCKS_PER_SEC * 1000);

    if (CurrentTime > 115) {
        printf("Reached 2.4\n");
}

	}

	void IntraRoute_Insertion(solution *TempSolution, int &intensity, int iteration, myList *ServedRequests_BoundaryMark)
	{
		int isGPU = Original_isGPU;
		int isFull = Original_isFull;
		int isRev = Original_isRev;
		int choice = 0;


		for (int p = 0; p < n; p++)
		{
			int *currvehicleList = new int[ExpectedScenarios];

			// 1) Remove one request from each scenario solution

			for (int s = 0; s < ExpectedScenarios; s++)
			{
				if (p < AvailableRequests[s].size)
					req_list[s][p] = AvailableRequests[s].list[p];
				else
					req_list[s][p] = AvailableRequests[s].list[0]; //by this way: we'll remove & reinsert atleast one req in all scenarios all the time.

				if (UnAvailableRequests[s].isFound(req_list[s][p]) == 1)
					printf("SOMETHING WRONG AT LINE 821??: disabled req is being removed!\n");

				if (existence[s][req_list[s][p]])
				{
					currvehicleList[s] = TempSolution[s].Request[req_list[s][p]].currvehicle;
					TempSolution[s].RemoveRequest(req_list[s][p], TempSolution[s].Request[req_list[s][p]].currvehicle);
					existence[s][req_list[s][p]] = 0;
				}

			}


			// 2) Reinsertion

			Explore_Neighborhood[0].Exploration_Decision_Maker(isGPU, isFull, isRev, p, TempSolution, Problem, req_list, Insertion, Coefficient, fsm, choice, currvehicleList, ServedRequests_BoundaryMark, 1);


			// 3) Update best solution (from within)

			for (int s = 0; s < ExpectedScenarios; s++)
			{
				TempSolution[s].OfflineEvaluateCost();
				Update_Existence(TempSolution, s);

				// Update the best solution
				Test_UpdateBestSolution(s, TempSolution, intensity, iteration, 4);
			}


			delete[] currvehicleList;
		}

//		for (int s = 0; s < ExpectedScenarios; s++)
	//		Test_UpdateBestSolution(s, TempSolution, intensity, iteration, 4);


/*
			for (int k = 0; k < m; k++)
				for (int j = 0; j < TempSolution[0].Vehicle[k].size; j++)
					if ( (TempSolution[0].Vehicle[k].path[j] - 1) < n )
					{
						int node = TempSolution[0].Vehicle[k].path[j] - 1;

						// Allocation of memory

						int intensity = 1;
						int *req, *found;

						if (intensity == 1)
						{
							req = new int;
							found = new int;
						}
						else
						{
							req = new int[intensity];
							found = new int[intensity];
						}

						for (int p = 0; p < intensity; p++)
						{
							req[p] = node;
							found[p] = 0;
						}

						// 1) Remove many requests together

						for (int p = 0; p < intensity; p++)
						{
							if (existence[0][req[p]])
							{
								TempSolution[0].RemoveRequest(req[p], TempSolution[0].Request[req[p]].currvehicle);
								existence[0][req[p]] = 0;
							}
						}

						// 2) Sequentially reinsert those requests

						///////////////////////////////////////////////////////////////////////////

						for (int p = 0; p < intensity; p++)
						{

							Explore_Neighborhood[0].Exploration_Decision_Maker(p, TempSolution, Problem, req, Coefficient, fsm, 0, 1);

						}
						TempSolution[0].OfflineEvaluateCost();
						//TempSolution[0].Offline_eight_step_evaluation();

						// Update existence
						Update_Existence(TempSolution, 0);
						///////////////////////////////////////////////////////////////////////////


						Test_UpdateBestSolution(0, TempSolution, intensity, iteration, 4);


						// Deallocation of memory

						if (intensity == 1)
						{
							delete req;
							delete found;
						}
						else
						{
							delete[] req;
							delete[] found;
						}

					}


*/

	}

	void Test_UpdateBestSolution(int s, solution *TempSolution, int &intensity, int iteration, int operation, int select = -1, int pick = -1, int grade = -1)
	{
		TempSolution[s].OfflineEvaluateCost();
		//TempSolution[s].Offline_eight_step_evaluation();

		//if (TempSolution[s].isComplete())
			if (TempSolution[s].Cost.getFeasibility())
			{

				float current_cost = Coefficient[s].cost_function(TempSolution[s].Cost);

				if (best_cost[s] == 0 || (current_cost < best_cost[s]))
				{

					best_cost[s] = current_cost;
					Best_Solution[s].Copy(TempSolution[s]);
					Best_Solution[s].OfflineEvaluateCost();
					//Best_Solution[s].Offline_eight_step_evaluation();

					if (Trigger_restart_mechanism)
						noImprovement[s] = 0;

					switch (operation)
					{
					case 1: //ALNS-case
					{
						printf("%0.0f	%0.01f	%0.01f	%d	%d	%d	%d	Sel:%d	Pick:%d	Grade:%d	%d	%s\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[s].cost_function(Best_Solution[s].Cost), Best_Solution[s].Cost.travel_cost, Best_Solution[s].total_req_served, fsm[0].FleetSize, iteration, intensity, select, pick, grade, 1, Best_Solution[s].Cost.isFeasible == 0 ? "NOT Feasible" : "");
						break;
					}
					case 2:
					{
						printf("%0.0f	%0.01f	%0.01f	%d	%d	%d	%d			%d	%s\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[s].cost_function(Best_Solution[s].Cost), Best_Solution[s].Cost.travel_cost, Best_Solution[s].total_req_served, fsm[0].FleetSize, iteration, intensity, 1, Best_Solution[s].Cost.isFeasible == 0 ? "NOT Feasible" : "");
						break;
					}
					case 3:
					{
						printf("%0.0f	%0.01f	%0.01f	%d	%d	%d	%d	ejection	%d	%s\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[s].cost_function(Best_Solution[s].Cost), Best_Solution[s].Cost.travel_cost, Best_Solution[s].total_req_served, fsm[0].FleetSize, iteration, intensity, 1, Best_Solution[s].Cost.isFeasible == 0 ? "NOT Feasible" : "");
						break;
					}
					case 4:
					{
						printf("%0.0f	%0.01f	%0.01f	%d	%d	%d	%d	intra		%d	%s\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[s].cost_function(Best_Solution[s].Cost), Best_Solution[s].Cost.travel_cost, Best_Solution[s].total_req_served, fsm[0].FleetSize, iteration, intensity, 1, Best_Solution[s].Cost.isFeasible == 0 ? "NOT Feasible" : "");
						break;
					}
					default:
					{
						break;
					}
					}

					//intensity = 1;

					//if (DO_intra == true)
						//IntraRoute_Insertion(TempSolution, intensity, iteration);

				}
			}




	}

	///// UTILITY FUNCTIONS

	void SortRequest(int& req_list)
	{
		int* list = &req_list;

		float **Array = new float*[n];
		for (int i = 0; i < n; i++)
			Array[i] = new float[2];

		for (int i = 0; i<n; i++)
		{
			Array[i][0] = i;
			Array[i][1] = Problem[0].Request[i].pickup.LatestTime + Problem[0].Request[i].dropoff.EarliestTime;
		}

		for (int i = 0; i<n; i++)
		{
			int biggest_j = 0;
			float biggest_value = 0;

			for (int j = 0; j<n - i; j++)
			{
				if (Array[j][1] >  biggest_value)
				{
					biggest_j = j;
					biggest_value = Array[j][1];
				}
			}

			float temp[2];

			temp[0] = Array[n - i - 1][0];
			temp[1] = Array[n - i - 1][1];

			Array[n - i - 1][0] = Array[biggest_j][0];
			Array[n - i - 1][1] = Array[biggest_j][1];

			Array[biggest_j][0] = temp[0];
			Array[biggest_j][1] = temp[1];
		}

		for (int i = 0; i<n; i++)
		{
			list[i] = Array[i][0];
		}

		for (int i = 0; i < n; i++)
			delete[] Array[i];

		delete[] Array;
	}

	void randomRequest(int& req_list)
	{
		int* list = &req_list;
		bool *isInList = new bool[n];

		for (int i = 0; i < n; i++)
			isInList[i] = 0;

		for (int i = 0; i<n; i++)
		{
			srand(time_now++);
			int count = rand() % (n - i);

			for (int j = 0; j<n; j++)
				if (isInList[j] == false)
				{
					if (count == 0)
					{
						isInList[j] = true;
						list[i] = j;
						break;
					}
					else
						count--;
				}
		}

		delete[] isInList;
	}

	void Initialize(int s, int reqID, solution *Solution, problem *Problem)
	{

		coefficient *Coefficient = new coefficient;

		sol *gpu = new sol[Solution[s].allinsertion];;
		cost *Full_Cost_Per_Route = new cost[TotalVehicles];
		solution *Temporary = new solution;
		cost *Temp = new cost[Solution[s].allinsertion];

		Coefficient[0].reset();

		Temporary[0].Copy(Solution[s]);
		Temporary[0].OfflineEvaluateCost();

		for (int k = 0; k < TotalVehicles; k++)
			Full_Cost_Per_Route[k] = Temporary[0].OfflineFetchRouteCost(k);


		// 1) SETUP
		int indices = 0, start = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			int size = Solution[s].Vehicle[k].size;

			for (int limiter = size; limiter >= 0; limiter--)
			{
				for (start = 0; start <= limiter; start++)
				{
					gpu[indices].start = start;
					gpu[indices].gap = size - limiter;
					gpu[indices].vehid = k;
					indices++;

					//if (indices >= Solution[s].allinsertion)
					//break;
				}

				//if (indices >= Solution[s].allinsertion)
				//break;
			}

			if (indices >= Solution[s].allinsertion)
				break;
		}

		// 2) INSERTION
		for (int idx = 0; idx < Solution[s].allinsertion; idx++)
		{
			gpu[idx].request = reqID - 1;
			gpu[idx].fixation(Solution[s].Vehicle[gpu[idx].vehid].size);
			gpu[idx].route[gpu[idx].start] = reqID;
			gpu[idx].route[gpu[idx].start + gpu[idx].gap + 1] = reqID + n;


			if (gpu[idx].routesize == 2)
			{
			}
			else
				for (int tid = 0; tid < gpu[idx].routesize; tid++)
				{
					bool boolone = tid >= (gpu[idx].start);
					bool booltwo = tid >= (gpu[idx].start + gpu[idx].gap);
					gpu[idx].route[tid + boolone + booltwo] = Solution[s].Vehicle[gpu[idx].vehid].path[tid];
				}

		}


		// 3) EVALUATION
		for (int idx = 0; idx < Solution[s].allinsertion; idx++)
		{
			gpu[idx].CPU_evaluateByRoute(Problem);
		}


		// 4) FEASIBILITY
		for (int idx = 0; idx < Solution[s].allinsertion; idx++)
		{
			Temp[idx].reset();

			for (int k = 0; k < TotalVehicles; k++)
			{

				if (k == gpu[idx].vehid)
				{
					Temp[idx].travel_cost += gpu[idx].Vehicle.Cost.travel_cost;
					Temp[idx].time_window += gpu[idx].Vehicle.Cost.time_window;
					Temp[idx].ride_time += gpu[idx].Vehicle.Cost.ride_time;
					Temp[idx].load += gpu[idx].Vehicle.Cost.load;
					Temp[idx].duration += gpu[idx].Vehicle.Cost.duration;
					Temp[idx].excessrideTime += gpu[idx].Vehicle.Cost.excessrideTime;
				}
				else
				{
					Temp[idx].travel_cost += Full_Cost_Per_Route[k].travel_cost;
					Temp[idx].time_window += Full_Cost_Per_Route[k].time_window;
					Temp[idx].ride_time += Full_Cost_Per_Route[k].ride_time;
					Temp[idx].load += Full_Cost_Per_Route[k].load;
					Temp[idx].duration += Full_Cost_Per_Route[k].duration;
					Temp[idx].excessrideTime += Full_Cost_Per_Route[k].excessrideTime;
				}

			}

			Temp[idx].isFeasible = Temp[idx].getFeasibility();
			gpu[idx].isFeasible = Temp[idx].getFeasibility();


		}



		// 5) REDUCTION
		int local_ID = 0;
		float best_cost = 0;
		float current_cost = 0;

		for (int idx = 0; idx < Solution[s].allinsertion; idx++)
		{
			//if (Temp[idx].isFeasible)
			if (Coefficient[0].cost_function(Temp[idx]) > 0)
				if ((best_cost == 0) || (Temp[idx].travel_cost < best_cost))
				{
					best_cost = Temp[idx].travel_cost;
					local_ID = idx;
				}
		}


		if (best_cost > 0)
		{
			int idx = local_ID;
			Solution[s].ReInsertRequest(gpu[idx].request, gpu[idx].vehid, gpu[idx].start, gpu[idx].start + gpu[idx].gap + 1);
		}

		//Solution[0].printinput();

		// 6) UPDATION
		Solution[s].ReBoot();



		delete[] gpu;
		delete[] Full_Cost_Per_Route;
		delete[] Temp;
		delete Temporary;
		delete Coefficient;
	}


	void Update_Existence(solution *TempSolution, int s)
	{
		//req on existence index starts from zero.
		// 1 - exist and 0 - not exist


		for (int i = 0; i < n; i++)
			existence[s][i] = 0;

		for (int k = 0; k < TotalVehicles; k++)
			for (int j = 0; j < TempSolution[s].Vehicle[k].size; j++)
				if (TempSolution[s].Vehicle[k].path[j] <= n)
					existence[s][TempSolution[s].Vehicle[k].path[j] - 1] = 1;

	}

	void Trigger_Restart_Mechanism(int s, float *currentbestcost, int iteration)
	{

		if (Trigger_restart_mechanism == true && best_cost[s] > 0)
		{
			// 1) update counter

			if (Coefficient[s].cost_function(Solution[s].Cost) > 1.1*currentbestcost[s])
				noImprovement[s]++;


			// 2) Trigger restart mechanism


			if (noImprovement[s] > noImprovementLimit[s])
			{
				Solution[s].Copy(Best_Solution[s]);
				Solution[s].OfflineEvaluateCost();

				Update_Existence(Solution, s);

				noImprovement[s] = 0;

				printf("reset\n");

			}
		}

	}

	void IH_Construction_heuristic(int s)
	{
		/////////////////////////////////////////
		int *request_list = new int[n];

		SortRequest(request_list[0]);
		Solution[s].Boot();

		for (int i = 0; i < n; i++)
		{
			int req = request_list[i];

			//this works only if feasible insertion found for all requests. (not possible for pr07 to p10).
			//Initialize(s, req + 1, Solution, &Solution[s].Problem);
		}

		Solution[s].IdentifyCurrvehicle();

		delete request_list;
		/////////////////////////////////////////
	}


};


