struct LNS
{

	clock_t TotalTime;

	bool DO_intra = true;
	bool Trigger_restart_mechanism = true; // CAUTION: THIS GIVES RUN-TIME ERROR FOR SOME INSTANCES!!!

	// parameters
	float *temperature;
	float *maxTemperature;
	float *minTemperature;
	float *coolingRate;
	int *noImprovement;
	int *noImprovementLimit;
	FSM *fsm;

	int *existence;
	float *best_cost;
	solution *Best_Solution;
	problem *Problem;
	solution *Solution;
	coefficient *Coefficient;

	explore_neighborhood *Explore_Neighborhood;

	reinsertion *Reinsertion;
	removal *Removal;

	insertion *Insertion;
	noise *Noise;

	LNS()
	{

	}

	~LNS()
	{

	}

	void Creation()
	{
		Reinsertion = new reinsertion;
		Removal = new removal;
		Insertion = new insertion[n];
		Noise = new noise;

		Explore_Neighborhood = new explore_neighborhood;
		Explore_Neighborhood[0].creation(); //cuda-mallocManaged is used here!! BEWARE!!
		Explore_Neighborhood[0].create_GPU_memories(); // heavily-time-consuming-task

		// parameters
		temperature = new float;
		maxTemperature = new float;
		minTemperature = new float;
		coolingRate = new float;
		fsm = new FSM;
		noImprovement = new int;
		noImprovementLimit = new int;

		existence = new int[n];
		for (int i = 0; i < n; i++)
			existence[i] = 0;
		best_cost = new float;
		Best_Solution = new solution;
		CHECK(cudaMallocHost((void**)&Problem, sizeof(problem)));
		CHECK(cudaMallocHost((void**)&Solution, sizeof(solution)));
		Coefficient = new coefficient;

		Parameterization();
		Import_Problem();
		Construct_Initial_Solution();

		///////////////////////////////////////////////////////
		/////// Just printing the total GPU memory usage //////
		///////////////////////////////////////////////////////

		size_t free, total;
		CHECK(cudaMemGetInfo(&free,&total));
		printf("GPU Memory usage info:\n Free = %0.02f GB,Total = %0.02f\n", ((float)free)/(1000000000), ((float)total/1000000000));
		printf("<<<BEWARE: if you explore neighborhood without using struct Explore_Neighborhood, there will error. \n Because vnd[0].gpu/temp has not been given enough memory!!!>>>\n");
	}

	void Destruction()
	{
		delete Reinsertion;
		delete Removal;
		delete[] Insertion;
		delete Noise;

		//below-three must be in this order.
		Explore_Neighborhood[0].destroy_GPU_memories(); // heavily-time-consuming-task
		Explore_Neighborhood[0].destroy(); //cuda-mallocManaged is used here!! BEWARE!!
		delete Explore_Neighborhood;

		delete temperature;
		delete maxTemperature;
		delete minTemperature;
		delete coolingRate;
		delete fsm;
		delete noImprovement;
		delete noImprovementLimit;

		delete[] existence;
		delete best_cost;
		delete Best_Solution;
		CHECK(cudaFreeHost(Problem));
		CHECK(cudaFreeHost(Solution));
		delete Coefficient;



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

				Coefficient[0].initialize();

				// SA-parameters
				coolingRate[0] = 0.99975;
				minTemperature[0] = 0.001;
				maxTemperature[0] = Defined_MaxTemperature; //20
				temperature[0] = maxTemperature[0];

				// restart parameter
				noImprovement[0] = 0;
				noImprovementLimit[0] = 200;

				// fsm parameter
				fsm[0].initialize();
				fsm[0].FSM_Enable = false;
	}

	void Import_Problem()
	{
		// 2) Import the problem instance

				BootProblem(Problem[0]);

				//Store Problem in GPU's Constant Memory
				Update_ConstantMemory(Problem);
				CHECK(cudaDeviceSynchronize());

	}

	void Construct_Initial_Solution()
	{
		// 2) Construct an initial solution

				Solution[0].Problem.Copy(Problem[0]);

				printf("Randomized_construction_heuristic is used!!\n");

				if (0)
					IH_Construction_heuristic();
				else
					Solution[0].RandomShuffle();

				Solution[0].IdentifyCurrvehicle();
				Solution[0].OfflineEvaluateCost();

				// Update existence
				Update_Existence(Solution);

				// Update
				int dummy;
				best_cost[0] = 0;
				Test_UpdateBestSolution(Solution, dummy, -1, -1);

		// noise mechanism (can be set only after booting problem)
				Noise[0].update_noise_parameters(Solution);

		// Fill candidate edges
				Explore_Neighborhood[0].fill_Candidate_Edges(Problem);

	}

	void Optimize(int num_of_iteration = 500000, int TimeLimit = 120)
	{

		TotalTime = clock();

		// Optimization begins

		solution *Temp_Solution = new solution;
		float currentbestcost = best_cost[0];

		int intensity = 1, reduction_choice = 0, iteration = 1;
		for (iteration = 1; iteration < 1000000000; iteration++)
		{

					Temp_Solution[0].Copy(Solution[0]);
					Temp_Solution[0].OfflineEvaluateCost();

			// utility: trigger exit, print itr time
					if ((float)((clock() - TotalTime) / CLOCKS_PER_SEC) > TimeLimit)
						goto displayResult;

					if (iteration % 500 == 0)
						printf("At Iter %d, TimeElapsed: %0.02f\n", iteration, (float)((clock() - TotalTime) / CLOCKS_PER_SEC));

			// 1) update parameters

					Coefficient[0].update(Solution[0].Cost);
					if (iteration % 20 == 0)
						Coefficient[0].reset();

					reduction_choice = 0;//rand() % 3; // 0 - reduce normal cost_function, 1 - reduce cost_function_with_batteryViols, 2 - reduce cost_function_with_batteryViols_2
					intensity = 1 + rand() % 5;
					if (intensity == 0 || intensity > n)
						printf("ERRORRR!!!! intensity is greater than n\n");

					int select = rand() % 3; //0-shaw, 1-random, 2-worst
					int pick = rand() % 3; //0-greedy, 1-regret, 2-sameOrder
					int grade = rand() % (m-1);

			// 2) explore the neighborhood

					if (0)//pick == 0 || pick == 1)//iteration % 10 == 0)
						Remove_and_Reinsert_Thoroughly(select, pick, grade, Solution, intensity, iteration, reduction_choice);
					else
					{
						for (int i=1; i<=2; i++)
							Remove_and_Reinsert_Loosely(i, Solution, intensity, iteration, reduction_choice);
					}


			// 3) restart mechanism

					Trigger_Restart_Mechanism(currentbestcost, iteration);

			// 4) SA: acceptance criterion

					if (Solution[0].Cost.getFeasibility())
					{
						if (Coefficient[0].cost_function(Solution[0].Cost) < currentbestcost)
							currentbestcost = Coefficient[0].cost_function(Solution[0].Cost);
						// 1) New best solution

						int intensity_dummy = 0;
						Test_UpdateBestSolution(Solution, intensity_dummy, iteration, 1);

					}
					else if ((Coefficient[0].cost_function(Solution[0].Cost) < Coefficient[0].cost_function(Temp_Solution[0].Cost)))
					{
						// 2) New improving solution
					}
					else
						{
							float p = exp((-1)*(Coefficient[0].cost_function(Solution[0].Cost) - Coefficient[0].cost_function(Temp_Solution[0].Cost)) / temperature[0]);

							float random_num = rand() % 10;
							random_num = random_num / 10;

							if (random_num < p)
							{
								// 3) New SA-accepted solution
							}
							else
							{
								Solution[0].Copy(Temp_Solution[0]);
								Solution[0].OfflineEvaluateCost();
							}
						}


			// SA: cooling down & temperature-check

					temperature[0] *= coolingRate[0];
					if (temperature[0] < minTemperature[0])
						temperature[0] = maxTemperature[0];

			/////////////////////////////////////////////
			//// FSM module: minimize the fleet size ////
			/////////////////////////////////////////////

			if (best_cost[0] > 0 && fsm[0].FSM_Enable == true)
			{
				//fsm[0].FSM_module_minimize_fleet_size(Solution, Best_Solution, best_cost[0], currentbestcost);
			}


		}

	displayResult:
		{

			// deleting Temp_Solution created for SA-acceptance
			delete Temp_Solution;

			if (best_cost[0] > 0)
			{
				Best_Solution[0].OfflineEvaluateCost();
				Best_Solution[0].print_with_battery_level();

				printf("Solution is %s\n", Best_Solution[0].Cost.isFeasible == true ? "Feasible" : " NOT Feasible");
				//if (Best_Solution[0].Cost.getFeasibility())
				//printf("-->> Cost = %f\n", Best_Solution[0].Cost.travel_cost);
			}
			else
				printf("CANNOT FIND FEASIBLE SOLUTION WITHIN GIVEN TIME!!\n");
			printf("Total Iterations: %d\n", iteration);
			printf("Execution time: %0.0f\n\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000);
		}
	}

	void Remove_and_Reinsert_Thoroughly(int select, int pick, int grade, solution *TempSolution, int &intensity, int iteration, int choice)
	{

		int isGPU = Original_isGPU;
		int isFull = Original_isFull;
		int isRev = Original_isRev;

		// Allocation of memory

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

		// 1) Remove many requests together
		// 1.1) choose operator
		Removal[0].operator_selection(select, intensity, TempSolution);
		for (int p = 0; p < intensity; p++)
		{
			req[p] = Removal[0].D_list[p];
			found[p] = 0;
		}

		// 1.2) remove requests
		for (int p = 0; p < intensity; p++)
		{
			if (existence[req[p]])
			{
				TempSolution[0].RemoveRequest(req[p], TempSolution[0].Request[req[p]].currvehicle);
				existence[req[p]] = 0;
			}
		}

		// 2) Sequentially reinsert those requests

		///////////////////////////////////////////////////////////////////////////
		for (int total_req_on_bank = intensity; total_req_on_bank > 0; total_req_on_bank--)
		{
			//printf("CheckPt 2.1\n");

			int *request = new int[total_req_on_bank];

			// Activate current bank si
			int index = 0;
			for (int j = 0; j < intensity; j++)
				if (existence[req[j]] == 0)
					request[index++] = req[j];


			// 2.1) Evaluate insertions of each request

			for (int p = 0; p < total_req_on_bank; p++)
			{
				// 1) explore neighborhood
				solution *NewSolution = new solution;
				NewSolution[0].Copy(TempSolution[0]);
				NewSolution[0].OfflineEvaluateCost();

				Explore_Neighborhood[0].Exploration_Decision_Maker(isGPU, isFull, isRev, p, NewSolution, Problem, request, Insertion, Coefficient, fsm, choice);

				delete NewSolution;
			}

			// 2.2) prepare the request_bank
			int noise_application = 0;// grade = rand() % (m-1);//0;
			Reinsertion[0].operator_selection(pick, Removal, Insertion, grade, noise_application, Noise, total_req_on_bank);

			// 2.3) re-insert only the best request from request bank

			int best_p = Reinsertion[0].sequence[0];
			TempSolution[0].ReInsertRequest(Insertion[best_p].request, Insertion[best_p].vehid, Insertion[best_p].start, Insertion[best_p].start + Insertion[best_p].gap + 1);
			TempSolution[0].OfflineEvaluateCost();

			delete[] request;


			// Update existence
			Update_Existence(TempSolution);


		}
		///////////////////////////////////////////////////////////////////////////

		// Update the best solution
		Test_UpdateBestSolution(TempSolution, intensity, iteration, 1, select, pick, grade);

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

	void Remove_and_Reinsert_Loosely(int removalOp, solution *TempSolution, int &intensity, int iteration, int choice)
	{
		int isGPU = Original_isGPU;
		int isFull = Original_isFull;
		int isRev = Original_isRev;

		// Allocation of memory

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

		//select ==> 0 shaw, 1 random, 2 worst
		int select = removalOp;

		// 1.1) choose operator
		Removal[0].operator_selection(select, intensity, TempSolution);
		for (int p = 0; p < intensity; p++)
		{
			req[p] = Removal[0].D_list[p];
			found[p] = 0;
		}

		// 1) Remove many requests together

		for (int p = 0; p < intensity; p++)
		{
			if (existence[req[p]])
			{
				TempSolution[0].RemoveRequest(req[p], TempSolution[0].Request[req[p]].currvehicle);
				existence[req[p]] = 0;
			}
		}

		// 2) Sequentially reinsert those requests

		///////////////////////////////////////////////////////////////////////////

		for (int p = 0; p < intensity; p++)
		{

			Explore_Neighborhood[0].Exploration_Decision_Maker(isGPU, isFull, isRev, p, TempSolution, Problem, req, Insertion, Coefficient, fsm, choice);

		}
		TempSolution[0].OfflineEvaluateCost();
		//TempSolution[0].Offline_eight_step_evaluation();

		// Update existence
		Update_Existence(TempSolution);
		///////////////////////////////////////////////////////////////////////////


		// Update the best solution
		Test_UpdateBestSolution(TempSolution, intensity, iteration, select+1);

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

	void IntraRoute_Insertion(solution *TempSolution, int &intensity, int iteration)
	{

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
							if (existence[req[p]])
							{
								TempSolution[0].RemoveRequest(req[p], TempSolution[0].Request[req[p]].currvehicle);
								existence[req[p]] = 0;
							}
						}

						// 2) Sequentially reinsert those requests

						///////////////////////////////////////////////////////////////////////////

						for (int p = 0; p < intensity; p++)
						{

							Explore_Neighborhood[0].Exploration_Decision_Maker(0, 1, 0, p, TempSolution, Problem, req, Insertion, Coefficient, fsm, 0, 1);

						}
						TempSolution[0].OfflineEvaluateCost();
						//TempSolution[0].Offline_eight_step_evaluation();

						// Update existence
						Update_Existence(TempSolution);
						///////////////////////////////////////////////////////////////////////////


						Test_UpdateBestSolution(TempSolution, intensity, iteration, 4);


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




	}

	void Test_UpdateBestSolution(solution *TempSolution, int &intensity, int iteration, int operation, int select = -1, int pick = -1, int grade = -1)
	{
		TempSolution[0].OfflineEvaluateCost();

		if (TempSolution[0].isComplete())
			if (TempSolution[0].Cost.getFeasibility())
			{

				float current_cost = Coefficient[0].cost_function(TempSolution[0].Cost);

				if (best_cost[0] == 0 || (current_cost < best_cost[0]))
				{

					best_cost[0] = current_cost;
					Best_Solution[0].Copy(TempSolution[0]);
					Best_Solution[0].OfflineEvaluateCost();

					if (Trigger_restart_mechanism)
						noImprovement[0] = 0;

					switch (operation)
					{
					case 1: //ALNS-case
					{
						printf("%0.0f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%d	%s	%d	%d	Sel:%d	Pick:%d	Grade:%d	%d\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[0].cost_function(Best_Solution[0].Cost), fsm[0].FleetSize, Best_Solution[0].Cost.travel_cost, Best_Solution[0].Cost.excessrideTime, Best_Solution[0].Cost.batteryViolation, Best_Solution[0].Cost.max_battery_Violation, Best_Solution[0].Cost.isFeasible == 0 ? "NOT Feasible" : "Feasible", iteration, intensity, select, pick, grade, 1);
						break;
					}
					case 2:
					{
						printf("%0.0f	%f	%f	%0.02f	%0.02f	%f	%d	%s	%d	%d			%d\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[0].cost_function(Best_Solution[0].Cost), Best_Solution[0].Cost.travel_cost, Best_Solution[0].Cost.excessrideTime, Best_Solution[0].Cost.batteryViolation, Best_Solution[0].Cost.max_battery_Violation, fsm[0].FleetSize, Best_Solution[0].Cost.isFeasible == 0 ? "NOT Feasible" : "Feasible", iteration, intensity, 1);
						break;
					}
					case 3:
					{
						printf("%0.0f	%f	%f	%0.02f	%0.02f	%f	%d	%s	%d	%d	ejection	%d\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[0].cost_function(Best_Solution[0].Cost), Best_Solution[0].Cost.travel_cost, Best_Solution[0].Cost.excessrideTime, Best_Solution[0].Cost.batteryViolation, Best_Solution[0].Cost.max_battery_Violation, fsm[0].FleetSize, Best_Solution[0].Cost.isFeasible == 0 ? "NOT Feasible" : "Feasible", iteration, intensity, 1);
						break;
					}
					case 4:
					{
						printf("%0.0f	%f	%f	%0.02f	%0.02f	%f	%d	%s	%d	%d	intra		%d\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[0].cost_function(Best_Solution[0].Cost), Best_Solution[0].Cost.travel_cost, Best_Solution[0].Cost.excessrideTime, Best_Solution[0].Cost.batteryViolation, Best_Solution[0].Cost.max_battery_Violation, fsm[0].FleetSize, Best_Solution[0].Cost.isFeasible == 0 ? "NOT Feasible" : "Feasible", iteration, intensity, 1);
						break;
					}
					default:
					{
						printf("%0.0f	%f	%f	%0.02f	%0.02f	%f	%d	%s	%d	%d	default		%d\n", (float)(clock() - TotalTime) / CLOCKS_PER_SEC * 1000, Coefficient[0].cost_function(Best_Solution[0].Cost), fsm[0].FleetSize, Best_Solution[0].Cost.travel_cost, Best_Solution[0].Cost.excessrideTime, Best_Solution[0].Cost.batteryViolation, Best_Solution[0].Cost.max_battery_Violation, Best_Solution[0].Cost.isFeasible == 0 ? "NOT Feasible" : "Feasible", iteration, intensity, 1);
						break;
					}
					}

					intensity = 1;

					if (DO_intra == true)
						IntraRoute_Insertion(TempSolution, intensity, iteration);

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

	void Update_Existence(solution *TempSolution)
	{
		//req on existence index starts from zero.
		// 1 - exist and 0 - not exist

		for (int k = 0; k < TotalVehicles; k++)
			for (int j = 0; j < TempSolution[0].Vehicle[k].size; j++)
				if (TempSolution[0].Vehicle[k].path[j] <= n)
					existence[TempSolution[0].Vehicle[k].path[j] - 1] = 1;

	}

	void Trigger_Restart_Mechanism(float currentbestcost, int iteration)
	{

		if (Trigger_restart_mechanism == true && best_cost[0] > 0)
		{
			// 1) update counter

			if (Coefficient[0].cost_function(Solution[0].Cost) > 1.1*currentbestcost)
				noImprovement[0]++;


			// 2) Trigger restart mechanism


			if (noImprovement[0] > noImprovementLimit[0])
			{
				Solution[0].Copy(Best_Solution[0]);
				Solution[0].OfflineEvaluateCost();

				Update_Existence(Solution);

				noImprovement[0] = 0;

				printf("reset\n");

			}
		}

	}

	void IH_Construction_heuristic()
	{
		/////////////////////////////////////////
		int *request_list = new int[n];

		SortRequest(request_list[0]);
		Solution[0].Boot();

		for (int i = 0; i < n; i++)
		{
			int req = request_list[i];

			//this works only if feasible insertion found for all requests. (not possible for pr07 to p10).
			Initialize(req + 1, Solution, &Solution[0].Problem);
		}

		Solution[0].IdentifyCurrvehicle();

		delete request_list;
		/////////////////////////////////////////
	}


};

