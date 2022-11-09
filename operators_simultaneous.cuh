
struct simple_request
{
	int req;
	float cost_effect;

	int vehicle;
};

struct simple_vehicle
{
	int route[ExpectedPath];
	int size;

	__host__ __device__
		simple_vehicle()
	{
		size = 0;
	}

	__host__ __device__
		~simple_vehicle() {}
};

struct simple_solution
{
	simple_request request[n];
	simple_vehicle vehicle[m];

	int *route_order;
	float *costfunction;

	float total_cost;

	__host__ __device__
		simple_solution()
	{
		route_order = new int[m];
		costfunction = new float[m];
		total_cost = 0;
	}

	__host__ __device__
		~simple_solution()
	{
		delete[] costfunction;
		delete[] route_order;
	}

	__host__ __device__
		void Copy(solution &inSolution)
	{
		// Copy route path for all vehicles.

		for (int k = 0; k < m; k++)
		{
			for (int j = 0; j < inSolution.Vehicle[k].size; j++)
			{
				int node = inSolution.Vehicle[k].path[j] - 1;
				if (node < n)
				{
					vehicle[k].route[vehicle[k].size] = node;
					vehicle[k].size++;
				}
			}
		}

		// Map vehicle label for all requests.

		for (int k = 0; k < m; k++)
		{
			for (int j = 0; j < vehicle[k].size; j++)
			{
				int node = vehicle[k].route[j];
				request[node].req = k;
			}
		}


		/*for (int i = 0; i < n; i++)
		{
		for (int k = 0; k < m; k++)
		{
		for (int j = 0; j < inSolution.Vehicle[k].size; j++)
		{
		int node = inSolution.Vehicle[k].path[j] - 1;
		if (node < n)
		{
		request[m].vehicle = k;
		}
		}
		}
		}*/
	}

};

struct probBound
{
	float lowerbound;
	float upperbound;
};

struct FSM
{
	int FSM_Enable;

	int FleetSize;
	int Minimum_Possible_FleetSize;

	int FSM_NotImproved;
	int FSM_NotImprovedLimit;

	int BarredList[TotalVehicles];
	int BarredCount;

	FSM() {}

	~FSM() {}

	void initialize()
	{
		FSM_Enable = false;

		FleetSize = m;
		Minimum_Possible_FleetSize = m - 1;

		FSM_NotImproved = 0;
		FSM_NotImprovedLimit = 50;

		BarredCount = 0;
		for (int k = 0; k < TotalVehicles; k++)
			BarredList[k] = 0;

		// Disable one vehicle
		//BarredList[m - 1] = 1;
	}

	void FSM_module_minimize_fleet_size(solution *Solution, solution *Best_Solution, float best_cost, float &currentbestcost)
	{
		// FSM code
		if (currentbestcost == 0 || best_cost < currentbestcost)
		{
			currentbestcost = best_cost;
			FSM_NotImproved = 0;
		}
		else
			FSM_NotImproved++;


		if (FleetSize > Minimum_Possible_FleetSize)
		{
			int EmptyVehicles = 0;
			for (int k = 0; k < m; k++)
				if (Best_Solution[0].Vehicle[k].size == 0)
					EmptyVehicles++;

			if (FSM_NotImproved > FSM_NotImprovedLimit)
			{
				// Reduce_Fleet_Size_by_One

				// Check whether Feasible Solution is there with required Fleet for further minimization.
				if (BarredCount == EmptyVehicles)
					if (Best_Solution[0].Cost.getFeasibility())
						//if (Best_Solution[0].isComplete())
						{
							printf("Before Turing:\n");
							Best_Solution[0].printinput();
							Solution[0].Copy(Best_Solution[0]);
							Solution[0].ReBoot();

							FSM_NotImproved = 0;
							best_cost = 10000;

							//Fleet Minimization (Bare one vehicle & remove all requests from them)
							Minimize_vehicle_size_by_one(Solution);

							//Increase Barred Count
							BarredCount++;
							FleetSize--;
							printf("Barred !!\n");

						}
			}


		}
	}

	void Minimize_vehicle_size_by_one(solution *Solution)
	{
		// 1) Find Target Vehicle (Vehicle with minimum number of requests)
		int min_k;
		int min_size;

		for (int k = 0; k < m; k++)
			if (BarredList[k] == 0)
			{
				min_k = k;
				min_size = Solution[0].Vehicle[k].size;
				break;
			}

		for (int k = 0; k < m; k++)
			if (BarredList[k] == 0)
				if (Solution[0].Vehicle[k].size < min_size)
				{
					min_size = Solution[0].Vehicle[k].size;
					min_k = k;
				}

		// 2) Bare it (Add it into Barred list)
		BarredList[min_k] = 1;

		// 3) Wipe out (Cleansing)
		for (int j = 0; j < Solution[0].Vehicle[min_k].size; j++)
			Solution[0].Vehicle[min_k].path[j] = 0;

		Solution[0].Vehicle[min_k].size = 0;

		Solution[0].ReBoot();


	}

};

struct noise
{
	// 0 - pure, 1 - noise

	// parameters
	float eta;
	float max_dij;
	float max_N;
	float min_N;

	float upper_probability_limit;

	float all_score[2];

	float score[2];
	float weight[2];

	float times_used[2];
	float all_times_used[2];

	noise()
	{
		for (int i = 0; i < 2; i++)
		{
			weight[i] = 1;
			times_used[i] = 0;
			all_times_used[i] = 0;
			score[i] = 0;

			all_score[i] = 0;
		}

		upper_probability_limit = weight[0] / (weight[0] + weight[1]);

		printf("UPL: %f\n", upper_probability_limit);
	}

	~noise()
	{

	}

	void update_noise_parameters(solution *Solution)
	{


		// 1) set eta

		eta = 0.025; // suggested by Ropke (2006)

					 // 2) update max_dij

		max_dij = 0;

		for (int i = 0; i<Var; i++)
			for (int j = 0; j < Var; j++)
			{
				if (max_dij == 0 || Solution[0].Problem.TravelTime[i][j] > max_dij)
					max_dij = Solution[0].Problem.TravelTime[i][j];
			}

		// 3) update noise intervals

		max_N = eta * max_dij;
		min_N = (-1) * max_N;

		printf("Noise-window:	[%f, %f]\n", min_N, max_N);
	}

	void adjust_weights(float r_factor)
	{
		// call this at end of every segment

		// 1) adjust weights
		for (int i = 0; i < 2; i++)
		{
			if (times_used[i] < 0.1)
				weight[i] = weight[i] * (1 - r_factor) + (0);
			else
				weight[i] = weight[i] * (1 - r_factor) + (r_factor * score[i] / times_used[i]);

		}

		// 2) reset score and times used
		for (int i = 0; i < 2; i++)
		{
			all_score[i] += score[i];

			score[i] = 0;
			times_used[i] = 0;
		}

		// 3) set upper_probability limit (< noise, >pure).

		upper_probability_limit = weight[0] / (weight[0] + weight[1]);

	}
};

struct operation
{
	int id;
	float score;
	float weight;

	// weight mechanism
	probBound *ProbBound;

	float times_used;
	float all_times_used;

	float all_score;

	operation()
	{
		weight = 1;
		times_used = 0;
		all_times_used = 0;
		score = 0;

		all_score = 0;

		ProbBound = new probBound;
	}

	~operation()
	{
		delete ProbBound;
	}

	void reset_score()
	{
		score = 0;
		//weight = 1 / total_used_operators;
	}

};

struct relatedness
{
	int *req_list;
	float *relatedness_degree;

	int size;

	float *alpha;
	float *beta;
	float *gamma;

	float sum_alpha;
	float sum_beta;
	float sum_gamma;

	relatedness()
	{

	}

	~relatedness()
	{

	}

	void construct(int inSize)
	{
		size = inSize;

		req_list = new int[size];
		relatedness_degree = new float[size];

		for (int i = 0; i < size; i++)
			relatedness_degree[i] = 0;

		alpha = new float[size];
		beta = new float[size];
		gamma = new float[size];

		sum_alpha = 0;
		sum_beta = 0;
		sum_gamma = 0;
	}

	void destruct()
	{
		delete[] req_list;
		delete[] relatedness_degree;

		delete[] alpha;
		delete[] beta;
		delete[] gamma;
	}

	void update_relatedness(int req, solution *TempSolution)
	{

		// 1) find relatedness parameters

		for (int i = 0; i < size; i++)
		{
			int couple_req = req_list[i];

			alpha[i] = ABS(TempSolution[0].Problem.TravelTime[couple_req] - TempSolution[0].Problem.TravelTime[req]) + ABS(TempSolution[0].Problem.TravelTime[couple_req + n] - TempSolution[0].Problem.TravelTime[req + n]);
			beta[i] = ABS(TempSolution[0].Request[couple_req].B - TempSolution[0].Request[req].B) + ABS(TempSolution[0].Request[couple_req + n].B - TempSolution[0].Request[req + n].B);
			gamma[i] = 0; // zero for homogeneous set of requests i.e. demand at each node is always 1 in cordeau instances.

			sum_alpha += alpha[i];
			sum_beta += beta[i];
			sum_gamma += gamma[i];

		}

		// 2) normalize parameters
		if (sum_alpha == 0 || sum_beta == 0)
		{
			printf("Stop: sum_alpha/beta became zero. WHY?\n");
			// system("PAUSE");
		}
		else
		{
			for (int i = 0; i < size; i++)
			{
				alpha[i] = (alpha[i] / sum_alpha);
				beta[i] = (beta[i] / sum_beta);
			}
		}

		// 3) find final_relatedness degree

		for (int i = 0; i < size; i++)
			relatedness_degree[i] = 2 * (alpha[i] + beta[i]) + gamma[i];

	}

	void sort_relatedness()
	{
		float temp_1;
		int temp_2;

		for (int i = 0; i < size; i++)
		{
			for (int j = i + 1; j < size; j++)
				if (relatedness_degree[i] > relatedness_degree[j]) //ascending sort: less the relatedness_degree => more the two requests are related.
				{
					temp_1 = relatedness_degree[i];
					relatedness_degree[i] = relatedness_degree[j];
					relatedness_degree[j] = temp_1;

					temp_2 = req_list[i];
					req_list[i] = req_list[j];
					req_list[j] = temp_2;
				}
		}

	}

	void setting(vector<int> vec)
	{
		for (int i = 0; i < size; i++)
			req_list[i] = vec[i];
	}

	float ABS(float value)
	{
		if (value < 0)
			value *= -1;
		else
			value *= 1;

		return value;
	}
};

struct worstcost
{
	int *req_list;
	float *cost_effect;
	int *currvehicle_list;

	int size;

	solution *Solution;
	coefficient *Coefficient;
	position *Position;

	worstcost()
	{

	}

	~worstcost()
	{

	}

	void construct(int inSize)
	{
		size = inSize;

		req_list = new int[size];
		cost_effect = new float[size];

		currvehicle_list = new int[size];

		Solution = new solution;
		Coefficient = new coefficient;
		Position = new position;

	}

	void destruct()
	{
		delete[] req_list;
		delete[] cost_effect;
		delete[] currvehicle_list;

		delete Solution;
		delete Coefficient;
		delete Position;
	}

	void update_costeffect(solution *TempSolution)
	{
		Coefficient[0].reset();

		// 1) Initialization

		Solution[0].Copy(TempSolution[0]);
		Solution[0].OfflineEvaluateCost();

		for (int i = 0; i < size; i++)
		{
			int node = req_list[i];
			currvehicle_list[i] = Solution[0].Request[node].currvehicle;
		}

		// 2) Main: update the cost_effect

		for (int i = 0; i < size; i++)
		{
			int node = req_list[i];
			int k = currvehicle_list[i];// Solution[0].Request[node].currvehicle;

			Solution[0].LocateRequest(Position, node);
			Solution[0].RemoveRequest(node, k);

			//NewSolution[0].OfflineEvaluateByRoute(k+1);
			Solution[0].OfflineEvaluateCost();

			cost_effect[node] = Coefficient[0].cost_function(Solution[0].Cost);

			Solution[0].ReInsertRequest(node, k, Position[0].start, Position[0].start + Position[0].gap);
		}

	}

	void sort_costeffect()
	{

		float temp_cost;
		int temp_req;

		for (int i = 0; i < size; i++)
		{
			for (int j = i + 1; j < size; j++)
				if (cost_effect[i] > cost_effect[j]) //ascending sort: less the cost_effect value => more the worstness of request placement.
				{
					temp_cost = cost_effect[i];
					cost_effect[i] = cost_effect[j];
					cost_effect[j] = temp_cost;

					temp_req = req_list[i];
					req_list[i] = req_list[j];
					req_list[j] = temp_req;
				}
		}

	}

	void setting(vector<int> vec)
	{
		for (int i = 0; i < size; i++)
			req_list[i] = vec[i];
	}

};

__global__ void update_cost_effect_Kernel(solution *Solution, solution *d_Solutions, float *d_cost_effect, int *d_reqList, int *d_size)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x; //threadIdx.x;
	int s = blockIdx.x;

	if (s < ExpectedScenarios)
	{
		//if (idx < d_size[s])
		if (idx + (n*s) < d_size[s])
		{
			//d_Solutions[idx + (n*s)] = (Solution[s]);
			d_Solutions[idx + (n*s)].Copy(Solution[s]);

			int req = d_reqList[idx + (n*s)];
			int k = d_Solutions[idx + (n*s)].Request[req].currvehicle;
			d_Solutions[idx + (n*s)].RemoveRequest(req, k);
			//d_Solutions[idx + (n*s)].RemoveRequest_Spl_Purpose(req, k);

			d_Solutions[idx + (n*s)].OfflineEvaluateCost();

			d_cost_effect[idx + (n*s)] = d_Solutions[idx + (n*s)].Cost.overall_ideal_cost;
		}
	}

}

__global__ void sort_cost_effect_Kernel(float *d_cost_effect, int *d_reqList, int *d_size)
{

		int s = blockIdx.x;

		__shared__ float temp_cost[ExpectedScenarios];
		__shared__ int temp_req[ExpectedScenarios];

		for (int i = 0; i < d_size[s]; i++)
			for (int j = i + 1; j < d_size[s]; j++)
				if (d_cost_effect[i + (n*s)] > d_cost_effect[j + (n*s)])
				{
					temp_cost[s] = d_cost_effect[i + (n*s)];
					d_cost_effect[i + (n*s)] = d_cost_effect[j + (n*s)];
					d_cost_effect[j + (n*s)] = temp_cost[s];

					temp_req[s] = d_reqList[i + (n*s)];
					d_reqList[i + (n*s)] = d_reqList[j + (n*s)];
					d_reqList[j + (n*s)] = temp_req[s];
				}


}

struct worstcostRemoval
{
	solution *Dev_Solution;

	solution *h_Solutions;
	solution *d_Solutions;

	float *h_cost_effect;
	float *d_cost_effect;

	int *h_reqList;
	int *d_reqList;

	int *h_size;
	int *d_size;

	worstcostRemoval() {}

	~worstcostRemoval() {}

	void create()
	{
		CHECK(cudaMalloc((void**)&Dev_Solution, ExpectedScenarios * sizeof(solution)));

		h_Solutions = new solution[n * ExpectedScenarios];
		CHECK(cudaMalloc((void**)&d_Solutions, n * ExpectedScenarios * sizeof(solution)));

		h_cost_effect = new float[n * ExpectedScenarios];
		CHECK(cudaMalloc((void**)&d_cost_effect, n * ExpectedScenarios * sizeof(float)));

		h_reqList = new int[n * ExpectedScenarios];
		CHECK(cudaMalloc((void**)&d_reqList, n * ExpectedScenarios * sizeof(int)));

		h_size = new int[ExpectedScenarios];
		CHECK(cudaMalloc((void**)&d_size, ExpectedScenarios * sizeof(int)));
	}

	void destroy()
	{
		CHECK(cudaFree(Dev_Solution));

		delete[] h_Solutions;
		CHECK(cudaFree(d_Solutions));

		delete[] h_cost_effect;
		CHECK(cudaFree(d_cost_effect));

		delete[] h_reqList;
		CHECK(cudaFree(d_reqList));

		delete[] h_size;
		CHECK(cudaFree(d_size));
	}

	void execute_in_host(solution *Solution, myList *AvailableRequests, coefficient *Coefficient)
	{

		// 1) Get availableRequests info

			for (int s = 0; s < ExpectedScenarios; s++)
			{
				h_size[s] = AvailableRequests[s].size;
				for (int i = 0; i < h_size[s]; i++)
					h_reqList[i + (n * s)] = AvailableRequests[s].list[i];
			}

			//printf("--Ckpoint 1\n--");

		// 2) CudaMemCpy

			CHECK(cudaMemcpy(d_reqList, h_reqList, n * ExpectedScenarios * sizeof(int), cudaMemcpyHostToDevice));
			CHECK(cudaMemcpy(d_size, h_size, ExpectedScenarios * sizeof(int), cudaMemcpyHostToDevice));

			//printf("--Ckpoint 2\n--");

		// 3) Launch kernels

			// 3.1) update cost effect
			for (int s = 0; s < ExpectedScenarios; s++)
				for (int idx = 0; idx < h_size[s]; idx++)
				{
					h_Solutions[idx + (n*s)].Copy(Solution[s]);

					int req = h_reqList[idx + (n*s)];
					int k = h_Solutions[idx + (n*s)].Request[req].currvehicle;
					h_Solutions[idx + (n*s)].RemoveRequest(req, k);

					h_Solutions[idx + (n*s)].OfflineEvaluateCost();

					h_cost_effect[idx + (n*s)] = Coefficient[s].cost_function(h_Solutions[idx + (n*s)].Cost);
				}


			//printf("--Ckpoint 3\n--");

			//  3.2) sort cost effect
			for (int s = 0; s < ExpectedScenarios; s++)
			{
				float temp_cost;
				int temp_req;

				for (int i = 0; i < h_size[s]; i++)
					for (int j = i + 1; j < h_size[s]; j++)
						if (h_cost_effect[i + (n*s)] > h_cost_effect[j + (n*s)])
						{
							temp_cost = h_cost_effect[i + (n*s)];
							h_cost_effect[i + (n*s)] = h_cost_effect[j + (n*s)];
							h_cost_effect[j + (n*s)] = temp_cost;

							temp_req = h_reqList[i + (n*s)];
							h_reqList[i + (n*s)] = h_reqList[j + (n*s)];
							h_reqList[j + (n*s)] = temp_req;
						}

			}

			//printf("--Ckpoint 4\n--");

		// 4 copy from reqList to D_List


	}

	void execute(solution *Solution, myList *AvailableRequests, coefficient *Coefficient)
	{

		// 1) Get availableRequests info
			clock_t Get_AvailReq_Time = clock();
			for (int s = 0; s < ExpectedScenarios; s++)
			{
				h_size[s] = AvailableRequests[s].size;
				for (int i = 0; i < h_size[s]; i++)
					h_reqList[i + (n * s)] = AvailableRequests[s].list[i];
			}
			//printf("Get_AvailReq_Time: %0.02f\n\n", (float)(clock() - Get_AvailReq_Time) / CLOCKS_PER_SEC * 1000);

			//printf("--Ckpoint 1\n--");

		// 2) CudaMemCpy

			clock_t MemCpyTime = clock();

			CHECK(cudaMemcpyAsync(d_reqList, h_reqList, n * ExpectedScenarios * sizeof(int), cudaMemcpyHostToDevice));
			CHECK(cudaMemcpyAsync(d_size, h_size, ExpectedScenarios * sizeof(int), cudaMemcpyHostToDevice));
			CHECK(cudaMemcpy(Dev_Solution, Solution, ExpectedScenarios * sizeof(solution), cudaMemcpyHostToDevice));
			//printf("MemCpyTime: %0.02f\n\n", (float)(clock() - MemCpyTime) / CLOCKS_PER_SEC * 1000);

			//printf("--Ckpoint 2\n--");

		// 3) Launch kernels

			clock_t update_cost_effect_time = clock();

			update_cost_effect_Kernel KERNEL_ARGS2(dim3(ExpectedScenarios), dim3(n))(Dev_Solution, d_Solutions, d_cost_effect, d_reqList, d_size);
			CHECK(cudaDeviceSynchronize());
			//printf("update_cost_effect_time: %0.02f\n\n", (float)(clock() - update_cost_effect_time) / CLOCKS_PER_SEC * 1000);


			// 3.1) update cost effect
			/*for (int s = 0; s < ExpectedScenarios; s++)
				for (int idx = 0; idx < h_size[s]; idx++)
				{
					h_Solutions[idx + (n*s)].Copy(Solution[s]);

					int req = h_reqList[idx + (n*s)];
					int k = h_Solutions[idx + (n*s)].Request[req].currvehicle;
					h_Solutions[idx + (n*s)].RemoveRequest(req, k);

					h_Solutions[idx + (n*s)].OfflineEvaluateCost();
					h_cost_effect[idx + (n*s)] = h_Solutions[idx + (n*s)].Cost.overall_ideal_cost;
				}
			printf("update_cost_effect_time: %0.02f\n\n", (float)(clock() - update_cost_effect_time) / CLOCKS_PER_SEC * 1000);
			 */

			//printf("--Ckpoint 3\n--");


			clock_t sort_cost_effect_kernelTime = clock();

			sort_cost_effect_Kernel KERNEL_ARGS4(dim3(ExpectedScenarios), dim3(1), ExpectedScenarios * (sizeof(int)+sizeof(float)), 0)(d_cost_effect, d_reqList, d_size);
			CHECK(cudaDeviceSynchronize());
			CHECK(cudaMemcpy(h_reqList, d_reqList, n * ExpectedScenarios * sizeof(int), cudaMemcpyDeviceToHost));
			//printf("sort_cost_effect_kernelTime: %0.02f\n\n", (float)(clock() - sort_cost_effect_kernelTime) / CLOCKS_PER_SEC * 1000);



			//  3.2) sort cost effect
/*
 			CHECK(cudaMemcpy(h_cost_effect, d_cost_effect, n * ExpectedScenarios * sizeof(float), cudaMemcpyDeviceToHost));
			for (int s = 0; s < ExpectedScenarios; s++)
			{
				float temp_cost;
				int temp_req;

				for (int i = 0; i < h_size[s]; i++)
					for (int j = i + 1; j < h_size[s]; j++)
						if (h_cost_effect[i + (n*s)] > h_cost_effect[j + (n*s)])
						{
							temp_cost = h_cost_effect[i + (n*s)];
							h_cost_effect[i + (n*s)] = h_cost_effect[j + (n*s)];
							h_cost_effect[j + (n*s)] = temp_cost;

							temp_req = h_reqList[i + (n*s)];
							h_reqList[i + (n*s)] = h_reqList[j + (n*s)];
							h_reqList[j + (n*s)] = temp_req;
						}

			}
			printf("sort_cost_effect_kernelTime: %0.02f\n\n", (float)(clock() - sort_cost_effect_kernelTime) / CLOCKS_PER_SEC * 1000);
*/


	}

};

struct removal
{

	int total = 3;
	int intensity;
	int target_req;

	float determinism_p = 6; // suggested by Ropke (2006)
	float determinism_p_worst = 3;

	// shaw, random, and worst removal
	operation *Operator;

	worstcost *WorstCost;
	relatedness *Related;

	float *Relatedness;

	int *req_list;
	int *disable_list; // 1-disable, 0-enable
	int *D_list; // resultant list - consist of requests in order which they need to be removed.
	int D_list_size;

	removal()
	{
		D_list_size = 0;

		Operator = new operation[total];

		WorstCost = new worstcost;
		Related = new relatedness;

		Relatedness = new float[n];

		disable_list = new int[n];
		D_list = new int[n];
		req_list = new int[n];

		// set initial probBounds
		set_probBounds();

	}

	~removal()
	{
		delete[] Operator;

		delete WorstCost;
		delete Related;

		delete[] Relatedness;

		delete[] disable_list;
		delete[] D_list;
		delete[] req_list;

	}

	// Funcs related to shaw removal

	void new_shaw_removal(solution *TempSolution)
	{
		D_list_size = 0;

		// 1) pick first req randomly // OR CAN INITIALIZE THAT AS THE WORST-COST REQ AND FIND RELATED OTHER REQUESTS

		vector<int> vec;
		for (int i = 0; i < n; i++)
			vec.push_back(i);

		//printf("intensity:	%d\n", intensity);

		for (int i = 0; i < intensity; i++)
		{
			if (i == 0)
			{
				D_list[0] = rand() % n;
				D_list_size++;

				vec.erase(remove(vec.begin(), vec.end(), D_list[i]), vec.end());
			}
			else
			{
				Related[0].construct(n - D_list_size);

				int index = rand() % D_list_size;
				int req = D_list[index];

				Related[0].setting(vec);
				Related[0].update_relatedness(req, TempSolution);
				Related[0].sort_relatedness();

				int choice_of_index = 0;

				int avoid_determinism = 1;
				if (avoid_determinism)
				{
					float y = RandomPick(0, 1);

					float product_1 = y;
					for (int times = 0; times < (determinism_p - 1); times++)
						product_1 *= y;

					float product_2 = (float)(D_list_size - 1);

					float product = product_1 * product_2;

					int integer_product = (int)(ceil)(product);

					choice_of_index = integer_product;

					/*printf("y = %f,	product_1 = %f,	product_2 = %f,	product = %f, integer_product = %d\n", y, product_1, product_2, product, integer_product);
					printf("	choice_of_index = %d\n", choice_of_index);*/

				}

				D_list[D_list_size] = Related[0].req_list[choice_of_index];

				vec.erase(remove(vec.begin(), vec.end(), D_list[D_list_size]), vec.end());
				D_list_size++;

				Related[0].destruct();


			}
		}

		//printf("-------------------------------\n");
	}

	void update_relatedness(int req, solution *TempSolution)
	{
		float *alpha = new float[n];
		float *beta = new float[n];
		float *gamma = new float[n];;

		float sum_alpha = 0;
		float sum_beta = 0;
		float sum_gamma = 0;

		// 1) to avoid req already in D_list. by simply disabling them

		/*for (int i = 0; i < n; i++)
		{
		if (i == req)
		disable_list[i] = 1;
		else
		{
		for (int j = 0; j < intensity; j++) // for req already in D_list.
		if (i == D_list[j])
		disable_list[i] = 1;
		}
		}*/

		for (int i = 0; i < D_list_size; i++)
			disable_list[D_list[i]] = 1;



		// 2) find relatedness parameters

		for (int i = 0; i < n; i++)
		{
			if (disable_list[i] == 0)
			{
				alpha[i] = ABS(TempSolution[0].Problem.TravelTime[i] - TempSolution[0].Problem.TravelTime[req]) + ABS(TempSolution[0].Problem.TravelTime[i + n] - TempSolution[0].Problem.TravelTime[req + n]);
				beta[i] = ABS(TempSolution[0].Request[i].B - TempSolution[0].Request[req].B) + ABS(TempSolution[0].Request[i + n].B - TempSolution[0].Request[req + n].B);
				gamma[i] = 0; // zero for homogeneous set of requests i.e. demand at each node is always 1 in cordeau instances.

				sum_alpha += alpha[i];
				sum_beta += beta[i];
				sum_gamma += gamma[i];
			}
		}

		// 3) find final_relatedness value

		for (int i = 0; i < n; i++)
		{
			if (disable_list[i] == 0)
			{
				alpha[i] = (alpha[i] / sum_alpha);
				beta[i] = (beta[i] / sum_beta);

				Relatedness[i] = 2 * (alpha[i] + beta[i]) + gamma[i];
			}
			else
			{
				Relatedness[i] = 10 ^ 4; // set high values to avoid disabled requests
			}

		}


		delete[] alpha;
		delete[] beta;
		delete[] gamma;
	};

	void sort_relatedness(int *req_list)
	{
		// ascending sort: lower value of "relatedness" corresponds to higher relation b/w two requests

		float temp_1;
		int temp_2;

		for (int i = 0; i<n; i++)
			for (int j = i + 1; j<n; j++)
				if (Relatedness[i] > Relatedness[j])
				{
					temp_1 = Relatedness[i];
					Relatedness[i] = Relatedness[j];
					Relatedness[j] = temp_1;

					temp_2 = req_list[i];
					req_list[i] = req_list[j];
					req_list[j] = temp_2;
				}


	}

	// Funcs related to random removal

	void random_removal(solution *TempSolution, myList *AvailableRequests)
	{
		D_list_size = 0;

		AvailableRequests[0].randomize();
		for (int i = 0; i < intensity; i++)
		{
			D_list[i] = AvailableRequests[0].list[i];
			D_list_size++;
		}
	}

	// Funcs related to worst removal

	void new_worst_removal(solution *TempSolution)
	{
		D_list_size = 0;

		solution *Entering_Solution = new solution;

		// 1) pick first req randomly // OR CAN INITIALIZE THAT AS THE WORST-COST REQ AND FIND RELATED OTHER REQUESTS

		vector<int> vec;
		for (int i = 0; i < n; i++)
			vec.push_back(i);


		//printf("------------ ENTRY --------------\n");

		// Initialization:
		Entering_Solution[0].Copy(TempSolution[0]);
		Entering_Solution[0].OfflineEvaluateCost();

		printf("----here 1\n");

		for (int i = 0; i < intensity; i++)
		{
			printf("----bef construct\n");

			WorstCost[0].construct(n - D_list_size);

			printf("----aft construct\n");

			WorstCost[0].setting(vec);
			WorstCost[0].update_costeffect(Entering_Solution);
			WorstCost[0].sort_costeffect();

			int choice_of_index = 0;

			int avoid_determinism = 0;
			if (avoid_determinism)
			{
				float y = RandomPick(0, 1);

				float product_1 = y;
				for (int times = 0; times < (determinism_p_worst - 1); times++)
					product_1 *= y;

				float product_2 = (float)(D_list_size - 1);

				float product = product_1 * product_2;

				int integer_product = (int)(ceil)(product);

				choice_of_index = integer_product;

				/*printf("y = %f,	product_1 = %f,	product_2 = %f,	product = %f, integer_product = %d\n", y, product_1, product_2, product, integer_product);
				printf("	choice_of_index = %d\n", choice_of_index);*/

			}

			D_list[D_list_size] = WorstCost[0].req_list[choice_of_index];

			// Immediate removal: as per Ropke (2006)
			int removable_req = D_list[D_list_size];
			Entering_Solution[0].RemoveRequest(removable_req, Entering_Solution[0].Request[removable_req].currvehicle);
			Entering_Solution[0].OfflineEvaluateCost();

			vec.erase(remove(vec.begin(), vec.end(), D_list[D_list_size]), vec.end());
			D_list_size++;

			printf("----bef destruct\n");

			WorstCost[0].destruct();

			printf("----aft destruct\n");

		}

		printf("-------------------------------\n");

		delete Entering_Solution;
	}

	void worst_removal(solution *TempSolution, myList *AvailableRequests)
	{
		D_list_size = 0;

		sort_worstcost(TempSolution, AvailableRequests);

		// avoid determinism: not done here

		for (int i = 0; i < n; i++)
		{
			D_list[i] = AvailableRequests[0].list[i];
			D_list_size++;
		}
	}

	void sort_worstcost(solution *TempSolution, myList *AvailableRequests)
	{
		// initialization

		solution *NewSolution = new solution;
		simple_solution *SimpleSolution = new simple_solution;
		coefficient *Coefficient = new coefficient;
		position *Position = new position;

		Coefficient[0].reset();

		// 1) Sort route_order descending in accordance with worst costfunction

		NewSolution[0].Copy(TempSolution[0]);
		NewSolution[0].OfflineEvaluateCost();
		SimpleSolution[0].Copy(NewSolution[0]);
		SimpleSolution[0].total_cost = 0;
		for (int k = 0; k < m; k++)
		{
			SimpleSolution[0].route_order[k] = k;
			SimpleSolution[0].costfunction[k] = Coefficient[0].cost_function(NewSolution[0].Vehicle[k].Cost);
			SimpleSolution[0].total_cost += SimpleSolution[0].costfunction[k];
		}
		Sorting(SimpleSolution[0].route_order, SimpleSolution[0].costfunction);


		for (int k = 0; k < m; k++)
		{
			for (int j = 0; j < NewSolution[0].Vehicle[k].size; j++)
			{
				int node = NewSolution[0].Vehicle[k].path[j] - 1;
				if (node < n)
				{
					NewSolution[0].LocateRequest(Position, node);
					NewSolution[0].RemoveRequest(node, k);

					//NewSolution[0].OfflineEvaluateByRoute(k+1);
					NewSolution[0].OfflineEvaluateCost();
					SimpleSolution[0].request[node].req = node;
					SimpleSolution[0].request[node].vehicle = k;
					SimpleSolution[0].request[node].cost_effect = Coefficient[0].cost_function(NewSolution[0].Cost);
					//printf("no problem!");
					NewSolution[0].ReInsertRequest(node, k, Position[0].start, Position[0].start + Position[0].gap);
				}
			}
		}

		// 2) For every route, sort requests descending in accordance with worst costfunction

		Route_Sorting(SimpleSolution);

		// 3) Update all requests descending into the

		Req_Sorting(SimpleSolution, AvailableRequests);

		delete Position;
		delete Coefficient;
		delete SimpleSolution;
		delete NewSolution;
	}

	void Sorting(int *route, float *costfunction)
	{
		float temp_route;
		float temp_costfunction;
		for (int i = 0; i<m; i++)
		{
			for (int j = i + 1; j<m; j++)
			{
				if (costfunction[i] > costfunction[j])
				{
					temp_costfunction = costfunction[i];
					costfunction[i] = costfunction[j];
					costfunction[j] = temp_costfunction;

					temp_route = route[i];
					route[i] = route[j];
					route[j] = temp_route;
				}
			}
		}

	}

	void Req_Sorting(simple_solution *SimpleSolution, myList *AvailableRequests)
	{
		int size = AvailableRequests[0].size;

		float *cost_effect = new float[size];

		for (int i = 0; i < size; i++)
		{
			int req = AvailableRequests[0].list[i];
			cost_effect[req] = SimpleSolution[0].request[req].cost_effect;
		}

		int temp_req_list;
		float temp_costeffect_list;
		for (int i=0; i < size; i++)
		{
			for (int j=i+1; j < size; j++)
			{
				if (cost_effect[i] > cost_effect[j])
				{
					temp_req_list = AvailableRequests[0].list[i];
					AvailableRequests[0].list[i] = AvailableRequests[0].list[j];
					AvailableRequests[0].list[j] = temp_req_list;

					temp_costeffect_list = cost_effect[i];
					cost_effect[i] = cost_effect[j];
					cost_effect[j] = temp_costeffect_list;

				}
			}
		}

		delete cost_effect;
	}

	void Route_Sorting(simple_solution *SimpleSolution)
	{

		for (int k = 0; k < m; k++)
		{
			int temp_req;
			int temp_cost;

			for (int i = 0; i < SimpleSolution[0].vehicle[k].size; i++)
			{
				for (int j = i + 1; j < SimpleSolution[0].vehicle[k].size; j++)
				{
					if (SimpleSolution[0].request[SimpleSolution[0].vehicle[k].route[i]].cost_effect > SimpleSolution[0].request[SimpleSolution[0].vehicle[k].route[j]].cost_effect)
					{
						temp_req = SimpleSolution[0].vehicle[k].route[i];
						SimpleSolution[0].vehicle[k].route[i] = SimpleSolution[0].vehicle[k].route[j];
						SimpleSolution[0].vehicle[k].route[j] = temp_req;
					}
				}
			}
		}

	}

	// Funcs - for utility

	void randomRequest(int& req_list)
	{
		int* list = &req_list;
		bool *isInList = new bool[n];

		for (int i = 0; i < n; i++)
			isInList[i] = 0;

		for (int i = 0; i < n; i++)
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

	void reset(int *req_list)
	{
		for (int i = 0; i < n; i++)
			req_list[i] = i;
	}

	float ABS(float value)
	{
		if (value < 0)
			value *= -1;
		else
			value *= 1;

		return value;
	}

	// Main func - for operator selection

	void operator_selection(int select, int original_intensity, solution *TempSolution, myList *AvailableRequests)
	{
		// 0 - shaw, 1 - random, 2 - worst
		intensity = original_intensity;

		int shaw = 0, random = 1, worst = 2;

		switch (select)
		{
			case 0:
			{
				new_shaw_removal(TempSolution); // later, update the score of each operator.
				Operator[shaw].times_used++;
				Operator[shaw].all_times_used++;
				break;
			}
			case 1:
			{
				random_removal(TempSolution, AvailableRequests);
				Operator[random].times_used++;
				Operator[random].all_times_used++;
				break;
			}
			case 2:
			{
				worst_removal(TempSolution, AvailableRequests);
				Operator[worst].times_used++;
				Operator[worst].all_times_used++;
				break;
			}
			default:
			{
				printf("SOMETHING is WRONG at line 5300's.\n");
				break;
			}
		}

	}

	// bound setting
	void set_probBounds()
	{
		// call this at end of every segment

		float *rho = new float[total];

		// pre-calculations
		float total_weight = 0;
		for (int i = 0; i < total; i++)
			total_weight += Operator[i].weight;

		for (int i = 0; i < total; i++)
		{
			rho[i] = Operator[i].weight / total_weight;
		}



		float lastUpperBound = 0;

		// set probBounds
		for (int i = 0; i < total; i++)
		{
			Operator[i].ProbBound[0].lowerbound = lastUpperBound;
			Operator[i].ProbBound[0].upperbound = lastUpperBound + rho[i];
			lastUpperBound = Operator[i].ProbBound[0].upperbound;
		}

		delete[] rho;
	}

	void adjust_weights(float r_factor)
	{
		// call this at end of every segment

		// 1) adjust weights
		for (int i = 0; i < total; i++)
		{
			if (Operator[i].times_used < 0.1)
				Operator[i].weight = Operator[i].weight * (1 - r_factor) + (0);
			else
				Operator[i].weight = Operator[i].weight * (1 - r_factor) + (r_factor * Operator[i].score / Operator[i].times_used);

		}

		// 2) reset score and times used
		for (int i = 0; i < total; i++)
		{
			Operator[i].all_score += Operator[i].score;

			Operator[i].score = 0;
			Operator[i].times_used = 0;
		}

	}
};

struct reinsertion
{
	operation *Operator;
	int total;

	int noise_application;

	int intensity;
	int grade; // vary between 0 to (m-1)
	//grade: // 0 - diff between min_best and second best, 1 - diff between min_best and third best, and so on for (m-1) times...

	int *RI_list;  // dummy list for processing
	int *sequence; // resultant list - consist of order in which reqs from D_list will be reinserted.

				   // 1) greedy, 2) regret-k

	reinsertion()
	{
		total = min(m + 1, 6);//m+1;

		Operator = new operation[total];

		RI_list = new int[n];
		sequence = new int[n];

		// set initial probBounds
		set_probBounds();
	}

	~reinsertion()
	{

		delete[] sequence;
		delete[] RI_list;

		delete[] Operator;
	}

	void greedy_reinsertion(insertion *Insertion, noise *Noise)
	{
		// greedy: min_best_cost_insertion

		for (int p = 0; p < intensity; p++)
			sequence[p] = p;



		float *bestcost_list;
		bestcost_list = new float[intensity];

		//printf("bestcost_list\n");
		for (int p = 0; p < intensity; p++)
		{
			bestcost_list[p] = Insertion[p].local_best[0];
			//printf("	%d	bestcost_list[%d] = %f\n", p, p, bestcost_list[p]);
		}
		//printf("continuing..\n");


		// noise-application
		if (noise_application)
		{
			for (int p = 0; p < intensity; p++)
			{
				float spin = RandomPick(Noise[0].min_N, Noise[0].max_N);
				//if (spin < 0.5)
				bestcost_list[p] = max(0, bestcost_list[p] + spin);

			}
		}

		// sort ascending order
		float temp_cost;
		int temp_sequence;
		for (int i = 0; i<intensity; i++)
		{
			for (int j = i + 1; j<intensity; j++)
			{
				if (bestcost_list[i] > bestcost_list[j]) // SORTED IN ascending ORDER: MINIMUM COST REQ WILL BE GREEDILY INSERTED FIRST....
				{
					temp_cost = bestcost_list[i];
					bestcost_list[i] = bestcost_list[j];
					bestcost_list[j] = temp_cost;

					temp_sequence = sequence[i];
					sequence[i] = sequence[j];
					sequence[j] = temp_sequence;
				}
			}
		}

		delete[] bestcost_list;

		/*for (int p = 0; p < intensity; p++)
		printf("	sequence[%d] = %d, cost = %f\n", p, sequence[p], bestcost_list[p]);

		printf("end");
		printf("----------\n\n");*/

	}

	void regret_reinsertion(removal *Removal, insertion *Insertion, int grade, noise *Noise)
	{
		// regret-k generic insertion

		if (grade > (m - 1))
		{
			printf("ERROR: grade exceeds limit. so reset back to upper limit.\n");
			grade = m - 1;
		}

		for (int p = 0; p < intensity; p++)
			sequence[p] = p;

		float *costdifference_list;
		costdifference_list = new float[intensity];

		//printf("costdifference_list\n");
		for (int p = 0; p < intensity; p++)
		{
			costdifference_list[p] = Insertion[p].cost_difference[grade];
			//printf("	%d	costdifference_list[%d] = %f\n", p, p, costdifference_list[p]);
		}
		//printf("continuing..\n");


		// noise-application
		if (noise_application)
		{
			for (int p = 0; p < intensity; p++)
			{
				float spin = RandomPick(Noise[0].min_N, Noise[0].max_N);
				//if (spin < 0.5)
				costdifference_list[p] = max(0, costdifference_list[p] + spin);
			}
		}

		// sort ascending order
		float temp_cost;
		int temp_sequence;
		for (int i = 0; i<intensity; i++)
		{
			for (int j = i + 1; j<intensity; j++)
			{
				if (costdifference_list[i] < costdifference_list[j]) // SORTED IN descending ORDER: HIGHEST COST REGRET REQ WILL BE INSERTED FIRST....
				{
					temp_cost = costdifference_list[i];
					costdifference_list[i] = costdifference_list[j];
					costdifference_list[j] = temp_cost;

					temp_sequence = sequence[i];
					sequence[i] = sequence[j];
					sequence[j] = temp_sequence;
				}
			}
		}

		delete[] costdifference_list;

		/*for (int p = 0; p < intensity; p++)
		printf("	sequence[%d] = %d, cost = %f\n", p, sequence[p], costdifference_list[p]);

		printf("end");
		printf("----------\n\n");*/

	}

	void operator_selection(int pick, removal *Removal, insertion *Insertion, int inGrade, int in_Noise_application, noise *Noise, int in_Intensity)
	{
		noise_application = in_Noise_application;

		grade = inGrade;
		intensity = in_Intensity;

		for (int i = 0; i < Removal[0].intensity; i++)
			RI_list[i] = Removal[0].D_list[i];

		int greedy = 0;

		switch (pick)
		{
		case 0:
		{
			greedy_reinsertion(Insertion, Noise); // later, update the score of each operator.
			Operator[greedy].times_used++;
			Operator[greedy].all_times_used++;
			break;
		}
		case 1:
		{
			regret_reinsertion(Removal, Insertion, grade, Noise);
			Operator[1 + grade].times_used++;
			Operator[1 + grade].all_times_used++;
			break;
		}
		default:
		{
			printf("SOMETHING is WRONG at line 5700's.\n");
			break;
		}
		}

	}

	// bound setting
	void set_probBounds()
	{
		// call this at end of every segment

		float *rho = new float[total];

		// pre-calculations
		float total_weight = 0;
		for (int i = 0; i < total; i++)
			total_weight += Operator[i].weight;

		for (int i = 0; i < total; i++)
		{
			rho[i] = Operator[i].weight / total_weight;
		}



		float lastUpperBound = 0;

		// set probBounds
		for (int i = 0; i < total; i++)
		{
			Operator[i].ProbBound[0].lowerbound = lastUpperBound;
			Operator[i].ProbBound[0].upperbound = lastUpperBound + rho[i];
			lastUpperBound = Operator[i].ProbBound[0].upperbound;
		}

		delete[] rho;
	}

	void adjust_weights(float r_factor)
	{
		// call this at end of every segment

		// 1) adjust weights
		for (int i = 0; i < total; i++)
		{
			if (Operator[i].times_used < 0.1)
				Operator[i].weight = Operator[i].weight * (1 - r_factor) + (0);
			else
				Operator[i].weight = Operator[i].weight * (1 - r_factor) + (r_factor * Operator[i].score / Operator[i].times_used);

		}

		// 2) reset score and times used
		for (int i = 0; i < total; i++)
		{
			Operator[i].all_score += Operator[i].score;

			Operator[i].score = 0;
			Operator[i].times_used = 0;
		}

	}

};
