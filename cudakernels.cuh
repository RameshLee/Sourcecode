__global__
void IH_parallelSETUP(sol *gpu, int reqID, solution *Solution)
{
	int indices = 0, start = 0;
	for (int k = 0; k < TotalVehicles; k++)
	{
		int size = Solution[0].Vehicle[k].size;

		for (int limiter = size; limiter >= 0; limiter--)
		{
			for (start = 0; start <= limiter; start++)
			{
				gpu[indices].start = start;
				gpu[indices].gap = size - limiter;
				gpu[indices].vehid = k;
				indices++;

				//if (indices >= Solution[0].allinsertion)
				//break;
			}

			//if (indices >= Solution[0].allinsertion)
			//break;
		}

		if (indices >= Solution[0].allinsertion)
			break;
	}
}

__global__
void IH_parallelSETUP_NEW(sol *gpu, container *gpu_Container, int reqID, solution *Solution)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < Solution[0].allinsertion)
	{
		gpu[idx].start = gpu_Container[idx].start;
		gpu[idx].gap = gpu_Container[idx].gap;
		gpu[idx].vehid = gpu_Container[idx].vehid;
	}
}

__global__
void IH_parallelINSERTION(sol *gpu, int reqID, solution *Solution)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < Solution[0].allinsertion)
	{
		gpu[idx].request = reqID - 1;
		gpu[idx].fixation(Solution[0].Vehicle[gpu[idx].vehid].size); //here, we update routesize.
		gpu[idx].route[gpu[idx].start] = reqID;
		gpu[idx].route[gpu[idx].start + gpu[idx].gap + 1] = reqID + n;


		/*if (gpu[idx].routesize == 2) {}

		else
			*/for (int tid = 0; tid < gpu[idx].routesize - 2; tid++)
			{
				bool boolone = tid >= (gpu[idx].start);
				bool booltwo = tid >= (gpu[idx].start + gpu[idx].gap);
				gpu[idx].route[tid + boolone + booltwo] = Solution[0].Vehicle[gpu[idx].vehid].path[tid];
			}
	}

}

__global__
void IH_parallelEVALUATION(sol *gpu, int reqID, solution *Solution, problem *Problem)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < Solution[0].allinsertion)
		gpu[idx].evaluateByRoute(Problem);

}

__global__
void IH_parallelUPDATE(sol *gpu, int reqID, solution *Solution, cost *Temp, cost *Full_Cost_Per_Route)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < Solution[0].allinsertion)
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



		Temp[idx].vehicle = gpu[idx].vehid;

		Temp[idx].isFeasible = Temp[idx].getFeasibility();
		gpu[idx].isFeasible = Temp[idx].getFeasibility();


		//vehicle breakdown code

		Temp[idx].vehid = gpu[idx].vehid;

		// Extract start value for haltIndex
		Temp[idx].start = gpu[idx].start;
		Temp[idx].vehicle = gpu[idx].vehid;


	}

}

//////////////////////////////////////////////
/////// PICKUP-specialized kernels //////////////////
//////////////////////////////////////////////

__global__
void IH_parallelSETUP_AnyNode_First(sol *gpu, int reqID, solution *Solution)
{
	int indices = 0, start = 0;
	for (int k = 0; k < TotalVehicles; k++)
	{
		int size = Solution[0].Vehicle[k].size;

		for (start = 0; start <= size; start++)
		{
			gpu[indices].request = reqID - 1;
			gpu[indices].start = start;
			gpu[indices].gap = 0; /// NEEEEEED TO REMOVEEEE THISSSS ONEEE!!!!!!!!!
			gpu[indices].vehid = k;
			indices++;
		}

		if (indices >= Solution[0].allinsertion)
			break;
	}
}

__global__
void IH_parallelINSERTION_AnyNode_First(sol *gpu, int reqID, solution *Solution)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < Solution[0].allinsertion)
	{
		gpu[idx].routesize = Solution[0].Vehicle[ gpu[idx].vehid ].size + 1;
		gpu[idx].route[gpu[idx].start] = reqID;

		if (gpu[idx].routesize == 2) {}
		else
			for (int tid = 0; tid < gpu[idx].routesize - 1; tid++)
			{
				bool boolone = tid >= (gpu[idx].start);
				gpu[idx].route[tid + boolone] = Solution[0].Vehicle[gpu[idx].vehid].path[tid];
			}
	}

}

//////////////////////////////////////////////
//////// DROPOFF-specialized kernels ////////////////
//////////////////////////////////////////////

__global__
void IH_parallelSETUP_AnyNode_Second(sol *gpu, int reqID, int start, int Curr_k, solution *Solution)
{
	for (int idx = 0; idx < Solution[0].allinsertion; idx++)
	{
		gpu[idx].gap = idx;

		gpu[idx].request = reqID - 1;
		gpu[idx].start = start;
		gpu[idx].vehid = Curr_k;
	}

}

__global__
void IH_parallelINSERTION_DropOff_Second(sol *gpu, int reqID, solution *Solution)
{

		int idx = threadIdx.x + blockDim.x * blockIdx.x;

		if (idx < Solution[0].allinsertion)
		{
			gpu[idx].routesize = Solution[0].Vehicle[ gpu[idx].vehid ].size + 1;

			gpu[idx].route[gpu[idx].start + gpu[idx].gap + 1] = reqID + n; // insert dropoff

			if (gpu[idx].routesize == 2) {}
			else
			{
				for (int tid = 0; tid < gpu[idx].routesize - 1; tid++)
				{
					bool booltwo = tid >= (gpu[idx].start + gpu[idx].gap + 1); // bypass
					gpu[idx].route[tid + booltwo] = Solution[0].Vehicle[gpu[idx].vehid].path[tid];

				}
			}
		}


}

__global__
void IH_parallelINSERTION_Pickup_Second(sol *gpu, int reqID, solution *Solution)
{

		int idx = threadIdx.x + blockDim.x * blockIdx.x;

		if (idx < Solution[0].allinsertion)
		{
			gpu[idx].routesize = Solution[0].Vehicle[ gpu[idx].vehid ].size + 1;
			gpu[idx].route[idx] = reqID; // reversal-active: insert pickup.

			if (gpu[idx].routesize == 2) {}
			else
			{
				for (int tid = 0; tid < gpu[idx].routesize - 1; tid++)
				{
					bool booltwo = tid >= idx; // reversal-active: bypass gap
					gpu[idx].route[tid + booltwo] = Solution[0].Vehicle[gpu[idx].vehid].path[tid];

				}
			}
		}


}


//////////////////////////////////////////////



__global__
void IH_CUDA_Dummy1(sol *gpu, container *gpu_Container, int reqID, solution *Solution, problem *Problem, cost *Temp, cost *Full_Cost_Per_Route)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < Solution[0].allinsertion)
	{
		bool boolone, booltwo;

		// 1) SETUP

		gpu[idx].start = gpu_Container[idx].start;
		gpu[idx].gap = gpu_Container[idx].gap;
		gpu[idx].vehid = gpu_Container[idx].vehid;

		// 2) INSERTION

		gpu[idx].request = reqID - 1;
		gpu[idx].fixation(Solution[0].Vehicle[gpu_Container[idx].vehid].size); //here, we update routesize.
		gpu[idx].route[gpu_Container[idx].start] = reqID;
		gpu[idx].route[gpu_Container[idx].start + gpu_Container[idx].gap + 1] = reqID + n;


		for (int tid = 0; tid < gpu[idx].routesize - 2; tid++)
		{
			boolone = tid >= (gpu_Container[idx].start);
			booltwo = tid >= (gpu_Container[idx].start + gpu_Container[idx].gap);
			gpu[idx].route[tid + boolone + booltwo] = Solution[0].Vehicle[gpu_Container[idx].vehid].path[tid];
		}
	}
}

__global__
void IH_CUDA_Dummy2(sol *gpu, container *gpu_Container, int reqID, solution *Solution, problem *Problem, cost *Temp, cost *Full_Cost_Per_Route)
{

		int idx = threadIdx.x + blockDim.x * blockIdx.x;

		if (idx < Solution[0].allinsertion)
		{
			// 3) EVALUATION

			gpu[idx].evaluateByRoute(Problem);
		}
}

__global__
void IH_CUDA_Dummy3(sol *gpu, container *gpu_Container, int reqID, solution *Solution, problem *Problem, cost *Temp, cost *Full_Cost_Per_Route)
{
		// 4) UPDATE

		int idx = threadIdx.x + blockDim.x * blockIdx.x;

		bool boolone, booltwo;
		Temp[idx].reset();
		for (int k = 0; k < TotalVehicles; k++)
		{
			boolone = (k == gpu[idx].vehid);
			booltwo = 1 - boolone;

			Temp[idx].travel_cost += (boolone * gpu[idx].Vehicle.Cost.travel_cost) + (booltwo * Full_Cost_Per_Route[k].travel_cost);
			Temp[idx].time_window += (boolone * gpu[idx].Vehicle.Cost.time_window) + (booltwo * Full_Cost_Per_Route[k].time_window);
			Temp[idx].ride_time += (boolone * gpu[idx].Vehicle.Cost.ride_time) + (booltwo * Full_Cost_Per_Route[k].ride_time);
			Temp[idx].load += (boolone * gpu[idx].Vehicle.Cost.load) + (booltwo * Full_Cost_Per_Route[k].load);
			Temp[idx].duration += (boolone * gpu[idx].Vehicle.Cost.duration) + (booltwo * Full_Cost_Per_Route[k].duration);
			Temp[idx].excessrideTime += (boolone * gpu[idx].Vehicle.Cost.excessrideTime) + (booltwo * Full_Cost_Per_Route[k].excessrideTime);


			/*if (k == gpu[idx].vehid)
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
			}*/

		}



		Temp[idx].vehicle = gpu[idx].vehid;

		Temp[idx].isFeasible = Temp[idx].getFeasibility();
		gpu[idx].isFeasible = Temp[idx].getFeasibility();


		//vehicle breakdown code

		Temp[idx].vehid = gpu[idx].vehid;

		// Extract start value for haltIndex
		Temp[idx].start = gpu[idx].start;
		Temp[idx].vehicle = gpu[idx].vehid;


}




__global__
void IH_parallel_OneInitialization(int s, sol *gpu, solution *Solution, holder *Dev_Holder)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
//#pragma unroll
	for (int k = 0; k < m; k++)
	{
		if (idx < Dev_Holder[Solution[s].Vehicle[k].size].total_insertions)
		{
			gpu[idx + Solution[s].Vehicle[k].cumulativeInsertions].start = Dev_Holder[Solution[s].Vehicle[k].size].start[idx];
			gpu[idx + Solution[s].Vehicle[k].cumulativeInsertions].gap = Dev_Holder[Solution[s].Vehicle[k].size].gap[idx];
			gpu[idx + Solution[s].Vehicle[k].cumulativeInsertions].vehid = k;
		}
	}

}


__global__
void IH_parallel_Setup(int s, int reqID, sol *gpu, cost *Temp, container *gpu_Container, cost *Full_Cost_Per_Route, solution *Solution, problem *Problem)//, float *Dev_Cost_array_flat, float *Dev_Cost_array, coefficient *Dev_Coefficient)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx < Solution[s].allinsertion)
	{
		bool boolone, booltwo;

		// 1) SETUP

		gpu[idx].start = gpu_Container[idx].start;
		gpu[idx].gap = gpu_Container[idx].gap;
		gpu[idx].vehid = gpu_Container[idx].vehid;
	}
}

__global__
void IH_parallel_OneKernel(int s, int reqID, sol *gpu, cost *Temp, container *gpu_Container, cost *Full_Cost_Per_Route, solution *Solution, problem *Problem, float *Dev_Cost_array_flat, float *Dev_Cost_array, coefficient *Dev_Coefficient, skeleton *d_SkeletonData)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx < Solution[s].allinsertion)
	{
		bool boolone, booltwo;

		// 1) SETUP

		/*gpu[idx].start = gpu_Container[idx].start;
		gpu[idx].gap = gpu_Container[idx].gap;
		gpu[idx].vehid = gpu_Container[idx].vehid;*/

		// 2) INSERTION

		gpu[idx].request = reqID - 1;
		gpu[idx].fixation(Solution[s].Vehicle[gpu[idx].vehid].size); //here, we update routesize.
		gpu[idx].route[gpu[idx].start] = reqID;
		gpu[idx].route[gpu[idx].start + gpu[idx].gap + 1] = reqID + n;

		for (int tid = 0; tid < gpu[idx].routesize - 2; tid++)
		{
			boolone = tid >= (gpu[idx].start);
			booltwo = tid >= (gpu[idx].start + gpu[idx].gap);
			gpu[idx].route[tid + boolone + booltwo] = Solution[s].Vehicle[gpu[idx].vehid].path[tid];
		}

		// 3) EVALUATION

		gpu[idx].evaluateByRoute(&Problem[s]);

		// 4) UPDATE

		Temp[idx].reset();
//#pragma unroll
		for (int k = 0; k < TotalVehicles; k++)
		{
			boolone = (k == gpu[idx].vehid);
			booltwo = 1 - boolone;

			Temp[idx].travel_cost += (boolone * gpu[idx].Vehicle.Cost.travel_cost) + (booltwo * Full_Cost_Per_Route[k].travel_cost);
			Temp[idx].time_window += (boolone * gpu[idx].Vehicle.Cost.time_window) + (booltwo * Full_Cost_Per_Route[k].time_window);
			Temp[idx].ride_time += (boolone * gpu[idx].Vehicle.Cost.ride_time) + (booltwo * Full_Cost_Per_Route[k].ride_time);
			Temp[idx].load += (boolone * gpu[idx].Vehicle.Cost.load) + (booltwo * Full_Cost_Per_Route[k].load);
			Temp[idx].duration += (boolone * gpu[idx].Vehicle.Cost.duration) + (booltwo * Full_Cost_Per_Route[k].duration);
			Temp[idx].excessrideTime += (boolone * gpu[idx].Vehicle.Cost.excessrideTime) + (booltwo * Full_Cost_Per_Route[k].excessrideTime);


			Temp[idx].overall_ideal_cost += (boolone * gpu[idx].Vehicle.Cost.overall_ideal_cost) + (booltwo * Full_Cost_Per_Route[k].overall_ideal_cost);

			/*if (k == gpu[idx].vehid)
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
			}*/

		}

		Temp[idx].vehicle = gpu[idx].vehid;
		Temp[idx].isFeasible = Temp[idx].getFeasibility();
		gpu[idx].isFeasible = Temp[idx].getFeasibility();
		Temp[idx].vehid = gpu[idx].vehid;
		Temp[idx].start = gpu[idx].start;

		Temp[idx].overall_current_cost = Dev_Coefficient[s].cost_function(Temp[idx]);
		Dev_Cost_array[idx] = Temp[idx].overall_current_cost;//Dev_Coefficient[s].cost_function(Temp[idx]);
		Dev_Cost_array_flat[idx +(s * Expectation)] = Dev_Coefficient[s].cost_function(Temp[idx]);

		d_SkeletonData[idx].req = gpu[idx].request;
		d_SkeletonData[idx].vehid = gpu[idx].vehid;
		d_SkeletonData[idx].start = gpu[idx].start;
		d_SkeletonData[idx].gap = gpu[idx].gap;

	}

}

__global__
void IH_parallel_OneKernel_Ver2(int s, int reqID, sol *gpu, cost *Temp, container *gpu_Container, cost *Full_Cost_Per_Route, solution *Solution, problem *Problem)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < Solution[0].allinsertion)
	{
		bool boolone, booltwo;

		// 1) SETUP

		gpu[idx].start = gpu_Container[idx].start;
		gpu[idx].gap = gpu_Container[idx].gap;
		gpu[idx].vehid = gpu_Container[idx].vehid;

		// 2) INSERTION

		gpu[idx].request = reqID - 1;
		gpu[idx].fixation(Solution[0].Vehicle[gpu_Container[idx].vehid].size); //here, we update routesize.
		gpu[idx].route[gpu_Container[idx].start] = reqID;
		gpu[idx].route[gpu_Container[idx].start + gpu_Container[idx].gap + 1] = reqID + n;


		for (int tid = 0; tid < gpu[idx].routesize - 2; tid++)
		{
			boolone = tid >= (gpu_Container[idx].start);
			booltwo = tid >= (gpu_Container[idx].start + gpu_Container[idx].gap);
			gpu[idx].route[tid + boolone + booltwo] = Solution[0].Vehicle[gpu_Container[idx].vehid].path[tid];
		}

		// 3) EVALUATION

		//gpu[idx].evaluateByRoute(&Problem[s]);
		gpu[idx].evaluateByRoute(&Problem[0]);

		// 4) UPDATE

		Temp[idx +(s * Expectation)].reset();
		for (int k = 0; k < TotalVehicles; k++)
		{
			boolone = (k == gpu[idx].vehid);
			booltwo = 1 - boolone;

			Temp[idx +(s * Expectation)].travel_cost += (boolone * gpu[idx].Vehicle.Cost.travel_cost) + (booltwo * Full_Cost_Per_Route[k].travel_cost);
			Temp[idx +(s * Expectation)].time_window += (boolone * gpu[idx].Vehicle.Cost.time_window) + (booltwo * Full_Cost_Per_Route[k].time_window);
			Temp[idx +(s * Expectation)].ride_time += (boolone * gpu[idx].Vehicle.Cost.ride_time) + (booltwo * Full_Cost_Per_Route[k].ride_time);
			Temp[idx +(s * Expectation)].load += (boolone * gpu[idx].Vehicle.Cost.load) + (booltwo * Full_Cost_Per_Route[k].load);
			Temp[idx +(s * Expectation)].duration += (boolone * gpu[idx].Vehicle.Cost.duration) + (booltwo * Full_Cost_Per_Route[k].duration);
			Temp[idx +(s * Expectation)].excessrideTime += (boolone * gpu[idx].Vehicle.Cost.excessrideTime) + (booltwo * Full_Cost_Per_Route[k].excessrideTime);

			Temp[idx +(s * Expectation)].overall_ideal_cost += (boolone * gpu[idx].Vehicle.Cost.overall_ideal_cost) + (booltwo * Full_Cost_Per_Route[k].overall_ideal_cost);


			/*if (k == gpu[idx].vehid)
			{
				Temp[idx +(s * Expectation)].travel_cost += gpu[idx].Vehicle.Cost.travel_cost;
				Temp[idx +(s * Expectation)].time_window += gpu[idx].Vehicle.Cost.time_window;
				Temp[idx +(s * Expectation)].ride_time += gpu[idx].Vehicle.Cost.ride_time;
				Temp[idx +(s * Expectation)].load += gpu[idx].Vehicle.Cost.load;
				Temp[idx +(s * Expectation)].duration += gpu[idx].Vehicle.Cost.duration;
				Temp[idx +(s * Expectation)].excessrideTime += gpu[idx].Vehicle.Cost.excessrideTime;
			}
			else
			{
				Temp[idx +(s * Expectation)].travel_cost += Full_Cost_Per_Route[k].travel_cost;
				Temp[idx +(s * Expectation)].time_window += Full_Cost_Per_Route[k].time_window;
				Temp[idx +(s * Expectation)].ride_time += Full_Cost_Per_Route[k].ride_time;
				Temp[idx +(s * Expectation)].load += Full_Cost_Per_Route[k].load;
				Temp[idx +(s * Expectation)].duration += Full_Cost_Per_Route[k].duration;
				Temp[idx +(s * Expectation)].excessrideTime += Full_Cost_Per_Route[k].excessrideTime;
			}*/

		}

		Temp[idx +(s * Expectation)].vehicle = gpu[idx].vehid;
		Temp[idx +(s * Expectation)].isFeasible = Temp[idx +(s * Expectation)].getFeasibility();
		Temp[idx +(s * Expectation)].vehid = gpu[idx].vehid;
		Temp[idx +(s * Expectation)].start = gpu[idx].start;
		Temp[idx +(s * Expectation)].vehicle = gpu[idx].vehid;
		gpu[idx].isFeasible = Temp[idx +(s * Expectation)].getFeasibility();


	}

}

__global__
void IH_parallel_OneKernel_New(int s, int reqID, sol *gpu, cost *Temp, container *gpu_Container, cost *New_Full_Cost_Per_Route, solution *Solution, problem *Problem, float *Dev_Cost_array_flat, float *Dev_Cost_array, coefficient *Dev_Coefficient, skeleton *d_SkeletonData)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx < Solution[s].allinsertion)
	{
		bool boolone, booltwo;

		// 1) SETUP

		/*gpu[idx].start = gpu_Container[idx].start;
		gpu[idx].gap = gpu_Container[idx].gap;
		gpu[idx].vehid = gpu_Container[idx].vehid;*/

		// 2) INSERTION

		gpu[idx].request = reqID - 1;
		gpu[idx].fixation(Solution[s].Vehicle[gpu[idx].vehid].size); //here, we update routesize.
		gpu[idx].route[gpu[idx].start] = reqID;
		gpu[idx].route[gpu[idx].start + gpu[idx].gap + 1] = reqID + n;

		for (int tid = 0; tid < gpu[idx].routesize - 2; tid++)
		{
			boolone = tid >= (gpu[idx].start);
			booltwo = tid >= (gpu[idx].start + gpu[idx].gap);
			gpu[idx].route[tid + boolone + booltwo] = Solution[s].Vehicle[gpu[idx].vehid].path[tid];
		}

		// 3) EVALUATION

		//gpu[idx].evaluateByRoute_standardDARP(&Problem[s]);
		gpu[idx].evaluateByRoute(&Problem[s]);

		// 4) UPDATE

		Temp[idx].reset();
//#pragma unroll
		for (int k = 0; k < TotalVehicles; k++)
		{
			boolone = (k == gpu[idx].vehid);
			booltwo = 1 - boolone;


			Temp[idx].on_board_user_wait_time += (boolone * gpu[idx].Vehicle.Cost.on_board_user_wait_time) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].on_board_user_wait_time);
			Temp[idx].route_duration += (boolone * gpu[idx].Vehicle.Cost.route_duration) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].route_duration);
			Temp[idx].early_arrival_time += (boolone * gpu[idx].Vehicle.Cost.early_arrival_time) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].early_arrival_time);

			Temp[idx].travel_cost += (boolone * gpu[idx].Vehicle.Cost.travel_cost) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].travel_cost);
			Temp[idx].time_window += (boolone * gpu[idx].Vehicle.Cost.time_window) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].time_window);
			Temp[idx].ride_time += (boolone * gpu[idx].Vehicle.Cost.ride_time) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].ride_time);
			Temp[idx].load += (boolone * gpu[idx].Vehicle.Cost.load) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].load);
			Temp[idx].duration += (boolone * gpu[idx].Vehicle.Cost.duration) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].duration);
			Temp[idx].excessrideTime += (boolone * gpu[idx].Vehicle.Cost.excessrideTime) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].excessrideTime);


			Temp[idx].overall_ideal_cost += (boolone * gpu[idx].Vehicle.Cost.overall_ideal_cost) + (booltwo * New_Full_Cost_Per_Route[k +(s * m)].overall_ideal_cost);

		}

		Temp[idx].vehicle = gpu[idx].vehid;
		Temp[idx].isFeasible = Temp[idx].getFeasibility();
		gpu[idx].isFeasible = Temp[idx].getFeasibility();
		Temp[idx].vehid = gpu[idx].vehid;
		Temp[idx].start = gpu[idx].start;

		Temp[idx].overall_current_cost = Dev_Coefficient[s].cost_function(Temp[idx]);
		Dev_Cost_array[idx] = Dev_Coefficient[s].cost_function(Temp[idx]);
		Dev_Cost_array_flat[idx +(s * Expectation)] = Dev_Coefficient[s].cost_function(Temp[idx]);

		d_SkeletonData[idx].req = gpu[idx].request;
		d_SkeletonData[idx].vehid = gpu[idx].vehid;
		d_SkeletonData[idx].start = gpu[idx].start;
		d_SkeletonData[idx].gap = gpu[idx].gap;

	}

}


__global__
void DUMMYKERNEL(int p, cost *Temp, insertion *Insertion, coefficient *Coefficient, solution *Solution, int s)
{


	int found = 0;
	Insertion[p].reset();
	for (int idx = 0; idx < Solution[s].allinsertion; idx++)
	{
		float currentcost = Coefficient[s].cost_function(Temp[idx]);
		if (Insertion[p].local_best[Temp[idx].vehicle] == 0 || currentcost < Insertion[p].local_best[Temp[idx].vehicle])
		{
				Insertion[p].local_best[Temp[idx].vehicle] = currentcost;
				Insertion[p].local_best_ID[Temp[idx].vehicle] = idx;
				found = 1;
		}
	}
	Insertion[p].sort_localbest();

	int min_cost_vehicle_k = Insertion[p].veh_sequence[0];
	int idx = Insertion[p].local_best_ID[min_cost_vehicle_k];

	for (int index = 1; index < m; index++) //beware - upto (m-1) only
		Insertion[p].cost_difference[index - 1] = Insertion[p].local_best[index] - Insertion[p].local_best[0];

	// here copy idx somewhere


}



__global__ void parallel_GPU_reduction_Kernel_first_initialization(int s, solution *Solution, float *Dev_Cost_array)
{

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx < Expectation)
		Dev_Cost_array[idx] = 3.40282e+038; //max value of float

	/*if (idx >= Solution[s].allinsertion)
		if (idx < (((Solution[s].allinsertion+32-1)/32) * 32) )
			Dev_Cost_array[idx] = 3.40282e+038; //max value of float*/

}

__global__
void parallel_GPU_reduction_Kernel(int s, int *Dev_Idx, float *Dev_Cost_array_flat, cost *Temp, float *Dev_Cost_array, solution *Solution)
{
	//thread_group block = this_thread_block();

	////////////////////
	int local_ID = 0;
	float localbest = 0;
	for (int idx = 0; idx < Solution[0].allinsertion; idx++)
	{
		float overall_current_cost = Temp[idx].overall_current_cost;
		//float overall_current_cost = Dev_Cost_array[idx];
		//float overall_current_cost = Dev_Cost_array_flat[idx +(s * Expectation)];
		if (localbest == 0 || overall_current_cost < localbest)
		{
			localbest = overall_current_cost;
			local_ID = idx;
		}
	}
	Dev_Idx[0] = local_ID;
	////////////////////

	//block.sync();

}

__global__ void parallel_GPU_reduction_Kernel_v2(float *g_idata, float *g_odata, int *g_i_index, int *g_o_index, int size)
{

	unsigned int tid = threadIdx.x;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	float *idata = g_idata + blockIdx.x * blockDim.x;
	int *i_index = g_i_index + blockIdx.x * blockDim.x;

	if (idx >= size) return;

	// in-place reduction in global memory
	for (int stride = blockDim.x / 2; stride > 0; stride >>= 1)
	{
		if (tid < stride)
		{
			idata[tid] = min(idata[tid], idata[tid + stride]);

			__syncthreads();

			if (idata[tid] == idata[tid + stride])
				i_index[tid] = i_index[tid + stride];

			bool Bool = (idata[tid] == idata[tid + stride]);
			i_index[tid] = Bool * i_index[tid + stride] + (1-Bool) * i_index[tid];

			__syncthreads();
		}

	}

	if (tid == 0)
	{
		g_o_index[blockIdx.x] = i_index[0];
		g_odata[blockIdx.x] = idata[0];
	}



}

__global__ void parallel_GPU_reduction_Kernel_v3(float *g_idata, float *g_odata, int *g_i_index, int *g_o_index, int size)
{

	unsigned int tid = threadIdx.x;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	float *idata = g_idata + blockIdx.x * blockDim.x;
	int *i_index = g_i_index + blockIdx.x * blockDim.x;

	if (idx >= size) return;

	for (int stride = 1; stride < blockDim.x; stride *= 2)
	{
		if ((tid % (2 * stride)) == 0)
		{
			idata[tid] = min(idata[tid], idata[tid + stride]);

			__syncthreads();

			bool Bool = (idata[tid] == idata[tid + stride]);
			i_index[tid] = Bool * i_index[tid + stride] + (1-Bool) * i_index[tid];

			__syncthreads();


			/*bool Bool = (idata[tid] > idata[tid + stride]);
			idata[tid] = Bool * idata[tid + stride] + (1-Bool) * idata[tid];
			__syncthreads();

			bool Bool2 = (idata[tid] == idata[tid + stride]);
			i_index[tid] = Bool2 * i_index[tid + stride] + (1-Bool2) * i_index[tid];
			__syncthreads();
   */


			idata[tid] = min(idata[tid], idata[tid + stride]);
			__syncthreads();

			if (idata[tid] == idata[tid + stride])
				i_index[tid] =i_index[tid + stride];
			else
				i_index[tid] += 0;
			__syncthreads();


		}
	}


	if (tid == 0)
	{
		g_o_index[blockIdx.x] = i_index[0];
		g_odata[blockIdx.x] = idata[0];
	}
}

__global__ void parallel_GPU_reduction_Kernel_v4(float *g_idata, float *g_odata, int *g_i_index, int *g_o_index, int size)
{

	// set thread ID
	unsigned int tid = threadIdx.x;
	float *idata = g_idata + blockIdx.x * blockDim.x;
	int *i_index = g_i_index + blockIdx.x * blockDim.x;
	// boundary check
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= size) return;

	// in-place reduction in global memory
	if (blockDim.x >= 1024 && tid < 512)
	{
		idata[tid] = min(idata[tid], idata[tid + 512]);
		__syncthreads();

		bool Bool = (idata[tid] == idata[tid + 512]);
		i_index[tid] = Bool * i_index[tid + 512] + (1-Bool) * i_index[tid];
		__syncthreads();
	}

	if (blockDim.x >= 512 && tid < 256)
	{
		idata[tid] = min(idata[tid], idata[tid + 256]);
		__syncthreads();

		bool Bool = (idata[tid] == idata[tid + 256]);
		i_index[tid] = Bool * i_index[tid + 256] + (1-Bool) * i_index[tid];
		__syncthreads();
	}

	if (blockDim.x >= 256 && tid < 128)
	{
		idata[tid] = min(idata[tid], idata[tid + 128]);
		__syncthreads();

		bool Bool = (idata[tid] == idata[tid + 128]);
		i_index[tid] = Bool * i_index[tid + 128] + (1-Bool) * i_index[tid];
		__syncthreads();
	}

	if (blockDim.x >= 128 && tid < 64)
	{
		idata[tid] = min(idata[tid], idata[tid + 64]);
		__syncthreads();

		bool Bool = (idata[tid] == idata[tid + 64]);
		i_index[tid] = Bool * i_index[tid + 64] + (1-Bool) * i_index[tid];
		__syncthreads();
	}

	// unrolling warp
	if (tid < 32)
	{
		volatile float *vsmem = idata;
		volatile int *vsmem_index = i_index;

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 32]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 32]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 16]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 16]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 8]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 8]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 4]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 4]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 2]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 2]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 1]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 1]);
		__syncthreads();
	}

	if (tid == 0)
	{
		g_o_index[blockIdx.x] = i_index[0];
		g_odata[blockIdx.x] = idata[0];
	}
}

__global__
void parallel_GPU_reduction_Kernel_v5(float *g_idata, float *g_odata, int *g_i_index, int *g_o_index, int size)
{
	__shared__ float smem[1024];
	__shared__ int smem_index[1024];

	// set thread ID
	unsigned int tid = threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	float *idata = g_idata + blockIdx.x * blockDim.x;
	int *i_index = g_i_index + blockIdx.x * blockDim.x;

	// set to smem by each threads
	smem[tid] = idata[tid];
	smem_index[tid] = i_index[tid];
	__syncthreads();

	if (idx >= size) return;

	// in-place reduction in shared memory
	if (blockDim.x >= 1024 && tid < 512)
	{
		smem[tid] = min(smem[tid], smem[tid + 512]);
		__syncthreads();

		bool Bool = (smem[tid] == smem[tid + 512]);
		smem_index[tid] = Bool * smem_index[tid + 512] + (1-Bool) * smem_index[tid];
		__syncthreads();
	}

	if (blockDim.x >= 512 && tid < 256)
	{	smem[tid] = min(smem[tid], smem[tid + 256]);
		__syncthreads();

		bool Bool = (smem[tid] == smem[tid + 256]);
		smem_index[tid] = Bool * smem_index[tid + 256] + (1-Bool) * smem_index[tid];
		__syncthreads();
	}

	if (blockDim.x >= 256 && tid < 128)
	{
		smem[tid] = min(smem[tid], smem[tid + 128]);
		__syncthreads();

		bool Bool = (smem[tid] == smem[tid + 128]);
		smem_index[tid] = Bool * smem_index[tid + 128] + (1-Bool) * smem_index[tid];
		__syncthreads();
	}

	if (blockDim.x >= 128 && tid < 64)
	{
		smem[tid] = min(smem[tid], smem[tid + 64]);
		__syncthreads();

		bool Bool = (smem[tid] == smem[tid + 64]);
		smem_index[tid] = Bool * smem_index[tid + 64] + (1-Bool) * smem_index[tid];
		__syncthreads();
	}
	// unrolling warp
	if (tid < 32)
	{
		volatile float *vsmem = smem;
		volatile int *vsmem_index = smem_index;

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 32]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 32]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 16]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 16]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 8]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 8]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 4]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 4]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 2]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 2]);
		__syncthreads();

		vsmem[tid] = min(vsmem[tid], vsmem[tid + 1]);
		__syncthreads();
		vsmem_index[tid] = min(vsmem_index[tid], vsmem_index[tid + 1]);
		__syncthreads();
	}

	if (tid == 0)
	{
		g_o_index[blockIdx.x] = smem_index[0];
		g_odata[blockIdx.x] = smem[0];
	}
}

__global__
void parallel_GPU_reduction_Kernel_final_reduction(int s, float *g_odata, int *g_o_index, int size, int *d_Final_indices)
{
	//for (int s = 0; s < ExpectedScenarios; s++)
	{
		__shared__ float min;
		__shared__ int minIndex;


		min = 0;
		minIndex = 0;

		for (int i = 0; i < size; i++)
		{
			if (min == 0 || g_odata[i] < min)
			{
				minIndex = g_o_index[i];
				min = g_odata[i];
			}
		}

		d_Final_indices[s] = minIndex;
	}

}


__global__ void parallelMemoization(solution *Solution, cost *New_Full_Cost_Per_Route)
{

		int s = threadIdx.x + blockDim.x * blockIdx.x;

		if (s < ExpectedScenarios)
		{
			Solution[s].OfflineEvaluateCost();
			for (int k = 0; k < TotalVehicles; k++)
				New_Full_Cost_Per_Route[k +(s * m)] = Solution[s].Vehicle[k].Cost;
		}

}
