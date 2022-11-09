
struct explore_neighborhood
{
	int isIntra;

	int Counter = 0, CounterSize = 5000;
	float NeighborhoodCounter[5000];
	bool DisplayAvgNeighborhood = true;

	int **Edge;

	solution *h_Solution;
	problem *h_Problem;
	sol *cpu;
	cost *h_Full_Cost_Per_Route;
	cost *h_Temp;

	problem *Dev_Problem;
	solution *Dev_Solution;
	sol *gpu;
	cost *Temp;
	cost *Full_Cost_Per_Route;

	container *cpu_Container;
	container *gpu_Container;

	//For constTime_Feasibility_Test
	float *F_slack;
	solution *Ref_Solution;
	int **maxPossible_Feasible_Pickups_next_to_this_index;
	float ***cumulativeWaitingTime;

	explore_neighborhood() {}

	~explore_neighborhood() {}

	//////////////////////////////////////////////////
	////// Memory allocation/deletion ///////////////
	//////////////////////////////////////////////////

	void creation()
	{
		// For constTime_Feasibility_Test
		F_slack = new float[2*n];
		Ref_Solution = new solution;

		maxPossible_Feasible_Pickups_next_to_this_index = new int*[m];
		for (int k = 0; k < m; k++)
			maxPossible_Feasible_Pickups_next_to_this_index[k] = new int[ExpectedPath];


		cumulativeWaitingTime = new float**[m];
		for (int k = 0; k < m; k++)
		{
			cumulativeWaitingTime[k] = new float*[ExpectedPath];

			for (int j = 0; j < ExpectedPath; j++)
			{
				cumulativeWaitingTime[k][j] = new float[ExpectedPath];
				for (int l = 0; l < ExpectedPath; l++) // unnecessary initialization
					cumulativeWaitingTime[k][j][l] = 0;
			}
		}


		// For Candidacy evaluation
		Edge = new int*[2*n];
		for(int i = 0; i < 2*n; i++)
			Edge[i] = new int[2*n];


		cpu_Container = new container[Expectation];
		h_Solution = new solution;
		h_Problem = new problem;
		cpu = new sol[Expectation];
		h_Full_Cost_Per_Route = new cost[TotalVehicles];
		h_Temp = new cost[Expectation];
	}

	void destroy()
		{
			//
			for (int k = 0; k < m; k++)
				delete[] maxPossible_Feasible_Pickups_next_to_this_index[k];
			delete[] maxPossible_Feasible_Pickups_next_to_this_index;

			for (int k = 0; k < m; k++)
			{
				for (int j = 0; j < ExpectedPath; j++)
					delete[] cumulativeWaitingTime[k][j];

				delete[] cumulativeWaitingTime[k];
			}
			delete[] cumulativeWaitingTime;

			delete[] F_slack;
			delete Ref_Solution;

			//
			for(int i = 0; i < (2*n); i++)
				delete[] Edge[i];
			delete[] Edge;


			delete[] cpu_Container;
			delete h_Solution;
			delete h_Problem;
			delete[] cpu;
			delete[] h_Full_Cost_Per_Route;
			delete[] h_Temp;
		}

	void create_GPU_memories()
	{

		CHECK(cudaMalloc((void**)&gpu_Container, Expectation * sizeof(container)));

		CHECK(cudaMalloc((void**)&Dev_Problem, sizeof(problem)));
		CHECK(cudaMalloc((void**)&Dev_Solution, sizeof(solution)));
		CHECK(cudaMalloc((void**)&gpu, Expectation * sizeof(sol)));
		CHECK(cudaMalloc((void**)&Temp, Expectation * sizeof(cost)));
		CHECK(cudaMalloc((void**)&Full_Cost_Per_Route, TotalVehicles * sizeof(cost)));

	}

	void destroy_GPU_memories()
	{
		CHECK(cudaFree(gpu_Container));

		CHECK(cudaFree(Dev_Problem));
		CHECK(cudaFree(Dev_Solution));
		CHECK(cudaFree(gpu));
		CHECK(cudaFree(Temp));
		CHECK(cudaFree(Full_Cost_Per_Route));
	}

	//////////////////////////////////////////////////////////////////////////
	////// Utility Functions - candiacy restriction & feasibility test ///////
	//////////////////////////////////////////////////////////////////////////

	float get_TW_Threshold_for_candidacy(problem *Problem)
	{
		float Final_Median;

		// allocate
		int size = 2*n*2*n;
		float *Unified_TW_Diff_Array = new float[size];
		int index = 0;
		for(int i = 0; i < 2*n; i++)
			for(int j = 0; j < 2*n; j++)
			{
				vertex Vertex_i = find_vertex(Problem[0], i);
				vertex Vertex_j = find_vertex(Problem[0], j);
				Unified_TW_Diff_Array[index] = abs((Vertex_i.EarliestTime+Vertex_i.LatestTime)/2 - (Vertex_j.EarliestTime+Vertex_j.LatestTime)/2);
				index++;
			}
		for(int i = 0; i < size; i++)
			for(int j = 0; j < size; j++)
			{
				if (Unified_TW_Diff_Array[i] > Unified_TW_Diff_Array[j])
				{
					int temp = Unified_TW_Diff_Array[i];
					Unified_TW_Diff_Array[i] = Unified_TW_Diff_Array[j];
					Unified_TW_Diff_Array[j] = temp;
				}
			}
		printf("Unified Median TW Threshold = %0.02f\n", Unified_TW_Diff_Array[(int)(size/2)]);
		delete[] Unified_TW_Diff_Array;
		///



		float *TW_Medians = new float[2*n];
		float **TW_Difference = new float*[2*n];
		for(int i = 0; i < 2*n; i++)
			TW_Difference[i] = new float[2*n];

		// 1) Loading TW_Difference Array
		for(int i = 0; i < 2*n; i++)
			for(int j = 0; j < 2*n; j++)
			{
				vertex Vertex_i = find_vertex(Problem[0], i);
				vertex Vertex_j = find_vertex(Problem[0], j);
				TW_Difference[i][j] = abs((Vertex_i.EarliestTime+Vertex_i.LatestTime)/2 - (Vertex_j.EarliestTime+Vertex_j.LatestTime)/2);
			}

		// 2) Sorting TW_Difference Array (2D-array, remember!)
		int x;
	    for (int k = 0; k < 2*n; k++)
	        for (int l = 0; l < 2*n; l++)
	        {
	            x = l+1;
	            for (int i = k; i < 2*n ; i++)
	            {
	                for (int j = x; j < 2*n; j++)
	                {
	                    if (TW_Difference[k][l] > TW_Difference[i][j])
	                        swap(TW_Difference[k][l], TW_Difference[i][j]);
	                }
	                x=0;
	            }
	        }

		/*for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 50; j++)
			{
				printf("TW_Difference[%d][%d]=%0.00f\n", i, j, TW_Difference[i][j]);
			}
			printf("----\n");
		}*/

	    // 3) Loading the Medians of each bucket of TW_Difference Array into TW_Medians Array
	    for (int i = 0; i < 2*n; i++)
	    {
	    	TW_Medians[i] = TW_Difference[i][n];
	    	printf("%0.0f ", TW_Medians[i]);
	    }

	    // 4) Sorting to find the Central Median value using the TW_Medians Arrray
		for(int i = 0; i < 2*n; i++)
			for(int j = i+1; j < 2*n; j++)
			{
				if (TW_Medians[i] > TW_Medians[j])
				{
					int temp = TW_Medians[i];
					TW_Medians[i] = TW_Medians[j];
					TW_Medians[j] = temp;
				}
			}
		Final_Median = TW_Medians[n];
		printf("\nMedian TW Threshold = %0.02f\n", Final_Median);

		// delete
		delete[] TW_Medians;
		for(int i = 0; i < (2*n); i++)
			delete[] TW_Difference[i];
		delete[] TW_Difference;

		return Final_Median;
	}

	float get_TravelTime_Threshold_for_candidacy(problem *Problem)
	{
		float Final_Median;


		// allocate
		int size = 2*n*2*n;
		float *Unified_TravelTime_Array = new float[size];
		int index = 0;
		for(int i = 0; i < 2*n; i++)
			for(int j = 0; j < 2*n; j++)
			{
				Unified_TravelTime_Array[index] = Problem[0].TravelTime[i][j];
				index++;
			}
		for(int i = 0; i < size; i++)
			for(int j = 0; j < size; j++)
			{
				if (Unified_TravelTime_Array[i] > Unified_TravelTime_Array[j])
				{
					int temp = Unified_TravelTime_Array[i];
					Unified_TravelTime_Array[i] = Unified_TravelTime_Array[j];
					Unified_TravelTime_Array[j] = temp;
				}
			}
		printf("Unified Median T.T. Threshold = %0.02f\n", Unified_TravelTime_Array[(int)(size/2)]);
		delete[] Unified_TravelTime_Array;
		///


		// allocate
		float *TravelTime_Medians = new float[2*n];
		float **TravelTime_Difference = new float*[2*n];
		for(int i = 0; i < 2*n; i++)
			TravelTime_Difference[i] = new float[2*n];

		// 1) Loading TW_Difference Array
		for(int i = 0; i < 2*n; i++)
			for(int j = 0; j < 2*n; j++)
			{
				TravelTime_Difference[i][j] = Problem[0].TravelTime[i][j];
			}

		// 2) Sorting TW_Difference Array (2D-array, remember!)
		int x;
	    for (int k = 0; k < 2*n; k++)
	        for (int l = 0; l < 2*n; l++)
	        {
	            x = l+1;
	            for (int i = k; i < 2*n ; i++)
	            {
	                for (int j = x; j < 2*n; j++)
	                {
	                    if (TravelTime_Difference[k][l] > TravelTime_Difference[i][j])
	                        swap(TravelTime_Difference[k][l], TravelTime_Difference[i][j]);
	                }
	                x=0;
	            }
	        }

		/*for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 300; j++)
			{
				printf("TravelTime_Difference[%d][%d]=%0.00f\n", i, j, TravelTime_Difference[i][j]);
			}
			printf("----\n");
		}*/

	    // 3) Loading the Medians of each bucket of TW_Difference Array into TW_Medians Array
	    for (int i = 0; i < 2*n; i++)
	    {
	    	TravelTime_Medians[i] = TravelTime_Difference[i][n];
	    	printf("%0.0f ", TravelTime_Medians[i]);
	    }

	    // 4) Sorting to find the Central Median value using the TW_Medians Arrray
		for(int i = 0; i < 2*n; i++)
			for(int j = i+1; j < 2*n; j++)
			{
				if (TravelTime_Medians[i] > TravelTime_Medians[j])
				{
					int temp = TravelTime_Medians[i];
					TravelTime_Medians[i] = TravelTime_Medians[j];
					TravelTime_Medians[j] = temp;
				}
			}
		Final_Median = TravelTime_Medians[n+(int)(0.9*n)];
		printf("\nMedian T.T. Threshold = %0.02f\n", Final_Median);

		// delete
		delete[] TravelTime_Medians;
		for(int i = 0; i < (2*n); i++)
			delete[] TravelTime_Difference[i];
		delete[] TravelTime_Difference;

		return Final_Median;
	}

	void fill_Candidate_Edges(problem *Problem)
	{
		float TravelTime_Difference_Threshold = get_TravelTime_Threshold_for_candidacy(Problem);
		float TW_Difference_Threshold = get_TW_Threshold_for_candidacy(Problem);

		int ActiveEdges = 0, InactiveEdges = 0;
		float LongestEdge = 0, ShortestEdge = 0;
		int TotalEdgesAboveThreshold = 0;
		for(int i = 0; i < 2*n; i++)
		{
			for(int j = 0; j < 2*n; j++)
			{
				///
				vertex Vertex_i = find_vertex(Problem[0], i);
				vertex Vertex_j = find_vertex(Problem[0], j);
				float TW_Difference = abs((Vertex_i.EarliestTime+Vertex_i.LatestTime)/2 - (Vertex_j.EarliestTime+Vertex_j.LatestTime)/2);

				float TravelTime = Problem[0].TravelTime[i][j];


				if (1)//TW_Difference < TW_Difference_Threshold && TravelTime < TravelTime_Difference_Threshold)
					Edge[i][j] = 1;
				else
					Edge[i][j] = 0;
				///

				////
				if (Edge[i][j])
					ActiveEdges++;
				else
					InactiveEdges++;

				if (LongestEdge == 0 || Problem[0].TravelTime[i][j] > LongestEdge)
					LongestEdge = Problem[0].TravelTime[i][j];
				if (ShortestEdge == 0 || Problem[0].TravelTime[i][j] < ShortestEdge)
					ShortestEdge = Problem[0].TravelTime[i][j];

				////
			}
		}
		printf("Total => ActiveEdges: %d, Inactive Edges: %d\n", ActiveEdges, InactiveEdges );
		printf("Total => LongestEdge: %0.02f, ShortestEdge: %0.02f\n", LongestEdge, ShortestEdge );
		printf("TotalEdgesAboveThreshold: %d\n", TotalEdgesAboveThreshold);


	}

	void Record_Neighborhood_Count()
	{
		if (Counter < CounterSize)
			NeighborhoodCounter[Counter++] = (float)h_Solution[0].allinsertion;
		else if (Counter == CounterSize && DisplayAvgNeighborhood)
		{
			DisplayAvgNeighborhood = false;
			float Average = 0;
			for (int i = 0; i < CounterSize; i++)
			{
				Average += NeighborhoodCounter[i];
				//printf("NeighborhoodCounter[%d] = %0.02f\n", i, NeighborhoodCounter[i]);
			}
			printf("Avg Neighborhood Explored: %0.02f\n", Average/CounterSize);
		}

	}

	bool isFeasibleInsertion(int req, int start, int gap, int k)
	{
		bool isFeasible = true;

		//Ref_Solution[0].print_with_battery_level();
		//Ref_Solution[0].printinput();

		// 1) Checking User TW feasibility

		for (int k=0; k<m; k++)
			for (int j=0;j<Ref_Solution[0].Vehicle[k].size;j++)
			{
				int node = Ref_Solution[0].Vehicle[k].path[j] - 1;
				//printf("F_slack[%d]=%0.02f\n", node, F_slack[node]);
			}

		int prev_node_to_pickup = Ref_Solution[0].Vehicle[k].path[start-1] - 1;
		int next_node_to_pickup = Ref_Solution[0].Vehicle[k].path[start] - 1;

		int prev_node_to_dropoff = Ref_Solution[0].Vehicle[k].path[start+gap-1] - 1;
		int next_node_to_dropoff = Ref_Solution[0].Vehicle[k].path[start+gap] - 1;

		// step-1: determine earliest possible service time at pickup node
		float B_pickup = max(Ref_Solution[0].Problem.Request[req].pickup.EarliestTime,
					(Ref_Solution[0].Request[prev_node_to_pickup].C + Ref_Solution[0].Problem.TravelTime[prev_node_to_pickup][req]));

		// step-2
		if (B_pickup > Ref_Solution[0].Problem.Request[req].pickup.LatestTime)
			isFeasible = false;

		// step-3: evaluate the time shift at next_node_to_pickup
		float delkPlus1 = ( B_pickup + Ref_Solution[0].Problem.Request[req].pickup.ServiceTime
					+ Ref_Solution[0].Problem.TravelTime[req][next_node_to_pickup] )
					- ( Ref_Solution[0].Request[next_node_to_pickup].B );

		// step-4
		if (delkPlus1 > F_slack[next_node_to_pickup])
			isFeasible = false;

		// step-5: determine earliest possible service time at dropoff node
		float DepartureTime_C_value_of_prev_node_to_dropoff = Ref_Solution[0].Request[prev_node_to_dropoff].C
													+ (delkPlus1 - cumulativeWaitingTime[k][next_node_to_pickup][prev_node_to_dropoff]); // replace it with ***cumulativeWaitingTime

		// step-6:
		float B_dropoff = max(Ref_Solution[0].Problem.Request[req].dropoff.EarliestTime,
						(DepartureTime_C_value_of_prev_node_to_dropoff + Ref_Solution[0].Problem.TravelTime[prev_node_to_dropoff][req+n]) );

		// step-7
		if (B_dropoff > Ref_Solution[0].Problem.Request[req].dropoff.LatestTime)
			isFeasible = false;

		// step-8: evaluate the time shift at next_node_to_dropoff
		float delLplus1 = ( B_dropoff + Ref_Solution[0].Problem.Request[req].dropoff.ServiceTime
					+ Ref_Solution[0].Problem.TravelTime[req+n][next_node_to_dropoff] )
					- ( Ref_Solution[0].Request[next_node_to_dropoff].B );

		// step-9
		if (delLplus1 > F_slack[next_node_to_dropoff])
			isFeasible = false;


		// 2) Checking Vehicle Load Feasibility

		if (maxPossible_Feasible_Pickups_next_to_this_index[k][start-1]  < 1) //load infeasibility check
			isFeasible = false;


		// 3) Checking Vehicle Route Duration Feasibility
		if (Ref_Solution[0].Vehicle[k].Cost.duration + delLplus1 > Ref_Solution[0].Problem.Vehicle[k].duration_constraint)
			isFeasible = false;

		// 4) Checking User ride time Feasibility

		return isFeasible;

	}

	void Setup_based_on_ConstantTime_FeasibilityTest(int req)
	{

		// 1) initialization

				Ref_Solution[0].Copy(h_Solution[0]);
				Ref_Solution[0].Offline_as_early_as_possible_evaluation(); //compute serve-as-early-as-possible schedule.
				//Ref_Solution[0].print_with_battery_level();

		// 2) Pre-computation of necessary variables

				// 2.1) Pre-compute seats_to_be_occupied_until_next_zeroSplitPoint (to ensure capacity constraint)

				for (int k = 0; k < m; k++)
				{
					for (int j = 0; j < Ref_Solution[0].Vehicle[k].size; j++)
					{
						int node = Ref_Solution[0].Vehicle[k].path[j] - 1;

						int maxLoad = 0;
						for (int l = j; l < Ref_Solution[0].Vehicle[k].size; l++)
						{
							int nextNode = Ref_Solution[0].Vehicle[k].path[l] - 1;
							if (l == j || Ref_Solution[0].Request[nextNode].p != 0) //zeroSplitPoint; l==0 => to inlcude the first point
								maxLoad = max(maxLoad, Ref_Solution[0].Request[nextNode].p);
							else
								break;
						}

						maxPossible_Feasible_Pickups_next_to_this_index[k][j] = max(0, Ref_Solution[0].Problem.Vehicle[k].capacity - maxLoad);

						//printf("maxPossible_Feasible_Pickups_next_to_this_index[%d][%d] = %d\n", k, j, maxPossible_Feasible_Pickups_next_to_this_index[k][j]);
					}
					//printf("---\n");
				}

				// 2.2) Pre-compute cumulativeWaitingTime
				for (int k = 0; k < m; k++)
					for (int j = 0; j < Ref_Solution[0].Vehicle[k].size; j++)
					{
						for (int l = 0; l < Ref_Solution[0].Vehicle[k].size; l++)
						{
							cumulativeWaitingTime[k][j][l] = 0;
							for (int index = j+1; index <= l; index++)
							{
								int node = Ref_Solution[0].Vehicle[k].path[index] - 1;
								cumulativeWaitingTime[k][j][l] += Ref_Solution[0].Request[node].W;
							}
							//printf("cumulativeWaitingTime[%d][%d][%d] = %0.02f\n", k, j, l, cumulativeWaitingTime[k][j][l]);

						}
						//printf("---");
					}

				// 2.3) Pre-compute all the min slacks
				for (int k = 0; k < m; k++)
					for (int index = 0; index < Ref_Solution[0].Vehicle[k].size; index++)
					{
						float cumulative_waiting_time = 0;
						int j = Ref_Solution[0].Vehicle[k].path[index] - 1;
						Ref_Solution[0].compute_slack_for_vertex(k, index, F_slack[j], cumulative_waiting_time);

						//printf("index: %d, VehicleSize: %d\n", index, Ref_Solution[0].Vehicle[k].size);
					}

				//for (int j = 0; j < 2*n; j++)
					//printf("k=%d, F_slack[%d] = %0.02f\n", Ref_Solution[0].Request[j].currvehicle, j, F_slack[j]);


		// 3) Setup operations for only feasible solutions

				h_Solution[0].allinsertion = 0;
				int indices = 0, start = 0;
				for (int k = 0; k < TotalVehicles; k++)
				{
					int size = h_Solution[0].Vehicle[k].size;
					for (int limiter = size; limiter >= 0; limiter--)
						for (start = 0; start <= limiter; start++)
						{
							cpu_Container[indices].start = start;
							cpu_Container[indices].gap = size - limiter;
							cpu_Container[indices].vehid = k;

							if (start > 0 && start < limiter && cpu_Container[indices].gap > 2)
							{
								if(isFeasibleInsertion(req, cpu_Container[indices].start, cpu_Container[indices].gap, cpu_Container[indices].vehid))
									indices++;

							}
							else
								indices++;

						}
				}
				h_Solution[0].allinsertion = indices;

				// module: to record the neighborhood size
				//Record_Neighborhood_Count();

				if (h_Solution[0].allinsertion == 0)
					printf("ERROR: allinsertion becomes zero due to candidacy issue!!!\n");


	}

	//////////////////////////////////////////////////
	////// SEQUENTIAL NEIGHBORHOOD EXPLORATION ///////
	//////////////////////////////////////////////////

	void explore_full_neighborhood(int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice)
	{

		// Copy in the solution

		h_Problem[0].Copy(TempProblem[0]);
		h_Solution[0].Copy(TempSolution[0]);
		h_Solution[0].OfflineEvaluateCost();

		// Do the stuff

		Memoization(h_Solution);

			//
			//clock_t KernelTime = clock();

		seq_SETUP(cpu, request[p] + 1, h_Solution);

			//printf("CPU Setup Execution time: %0.06f\n\n", (float)(clock() - KernelTime) / CLOCKS_PER_SEC * 1000);
			//

		seq_INSERTION(cpu, request[p] + 1, h_Solution);
		seq_EVALUATION(cpu, request[p] + 1, h_Solution, h_Problem);
		seq_UPDATE(cpu, request[p] + 1, h_Solution, h_Temp, h_Full_Cost_Per_Route);

		//inbuilt_Reduction(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);
		inbuilt_Reduction_v2(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);

		h_Solution[0].ReInsertRequest(Insertion[p].request, Insertion[p].vehid, Insertion[p].start, Insertion[p].start + Insertion[p].gap + 1);
		h_Solution[0].OfflineEvaluateCost();


		////// Copy back the solution

		TempSolution[0].Copy(h_Solution[0]);
		TempSolution[0].OfflineEvaluateCost();


	}

	void explore_reduced_neighborhood(int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice)
	{


		h_Problem[0].Copy(TempProblem[0]);

		///////////////////////////////////
		//// 1) INSERT PICK UP POINT
		///////////////////////////////////

			h_Solution[0].Copy(TempSolution[0]);
			h_Solution[0].OfflineEvaluateCost();

			// 2.1.1) update allinsertion variable
			h_Solution[0].allinsertion = 0;
			for (int k = 0; k < m; k++)
				h_Solution[0].allinsertion += h_Solution[0].Vehicle[k].size + 1;

			// Launch kernels

			Memoization(h_Solution);

			seq_SETUP_PickUp_First(cpu, request[p] + 1, h_Solution);
			seq_INSERTION_PickUp_First(cpu, request[p] + 1, h_Solution);
			seq_EVALUATION(cpu, request[p] + 1, h_Solution, h_Problem);
			seq_UPDATE(cpu, request[p] + 1, h_Solution, h_Temp, h_Full_Cost_Per_Route);

			inbuilt_Reduction(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);
			//inbuilt_Reduction_v2(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);


			// Print
/*			printf("PICKUP::");
			for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
			{
				printf("cpu[%d].gap=%d, cpu[idx].routesize=%d, ", idx, cpu[idx].gap, cpu[idx].routesize);
				if (choice == 0)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function(h_Temp[idx]));
				else if (choice == 1)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]));
				else
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]));
			}
			printf("==> PICKUP: bestvehicle: %d, start: %d, gap: %d,", Insertion[p].vehid, Insertion[p].start, Insertion[p].gap);
			if (choice == 0)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
			else if (choice == 1)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
			else
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
			printf("-----\n");*/

			// Re-insert pickup point

			h_Solution[0].ReInsertRequest_PickUp_Point(Insertion[p].request, Insertion[p].vehid, Insertion[p].start);
			h_Solution[0].OfflineEvaluateCost();

//			printf("---");
//			h_Solution[0].Cost.print();

		///////////////////////////////////
		//// 2) INSERT DROP OFF POINT
		///////////////////////////////////

			int gpuStart = Insertion[p].start;
			int gpuCurr_k = Insertion[p].vehid;

			h_Solution[0].allinsertion = h_Solution[0].Vehicle[gpuCurr_k].size - Insertion[p].start;


			int Prev_allInsertion = -1;
			for (int k = 0; k < m; k++)
				Prev_allInsertion += h_Solution[0].Vehicle[k].size + 1;

			//printf("At Iter, allinsertion: %d\n", Prev_allInsertion + h_Solution[0].allinsertion);



			/*if (Insertion[p].vehid > 1)
			{
				printf("Hit on BreakPoint\n");
			}*/
			//printf("gpuStart: %d, gpuCurr_k: %d\n", gpuStart, gpuCurr_k);

			// Launch Kernels

			Memoization(h_Solution);

			seq_SETUP_DropOff_Second(cpu, request[p] + 1, gpuStart, gpuCurr_k, h_Solution);
			seq_INSERTION_DropOff_Second(cpu, request[p] + 1, h_Solution);
			seq_EVALUATION(cpu, request[p] + 1, h_Solution, h_Problem);
			seq_UPDATE(cpu, request[p] + 1, h_Solution, h_Temp, h_Full_Cost_Per_Route);

			//inbuilt_Reduction(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);
			inbuilt_Reduction_v2(p, h_Solution, Insertion, Coefficient, fsm, choice, 0);

			// Re-insert dropoff point

			h_Solution[0].ReInsertRequest_DropOff_Point(Insertion[p].request, Insertion[p].vehid, Insertion[p].start + Insertion[p].gap + 1);
			h_Solution[0].OfflineEvaluateCost();


//			printf("---");
//			h_Solution[0].Cost.print();


			/*if ((h_Solution[0].Vehicle[0].size + h_Solution[0].Vehicle[1].size) > 32)
			{
				printf("HIT the breakpoint here\n");
			}*/

			// Print
/*			printf("DROPOFF::");
			for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
			{
				printf("cpu[%d].gap=%d, cpu[idx].routesize=%d, ", idx, cpu[idx].gap, cpu[idx].routesize);
				if (choice == 0)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function(h_Temp[idx]));
				else if (choice == 1)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]));
				else
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]));
			}
			printf("==> DROPOFF: bestvehicle: %d, start: %d, gap: %d, ", Insertion[p].vehid, Insertion[p].start, Insertion[p].gap);
			if (choice == 0)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
			else if (choice == 1)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
			else
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
			printf("-----\n");
			printf("-----\n");*/


			////// Copy back the solution

			TempSolution[0].Copy(h_Solution[0]);
			TempSolution[0].OfflineEvaluateCost();

	}

	void explore_reduced_neighborhood_reverse(int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice)
	{

		h_Problem[0].Copy(TempProblem[0]);

		///////////////////////////////////
		//// 1) INSERT PICK UP POINT
		///////////////////////////////////

			h_Solution[0].Copy(TempSolution[0]);
			h_Solution[0].OfflineEvaluateCost();

			// 2.1.1) update allinsertion variable
			h_Solution[0].allinsertion = 0;
			for (int k = 0; k < m; k++)
				h_Solution[0].allinsertion += h_Solution[0].Vehicle[k].size + 1;

			// Launch kernels

			Memoization(h_Solution);

			seq_SETUP_PickUp_First(cpu, request[p] + 1, h_Solution);
			seq_INSERTION_PickUp_First(cpu, request[p] + 1 + n, h_Solution);
			seq_EVALUATION(cpu, request[p] + 1, h_Solution, h_Problem);
			seq_UPDATE(cpu, request[p] + 1, h_Solution, h_Temp, h_Full_Cost_Per_Route);

			inbuilt_Reduction(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);
			//inbuilt_Reduction_v2(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);


			// Print
/*			printf("PICKUP::");
			for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
			{
				printf("cpu[%d].gap=%d, cpu[idx].routesize=%d, ", idx, cpu[idx].gap, cpu[idx].routesize);
				if (choice == 0)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function(h_Temp[idx]));
				else if (choice == 1)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]));
				else
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]));
			}
			printf("==> PICKUP: bestvehicle: %d, start: %d, gap: %d,", Insertion[p].vehid, Insertion[p].start, Insertion[p].gap);
			if (choice == 0)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
			else if (choice == 1)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
			else
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
			printf("-----\n");*/

			// Re-insert pickup point

			h_Solution[0].ReInsertRequest_PickUp_Point_Reverse(Insertion[p].request, Insertion[p].vehid, Insertion[p].start);
			h_Solution[0].OfflineEvaluateCost();


			//printf("AFTER INSERTION OF PICKUP:\n");
			//h_Solution[0].printinput();

//			printf("---");
//			h_Solution[0].Cost.print();

		///////////////////////////////////
		//// 2) INSERT DROP OFF POINT
		///////////////////////////////////

			int gpuStart = Insertion[p].start;
			int gpuCurr_k = Insertion[p].vehid;

			h_Solution[0].allinsertion = Insertion[p].start + 1;

			/*if (Insertion[p].vehid > 1)
			{
				printf("Hit on BreakPoint\n");
			}*/
			//printf("gpuStart: %d, gpuCurr_k: %d\n", gpuStart, gpuCurr_k);

			// Launch Kernels

			Memoization(h_Solution);

			seq_SETUP_DropOff_Second(cpu, request[p] + 1, gpuStart, gpuCurr_k, h_Solution);
			seq_INSERTION_DropOff_Second(cpu, request[p] + 1 - n, h_Solution, 1);
			seq_EVALUATION(cpu, request[p] + 1, h_Solution, h_Problem);
			seq_UPDATE(cpu, request[p] + 1, h_Solution, h_Temp, h_Full_Cost_Per_Route);

			//inbuilt_Reduction(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);
			inbuilt_Reduction_v2(p, h_Solution, Insertion, Coefficient, fsm, choice, 0);

			// Re-insert dropoff point

			h_Solution[0].ReInsertRequest_DropOff_Point_Reverse(Insertion[p].request, Insertion[p].vehid, Insertion[p].gap);

			h_Solution[0].OfflineEvaluateCost();

			//printf("AFTER INSERTION OF DROPOFF:\n");
			//h_Solution[0].printinput();


//			printf("---");
//			h_Solution[0].Cost.print();


			/*if ((h_Solution[0].Vehicle[0].size + h_Solution[0].Vehicle[1].size) > 32)
			{
				printf("HIT the breakpoint here\n");
			}*/

			// Print
/*			printf("DROPOFF::");
			for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
			{
				printf("cpu[%d].gap=%d, cpu[idx].routesize=%d, ", idx, cpu[idx].gap, cpu[idx].routesize);
				if (choice == 0)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function(h_Temp[idx]));
				else if (choice == 1)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]));
				else
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]));
			}
			printf("==> DROPOFF: bestvehicle: %d, start: %d, gap: %d, ", Insertion[p].vehid, Insertion[p].start, Insertion[p].gap);
			if (choice == 0)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
			else if (choice == 1)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
			else
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
			printf("-----\n");
			printf("-----\n");*/


			////// Copy back the solution

			TempSolution[0].Copy(h_Solution[0]);
			TempSolution[0].OfflineEvaluateCost();


			////////////////////////////////////////////////////////////////
			// FINAL LAYER: Updating gap value for generic reinsertion.////
			////////////////////////////////////////////////////////////////
			//printf("BEFORE FINAL LAYER: start=%d, gap=%d\n", Insertion[p].start, Insertion[p].gap);
			int Tmp = Insertion[p].start;
			Insertion[p].start = Insertion[p].gap;
			Insertion[p].gap = Tmp - Insertion[p].gap;
			//printf("AFTER FINAL LAYER: start=%d, gap=%d\n", Insertion[p].start, Insertion[p].gap);


	}

	void inbuilt_Reduction(int p, solution *TempSolution, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp)
	{
		// 1) Reduction
		int found = 0;

		// Alternative one

		Insertion[p].reset();

		for (int idx = 0; idx < TempSolution[0].allinsertion; idx++)
		{
			float currentcost;

				if (choice == 0)
					currentcost = Coefficient[0].cost_function(h_Temp[idx]);
				else if (choice == 1)
					currentcost = Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]);
				else
					currentcost = Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]);

			//if (h_h_Temp[idx].isFeasible)
			{
				if (Insertion[p].local_best[h_Temp[idx].vehicle] == 0 || currentcost < Insertion[p].local_best[h_Temp[idx].vehicle])
				{
					/////////////////////////////////////////////////////////

					if (!fsm[0].BarredList[h_Temp[idx].vehid])
					{
						// if haltindex is the depot on return visit, then last node is nodeBeforeHaltIndex

						Insertion[p].local_best[h_Temp[idx].vehicle] = currentcost;
						Insertion[p].local_best_ID[h_Temp[idx].vehicle] = idx;
						found = 1;

						//printf("	CPU Entry: %d,	vehicle: %d\n", idx, h_Temp[idx].vehicle);

					}
					else
						Insertion[p].local_best[h_Temp[idx].vehid] = 100000; // setting high-cost just to avoid barred vehicles
				}
			}

			/*if (PickUp)
				printf("PICKUP: h_Temp[%d]: currentcost = %0.02f\n", idx, currentcost);
			else
				printf("DROPOFF: h_Temp[%d]: currentcost = %0.02f\n", idx, currentcost);*/

		}

		if (PickUp)
			Insertion[p].sort_localbest();

		int min_cost_vehicle_k = Insertion[p].veh_sequence[0];
		int idx = Insertion[p].local_best_ID[min_cost_vehicle_k];

		// calculating gradings

		for (int index = 1; index < m; index++) //beware - upto (m-1) only
			Insertion[p].cost_difference[index - 1] = Insertion[p].local_best[index] - Insertion[p].local_best[0];

		// 2) Extraction
		//int idx = local_ID;

		Insertion[p].Extract(cpu[idx]);

		/*printf("best idx: %d; best vehicle: %d(%d); ", idx, min_cost_vehicle_k, Insertion[p].vehid);
		if (choice == 0)
			printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
		else if (choice == 1)
			printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
		else
			printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
*/
		/*printf("-------------------------------------------------");
		for (int k = 0; k < m; k++)
		{
		printf("For vehicle: %d\n", k);
		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		if (h_Temp[idx].vehid == k)
		printf("	%d	%0.0f	%s\n", idx, h_Temp[idx].travel_cost, h_Temp[idx].isFeasible == 0 ? "NOT Feasible" : "Feasible");

		}*/


		/*printf("CPU_ veh_Sequence: \n");
		for (int k = 0; k < m; k++)
		{
		printf("	%d	%f\n", Insertion[p].veh_sequence[k], Insertion[p].local_best[k]);
		}
		printf("\n");


		printf("\nCPU BEST SELECTED: %d\n", idx);
		h_Temp[idx].print();
		Insertion[p].print();
		// system("PAUSE");*/

	}

	void inbuilt_Reduction_v2(int p, solution *TempSolution, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp)
	{
		// 1) Reduction
		int found = 0;
		int local_ID = 0; //overall best Insertion
		float localbest = 0;

		// Alternative one

		Insertion[p].reset();

		for (int idx = 0; idx < TempSolution[0].allinsertion; idx++)
		{
			float currentcost;

			if (choice == 0)
				currentcost = Coefficient[0].cost_function(h_Temp[idx]);
			else if (choice == 1)
				currentcost = Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]);
			else
				currentcost = Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]);

			//if (h_h_Temp[idx].isFeasible)
			if (localbest == 0 || currentcost < localbest)
			{
				if (!fsm[0].BarredList[h_Temp[idx].vehid])
				{
					localbest = currentcost;
					local_ID = idx;
					found = 1;
				}
			}

		}


		int idx = local_ID;
		Insertion[p].Extract(cpu[idx]);

	}

	void Memoization(solution *Solution)
	{
		solution *Dummy_Solution = new solution;
		Dummy_Solution[0].Copy(Solution[0]);

		for (int k = 0; k < TotalVehicles; k++)
			h_Full_Cost_Per_Route[k] = Dummy_Solution[0].OfflineFetchRouteCost(k);

		delete Dummy_Solution;
	}

	void seq_SETUP(sol *cpu, int reqID, solution *h_Solution)
	{
		/*int indices = 0, start = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			int size = h_Solution[0].Vehicle[k].size;

			for (int limiter = size; limiter >= 0; limiter--)
			{
				for (start = 0; start <= limiter; start++)
				{
					cpu[indices].start = start;
					cpu[indices].gap = size - limiter;
					cpu[indices].vehid = k;
					indices++;

					//if (indices >= Solution[0].allinsertion)
					//break;
				}

				//if (indices >= Solution[0].allinsertion)
				//break;
			}

			if (indices >= h_Solution[0].allinsertion)
				break;
		}
		*/


		// Importing ConstTime_TW_feasibility test
		Setup_based_on_ConstantTime_FeasibilityTest(reqID-1);
		for (int i = 0; i < h_Solution[0].allinsertion; i++)
		{
			cpu[i].start = cpu_Container[i].start;
			cpu[i].gap = cpu_Container[i].gap;
			cpu[i].vehid = cpu_Container[i].vehid;
		}


		//// check-for-intra
		if (isIntra)
		{
			int index = 0;
			for (int i = 0; i < h_Solution[0].allinsertion; i++)
				if (h_Solution[0].Request[reqID-1].currvehicle == cpu_Container[i].vehid)
				{
					cpu[index].start = cpu_Container[i].start;
					cpu[index].gap = cpu_Container[i].gap;
					cpu[index].vehid = cpu_Container[i].vehid;
					index++;
				}
			h_Solution[0].allinsertion = index;
		}
	}

	void seq_INSERTION(sol *cpu, int reqID, solution *h_Solution)
	{

		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		{
			cpu[idx].request = reqID - 1;
			cpu[idx].fixation(h_Solution[0].Vehicle[cpu[idx].vehid].size);
			cpu[idx].route[cpu[idx].start] = reqID;
			cpu[idx].route[cpu[idx].start + cpu[idx].gap + 1] = reqID + n;


			if (cpu[idx].routesize == 2) {}

			else
				for (int tid = 0; tid < cpu[idx].routesize; tid++)
				{
					bool boolone = tid >= (cpu[idx].start);
					bool booltwo = tid >= (cpu[idx].start + cpu[idx].gap);
					cpu[idx].route[tid + boolone + booltwo] = h_Solution[0].Vehicle[cpu[idx].vehid].path[tid];
				}
		}

	}

	void seq_EVALUATION(sol *cpu, int reqID, solution *h_Solution, problem *h_Problem)
	{
		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		{

			cpu[idx].CPU_evaluateByRoute(&h_Solution[0].Problem);
			//cpu[idx].CPU_evaluateByRoute(h_Problem);
			//h_Solution[0].OfflineEvaluateCost();
		}

	}

	void seq_UPDATE(sol *cpu, int reqID, solution *h_Solution, cost *h_Temp, cost *h_Full_Cost_Per_Route)
	{

		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		{
			h_Temp[idx].reset();

			for (int k = 0; k < TotalVehicles; k++)
			{
				if (k == cpu[idx].vehid)
				{
					cpu[idx].CPU_evaluateByRoute(h_Problem);

					h_Temp[idx].travel_cost += cpu[idx].Vehicle.Cost.travel_cost;
					h_Temp[idx].time_window += cpu[idx].Vehicle.Cost.time_window;
					h_Temp[idx].ride_time += cpu[idx].Vehicle.Cost.ride_time;
					h_Temp[idx].load += cpu[idx].Vehicle.Cost.load;
					h_Temp[idx].duration += cpu[idx].Vehicle.Cost.duration;
					h_Temp[idx].excessrideTime += cpu[idx].Vehicle.Cost.excessrideTime;
				}
				else
				{
					//h_Solution[0].OfflineEvaluateByRoute(k);

					h_Temp[idx].travel_cost += h_Full_Cost_Per_Route[k].travel_cost;
					h_Temp[idx].time_window += h_Full_Cost_Per_Route[k].time_window;
					h_Temp[idx].ride_time += h_Full_Cost_Per_Route[k].ride_time;
					h_Temp[idx].load += h_Full_Cost_Per_Route[k].load;
					h_Temp[idx].duration += h_Full_Cost_Per_Route[k].duration;
					h_Temp[idx].excessrideTime += h_Full_Cost_Per_Route[k].excessrideTime;
				}

			}



			h_Temp[idx].vehicle = cpu[idx].vehid;

			h_Temp[idx].isFeasible = h_Temp[idx].getFeasibility();
			cpu[idx].isFeasible = h_Temp[idx].getFeasibility();




			h_Temp[idx].vehid = cpu[idx].vehid;

			// Extract start value for haltIndex
			h_Temp[idx].start = cpu[idx].start;
			h_Temp[idx].vehicle = cpu[idx].vehid;


			//printf("h_Temp[%d].vehicle = %d\n", idx, h_Temp[idx].vehicle);

		}


	}

	void seq_SETUP_PickUp_First(sol *cpu, int reqID, solution *Solution)
	{
		int indices = 0, start = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			int size = Solution[0].Vehicle[k].size;

			for (start = 0; start <= size; start++)
			{
				cpu[indices].request = reqID - 1;
				cpu[indices].start = start;
				cpu[indices].gap = 0; /// NEEEEEED TO REMOVEEEE THISSSS ONEEE!!!!!!!!!
				cpu[indices].vehid = k;
				indices++;
			}

			if (indices >= Solution[0].allinsertion)
				break;
		}
	}

	void seq_INSERTION_PickUp_First(sol *cpu, int reqID, solution *Solution)
	{

		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		{
			cpu[idx].routesize = Solution[0].Vehicle[ cpu[idx].vehid ].size + 1;
			cpu[idx].route[cpu[idx].start] = reqID;

			if (cpu[idx].routesize == 2) {}
			else
				for (int tid = 0; tid < cpu[idx].routesize - 1; tid++)
				{
					bool boolone = tid >= (cpu[idx].start);
					cpu[idx].route[tid + boolone] = Solution[0].Vehicle[cpu[idx].vehid].path[tid];
				}
		}

		/*printf("PICKUP:\n");
		for (int idx = 0; idx < Solution[0].allinsertion; idx++)
		{
			printf("cpu[%d].route= ", idx);
			for (int tid = 0; tid < cpu[idx].routesize; tid++)
				printf(" %d", cpu[idx].route[tid]);
			printf("\n");
		}
		printf("-----\n");*/
	}

	void seq_SETUP_DropOff_Second(sol *cpu, int reqID, int start, int Curr_k, solution *Solution)
	{

		for (int idx = 0; idx < Solution[0].allinsertion; idx++)
		{
			cpu[idx].gap = idx;

			cpu[idx].request = reqID - 1;
			cpu[idx].start = start;
			cpu[idx].vehid = Curr_k;

			//printf("DROPOFF: cpu[%d].vehid = %d\n", idx, cpu[idx].vehid);

		}

	}

	void seq_INSERTION_DropOff_Second(sol *cpu, int reqID, solution *Solution, int Reversal = 0)
	{

		for (int idx = 0; idx < Solution[0].allinsertion; idx++)
		{
			cpu[idx].routesize = Solution[0].Vehicle[ cpu[idx].vehid ].size + 1;


			if (Reversal == 0)
				cpu[idx].route[cpu[idx].start + cpu[idx].gap + 1] = reqID + n;
			else
				cpu[idx].route[idx] = reqID + n; // reversal-active

			//printf("JUST in case: routesize, vehid => %d %d\n", cpu[idx].routesize, cpu[idx].vehid);

			if (cpu[idx].routesize == 2) {}
			else
			{
				for (int tid = 0; tid < cpu[idx].routesize - 1; tid++)
				{
					bool booltwo;

					if (Reversal == 0)
						booltwo = tid >= (cpu[idx].start + cpu[idx].gap + 1);
					else
						booltwo = tid >= idx;

					cpu[idx].route[tid + booltwo] = Solution[0].Vehicle[cpu[idx].vehid].path[tid];

					//printf("	--> tid=%d, cpu[%d].routesize=%d, cpu[idx].vehid=%d, cpu[idx].route[%d]=%d\n",
						//	tid, idx, cpu[idx].routesize, cpu[idx].vehid, tid+booltwo, cpu[idx].route[tid + booltwo]);

				}
			}
		}

		/*printf("DROPOFF:\n");
		for (int idx = 0; idx < Solution[0].allinsertion; idx++)
		{
			printf("cpu[%d].route= ", idx);
			for (int tid = 0; tid < cpu[idx].routesize; tid++)
				printf(" %d", cpu[idx].route[tid]);
			printf("\n");
		}
		printf("-----\n");*/

	}

	////////////////////////////////////////////////
	////// PARALLEL NEIGHBORHOOD EXPLORATION ///////
	////////////////////////////////////////////////

	void explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest(int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice)
	{

		// Copy in the solution

		h_Problem[0].Copy(TempProblem[0]);
		h_Solution[0].Copy(TempSolution[0]);
		h_Solution[0].OfflineEvaluateCost();

		// Do the stuff

		CHECK(cudaMemcpyAsync(Dev_Problem, h_Problem, sizeof(problem), cudaMemcpyHostToDevice));

		Memoization(h_Solution);
		CHECK(cudaMemcpy(Full_Cost_Per_Route, h_Full_Cost_Per_Route, TotalVehicles * sizeof(cost), cudaMemcpyHostToDevice));

	// 1) SETUP container

		Setup_based_on_ConstantTime_FeasibilityTest(request[p]);

		CHECK(cudaMemcpy(Dev_Solution, h_Solution, sizeof(solution), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(gpu_Container, cpu_Container, h_Solution[0].allinsertion * sizeof(container), cudaMemcpyHostToDevice));
		//IH_parallel_OneKernel KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, gpu_Container, request[p] + 1, Dev_Solution, Dev_Problem, Temp, Full_Cost_Per_Route);


		IH_CUDA_Dummy1 KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, gpu_Container, request[p] + 1, Dev_Solution, Dev_Problem, Temp, Full_Cost_Per_Route);

		int TotalScenarios = 1;
		for (int i = 0; i < TotalScenarios; i++)
			IH_CUDA_Dummy2 KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, gpu_Container, request[p] + 1, Dev_Solution, Dev_Problem, Temp, Full_Cost_Per_Route);

		IH_CUDA_Dummy3 KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, gpu_Container, request[p] + 1, Dev_Solution, Dev_Problem, Temp, Full_Cost_Per_Route);

		CHECK(cudaMemcpy(h_Temp, Temp, h_Solution[0].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));

		inbuilt_GPU_Reduction(p, Insertion, Coefficient, fsm, choice, 1);
		//inbuilt_GPU_Reduction_v2(p, Insertion, Coefficient, fsm, choice, 1);

		h_Solution[0].ReInsertRequest(Insertion[p].request, Insertion[p].vehid, Insertion[p].start, Insertion[p].start + Insertion[p].gap + 1);
		h_Solution[0].OfflineEvaluateCost();

		////// Copy back the solution

		TempSolution[0].Copy(h_Solution[0]);
		TempSolution[0].OfflineEvaluateCost();


	}

	void explore_full_neighborhood_in_GPU_with_candidacy(int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice)
	{

		// Copy in the solution

		h_Problem[0].Copy(TempProblem[0]);
		h_Solution[0].Copy(TempSolution[0]);
		h_Solution[0].OfflineEvaluateCost();

		// Do the stuff

		CHECK(cudaMemcpyAsync(Dev_Problem, h_Problem, sizeof(problem), cudaMemcpyHostToDevice));

		Memoization(h_Solution);
		CHECK(cudaMemcpy(Full_Cost_Per_Route, h_Full_Cost_Per_Route, TotalVehicles * sizeof(cost), cudaMemcpyHostToDevice));

	// 1) SETUP container

		h_Solution[0].allinsertion = 0;
		int indices = 0, start = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			int size = h_Solution[0].Vehicle[k].size;
			for (int limiter = size; limiter >= 0; limiter--)
				for (start = 0; start <= limiter; start++)
				{
					cpu_Container[indices].start = start;
					cpu_Container[indices].gap = size - limiter;
					cpu_Container[indices].vehid = k;

					int befPickup = 0, aftPickup = 0, A = 0, B = 0;

					if (start > 0 && start < limiter)
					{
						befPickup = h_Solution[0].Vehicle[k].path[start-1] - 1;
						aftPickup = h_Solution[0].Vehicle[k].path[start+1] - 1;

						A = Edge[befPickup][request[p]];
						B = Edge[request[p]][aftPickup];

						if (A && B)
							indices++;
					}
					else
						indices++;

				}
		}
		h_Solution[0].allinsertion = indices;
		//printf("allinsertion: %d\n", h_Solution[0].allinsertion);
		if (indices == 0)
			printf("ERROR: allinsertion becomes zero due to candidacy issue!!!\n");
		CHECK(cudaMemcpy(Dev_Solution, h_Solution, sizeof(solution), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(gpu_Container, cpu_Container, h_Solution[0].allinsertion * sizeof(container), cudaMemcpyHostToDevice));
		IH_parallel_OneKernel KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, gpu_Container, request[p] + 1, Dev_Solution, Dev_Problem, Temp, Full_Cost_Per_Route);

		CHECK(cudaMemcpy(h_Temp, Temp, h_Solution[0].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));

		inbuilt_GPU_Reduction(p, Insertion, Coefficient, fsm, choice, 1);
		//inbuilt_GPU_Reduction_v2(p, Insertion, Coefficient, fsm, choice, 1);

		h_Solution[0].ReInsertRequest(Insertion[p].request, Insertion[p].vehid, Insertion[p].start, Insertion[p].start + Insertion[p].gap + 1);
		h_Solution[0].OfflineEvaluateCost();

		////// Copy back the solution

		TempSolution[0].Copy(h_Solution[0]);
		TempSolution[0].OfflineEvaluateCost();


	}

	void explore_full_neighborhood_in_GPU(int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice)
	{

		// Copy in the solution

		h_Problem[0].Copy(TempProblem[0]);
		h_Solution[0].Copy(TempSolution[0]);
		h_Solution[0].OfflineEvaluateCost();

		// Do the stuff

		CHECK(cudaMemcpyAsync(Dev_Problem, h_Problem, sizeof(problem), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpyAsync(Dev_Solution, h_Solution, sizeof(solution), cudaMemcpyHostToDevice));

		Memoization(h_Solution);
		CHECK(cudaMemcpy(Full_Cost_Per_Route, h_Full_Cost_Per_Route, TotalVehicles * sizeof(cost), cudaMemcpyHostToDevice));


				//////
				//clock_t KernelTime = clock();
				//seq_SETUP(cpu, request[p] + 1, h_Solution);
				//CHECK(cudaMemcpy(gpu, cpu, sizeof(h_Solution[0].allinsertion), cudaMemcpyHostToDevice));

				//printf("cpu-gpu MemCpy Execution time: %0.06f\n\n", (float)(clock() - KernelTime) / CLOCKS_PER_SEC * 1000);
								//////

				//////
				//clock_t KernelTime = clock();


	// 1) SETUP container

		int indices = 0, start = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			int size = h_Solution[0].Vehicle[k].size;
			for (int limiter = size; limiter >= 0; limiter--)
				for (start = 0; start <= limiter; start++)
				{
					cpu_Container[indices].start = start;
					cpu_Container[indices].gap = size - limiter;
					cpu_Container[indices].vehid = k;
					indices++;
				}
			if (indices >= h_Solution[0].allinsertion)
				break;
		}

		// module: to record the neighborhood size
		//Record_Neighborhood_Count();


		//printf("allinsertion: %d\n", h_Solution[0].allinsertion);
		CHECK(cudaMemcpy(gpu_Container, cpu_Container, h_Solution[0].allinsertion * sizeof(container), cudaMemcpyHostToDevice));


		//IH_parallel_OneKernel KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, gpu_Container, request[p] + 1, Dev_Solution, Dev_Problem, Temp, Full_Cost_Per_Route);

		IH_parallelSETUP KERNEL_ARGS2(dim3(1), dim3(1))(gpu, request[p] + 1, Dev_Solution);
		IH_parallelINSERTION KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution);

		IH_parallelEVALUATION KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Dev_Problem);


		IH_parallelUPDATE KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Temp, Full_Cost_Per_Route);
		CHECK(cudaDeviceSynchronize());

		CHECK(cudaMemcpy(h_Temp, Temp, h_Solution[0].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));

		inbuilt_GPU_Reduction(p, Insertion, Coefficient, fsm, choice, 1);
		//inbuilt_GPU_Reduction_v2(p, Insertion, Coefficient, fsm, choice, 1);

		h_Solution[0].ReInsertRequest(Insertion[p].request, Insertion[p].vehid, Insertion[p].start, Insertion[p].start + Insertion[p].gap + 1);
		h_Solution[0].OfflineEvaluateCost();

		////// Copy back the solution

		TempSolution[0].Copy(h_Solution[0]);
		TempSolution[0].OfflineEvaluateCost();


	}

	void explore_reduced_neighborhood_in_GPU(int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice)
	{


		h_Problem[0].Copy(TempProblem[0]);

		///////////////////////////////////
		//// 1) INSERT PICK UP POINT
		///////////////////////////////////

			h_Solution[0].Copy(TempSolution[0]);
			h_Solution[0].OfflineEvaluateCost();

			// 2.1.1) update allinsertion variable
			h_Solution[0].allinsertion = 0;
			for (int k = 0; k < m; k++)
				h_Solution[0].allinsertion += h_Solution[0].Vehicle[k].size + 1;
			printf("reduced allInsertion: %d\n", h_Solution[0].allinsertion);
			// Launch kernels


		CHECK(cudaMemcpyAsync(Dev_Problem, h_Problem, sizeof(problem), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpyAsync(Dev_Solution, h_Solution, sizeof(solution), cudaMemcpyHostToDevice));

		Memoization(h_Solution);
		CHECK(cudaMemcpy(Full_Cost_Per_Route, h_Full_Cost_Per_Route, TotalVehicles * sizeof(cost), cudaMemcpyHostToDevice));

		IH_parallelSETUP_AnyNode_First KERNEL_ARGS2(dim3(1), dim3(1))(gpu, request[p] + 1, Dev_Solution);
		IH_parallelINSERTION_AnyNode_First KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution);
		IH_parallelEVALUATION KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Dev_Problem);
		IH_parallelUPDATE KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Temp, Full_Cost_Per_Route);
		CHECK(cudaDeviceSynchronize());

		CHECK(cudaMemcpy(h_Temp, Temp, h_Solution[0].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));

		inbuilt_GPU_Reduction(p, Insertion, Coefficient, fsm, choice, 1);
		//inbuilt_GPU_Reduction_v2(p, Insertion, Coefficient, fsm, choice, 1);


			/*	printf("PICKUP:\n");
				for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
				{
					printf("gpu[%d].route= ", idx);
					for (int tid = 0; tid < gpu[idx].routesize; tid++)
						printf(" %d", gpu[idx].route[tid]);
					printf("\n");
				}
				printf("-----\n");

			printf("PICKUP::");
			for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
			{
				printf("gpu[%d].gap=%d, gpu[idx].routesize=%d, ", idx, gpu[idx].gap, gpu[idx].routesize);
				if (choice == 0)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function(h_Temp[idx]));
				else if (choice == 1)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]));
				else
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]));
			}
			printf("==> PICKUP: bestvehicle: %d, start: %d, gap: %d,", Insertion[p].vehid, Insertion[p].start, Insertion[p].gap);
			if (choice == 0)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
			else if (choice == 1)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
			else
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
			printf("-----\n");*/



		h_Solution[0].ReInsertRequest_PickUp_Point(Insertion[p].request, Insertion[p].vehid, Insertion[p].start);
		h_Solution[0].OfflineEvaluateCost();


		///////////////////////////////////
		//// 2) INSERT DROP OFF POINT
		///////////////////////////////////

			int gpuStart = Insertion[p].start;
			int gpuCurr_k = Insertion[p].vehid;

			h_Solution[0].allinsertion = h_Solution[0].Vehicle[gpuCurr_k].size - Insertion[p].start;


			int Prev_allInsertion = -1;
			for (int k = 0; k < m; k++)
				Prev_allInsertion += h_Solution[0].Vehicle[k].size + 1;

			// Launch Kernels

			if (SecondInsertionInCPU)
			{
				Memoization(h_Solution);

				seq_SETUP_DropOff_Second(cpu, request[p] + 1, gpuStart, gpuCurr_k, h_Solution);
				seq_INSERTION_DropOff_Second(cpu, request[p] + 1, h_Solution);
				seq_EVALUATION(cpu, request[p] + 1, h_Solution, h_Problem);
				seq_UPDATE(cpu, request[p] + 1, h_Solution, h_Temp, h_Full_Cost_Per_Route);

				//inbuilt_Reduction(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);
				inbuilt_Reduction_v2(p, h_Solution, Insertion, Coefficient, fsm, choice, 0);

			}
			else
			{
				CHECK(cudaMemcpyAsync(Dev_Problem, h_Problem, sizeof(problem), cudaMemcpyHostToDevice));
				CHECK(cudaMemcpyAsync(Dev_Solution, h_Solution, sizeof(solution), cudaMemcpyHostToDevice));

				Memoization(h_Solution);
				CHECK(cudaMemcpy(Full_Cost_Per_Route, h_Full_Cost_Per_Route, TotalVehicles * sizeof(cost), cudaMemcpyHostToDevice));

				IH_parallelSETUP_AnyNode_Second KERNEL_ARGS2(dim3(1), dim3(1))(gpu, request[p] + 1, gpuStart, gpuCurr_k, Dev_Solution);
				IH_parallelINSERTION_DropOff_Second KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution);
				IH_parallelEVALUATION KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Dev_Problem);
				IH_parallelUPDATE KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Temp, Full_Cost_Per_Route);
				CHECK(cudaDeviceSynchronize());

				CHECK(cudaMemcpy(h_Temp, Temp, h_Solution[0].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));

				//inbuilt_GPU_Reduction(p, Insertion, Coefficient, fsm, choice, 0);
				inbuilt_GPU_Reduction_v2(p, Insertion, Coefficient, fsm, choice, 0);

			}

			// Print

			/*printf("DROPOFF INSERTIONS:\n");
				for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
				{
					printf("gpu[%d].route= ", idx);
					for (int tid = 0; tid < gpu[idx].routesize; tid++)
						printf(" %d", gpu[idx].route[tid]);
					printf("\n");
				}
				printf("-----\n");

		printf("DROPOFF::");
		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		{
			printf("gpu[%d].gap=%d, gpu[idx].routesize=%d, ", idx, gpu[idx].gap, gpu[idx].routesize);
			if (choice == 0)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function(h_Temp[idx]));
			else if (choice == 1)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]));
			else
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]));
		}
		printf("==> DROPOFF: bestvehicle: %d, start: %d, gap: %d, ", Insertion[p].vehid, Insertion[p].start, Insertion[p].gap);
		if (choice == 0)
			printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
		else if (choice == 1)
			printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
		else
			printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
		printf("-----\n");
		printf("-----\n");*/

			// Re-insert dropoff point

			h_Solution[0].ReInsertRequest_DropOff_Point(Insertion[p].request, Insertion[p].vehid, Insertion[p].start + Insertion[p].gap + 1);
			h_Solution[0].OfflineEvaluateCost();


			////// Copy back the solution

			TempSolution[0].Copy(h_Solution[0]);
			TempSolution[0].OfflineEvaluateCost();

	}

	void explore_reduced_neighborhood_reverse_in_GPU(int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice)
	{


		h_Problem[0].Copy(TempProblem[0]);

		///////////////////////////////////
		//// 1) INSERT PICK UP POINT
		///////////////////////////////////

			h_Solution[0].Copy(TempSolution[0]);
			h_Solution[0].OfflineEvaluateCost();

			// 2.1.1) update allinsertion variable
			h_Solution[0].allinsertion = 0;
			for (int k = 0; k < m; k++)
				h_Solution[0].allinsertion += h_Solution[0].Vehicle[k].size + 1;

			// Launch kernels

		CHECK(cudaMemcpyAsync(Dev_Problem, h_Problem, sizeof(problem), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpyAsync(Dev_Solution, h_Solution, sizeof(solution), cudaMemcpyHostToDevice));

		Memoization(h_Solution);
		CHECK(cudaMemcpy(Full_Cost_Per_Route, h_Full_Cost_Per_Route, TotalVehicles * sizeof(cost), cudaMemcpyHostToDevice));

		IH_parallelSETUP_AnyNode_First KERNEL_ARGS2(dim3(1), dim3(1))(gpu, request[p] + 1, Dev_Solution);
		IH_parallelINSERTION_AnyNode_First KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1 + n, Dev_Solution);
		IH_parallelEVALUATION KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Dev_Problem);
		IH_parallelUPDATE KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Temp, Full_Cost_Per_Route);
		CHECK(cudaDeviceSynchronize());

		CHECK(cudaMemcpy(h_Temp, Temp, h_Solution[0].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));

		inbuilt_GPU_Reduction(p, Insertion, Coefficient, fsm, choice, 1);
		//inbuilt_GPU_Reduction_v2(p, Insertion, Coefficient, fsm, choice, 1);


		/*		printf("PICKUP:\n");
				for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
				{
					printf("gpu[%d].route= ", idx);
					for (int tid = 0; tid < gpu[idx].routesize; tid++)
						printf(" %d", gpu[idx].route[tid]);
					printf("\n");
				}
				printf("-----\n");

			printf("PICKUP::");
			for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
			{
				printf("gpu[%d].gap=%d, gpu[idx].routesize=%d, ", idx, gpu[idx].gap, gpu[idx].routesize);
				if (choice == 0)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function(h_Temp[idx]));
				else if (choice == 1)
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]));
				else
					printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]));
			}
			printf("==> PICKUP: bestvehicle: %d, start: %d, gap: %d,", Insertion[p].vehid, Insertion[p].start, Insertion[p].gap);
			if (choice == 0)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
			else if (choice == 1)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
			else
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
			printf("-----\n");*/



		h_Solution[0].ReInsertRequest_PickUp_Point_Reverse(Insertion[p].request, Insertion[p].vehid, Insertion[p].start);
		h_Solution[0].OfflineEvaluateCost();

		/*printf("PARTIAL SOLUTION:\n");
		h_Solution[0].printinput();*/


		///////////////////////////////////
		//// 2) INSERT DROP OFF POINT
		///////////////////////////////////

			int gpuStart = Insertion[p].start;
			int gpuCurr_k = Insertion[p].vehid;

			h_Solution[0].allinsertion = Insertion[p].start + 1;

			//printf("gpuStart=%d, gpuCurr_k=%d, allinsertion=%d\n", gpuStart, gpuCurr_k, h_Solution[0].allinsertion);

			int Prev_allInsertion = -1;
			for (int k = 0; k < m; k++)
				Prev_allInsertion += h_Solution[0].Vehicle[k].size + 1;

			// Launch Kernels

			if (SecondInsertionInCPU)
			{
				Memoization(h_Solution);

				seq_SETUP_DropOff_Second(cpu, request[p] + 1, gpuStart, gpuCurr_k, h_Solution);
				seq_INSERTION_DropOff_Second(cpu, request[p] + 1 - n, h_Solution, 1);
				seq_EVALUATION(cpu, request[p] + 1, h_Solution, h_Problem);
				seq_UPDATE(cpu, request[p] + 1, h_Solution, h_Temp, h_Full_Cost_Per_Route);

				//inbuilt_Reduction(p, h_Solution, Insertion, Coefficient, fsm, choice, 1);
				inbuilt_Reduction_v2(p, h_Solution, Insertion, Coefficient, fsm, choice, 0);
			}
			else
			{
				CHECK(cudaMemcpyAsync(Dev_Problem, h_Problem, sizeof(problem), cudaMemcpyHostToDevice));
				CHECK(cudaMemcpyAsync(Dev_Solution, h_Solution, sizeof(solution), cudaMemcpyHostToDevice));

				Memoization(h_Solution);
				CHECK(cudaMemcpy(Full_Cost_Per_Route, h_Full_Cost_Per_Route, TotalVehicles * sizeof(cost), cudaMemcpyHostToDevice));

				IH_parallelSETUP_AnyNode_Second KERNEL_ARGS2(dim3(1), dim3(1))(gpu, request[p] + 1, gpuStart, gpuCurr_k, Dev_Solution);
				IH_parallelINSERTION_Pickup_Second KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution);
				IH_parallelEVALUATION KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Dev_Problem);
				IH_parallelUPDATE KERNEL_ARGS2(dim3(round((h_Solution[0].allinsertion) / (32)) + 1), dim3(32))(gpu, request[p] + 1, Dev_Solution, Temp, Full_Cost_Per_Route);
				CHECK(cudaDeviceSynchronize());

				CHECK(cudaMemcpy(h_Temp, Temp, h_Solution[0].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));

				//inbuilt_GPU_Reduction(p, Insertion, Coefficient, fsm, choice, 0);
				inbuilt_GPU_Reduction_v2(p, Insertion, Coefficient, fsm, choice, 0);
			}


			// Print

			/*printf("DROPOFF INSERTIONS:\n");
				for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
				{
					printf("gpu[%d].route= ", idx);
					for (int tid = 0; tid < gpu[idx].routesize; tid++)
						printf(" %d", gpu[idx].route[tid]);
					printf("\n");
				}
				printf("-----\n");

		printf("DROPOFF::");
		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		{
			printf("gpu[%d].gap=%d, gpu[idx].routesize=%d, ", idx, gpu[idx].gap, gpu[idx].routesize);
			if (choice == 0)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function(h_Temp[idx]));
			else if (choice == 1)
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]));
			else
				printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]));
		}
		printf("==> DROPOFF: bestvehicle: %d, start: %d, gap: %d, ", Insertion[p].vehid, Insertion[p].start, Insertion[p].gap);
		if (choice == 0)
			printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
		else if (choice == 1)
			printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
		else
			printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
		printf("-----\n");
		printf("-----\n");*/

			// Re-insert dropoff point

			h_Solution[0].ReInsertRequest_DropOff_Point_Reverse(Insertion[p].request, Insertion[p].vehid, Insertion[p].gap);
			h_Solution[0].OfflineEvaluateCost();


			////// Copy back the solution

			TempSolution[0].Copy(h_Solution[0]);
			TempSolution[0].OfflineEvaluateCost();


			////////////////////////////////////////////////////////////////
			// FINAL LAYER: Updating gap value for generic reinsertion.////
			////////////////////////////////////////////////////////////////
			//printf("BEFORE FINAL LAYER: start=%d, gap=%d\n", Insertion[p].start, Insertion[p].gap);
			int Tmp = Insertion[p].start;
			Insertion[p].start = Insertion[p].gap;
			Insertion[p].gap = Tmp - Insertion[p].gap;
			//printf("AFTER FINAL LAYER: start=%d, gap=%d\n", Insertion[p].start, Insertion[p].gap);

	}

	void inbuilt_GPU_Reduction(int p, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp)
	{

		// 1) Reduction
		int found = 0;

		// Alternative one

		Insertion[p].reset();

		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		{
			float currentcost;

				if (choice == 0)
					currentcost = Coefficient[0].cost_function(h_Temp[idx]);
				else if (choice == 1)
					currentcost = Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]);
				else
					currentcost = Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]);

			//if (h_h_Temp[idx].isFeasible)
			{
				if (Insertion[p].local_best[h_Temp[idx].vehicle] == 0 || currentcost < Insertion[p].local_best[h_Temp[idx].vehicle])
				{
					/////////////////////////////////////////////////////////

					if (!fsm[0].BarredList[h_Temp[idx].vehid])
					{
						// if haltindex is the depot on return visit, then last node is nodeBeforeHaltIndex

						Insertion[p].local_best[h_Temp[idx].vehicle] = currentcost;
						Insertion[p].local_best_ID[h_Temp[idx].vehicle] = idx;
						found = 1;

						//printf("	CPU Entry: %d,	vehicle: %d\n", idx, h_Temp[idx].vehicle);

					}
					else
						Insertion[p].local_best[h_Temp[idx].vehid] = 100000; // setting high-cost just to avoid barred vehicles
				}
			}

			/*if (PickUp)
				printf("PICKUP: h_Temp[%d]: currentcost = %0.02f\n", idx, currentcost);
			else
				printf("DROPOFF: h_Temp[%d]: currentcost = %0.02f\n", idx, currentcost);*/

		}

		if (PickUp)
			Insertion[p].sort_localbest();

		int min_cost_vehicle_k = Insertion[p].veh_sequence[0];
		int idx = Insertion[p].local_best_ID[min_cost_vehicle_k];

		// calculating gradings

		for (int index = 1; index < m; index++) //beware - upto (m-1) only
			Insertion[p].cost_difference[index - 1] = Insertion[p].local_best[index] - Insertion[p].local_best[0];

		// 2) Extraction
		//int idx = local_ID;

		sol *bestSol = new sol;
		CHECK(cudaMemcpy(bestSol, &gpu[idx], sizeof(sol), cudaMemcpyDeviceToHost));
		Insertion[p].Extract(bestSol[0]);
		delete bestSol;

	}

	void inbuilt_GPU_Reduction_v2(int p, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp)
	{
		// 1) Reduction
		int found = 0;
		int local_ID = 0; //overall best Insertion
		float localbest = 0;

		// Alternative one

		Insertion[p].reset();

		for (int idx = 0; idx < h_Solution[0].allinsertion; idx++)
		{
			float currentcost;

			if (choice == 0)
				currentcost = Coefficient[0].cost_function(h_Temp[idx]);
			else if (choice == 1)
				currentcost = Coefficient[0].cost_function_with_batteryViolation(h_Temp[idx]);
			else
				currentcost = Coefficient[0].cost_function_with_batteryViolation_2(h_Temp[idx]);

			//if (h_h_Temp[idx].isFeasible)
			if (localbest == 0 || currentcost < localbest)
			{
				/////////////////////////////////////////////////////////

				if (!fsm[0].BarredList[h_Temp[idx].vehid])
				{
					localbest = currentcost;
					local_ID = idx;
					found = 1;
				}
			}


			// if (PickUp)
			// 	printf("PICKUP: h_Temp[%d]: currentcost = %0.02f\n", idx, currentcost);
			// else
			// 	printf("DROPOFF: h_Temp[%d]: currentcost = %0.02f\n", idx, currentcost);


		}


		int idx = local_ID;

		sol *bestSol = new sol;
		CHECK(cudaMemcpy(bestSol, &gpu[idx], sizeof(sol), cudaMemcpyDeviceToHost));
		Insertion[p].Extract(bestSol[0]);
		delete bestSol;


/*		if (PickUp)
			printf("Pickup is reinserted into vehid %d with start %d\n", cpu[idx].vehid, cpu[idx].start);
		else
			printf("Dropoff is reinserted into vehid %d with gap %d\n", cpu[idx].vehid, cpu[idx].gap);
*/

/*			printf("best idx: %d; best vehicle: %d(%d); ", idx, cpu[idx].vehid, Insertion[p].vehid);
		if (choice == 0)
			printf("Cost: %0.02f\n", Coefficient[0].cost_function(Insertion[p].Cost));
		else if (choice == 1)
			printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation(Insertion[p].Cost));
		else
			printf("Cost: %0.02f\n", Coefficient[0].cost_function_with_batteryViolation_2(Insertion[p].Cost));
*/
		//printf("localbest: %0.02f\n", localbest);

		/*h_Temp[idx].print();
		printf("----");
		Insertion[p].Cost.print();
		printf("-----");*/
	}

	//////////////////////////////////////////////////
	////// EXPLORATION DECISION MAKER FUNCTION ///////
	//////////////////////////////////////////////////

	void Exploration_Decision_Maker(int isGPU, int isFull, int isRev, int &p, solution *TempSolution, problem *TempProblem, int *request, insertion *Insertion, coefficient *Coefficient, FSM *fsm, int choice, int intraRouteInsertion = 0)
	{
		isIntra = 0;
		if (intraRouteInsertion)
		{
			isIntra = 1;
			explore_full_neighborhood(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
		}
		else if (isGPU) //Parallel
		{
			if (isFull)
			{
				//explore_full_neighborhood_in_GPU(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
				//explore_full_neighborhood_in_GPU_with_candidacy(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
				explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
			}
			else
			{
				if (isRev)
					explore_reduced_neighborhood_in_GPU(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
				else
					explore_reduced_neighborhood_reverse_in_GPU(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
			}
		}
		else //Sequential
		{
			if (isFull)
			{
				explore_full_neighborhood(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
			}
			else
			{
				if (isRev)
					explore_reduced_neighborhood(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
				else
					explore_reduced_neighborhood_reverse(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice);
			}
		}
	}

};
