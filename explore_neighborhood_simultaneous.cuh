struct parallelReduction
{
	int T;
	int Expected_B;

	int *B;
	int *h_i_index[ExpectedScenarios];
	int *h_o_index[ExpectedScenarios];

	int *d_i_index[ExpectedScenarios];
	int *d_o_index[ExpectedScenarios];

	float *h_o_data[ExpectedScenarios];
	float *d_o_data[ExpectedScenarios];

	int *h_Final_indices;
	int *d_Final_indices;

	sol *bestSol;

	parallelReduction() {}

	~parallelReduction() {}

	void creation()
	{
		h_Final_indices = new int[ExpectedScenarios];
		CHECK(cudaMalloc((void**)&d_Final_indices, ExpectedScenarios * sizeof(int)));

		T = TotalThreads;//max(size, 1024);
		Expected_B = Expected_Blocks;

		bestSol = new sol[ExpectedScenarios];
		B = new int[ExpectedScenarios];
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			h_i_index[s] = new int[Expectation];
			h_o_index[s] = new int[Expected_B];

			CHECK(cudaMalloc((void**)&d_i_index[s], Expectation * sizeof(int)));
			CHECK(cudaMalloc((void**)&d_o_index[s], Expected_B * sizeof(int)));

			h_o_data[s] = new float[Expected_B];
			CHECK(cudaMalloc((void**)&d_o_data[s], Expected_B * sizeof(float)));

		}
	}

	void destroy()
	{
		delete[] h_Final_indices;
		CHECK(cudaFree(d_Final_indices));

		delete[] bestSol;
		delete[] B;
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			delete[] h_i_index[s];
			delete[] h_o_index[s];

			CHECK(cudaFree(d_i_index[s]));
			CHECK(cudaFree(d_o_index[s]));

			delete[] h_o_data[s];
			CHECK(cudaFree(d_o_data[s]));
		}
	}

};

struct explore_neighborhood
{
	int isIntra;

	//for candidacy test
	int Counter = 0, CounterSize = 5000;
	float NeighborhoodCounter[5000];
	bool DisplayAvgNeighborhood = true;
	int **Edge;

	//for constTime_Feasibility_Test
	float *F_slack;
	solution *Ref_Solution;
	int **maxPossible_Feasible_Pickups_next_to_this_index;
	float ***cumulativeWaitingTime;


	solution *h_Solution; //multiple
	problem *h_Problem; //multiple
	sol *cpu[ExpectedScenarios]; //multiple
	cost *h_Temp[ExpectedScenarios];
	cost *h_Full_Cost_Per_Route[ExpectedScenarios];
	container *cpu_Container[ExpectedScenarios];

	problem *Dev_Problem;
	solution *Dev_Solution;
	sol *gpu[ExpectedScenarios];
	cost *Temp[ExpectedScenarios];
	cost *Full_Cost_Per_Route[ExpectedScenarios];
	container *gpu_Container[ExpectedScenarios];

	coefficient *Dev_Coefficient;
	insertion *Dev_Insertion[ExpectedScenarios];

	cudaStream_t streams[ExpectedScenarios];

	cost *h_Temp_flat;
	cost *Temp_flat;

	float *h_Cost_array_flat;
	float *Dev_Cost_array_flat;

	float *h_Cost_array[ExpectedScenarios];
	float *Dev_Cost_array[ExpectedScenarios];

	parallelReduction *Reduction;

	int TotalHolders;
	holder *h_Holder;
	holder *Dev_Holder;

	skeleton *h_SkeletonData[ExpectedScenarios];
	skeleton *d_SkeletonData[ExpectedScenarios];

	container *Dummy_Container[ExpectedScenarios];

	cost *New_Full_Cost_Per_Route;

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



		h_Solution = new solution[ExpectedScenarios];
		h_Problem = new problem[ExpectedScenarios];
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			cpu[s] = new sol[Expectation];
			h_Temp[s] = new cost[Expectation];
			h_Full_Cost_Per_Route[s] = new cost[TotalVehicles];
			cpu_Container[s] = new container[Expectation];
			Dummy_Container[s] = new container[Expectation];
		}


		CHECK(cudaMalloc((void**)&Dev_Problem, ExpectedScenarios * sizeof(problem)));
		CHECK(cudaMalloc((void**)&Dev_Solution, ExpectedScenarios * sizeof(solution)));
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			CHECK(cudaMalloc((void**)&gpu[s], Expectation * sizeof(sol)));
			CHECK(cudaMalloc((void**)&Temp[s], Expectation * sizeof(cost)));
			CHECK(cudaMalloc((void**)&Full_Cost_Per_Route[s], TotalVehicles * sizeof(cost)));
			CHECK(cudaMalloc((void**)&gpu_Container[s], Expectation * sizeof(container)));
			CHECK(cudaStreamCreate(&streams[s]));
		}

		CHECK(cudaMalloc((void**)&Dev_Coefficient, ExpectedScenarios * sizeof(coefficient)));
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMalloc((void**)&Dev_Insertion[s], TotalRequests * sizeof(insertion)));


		h_Temp_flat = new cost[Expectation * ExpectedScenarios];
		CHECK(cudaMalloc((void**)&Temp_flat, Expectation * ExpectedScenarios * sizeof(cost)));

		h_Cost_array_flat = new float[Expectation * ExpectedScenarios];
		CHECK(cudaMalloc((void**)&Dev_Cost_array_flat, Expectation * ExpectedScenarios * sizeof(float)));

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			h_Cost_array[s] = new float[Expectation];
			CHECK(cudaMalloc((void**)&Dev_Cost_array[s], Expectation * sizeof(float)));
		}

		///////////////////////
		//some initializaiton
		for (int i=0; i<Expectation*ExpectedScenarios; i++)
		{
			h_Temp_flat[i].reset();
			h_Cost_array_flat[i] = 0;
		}
		CHECK(cudaMemcpy(Temp_flat, h_Temp_flat, Expectation * ExpectedScenarios * sizeof(cost), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(Dev_Cost_array_flat, h_Cost_array_flat, Expectation * ExpectedScenarios * sizeof(float), cudaMemcpyHostToDevice));

		for (int s = 0; s < ExpectedScenarios; s++)
			for (int i = 0; i < Expectation; i++)
				h_Cost_array[s][i] = 0;
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpy(Dev_Cost_array[s], h_Cost_array[s], Expectation * sizeof(float), cudaMemcpyHostToDevice));


		Reduction = new parallelReduction;
		Reduction[0].creation();
		// initialization
		for (int s = 0; s < ExpectedScenarios; s++)
			for (int i = 0; i < Expectation; i++)
				Reduction[0].h_i_index[s][i] = i;

		h_Holder = new holder[ExpectedPath];
		CHECK(cudaMalloc((void**)&Dev_Holder, ExpectedPath * sizeof(holder)));

		TotalHolders = 0;
		for (int vehid_size = 0; vehid_size < ExpectedPath; vehid_size++)
		{
			int total_insertions = 0;
			for (int l = 0; l <= vehid_size+1; l++)
				total_insertions += l;
			if (total_insertions >= Expectation)
				break;

			h_Holder[vehid_size].initialize(vehid_size, total_insertions);
			TotalHolders++;


		}
		printf("TotalHolders: %d\n", TotalHolders);
		CHECK(cudaMemcpy(Dev_Holder, h_Holder, TotalHolders * sizeof(holder), cudaMemcpyHostToDevice));


		/////// CHECKING SETUP TIME ///////
		/*for (int i = 0; i < TotalHolders; i++)
			printf("VehSize: %d --> TotaIns: %d\n", h_Holder[i].vehid_size, h_Holder[i].total_insertions);

		int s = 0;
		for (int k = 0; k < m; k++)
			h_Solution[0].Vehicle[k].size = 2+k;//(int)(ExpectedPath/m);//2 + k;
		h_Solution[s].ReBoot();
		CHECK(cudaMemcpy(Dev_Solution, h_Solution, ExpectedScenarios * sizeof(solution), cudaMemcpyHostToDevice));

		clock_t SetupTime = clock();
		for (int i = 0; i < ExpectedScenarios; i++)
			IH_parallel_OneInitialization KERNEL_ARGS4(dim3((Expectation/32) + 1), dim3(32), 0, streams[i])(s, gpu[s], Dev_Solution, Dev_Holder);
		CHECK(cudaDeviceSynchronize());
		printf("CHECKING SETUP TIME: %0.04f ms ", (float)(clock() - SetupTime) / CLOCKS_PER_SEC * 1000);

		CHECK(cudaMemcpy(h_Solution, Dev_Solution, ExpectedScenarios * sizeof(solution), cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(cpu[s], gpu[s], Expectation * sizeof(sol), cudaMemcpyDeviceToHost));

		printf("----------------\n");
		printf("allInsertion: %d\n", h_Solution[s].allinsertion);
		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			if (idx > 1)
				if (cpu[s][idx].vehid != cpu[s][idx-1].vehid)
					printf("xxxx\n");

			printf("cpu[%d]: %d %d %d\n", idx, cpu[s][idx].start, cpu[s][idx].gap, cpu[s][idx].vehid);
		}
		printf("----------------\n");*/

		///////////////////////////////

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			h_SkeletonData[s] = new skeleton[Expectation];
			CHECK(cudaMalloc((void**)&d_SkeletonData[s], Expectation * sizeof(skeleton)));
		}


		CHECK(cudaMalloc((void**)&New_Full_Cost_Per_Route, ExpectedScenarios * m * sizeof(cost)));
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


		delete[] h_Solution;
		delete[] h_Problem;
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			delete[] cpu[s];
			delete[] h_Temp[s];
			delete[] h_Full_Cost_Per_Route[s];
			delete[] cpu_Container[s];
		}


		CHECK(cudaFree(Dev_Problem));
		CHECK(cudaFree(Dev_Solution));
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			CHECK(cudaFree(gpu[s]));
			CHECK(cudaFree(Temp[s]));
			CHECK(cudaFree(Full_Cost_Per_Route[s]));
			CHECK(cudaFree(gpu_Container[s]));;
			CHECK(cudaStreamDestroy(streams[s]));
		}

		CHECK(cudaFree(Dev_Coefficient));
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaFree(Dev_Insertion[s]));


		delete[] h_Temp_flat;
		CHECK(cudaFree(Temp_flat));

		delete[] h_Cost_array_flat;
		CHECK(cudaFree(Dev_Cost_array_flat));

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			delete[] h_Cost_array[s];
			CHECK(cudaFree(Dev_Cost_array[s]))
		}

		Reduction[0].destroy();
		delete Reduction;

		delete[] h_Holder;
		CHECK(cudaFree(Dev_Holder));

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			delete[] h_SkeletonData[s];
			CHECK(cudaFree(d_SkeletonData[s]));
		}

		CHECK(cudaFree(New_Full_Cost_Per_Route));
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

	void Setup_based_on_ConstantTime_FeasibilityTest(int req, int s)
	{

		// 1) initialization

				Ref_Solution[0].Copy(h_Solution[s]);
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

				////////////////////////
				int choice = rand() % 2;
				if (0)
				{

					int indices = 0, start = 0;
					h_Solution[s].ReBoot();
					int terminate = 0;
					for (int k = 0; k < m; k++)
					{
						for (int idx = 0; idx < h_Holder[h_Solution[s].Vehicle[k].size].total_insertions; idx++, indices++)
						{
							cpu_Container[s][idx + h_Solution[s].Vehicle[k].cumulativeInsertions].start = h_Holder[h_Solution[s].Vehicle[k].size].start[idx];
							cpu_Container[s][idx + h_Solution[s].Vehicle[k].cumulativeInsertions].gap = h_Holder[h_Solution[s].Vehicle[k].size].gap[idx];
							cpu_Container[s][idx + h_Solution[s].Vehicle[k].cumulativeInsertions].vehid = k;

							if (indices >= Expectation)
							{
								printf("ERROR INSIDE void SetupCTFeasibilityTest: indices exceeds the Expectation value!!!\n");
								terminate = 1;
								goto next_statement0;
							}
						}

						next_statement0:
							if (terminate == 1)
								break;
					}
					//h_Solution[s].allinsertion = indices;
					/*printf("For Scenario s=%d\n", s);
					printf("comparison: %d %d\n", h_Solution[s].allinsertion, Dev_Solution[s].allinsertion);

					for (int i = 0; i < h_Solution[s].allinsertion; i++)
					{
						printf("In-depth comparison s=%d: %d(%d), %d(%d), %d(%d)\n", s, cpu_Container[s][i].start, Dummy_Container[s][i].start, cpu_Container[s][i].gap, Dummy_Container[s][i].gap, cpu_Container[s][i].vehid, Dummy_Container[s][i].vehid);
					}
					printf("xxxxxx\n");*/


				}
				else if(1) //elegant-setup-CPU
				{
					h_Solution[s].allinsertion = 0;
					int indices = 0, start = 0;

					h_Solution[s].ReBoot();
					int terminate = 0;

					for (int k = 0; k < TotalVehicles; k++)
					{
						int i_start = 0;
						int size = h_Solution[s].Vehicle[k].size;
						for (int bound = size+1; bound>0; bound--, i_start++)
						{
							for (int i_gap = 0; i_gap < bound; i_gap++, indices++)
							{
								cpu_Container[s][indices].start = i_start;
								cpu_Container[s][indices].gap = i_gap;
								cpu_Container[s][indices].vehid = k;

								if (indices >= Expectation)
								{
									printf("ERROR INSIDE void SetupCTFeasibilityTest: indices exceeds the Expectation value!!!\n");
									terminate = 1;
									goto next_statement;
								}
							}
						}

						next_statement:
							if (terminate == 1)
								break;
					}
					h_Solution[s].allinsertion = indices;
				}
				else
				{
					h_Solution[s].allinsertion = 0;
					int indices = 0, start = 0;
					int terminate = 0;
					for (int k = 0; k < TotalVehicles; k++)
					{
						int size = h_Solution[s].Vehicle[k].size;
						for (int limiter = size; limiter >= 0; limiter--)
							for (start = 0; start <= limiter; start++)
							{

								if (indices >= Expectation)
								{
									printf("ERROR INSIDE void SetupCTFeasibilityTest: indices exceeds the Expectation value!!!\n");
									terminate = 1;
									goto next_statement1;
								}

								cpu_Container[s][indices].start = start;
								cpu_Container[s][indices].gap = size - limiter;
								cpu_Container[s][indices].vehid = k;

								if (0)//start > 0 && start < limiter && cpu_Container[s][indices].gap > 2)
								{
									if(isFeasibleInsertion(req, cpu_Container[s][indices].start, cpu_Container[s][indices].gap, cpu_Container[s][indices].vehid))
										indices++;

								}
								else
									indices++;

							}

						next_statement1:
							if (terminate == 1)
								break;
					}
					h_Solution[s].allinsertion = indices;
				}

				//////////////////////////////


				//for (int i = 0; i < h_Solution[s].allinsertion; i++)
					//printf("s:%d sz=%d, start=%d, gap=%d\n",s, h_Solution[s].Vehicle[cpu_Container[s][indices].vehid].size, cpu_Container[s][indices].start, cpu_Container[s][indices].gap);


				// module: to record the neighborhood size
				//Record_Neighborhood_Count();

				if (h_Solution[s].allinsertion == 0 || h_Solution[s].allinsertion > Expectation)
					printf("ERROR: allinsertion becomes zero due to candidacy issue or exceeds Expectation Value!!!\n");


	}

	//////////////////////////////////////////////////
	////// SEQUENTIAL NEIGHBORHOOD EXPLORATION ///////
	//////////////////////////////////////////////////

	void explore_full_neighborhood(int &p, solution *TempSolution, problem *TempProblem, int **request, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, int *currvehicleList, int isIntra = 0)
	{

		// Copy in the solution
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			h_Problem[s].Copy(TempProblem[s]);
			h_Solution[s].Copy(TempSolution[s]);
			h_Solution[s].OfflineEvaluateCost();

			Memoization(s, h_Full_Cost_Per_Route[s]);


			seq_SETUP(cpu[s], request[s][p] + 1, h_Solution, s);
			seq_OneKernel(cpu[s], request[s][p] + 1, h_Temp[s], s, h_Full_Cost_Per_Route[s]);
			inbuilt_Reduction_v2(p, h_Solution, Insertion, Coefficient, fsm, choice, 1, h_Temp[s], cpu[s], s, currvehicleList, isIntra);

			h_Solution[s].ReInsertRequest(Insertion[s][p].request, Insertion[s][p].vehid, Insertion[s][p].start, Insertion[s][p].start + Insertion[s][p].gap + 1);
			h_Solution[s].OfflineEvaluateCost();

			// copy back
			TempSolution[s].Copy(h_Solution[s]);
			TempSolution[s].OfflineEvaluateCost();
		}

	}

	void inbuilt_Reduction_v2(int p, solution *TempSolution, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp, cost *h_Temp, sol *cpu, int s, int *currvehicleList, int isIntra)
	{
		// 1) Reduction
		int found = 0;
		int local_ID = 0; //overall best Insertion
		float localbest = 0;

		// Alternative one

		Insertion[s][p].reset();

		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			float currentcost;

			currentcost = Coefficient[s].cost_function(h_Temp[idx]);

			//if (h_Temp[idx].isFeasible)
			if (isIntra)
			{
				if (h_Temp[idx].vehid == currvehicleList[s])
				{
					if (localbest == 0 || currentcost < localbest)
					{
						localbest = currentcost;
						local_ID = idx;
						found = 1;
					}
				}
			}
			else
			{
				if (localbest == 0 || currentcost < localbest)
				{
					localbest = currentcost;
					local_ID = idx;
					found = 1;
				}
			}


		}


		int idx = local_ID;
		Insertion[s][p].Extract(cpu[idx]);

	}

	void Memoization(int s, cost *h_Full_Cost_Per_Route)
	{
		solution *Dummy_Solution = new solution;
		Dummy_Solution[0].Copy(h_Solution[s]);

		for (int k = 0; k < TotalVehicles; k++)
		{
			h_Full_Cost_Per_Route[k] = Dummy_Solution[0].OfflineFetchRouteCost(k);

			//Dummy_Solution[0].OfflineEvaluateCost();
			//h_Full_Cost_Per_Route[k] = Dummy_Solution[0].Vehicle[k].Cost;

		}
		delete Dummy_Solution;
	}

	void seq_SETUP(sol *cpu, int reqID, solution *h_Solution, int s)
	{

		// Importing ConstTime_TW_feasibility test
		Setup_based_on_ConstantTime_FeasibilityTest(reqID-1, s);
		for (int i = 0; i < h_Solution[s].allinsertion; i++)
		{
			cpu[i].start = cpu_Container[s][i].start;
			cpu[i].gap = cpu_Container[s][i].gap;
			cpu[i].vehid = cpu_Container[s][i].vehid;
		}


		//// check-for-intra
		if (isIntra)
		{
			int index = 0;
			for (int i = 0; i < h_Solution[s].allinsertion; i++)
				if (h_Solution[s].Request[reqID-1].currvehicle == cpu_Container[s][i].vehid)
				{
					cpu[index].start = cpu_Container[s][i].start;
					cpu[index].gap = cpu_Container[s][i].gap;
					cpu[index].vehid = cpu_Container[s][i].vehid;
					index++;
				}
			h_Solution[s].allinsertion = index;
		}
	}

	void seq_OneKernel(sol *cpu, int reqID, cost *h_Temp, int s, cost *h_Full_Cost_Per_Route)
	{

		// 1) Insertion

		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			cpu[idx].request = reqID - 1;
			cpu[idx].fixation(h_Solution[s].Vehicle[cpu[idx].vehid].size);
			cpu[idx].route[cpu[idx].start] = reqID;
			cpu[idx].route[cpu[idx].start + cpu[idx].gap + 1] = reqID + n;

			for (int tid = 0; tid < cpu[idx].routesize; tid++)
			{
				bool boolone = tid >= (cpu[idx].start);
				bool booltwo = tid >= (cpu[idx].start + cpu[idx].gap);
				cpu[idx].route[tid + boolone + booltwo] = h_Solution[s].Vehicle[cpu[idx].vehid].path[tid];
			}
		}

		// 2) Evaluation

		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			cpu[idx].CPU_evaluateByRoute(&h_Solution[s].Problem);
			//cpu[idx].CPU_evaluateByRoute(h_Problem);
			//h_Solution[s].OfflineEvaluateCost();
		}

		// 3) Calculate objective function

		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			h_Temp[idx].reset();

			for (int k = 0; k < TotalVehicles; k++)
			{
				if (k == cpu[idx].vehid)
				{
					cpu[idx].CPU_evaluateByRoute(&h_Solution[s].Problem);

					h_Temp[idx].travel_cost += cpu[idx].Vehicle.Cost.travel_cost;
					h_Temp[idx].time_window += cpu[idx].Vehicle.Cost.time_window;
					h_Temp[idx].ride_time += cpu[idx].Vehicle.Cost.ride_time;
					h_Temp[idx].load += cpu[idx].Vehicle.Cost.load;
					h_Temp[idx].duration += cpu[idx].Vehicle.Cost.duration;
					h_Temp[idx].excessrideTime += cpu[idx].Vehicle.Cost.excessrideTime;
				}
				else
				{
					//h_Solution[s].OfflineEvaluateByRoute(k);

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
		}

	}

	////////////////////////////////////////////////
	////// PARALLEL NEIGHBORHOOD EXPLORATION ///////
	////////////////////////////////////////////////

	void explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest(int &p, solution *TempSolution, problem *TempProblem, int **request, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, myList *ServedRequests_BoundaryMark)
	{

		// Copy in the solution
		clock_t InitialSolutionCopy = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{

			h_Problem[s] = TempProblem[s];
			h_Solution[s] = TempSolution[s];
		}
		//printf("InitialSolutionCopy: %0.04f ms ", (float)(clock() - InitialSolutionCopy) / CLOCKS_PER_SEC * 1000);

if (CurrentTime > 115) {
    printf("Reached 2.2.1: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

		clock_t PblmCopyTime = clock();
		CHECK(cudaMemcpy(Dev_Problem, h_Problem, ExpectedScenarios * sizeof(problem), cudaMemcpyHostToDevice));
		//printf("PblmCopyTime: %0.04f ms ", (float)(clock() - PblmCopyTime) / CLOCKS_PER_SEC * 1000);


if (CurrentTime > 115) {
    printf("Reached 2.2.2: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

        /*
		clock_t MemoizationTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			Memoization(s, h_Full_Cost_Per_Route[s]);
		//printf("Memoization: %0.04f ms ", (float)(clock() - MemoizationTime) / CLOCKS_PER_SEC * 1000);

		clock_t FirstMemCpyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpy(Full_Cost_Per_Route[s], h_Full_Cost_Per_Route[s], TotalVehicles * sizeof(cost), cudaMemcpyHostToDevice));
		//printf("FirstMemCpyTime: %0.04f ms \n", (float)(clock() - FirstMemCpyTime) / CLOCKS_PER_SEC * 1000);
        */

	// 1) Simultaneous exploration


/*
		clock_t SetupTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{
						//printf("Enter 3.1\n");

			Setup_based_on_ConstantTime_FeasibilityTest(request[s][p], s);
			//printf("allinsertion: %d\n", h_Solution[s].allinsertion);

						//printf("Enter 3.2\n");

		}
		//printf("SetupTime: %0.04f ms ", (float)(clock() - SetupTime) / CLOCKS_PER_SEC * 1000);

		clock_t ContainerCopyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpy(gpu_Container[s], cpu_Container[s], h_Solution[s].allinsertion * sizeof(container), cudaMemcpyHostToDevice));

*/

		//printf("ContainerCopyTime: %0.04f ms ", (float)(clock() - ContainerCopyTime) / CLOCKS_PER_SEC * 1000);
		//printf("False SetupTime: %0.04f ms ", (float)(clock() - SetupTime) / CLOCKS_PER_SEC * 1000);

		clock_t FwCopyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			h_Solution[s].ReBoot();
		CHECK(cudaMemcpy(Dev_Solution, h_Solution, ExpectedScenarios * sizeof(solution), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(Dev_Coefficient, Coefficient, ExpectedScenarios * sizeof(coefficient), cudaMemcpyHostToDevice));
		//printf("FwCpyTime: %0.04f ms ", (float)(clock() - FwCopyTime) / CLOCKS_PER_SEC * 1000);
		CHECK(cudaDeviceSynchronize());

if (CurrentTime > 115) {
    printf("Reached 2.2.3: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

		clock_t PLLInitializationTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			parallel_GPU_reduction_Kernel_first_initialization KERNEL_ARGS4(dim3(round((Expectation) / (32)) + 1), dim3(32), 0, streams[s])(s, Dev_Solution, Dev_Cost_array[s]);
		CHECK(cudaDeviceSynchronize());
		//printf("PLLInitializationTime: %0.04f ms ", (float)(clock() - PLLInitializationTime) / CLOCKS_PER_SEC * 1000);

if (CurrentTime > 115) {
    printf("Reached 2.2.4: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

		clock_t NewMemoizationTime = clock();
		parallelMemoization KERNEL_ARGS2(dim3(1), dim3(ExpectedScenarios))(Dev_Solution, New_Full_Cost_Per_Route);
		CHECK(cudaDeviceSynchronize());
		//printf("NewMemoizationTime: %0.04f ms \n", (float)(clock() - NewMemoizationTime) / CLOCKS_PER_SEC * 1000);


if (CurrentTime > 115) {
    printf("Reached 2.2.5: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

/////
		/*clock_t KernelStTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{	IH_parallel_Setup KERNEL_ARGS4(dim3(round((h_Solution[s].allinsertion) / (32)) + 1), dim3(32), 0, streams[0])(s, request[s][p] + 1, gpu[s], Temp[s], gpu_Container[s], Full_Cost_Per_Route[s], Dev_Solution, Dev_Problem);//, Dev_Cost_array_flat, Dev_Cost_array[s], Dev_Coefficient);
			//CHECK(cudaDeviceSynchronize());
		}
		CHECK(cudaDeviceSynchronize());*/
		//printf("KernelSetupTime: %0.04f ms ", (float)(clock() - KernelStTime) / CLOCKS_PER_SEC * 1000);

		clock_t New_SetupTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			IH_parallel_OneInitialization KERNEL_ARGS4(dim3(round((h_Solution[s].allinsertion) / (32)) + 1), dim3(32), 0, streams[s])(s, gpu[s], Dev_Solution, Dev_Holder);
		CHECK(cudaDeviceSynchronize());
		//printf("New_SetupTime: %0.04f ms ", (float)(clock() - New_SetupTime) / CLOCKS_PER_SEC * 1000);


if (CurrentTime > 115) {
    printf("Reached 2.2.6: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

		///// ROUTINE TO TEST SETUP /////
		/*for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpy(cpu[s], gpu[s], h_Solution[s].allinsertion * sizeof(sol), cudaMemcpyDeviceToHost));

		for (int s = 0; s < ExpectedScenarios; s++)
			IH_parallel_OneInitialization KERNEL_ARGS4(dim3((Expectation/32) + 1), dim3(32), 0, streams[s])(s, gpu[s], Dev_Solution, Dev_Holder);
		CHECK(cudaDeviceSynchronize());

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			printf("For Scenario s=%d\n", s);
			printf("comparison: %d %d\n", h_Solution[s].allinsertion, Dev_Solution[s].allinsertion);

			for (int i = 0; i < h_Solution[s].allinsertion; i++)
			{
				printf("In-depth comparison s=%d: %d(%d), %d(%d), %d(%d)\n", s, cpu[s][i].start, gpu[s][i].start, cpu[s][i].gap, gpu[s][i].gap, cpu[s][i].vehid, gpu[s][i].vehid);
			}
			printf("xxxxxx\n");
		}*/

		/////////////////////////////////

		clock_t KernelTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			//IH_parallel_OneKernel KERNEL_ARGS4(dim3(round((h_Solution[s].allinsertion) / (32)) + 1), dim3(32), 0, streams[0])(s, request[s][p] + 1, gpu[s], Temp[s], gpu_Container[s], Full_Cost_Per_Route[s], Dev_Solution, Dev_Problem, Dev_Cost_array_flat, Dev_Cost_array[s], Dev_Coefficient, d_SkeletonData[s]);
			IH_parallel_OneKernel_New KERNEL_ARGS4(dim3(round((h_Solution[s].allinsertion) / (32)) + 1), dim3(32), 0, streams[s])(s, request[s][p] + 1, gpu[s], Temp[s], gpu_Container[s], New_Full_Cost_Per_Route, Dev_Solution, Dev_Problem, Dev_Cost_array_flat, Dev_Cost_array[s], Dev_Coefficient, d_SkeletonData[s]);

			//IH_parallel_OneKernel_Ver2 KERNEL_ARGS4(dim3(round((h_Solution[s].allinsertion) / (32)) + 1), dim3(32), 0, streams[s])(s, request[s][p] + 1, gpu[s], Temp_flat, gpu_Container[s], Full_Cost_Per_Route[s], &Dev_Solution[s], &Dev_Problem[s]);
			//CHECK(cudaDeviceSynchronize());
		}
		CHECK(cudaDeviceSynchronize());
		//printf("KernelTime: %0.04f ms \n", (float)(clock() - KernelTime) / CLOCKS_PER_SEC * 1000);

		//clock_t Dumm = clock();

		//for (int s = 0; s < ExpectedScenarios; s++)
			//CHECK(cudaMemcpyAsync(Dev_Insertion[s], Insertion[s], n * sizeof(insertion), cudaMemcpyHostToDevice));
		//CHECK(cudaMemcpy(Dev_Insertion[0], Insertion[0], n * sizeof(insertion), cudaMemcpyHostToDevice));
		//printf("%0.0f ms \n", (float)(clock() - Dumm) / CLOCKS_PER_SEC * 1000);
/*		for (int s = 0; s < ExpectedScenarios; s++)
		{
						//printf("Enter 3.3\n");
			DUMMYKERNEL KERNEL_ARGS4(dim3(round((h_Solution[s].allinsertion) / (32)) + 1), dim3(32), 0, streams[s])(p, Temp[s], &Dev_Insertion[s], &Dev_Coefficient[s], &Dev_Solution[s], s);
		}
		CHECK(cudaDeviceSynchronize());
		printf("%0.0f ms ", (float)(clock() - Dumm) / CLOCKS_PER_SEC * 1000);*/




if (CurrentTime > 115) {
    printf("Reached 2.2.7: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}



		clock_t Flat_CopyTime = clock();
		CHECK(cudaMemcpy(h_Temp_flat, Temp_flat, Expectation * ExpectedScenarios * sizeof(cost), cudaMemcpyDeviceToHost));
		//printf("FlatCpy: %0.0f ms \n", (float)(clock() - Flat_CopyTime) / CLOCKS_PER_SEC * 1000);

		clock_t CopyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpy(h_Temp[s], Temp[s], h_Solution[s].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));
		//printf("Cpy: %0.06f ms \n", (float)(clock() - CopyTime) / CLOCKS_PER_SEC * 1000);


if (CurrentTime > 115) {
    printf("Reached 2.2.8: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

/*
		clock_t Copy2Time = clock();
		//for (int s = 0; s < ExpectedScenarios; s++)
			//CHECK(cudaMemcpy(h_Cost_array[s], Dev_Cost_array[s], h_Solution[s].allinsertion * sizeof(float), cudaMemcpyDeviceToHost));
		//printf("Cpy2: %0.06f ms \n\n", (float)(clock() - Copy2Time) / CLOCKS_PER_SEC * 1000);

		clock_t Copy3Time = clock();
		//for (int s = 0; s < ExpectedScenarios; s++)
			//CHECK(cudaMemcpy(h_Cost_array[s], Dev_Cost_array[s], (h_Solution[s].allinsertion/32) + 1 * sizeof(float), cudaMemcpyDeviceToHost));
		//printf("Cpy3: %0.06f ms \n\n", (float)(clock() - Copy3Time) / CLOCKS_PER_SEC * 1000);

		//clock_t Copy4Time = clock();
	//	CHECK(cudaMemcpy(h_Cost_array_flat, Dev_Cost_array_flat, Expectation * ExpectedScenarios * sizeof(float), cudaMemcpyDeviceToHost));
		//printf("Cpy4: %0.06f ms \n\n", (float)(clock() - Copy4Time) / CLOCKS_PER_SEC * 1000);


		clock_t Copy5Time = clock();
		//CHECK(cudaMemcpy(h_Cost_array_flat, Dev_Cost_array_flat, (Expectation/32)+1 * ExpectedScenarios * sizeof(float), cudaMemcpyDeviceToHost));
		//printf("Cpy5: %0.06f ms \n", (float)(clock() - Copy5Time) / CLOCKS_PER_SEC * 1000);

*/




/*
		printf("%0.0f ms \n\n", (float)(clock() - CopyTime) / CLOCKS_PER_SEC * 1000);

		clock_t Copy2Time = clock();

		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpy(h_Temp[s], Temp[s], ((int)(h_Solution[s].allinsertion/32)+1) * sizeof(cost), cudaMemcpyDeviceToHost));

		printf("%0.0f ms \n\n", (float)(clock() - Copy2Time) / CLOCKS_PER_SEC * 1000);
*/

		/*for (int s = 0; s < ExpectedScenarios; s++)
		{
			for (int i = 0; i < Expectation; i++)
				printf("%d %d	%0.02f	%0.02f\n", i, i + (s*Expectation), h_Temp_flat[i + (s*Expectation)].travel_cost, h_Temp[s][i].travel_cost);
			printf("---------\n");
		}
		printf("---------\n");
		printf("---------\n");
		printf("---------\n");
		printf("---------\n");*/

		clock_t ReductionTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{

						//printf("Enter 3.4\n");
			inbuilt_GPU_Reduction(p, Insertion, Coefficient, fsm, choice, 1, gpu[s], h_Temp[s], s, ServedRequests_BoundaryMark);
			//inbuilt_GPU_Reduction_Ver2(p, Insertion, Coefficient, fsm, gpu[s], s);
			//inbuilt_GPU_Reduction_BlockBased_FinalRun(p, Insertion, Coefficient, fsm, choice, 1, gpu[s], h_Temp[s], h_Cost_array[s], s);
			//inbuilt_parallel_GPU_Reduction_SEQ(p, Insertion, gpu[s], s);
		}
		//inbuilt_parallel_GPU_Reduction_PLL(p, Insertion);
	//inbuilt_parallel_GPU_Reduction_PLL_2(p, Insertion, Coefficient);
		//Compare_SEQ_and_PLL_reduction(p, Insertion, Coefficient, fsm, choice, 1);
		//printf("ReductionTime: %0.04f ms \n", (float)(clock() - ReductionTime) / CLOCKS_PER_SEC * 1000);



if (CurrentTime > 115) {
    printf("Reached 2.2.9: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

		clock_t ReInsertionTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{
						//printf("Enter 3.5\n");
						//h_Solution[s].printinput();
						//printf("s:%d, req:%d, k=%d, start=%d, gap=%d\n", s, Insertion[s][p].request, Insertion[s][p].vehid, Insertion[s][p].start, Insertion[s][p].gap);
			h_Solution[s].ReInsertRequest(Insertion[s][p].request, Insertion[s][p].vehid, Insertion[s][p].start, Insertion[s][p].start + Insertion[s][p].gap + 1);
						//h_Solution[s].printinput();

						//printf("Enter 3.6\n");
			h_Solution[s].OfflineEvaluateCost();

						//printf("---\n");
		}
		//printf("ReInsertionTime: %0.04f ms ", (float)(clock() - ReInsertionTime) / CLOCKS_PER_SEC * 1000);
						//printf("Exit 3\n");


if (CurrentTime > 115) {
    printf("Reached 2.2.10: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

		////// Copy back the solution

		clock_t FinalSolutionCopyBackTime = clock();
						//printf("Enter 4\n");
		for (int s = 0; s < ExpectedScenarios; s++)
		{
						//printf("Enter 4.1\n");
			TempSolution[s] = h_Solution[s];
			/*TempSolution[s].Copy(h_Solution[s]);
			TempSolution[s].OfflineEvaluateCost();*/
		}
						//printf("Exit 4\n");

		//printf("FinalSolutionCopyBackTime: %0.04f ms\n-----\n", (float)(clock() - FinalSolutionCopyBackTime) / CLOCKS_PER_SEC * 1000);

if (CurrentTime > 115) {
    printf("Reached 2.2.11: explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest\n");
}

	}

	void explore_full_neighborhood_in_GPU_for_Intra(int &p, solution *TempSolution, problem *TempProblem, int **request, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, int *currvehicleList, int isIntra = 0)
	{


		clock_t InitialSolutionCopy = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{

			h_Problem[s] = TempProblem[s];
			h_Solution[s] = TempSolution[s];
		}
		CHECK(cudaMemcpy(Dev_Problem, h_Problem, ExpectedScenarios * sizeof(problem), cudaMemcpyHostToDevice));


	// 1) Simultaneous exploration

/*
		clock_t SetupTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{
						//printf("Enter 3.1\n");

			Setup_based_on_ConstantTime_FeasibilityTest(request[s][p], s);
			//printf("allinsertion: %d\n", h_Solution[s].allinsertion);

						//printf("Enter 3.2\n");

		}
		//printf("SetupTime: %0.04f ms ", (float)(clock() - SetupTime) / CLOCKS_PER_SEC * 1000);

		clock_t ContainerCopyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpy(gpu_Container[s], cpu_Container[s], h_Solution[s].allinsertion * sizeof(container), cudaMemcpyHostToDevice));

*/

		for (int s = 0; s < ExpectedScenarios; s++)
			h_Solution[s].ReBoot();
		CHECK(cudaMemcpy(Dev_Solution, h_Solution, ExpectedScenarios * sizeof(solution), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(Dev_Coefficient, Coefficient, ExpectedScenarios * sizeof(coefficient), cudaMemcpyHostToDevice));
		CHECK(cudaDeviceSynchronize());

		for (int s = 0; s < ExpectedScenarios; s++)
			parallel_GPU_reduction_Kernel_first_initialization KERNEL_ARGS4(dim3(round((Expectation) / (32)) + 1), dim3(32), 0, streams[s])(s, Dev_Solution, Dev_Cost_array[s]);
		CHECK(cudaDeviceSynchronize());

		parallelMemoization KERNEL_ARGS2(dim3(1), dim3(ExpectedScenarios))(Dev_Solution, New_Full_Cost_Per_Route);
		CHECK(cudaDeviceSynchronize());

		for (int s = 0; s < ExpectedScenarios; s++)
			IH_parallel_OneInitialization KERNEL_ARGS4(dim3(round((h_Solution[s].allinsertion) / (32)) + 1), dim3(32), 0, streams[s])(s, gpu[s], Dev_Solution, Dev_Holder);
		CHECK(cudaDeviceSynchronize());

		for (int s = 0; s < ExpectedScenarios; s++)
			IH_parallel_OneKernel_New KERNEL_ARGS4(dim3(round((h_Solution[s].allinsertion) / (32)) + 1), dim3(32), 0, streams[s])(s, request[s][p] + 1, gpu[s], Temp[s], gpu_Container[s], New_Full_Cost_Per_Route, Dev_Solution, Dev_Problem, Dev_Cost_array_flat, Dev_Cost_array[s], Dev_Coefficient, d_SkeletonData[s]);
		CHECK(cudaDeviceSynchronize());

		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpy(h_Temp[s], Temp[s], h_Solution[s].allinsertion * sizeof(cost), cudaMemcpyDeviceToHost));


		// 3) Reduction

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			intra_reduction_GPU(p, Insertion, Coefficient, fsm, choice, 1, gpu[s], h_Temp[s], s, currvehicleList);
		}

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			h_Solution[s].ReInsertRequest(Insertion[s][p].request, Insertion[s][p].vehid, Insertion[s][p].start, Insertion[s][p].start + Insertion[s][p].gap + 1);
			h_Solution[s].OfflineEvaluateCost();
		}

		////// Copy back the solution
		for (int s = 0; s < ExpectedScenarios; s++)
			TempSolution[s] = h_Solution[s];

	}

	void intra_reduction_GPU(int p, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp, sol *gpu, cost *h_Temp, int s, int *currvehicleList)
	{
		float local_best = 0;
		int local_ID = 0;

		Insertion[s][p].reset();
		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			float currentcost = Coefficient[s].cost_function(h_Temp[idx]);

			if (local_best == 0 || currentcost < local_best)
			{
				local_best = currentcost;
				local_ID = idx;
			}

		}

		sol *bestSol = new sol;
		CHECK(cudaMemcpy(bestSol, &gpu[local_ID], sizeof(sol), cudaMemcpyDeviceToHost));
		Insertion[s][p].Extract(bestSol[0]);
		//printf("[%d] bestID 1: %d, cost=%0.02f\n", s, idx, Coefficient[s].cost_function(bestSol[0].Vehicle.Cost));
		delete bestSol;


	}

	void inbuilt_GPU_Reduction(int p, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp, sol *gpu, cost *h_Temp, int s, myList *ServedRequests_BoundaryMark)
	{

		int found = 0;
		Insertion[s][p].reset();
		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			float currentcost = Coefficient[s].cost_function(h_Temp[idx]);

			if (Insertion[s][p].local_best[h_Temp[idx].vehicle] == 0 || currentcost < Insertion[s][p].local_best[h_Temp[idx].vehicle])
			{
				if (!fsm[0].BarredList[h_Temp[idx].vehid])
				{
					// ENSURES SERVED REQUEST BEING NOT REMOVED!!!
					if (h_Temp[idx].start > (ServedRequests_BoundaryMark[0].list[h_Temp[idx].vehicle] - 1) )
					{
						Insertion[s][p].local_best[h_Temp[idx].vehicle] = currentcost;
						Insertion[s][p].local_best_ID[h_Temp[idx].vehicle] = idx;
						found = 1;
					}
				}
				else
					Insertion[s][p].local_best[h_Temp[idx].vehid] = 100000; // setting high-cost just to avoid barred vehicles*/
			}

		}


if (CurrentTime > 115) {
    printf("Reached 2.2.8.1: inbuilt_GPU_Reduction()\n");
}

		Insertion[s][p].sort_localbest();
		int min_cost_vehicle_k = Insertion[s][p].veh_sequence[0];
		int idx = Insertion[s][p].local_best_ID[min_cost_vehicle_k];
		// calculating gradings
		for (int index = 1; index < m; index++) //beware - upto (m-1) only
			Insertion[s][p].cost_difference[index - 1] = Insertion[s][p].local_best[index] - Insertion[s][p].local_best[0];

if (CurrentTime > 115) {
    printf("Reached 2.2.8.2: inbuilt_GPU_Reduction()\n");
}

		sol *bestSol = new sol;
		CHECK(cudaMemcpy(bestSol, &gpu[idx], sizeof(sol), cudaMemcpyDeviceToHost));
		Insertion[s][p].Extract(bestSol[0]);
		//printf("[%d] bestID 1: %d, cost=%0.02f\n", s, idx, Coefficient[s].cost_function(bestSol[0].Vehicle.Cost));
		delete bestSol;

if (CurrentTime > 115) {
    printf("Reached 2.2.8.3: inbuilt_GPU_Reduction()\n");
}

	}

	void Compare_SEQ_and_PLL_reduction(int p, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp)
	{

		int T = Reduction[0].T;//max(size, 1024);
		int Expected_B = Reduction[0].Expected_B;

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			Reduction[0].B[s] = ( ((h_Solution[s].allinsertion+T-1)/T)); //BEWARE: last block is ignored!!
			//printf("%d %d %d %d\n", Reduction[0].B[s], T, T*Reduction[0].B[s], h_Solution[s].allinsertion);
		}

		clock_t InitialCopyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpyAsync(Reduction[0].d_i_index[s], Reduction[0].h_i_index[s], Expectation * sizeof(int), cudaMemcpyHostToDevice));


		sol *bestSol = new sol;
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			CHECK(cudaMemcpy(h_Cost_array[s], Dev_Cost_array[s], Expectation * sizeof(float), cudaMemcpyDeviceToHost));

			int idx;
			if (0) // 1) SEQ REDUCTION
			{
/*
					int found = 0;
					Insertion[s][p].reset();
					for (int i = 0; i < h_Solution[s].allinsertion; i++)
					{
						float currentcost = Coefficient[s].cost_function(h_Temp[s][i]);
						if (Insertion[s][p].local_best[h_Temp[s][i].vehicle] == 0 || currentcost < Insertion[s][p].local_best[h_Temp[s][i].vehicle])
						{
							if (!fsm[0].BarredList[h_Temp[s][i].vehid])
							{
								Insertion[s][p].local_best[h_Temp[s][i].vehicle] = currentcost;
								Insertion[s][p].local_best_ID[h_Temp[s][i].vehicle] = i;
								found = 1;
							}
							else
								Insertion[s][p].local_best[h_Temp[s][i].vehid] = 100000; // setting high-cost just to avoid barred vehicles
						}
					}

					Insertion[s][p].sort_localbest();
					int min_cost_vehicle_k = Insertion[s][p].veh_sequence[0];
					idx = Insertion[s][p].local_best_ID[min_cost_vehicle_k];
					// calculating gradings
					for (int index = 1; index < m; index++) //beware - upto (m-1) only
						Insertion[s][p].cost_difference[index - 1] = Insertion[s][p].local_best[index] - Insertion[s][p].local_best[0];
					CHECK(cudaMemcpy(bestSol, &gpu[s][idx], sizeof(sol), cudaMemcpyDeviceToHost));
					Insertion[s][p].Extract(bestSol[0]);
					//printf("[%d] bestID 1: %d, cost=%0.02f\n", s, idx, Coefficient[s].cost_function(bestSol[0].Vehicle.Cost));
*/

				float bestcostnow = 0;
				Insertion[s][p].reset();
				for (int i = 0; i < h_Solution[s].allinsertion; i++)
				{
					float currentcost = h_Temp[s][i].overall_current_cost;// h_Cost_array[s][i];//Coefficient[s].cost_function(h_Temp[s][i]);
					if (bestcostnow == 0 || currentcost <= bestcostnow)
					{
							bestcostnow = currentcost;
							idx = i;
					}
				}
				CHECK(cudaMemcpy(bestSol, &gpu[s][idx], sizeof(sol), cudaMemcpyDeviceToHost));
				Insertion[s][p].Extract(bestSol[0]);
				//printf("[%d] bestID 1: %d, cost=%0.02f\n", s, idx, Coefficient[s].cost_function(bestSol[0].Vehicle.Cost));





				/*for (int i = 0; i < h_Solution[s].allinsertion; i++)
				{
					if (h_Temp[s][i].overall_current_cost != h_Cost_array[s][i])
						printf("	index[%d] = %0.0f (%0.0f)\n", i, h_Temp[s][i].overall_current_cost, h_Cost_array[s][i]);
				}
*/

			}
			if (1) // 2) PPLL REDUCTION
			{

					//parallel_GPU_reduction_Kernel_v2 KERNEL_ARGS4(dim3(Reduction[0].B[s]), dim3(T), 0, streams[s])(Dev_Cost_array[s], Reduction[0].d_o_data[s], Reduction[0].d_i_index[s], Reduction[0].d_o_index[s], h_Solution[s].allinsertion);
					parallel_GPU_reduction_Kernel_v2 KERNEL_ARGS4(dim3(Reduction[0].B[s]), dim3(T), 0, streams[s])(Dev_Cost_array[s], Reduction[0].d_o_data[s], Reduction[0].d_i_index[s], Reduction[0].d_o_index[s], Reduction[0].B[s]*T);
					CHECK(cudaDeviceSynchronize());

					/*parallel_GPU_reduction_Kernel_final_reduction KERNEL_ARGS4(dim3(1), dim3(1), sizeof(int)+sizeof(float), streams[s])(s, Reduction[0].d_o_data[s], Reduction[0].d_o_index[s], Reduction[0].B[s], Reduction[0].d_Final_indices);
					CHECK(cudaDeviceSynchronize());*/

					CHECK(cudaMemcpy(Reduction[0].h_Final_indices, Reduction[0].d_Final_indices, ExpectedScenarios * sizeof(int), cudaMemcpyDeviceToHost));




					CHECK(cudaMemcpy(Reduction[0].h_o_data[s], Reduction[0].d_o_data[s], Reduction[0].B[s] * sizeof(float), cudaMemcpyDeviceToHost));
					CHECK(cudaMemcpy(Reduction[0].h_o_index[s], Reduction[0].d_o_index[s], Reduction[0].B[s] * sizeof(int), cudaMemcpyDeviceToHost));
					float min = 0;
					int minIndex = 0;
					for (int i = 0; i < Reduction[0].B[s]; i++)
					{
						if (min == 0 ||  Reduction[0].h_o_data[s][i] < min)
						{
							minIndex =  Reduction[0].h_o_index[s][i];
							min = Reduction[0].h_o_data[s][i];
						}
					}
					Reduction[0].h_Final_indices[s] = minIndex;





					int local_ID = Reduction[0].h_Final_indices[s];
					CHECK(cudaMemcpy(&Reduction[0].bestSol[s], &gpu[s][local_ID], sizeof(sol), cudaMemcpyDeviceToHost));
					Insertion[s][p].reset();
					Insertion[s][p].Extract(Reduction[0].bestSol[s]);
					//printf("Final Reduction value: %0.02f\n", Coefficient[s].cost_function(Reduction[0].bestSol[s].Vehicle.Cost));


					/*if (h_Cost_array[s][idx] != h_Cost_array[s][local_ID])
					{
						printf("%d, %d==%d, %0.0f!=%0.0f\n", s, idx, local_ID, h_Cost_array[s][idx], h_Cost_array[s][local_ID]);
						printf("MISMATCH has happened!!!\n");
					}*/

					/*if (idx != local_ID)
					{
						printf("[%d] %d(%d) (%0.0f(%0.0f)==%0.0f(%0.0f)) (%d > %d)\n", s, idx, local_ID, h_Temp[s][idx].overall_current_cost, h_Cost_array[s][idx], h_Temp[s][local_ID].overall_current_cost, h_Cost_array[s][local_ID], Reduction[0].B[s]*T, h_Solution[s].allinsertion);



						for (int i = 0; i < Reduction[0].B[s]; i++)
							printf("	At index %d, data %0.0f (%0.0f)\n", Reduction[0].h_o_index[s][i], Reduction[0].h_o_data[s][i], h_Cost_array[s][Reduction[0].h_o_index[s][i]]);


						for (int i = 0; i < h_Solution[s].allinsertion; i++)
						{
							//if (h_Temp[s][i].overall_current_cost != h_Cost_array[s][i])
								printf("	index[%d] = %0.0f (%0.0f)\n", i, h_Temp[s][i].overall_current_cost, h_Cost_array[s][i]);
						}


						printf("	--MinIndex %d, min %0.0f\n", minIndex, min);
						printf("-------------\n");

						for (int i = 0; i < Reduction[0].B[s]; i++)
							printf("--inner: index[%d] = %0.0f\n", Reduction[0].h_o_index[s][i], Reduction[0].h_o_data[s][i]);
						printf("---xxxxxxxxxxx------\n");
					}*/
			}

		}
		delete bestSol;

	}

	void inbuilt_GPU_Reduction_Ver2(int p, insertion **Insertion, coefficient *Coefficient, FSM *fsm, sol *gpu, int s)
	{

		int found = 0;
		Insertion[s][p].reset();
		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			float currentcost = Coefficient[s].cost_function(h_Temp_flat[idx + (s*Expectation)]);

			if (Insertion[s][p].local_best[h_Temp_flat[idx + (s*Expectation)].vehicle] == 0 || currentcost < Insertion[s][p].local_best[h_Temp_flat[idx + (s*Expectation)].vehicle])
			{
				if (!fsm[0].BarredList[h_Temp_flat[idx + (s*Expectation)].vehid])
				{
					Insertion[s][p].local_best[h_Temp_flat[idx + (s*Expectation)].vehicle] = currentcost;
					Insertion[s][p].local_best_ID[h_Temp_flat[idx + (s*Expectation)].vehicle] = idx;
					found = 1;
				}
				else
					Insertion[s][p].local_best[h_Temp_flat[idx + (s*Expectation)].vehid] = 100000; // setting high-cost just to avoid barred vehicles
			}
		}

		Insertion[s][p].sort_localbest();

		int min_cost_vehicle_k = Insertion[s][p].veh_sequence[0];
		int idx = Insertion[s][p].local_best_ID[min_cost_vehicle_k];
		for (int index = 1; index < m; index++) //beware - upto (m-1) only
			Insertion[s][p].cost_difference[index - 1] = Insertion[s][p].local_best[index] - Insertion[s][p].local_best[0];

		// 2) Extraction

		sol *bestSol = new sol;
		CHECK(cudaMemcpy(bestSol, &gpu[idx], sizeof(sol), cudaMemcpyDeviceToHost));
		Insertion[s][p].Extract(bestSol[0]);
		delete bestSol;


	}

	void inbuilt_GPU_Reduction_BlockBased_FinalRun(int p, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, int PickUp, sol *gpu, cost *h_Temp, float *h_Cost_array, int s)
	{
		int local_ID = 0;
		float localbest = 0;

		Insertion[s][p].reset();

		for (int idx = 0; idx < h_Solution[s].allinsertion; idx++)
		{
			//float overall_current_cost = h_Temp[idx].overall_current_cost;
			//float overall_current_cost = h_Cost_array[idx];
			float overall_current_cost = h_Cost_array_flat[idx +(s * Expectation)];

			//printf("%0.0f %0.0f %0.0f\n", h_Temp[idx].overall_current_cost, h_Cost_array[idx], h_Cost_array_flat[idx +(s * Expectation)]);

			if (abs(h_Temp[idx].overall_current_cost - h_Cost_array[idx]) > 1)
				printf("s=%d DOES NOT MATCH!!!!!\n", s);

			if (abs(h_Temp[idx].overall_current_cost - h_Cost_array_flat[idx +(s * Expectation)]) > 1)
				printf("s=%d DOES NOT MATCH!!!!!\n", s);

			if (localbest == 0 || overall_current_cost < localbest)
			{
				localbest = overall_current_cost;
				local_ID = idx;
			}
		}

		sol *bestSol = new sol;
		CHECK(cudaMemcpy(bestSol, &gpu[local_ID], sizeof(sol), cudaMemcpyDeviceToHost));
		Insertion[s][p].Extract(bestSol[0]);
		delete bestSol;


	}

	void inbuilt_parallel_GPU_Reduction_SEQ(int p, insertion **Insertion, sol *gpu, int s)
	{

		sol *bestSol;
		int *h_Idx, *Dev_Idx;

		bestSol = new sol;
		h_Idx = new int;
		CHECK(cudaMalloc((void**)&Dev_Idx, sizeof(int)));

		parallel_GPU_reduction_Kernel KERNEL_ARGS4(dim3(1), dim3(1), 0, streams[s])(s, Dev_Idx, Dev_Cost_array_flat, Temp[s], Dev_Cost_array[s], &Dev_Solution[s]);
		CHECK(cudaDeviceSynchronize());

		CHECK(cudaMemcpy(h_Idx, Dev_Idx, sizeof(int), cudaMemcpyDeviceToHost));

		int local_ID = h_Idx[0];
		CHECK(cudaMemcpy(bestSol, &gpu[local_ID], sizeof(sol), cudaMemcpyDeviceToHost));
		Insertion[s][p].reset();
		Insertion[s][p].Extract(bestSol[0]);

		delete bestSol;
		delete h_Idx;
		CHECK(cudaFree(Dev_Idx));
	}

	void inbuilt_parallel_GPU_Reduction_PLL(int p, insertion **Insertion)
	{

		int T = Reduction[0].T;//max(size, 1024);
		int Expected_B = Reduction[0].Expected_B;



		for (int s = 0; s < ExpectedScenarios; s++)
			Reduction[0].B[s] = ((h_Solution[s].allinsertion+T-1)/T) - 1; //BEWARE: last block is ignored!!

		clock_t InitialCopyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpyAsync(Reduction[0].d_i_index[s], Reduction[0].h_i_index[s], Expectation * sizeof(int), cudaMemcpyHostToDevice));
		//printf("InitialCopyTime: %0.04f ms \n", (float)(clock() - InitialCopyTime) / CLOCKS_PER_SEC * 1000);


		/*for (int i = 0; i < h_Solution[0].allinsertion; i++)
		{
			printf("h_i_index[%d]=%d, CostArray=%0.01f\n", i , Reduction[0].h_i_index[0][i], h_Cost_array[0][i]);
			if (i > 0 && i % T == 0)
				printf("----\n");
		}*/

		clock_t InnerKernelTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			parallel_GPU_reduction_Kernel_v2 KERNEL_ARGS4(dim3(Reduction[0].B[s]), dim3(T), 0, streams[s])(Dev_Cost_array[s], Reduction[0].d_o_data[s], Reduction[0].d_i_index[s], Reduction[0].d_o_index[s], h_Solution[s].allinsertion);
		CHECK(cudaDeviceSynchronize());
		//printf("InnerKernelTime: %0.04f ms \n", (float)(clock() - InnerKernelTime) / CLOCKS_PER_SEC * 1000);

		clock_t InnerCopyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			CHECK(cudaMemcpy(Reduction[0].h_o_data[s], Reduction[0].d_o_data[s], Reduction[0].B[s] * sizeof(float), cudaMemcpyDeviceToHost));
			CHECK(cudaMemcpy(Reduction[0].h_o_index[s], Reduction[0].d_o_index[s], Reduction[0].B[s] * sizeof(int), cudaMemcpyDeviceToHost));
		}
		//printf("InnerCopyTime: %0.04f ms \n", (float)(clock() - InnerCopyTime) / CLOCKS_PER_SEC * 1000);


		/*
		for (int i = 0; i < Reduction[0].B[0]; i++)
			printf("h_o_index[%d]=%d, h_o_data=%0.01f\n", i , Reduction[0].h_o_index[0][i], Reduction[0].h_o_data[0][i]);
		*/

		clock_t FinalRunTime = clock();
		float *min = new float[ExpectedScenarios];
		int *minIndex = new int[ExpectedScenarios];
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			min[s] = 0; minIndex[s] = 0;
			for (int i = 0; i < Reduction[0].B[s]; i++)
			{
				if (min[s] == 0 || Reduction[0].h_o_data[s][i] < min[s])
				{
					minIndex[s] = Reduction[0].h_o_index[s][i];
					min[s] = Reduction[0].h_o_data[s][i];
				}
			}

			//printf("min: %0.01f at minIndex %d\n", min, minIndex);
		}
//		printf("FinalRunTime: %0.04f ms \n", (float)(clock() - FinalRunTime) / CLOCKS_PER_SEC * 1000);
		for (int s = 0; s < ExpectedScenarios; s++)
		{
			int local_ID = minIndex[s];
			//printf("local_ID: %d\n", local_ID);
			CHECK(cudaMemcpy(Reduction[0].bestSol, &gpu[s][local_ID], sizeof(sol), cudaMemcpyDeviceToHost));
			Insertion[s][p].reset();
			Insertion[s][p].Extract(Reduction[0].bestSol[0]);
		}
		delete[] min;
		delete[] minIndex;
//		printf("FinalRunTime: %0.04f ms \n", (float)(clock() - FinalRunTime) / CLOCKS_PER_SEC * 1000);


	}

	void inbuilt_parallel_GPU_Reduction_PLL_2(int p, insertion **Insertion, coefficient *Coefficient)
	{

		int T = Reduction[0].T;//max(size, 1024);
		int Expected_B = Reduction[0].Expected_B;

		for (int s = 0; s < ExpectedScenarios; s++)
			Reduction[0].B[s] = ((h_Solution[s].allinsertion+T-1)/T); //BEWARE: last block is ignored!!

		clock_t InitialCopyTime = clock();
		for (int s = 0; s < ExpectedScenarios; s++)
			CHECK(cudaMemcpyAsync(Reduction[0].d_i_index[s], Reduction[0].h_i_index[s], Expectation * sizeof(int), cudaMemcpyHostToDevice));


		clock_t RedInnerTime = clock();


		/*for (int s = 0; s < ExpectedScenarios; s++)
			parallel_GPU_reduction_Kernel_first_initialization KERNEL_ARGS4(dim3(Reduction[0].B[s]), dim3(T), 0, streams[s])(Reduction[0].B[s], T, Dev_Cost_array[s], h_Solution[s].allinsertion);
		CHECK(cudaDeviceSynchronize());*/


		for (int s = 0; s < ExpectedScenarios; s++)
		{
			parallel_GPU_reduction_Kernel_v2 KERNEL_ARGS4(dim3(Reduction[0].B[s]), dim3(T), 0, streams[s])(Dev_Cost_array[s], Reduction[0].d_o_data[s], Reduction[0].d_i_index[s], Reduction[0].d_o_index[s], h_Solution[s].allinsertion);
			//parallel_GPU_reduction_Kernel_v3 KERNEL_ARGS4(dim3(Reduction[0].B[s]), dim3(T), 0, streams[s])(Dev_Cost_array[s], Reduction[0].d_o_data[s], Reduction[0].d_i_index[s], Reduction[0].d_o_index[s], h_Solution[s].allinsertion);
			//parallel_GPU_reduction_Kernel_v4 KERNEL_ARGS4(dim3(Reduction[0].B[s]), dim3(T), 0, streams[s])(Dev_Cost_array[s], Reduction[0].d_o_data[s], Reduction[0].d_i_index[s], Reduction[0].d_o_index[s], h_Solution[s].allinsertion);
			//parallel_GPU_reduction_Kernel_v5 KERNEL_ARGS4(dim3(Reduction[0].B[s]), dim3(T), (T*sizeof(int))+(T*sizeof(float)), streams[s])(Dev_Cost_array[s], Reduction[0].d_o_data[s], Reduction[0].d_i_index[s], Reduction[0].d_o_index[s], h_Solution[s].allinsertion);
		}
		CHECK(cudaDeviceSynchronize());
		//printf("RedInnerTime: %0.04f ms \n", (float)(clock() - RedInnerTime) / CLOCKS_PER_SEC * 1000);


		for (int s = 0; s < ExpectedScenarios; s++)
		{

			//CHECK(cudaDeviceSynchronize());

			parallel_GPU_reduction_Kernel_final_reduction KERNEL_ARGS4(dim3(1), dim3(1), sizeof(int)+sizeof(float), streams[s])(s, Reduction[0].d_o_data[s], Reduction[0].d_o_index[s], Reduction[0].B[s], Reduction[0].d_Final_indices);
		}
		CHECK(cudaDeviceSynchronize());


		CHECK(cudaMemcpy(Reduction[0].h_Final_indices, Reduction[0].d_Final_indices, ExpectedScenarios * sizeof(int), cudaMemcpyDeviceToHost));

		/*for (int s = 0; s < ExpectedScenarios; s++)
		{
			int local_ID = Reduction[0].h_Final_indices[s];
			printf("[%d] local_ID: %d ", s, local_ID);
			CHECK(cudaMemcpy(&Reduction[0].bestSol[s], &gpu[s][local_ID], sizeof(sol), cudaMemcpyDeviceToHost));

			Insertion[s][p].reset();
			Insertion[s][p].Extract(Reduction[0].bestSol[s]);
			printf("Final Reduction value: %0.02f\n", Coefficient[s].cost_function(Reduction[0].bestSol[s].Vehicle.Cost));
		}*/

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			int local_ID = Reduction[0].h_Final_indices[s];
			//printf("local_ID: %d\n---\n", local_ID);
			CHECK(cudaMemcpy(&h_SkeletonData[s][local_ID], &d_SkeletonData[s][local_ID], sizeof(skeleton), cudaMemcpyDeviceToHost));
			Insertion[s][p].Extract_2(h_SkeletonData[s][local_ID]);
		}



	/*	for (int s = 0; s < ExpectedScenarios; s++)
			printf("[%d] (B,T)=(%d,%d) => %d, allinsertion: %d\n---\n", s, Reduction[0].B[s], T, Reduction[0].B[s]*T, h_Solution[s].allinsertion);

		for (int s = 0; s < ExpectedScenarios; s++)
		{
			for (int i = 0; i < Reduction[0].B[s]; i++)
			{
				printf("data[%d]=%0.02f, index=%d\n", i, Reduction[0].d_o_data[s][i], Reduction[0].d_o_index[s][i]);
			}
			printf("Final indice: %d (%d)\n", Reduction[0].d_Final_indices[s], Reduction[0].h_Final_indices[s]);
		}
		printf("xxxxxxxxxxxxxxxxxxxxxxxxx\n");*/

	}

	//////////////////////////////////////////////////
	////// EXPLORATION DECISION MAKER FUNCTION ///////
	//////////////////////////////////////////////////

	void Exploration_Decision_Maker(int isGPU, int isFull, int isRev, int &p, solution *TempSolution, problem *TempProblem, int **request, insertion **Insertion, coefficient *Coefficient, FSM *fsm, int choice, int *currvehicleList, myList *ServedRequests_BoundaryMark, int intraRouteInsertion = 0)
	{
		isIntra = 0;
		if (intraRouteInsertion)
		{
			//clock_t IntraExplorationTime = clock();
			isIntra = 1;
			if (!isGPU)
				explore_full_neighborhood(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice, currvehicleList, 1);
			else
				explore_full_neighborhood_in_GPU_for_Intra(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice, currvehicleList, 1);
			//printf("IntraExplorationTime: %0.02f ms\n-----\n", (float)(clock() - IntraExplorationTime) / CLOCKS_PER_SEC * 1000);

		}
		else if (isGPU) //Parallel
		{
			explore_full_neighborhood_in_GPU_with_ConstTime_feasibilityTest(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice, ServedRequests_BoundaryMark);
		}
		else //Sequential
		{
			explore_full_neighborhood(p, TempSolution, TempProblem, request, Insertion, Coefficient, fsm, choice, currvehicleList);

		}
	}

};
