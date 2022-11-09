
struct insertion
{
	// for regret insertion
	float *cost_difference; // 0 - diff between min_best and second best, 1 - diff between min_best and third best, and so on for (m-1) times...

	float *local_best; // consist of best insertion cost for each vehicle for a req i
	int *veh_sequence; // consist of vehicle sequence (like indexing) for *local_best

	int *local_best_ID; // consist of local_ID for min cost for each vehicle for a req i

						// others
	int local_ID;

	int solID;
	cost Cost;
	float current_cost;

	bool isFeasible;

	int request;
	int vehid;

	int start;
	int gap;

	__host__ __device__ insertion()
	{
		cost_difference = new float[m - 1]; //BEWARE, upto m-1 only

		local_best_ID = new int[m];
		local_best = new float[m];
		veh_sequence = new int[m];

		for (int k = 0; k < m; k++)
			local_best[k] = 0;

		solID = 0;
		Cost.reset();
		current_cost = 0;

		request = 0;
		vehid = 0;

		start = 0;
		gap = 0;
	}

	__host__ __device__ ~insertion()
	{
		delete[] cost_difference;

		delete[] local_best_ID;
		delete[] local_best;
		delete[] veh_sequence;
	}

	__host__ __device__
		void reset()
	{
		for (int k = 0; k < m; k++)
		{
			local_best[k] = 0;
		}

		solID = 0;
		Cost.reset();
		current_cost = 0;

		request = 0;
		vehid = 0;

		start = 0;
		gap = 0;
	}

	__host__ __device__
		void Copy(insertion &inInsertion)
	{
		solID = inInsertion.solID;
		Cost.Copy(inInsertion.Cost);
		current_cost = inInsertion.current_cost;

		isFeasible = inInsertion.isFeasible;

		request = inInsertion.request;
		vehid = inInsertion.vehid;

		start = inInsertion.start;
		gap = inInsertion.gap;
	}

	__host__ __device__
		void print()
	{
		printf("	request, vehid, start, gap: %d	%d	%d	%d\n", request, vehid, start, gap);
		printf("	cost: %0.02f	%s\n", Cost.travel_cost, Cost.isFeasible == 0 ? "NOT Feasible" : "Feasible");
		Cost.print();
		printf("	current_cost: %f\n", current_cost);
		printf("----------------\n");
	}

	__host__ __device__
		void print_short()
	{
		printf("	request, vehid, start, gap: %d	%d	%d	%d\n", request, vehid, start, gap);
	}

	__host__ __device__
		void Extract(sol &inSol)
	{
		// selective copies. ///

		solID = inSol.solID;
		Cost.Copy(inSol.Vehicle.Cost);
		current_cost = inSol.current_cost;

		isFeasible = inSol.isFeasible;

		request = inSol.request;
		vehid = inSol.vehid;

		start = inSol.start;
		gap = inSol.gap;
	}

	__host__ __device__
		void Extract_2(skeleton &inSkeletonData)
	{
		request = inSkeletonData.req;
		vehid = inSkeletonData.vehid;

		start = inSkeletonData.start;
		gap = inSkeletonData.gap;
	}

	__host__ __device__
		void sort_localbest()
	{
		// ascending sosrt

		for (int k = 0; k < m; k++)
			veh_sequence[k] = k;

		float temp_cost;
		float temp_veh;
		for (int i = 0; i<m; i++)
		{
			for (int j = i + 1; j<m; j++)
			{
				if (local_best[i] > local_best[j])
				{
					temp_cost = local_best[i];
					local_best[i] = local_best[j];
					local_best[j] = temp_cost;

					temp_veh = veh_sequence[i];
					veh_sequence[i] = veh_sequence[j];
					veh_sequence[j] = temp_veh;
				}
			}
		}
	}
};

struct solution
{

	vehicle_detail Vehicle[TotalVehicles];
	request_detail Request[2 * TotalRequests];
	problem Problem;
	cost Cost;
	cost Temp;
	int allinsertion;

	int Feasible_Req_Count;
	int total_req_served;

	__host__ __device__
		solution() {}

	__host__ __device__
		void print()
	{
		printf("SOLUTION CONSIST OF:\n");

		float TotalTravelCost = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			Vehicle[k].StartingPoint = 2 * n + k;
			Vehicle[k].EndingPoint = 2 * (n)+m + k;

			printf("VEHICLE %d:  Startime %f:	Endtime: %f:	RouteDuration: %f	TravelCost: %f\n", k, Vehicle[k].start_time, Vehicle[k].end_time, Vehicle[k].end_time - Vehicle[k].start_time, Vehicle[k].Cost.travel_cost);
			printf("Request	Load	ArrivalTime	StartofService	EndofService	DepartureTime	WaitingTime	RideTime	TravelCost\n");
			printf("%d	%d	%f	%f	%f	%f	%f	%f	%f\n", Vehicle[k].StartingPoint, 0, 0, 0, 0, 0, 0, 0, 0);

			//printf("TravelTime[%d][%d] = %f\n", Vehicle[k].StartingPoint, Vehicle[k].path[0]-1, Problem.TravelTime[Vehicle[k].StartingPoint][Vehicle[k].path[0] - 1]);
			float TravelCost = 0;

			int depart_location = 2 * n + k;
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				//if (Vehicle[k].path[j] <= n / 2)
				{
					int reqID = Vehicle[k].path[j];
					int node = reqID - 1;
					TravelCost += Problem.TravelTime[depart_location][node];
					Request[reqID - 1].D = Request[reqID - 1].C;
					printf("%d	%d	%f	%f	%f	%f	%f	%f	%f\n", reqID, Request[reqID - 1].p, Request[reqID - 1].A, Request[reqID - 1].B, Request[reqID - 1].C, Request[reqID - 1].D, Request[reqID - 1].W, Request[reqID - 1].L, TravelCost);// Problem.TravelTime[depart_location][node]);
																																																												//printf("TravelTime[%d][%d] = %f\n", depart_location, node, Problem.TravelTime[depart_location][node]);
					depart_location = node;


				}
			}
			TravelCost += Problem.TravelTime[depart_location][Vehicle[k].EndingPoint];
			printf("%d	%d	%f	%f	%f	%f	%f	%f	%f\n", Vehicle[k].EndingPoint, 0, 0, 0, 0, 0, 0, 0, TravelCost);// Problem.TravelTime[depart_location][Vehicle[k].EndingPoint]);

																													//printf("TravelTime[%d][%d] = %f\n", depart_location, 2*n+m+k, Problem.TravelTime[depart_location][2 * n + m + k]);

			printf("Travelcost:	%f\n", TravelCost);
			printf("\n\n");

			if (Vehicle[k].size == 0)
			{
				TotalTravelCost += 0;
			}
			else
				TotalTravelCost += Vehicle[k].Cost.travel_cost;
		}

		printf("Total Travel Cost: %f\n", TotalTravelCost);

		// Find Waiting Time
		float TotalWaitingTime = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			for (int j = 0; j < Vehicle[k].size; j++)
				TotalWaitingTime += Request[Vehicle[k].path[j] - 1].W;
		}

		printf("Total Vehicle Waiting Time: %f\n", TotalWaitingTime);

		// Find Ride Time
		float TotalRideTime = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			for (int j = 0; j < Vehicle[k].size; j++)
				TotalRideTime += Request[Vehicle[k].path[j] - 1].L;
		}

		printf("Total Passenger Ride Time: %f\n", TotalRideTime);
		printf("Total Feasible Req_Count: %d\n", Feasible_Req_Count);

		int total_requests = 0;
		for (int k = 0; k < m; k++)
			total_requests += (Vehicle[k].size / 2);
		printf("Total Req_Count: %d\n", total_requests);


		printf("\n\n");
	}

	__host__ __device__
		void print_with_battery_level()
	{


		for (int k = 0; k < m; k++)
		{
			//	Offline_battery_evaluation(k);
		}

		printf("SOLUTION CONSIST OF:\n");
		printf("Note:all nodes including depots have +1 added.\n");
		float TotalTravelCost = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{

			printf("VEHICLE %d:  Startime %0.02f:	Endtime: %0.02f:	RouteDuration: %0.02f	TravelCost: %0.02f\n",
				k, Vehicle[k].start_time, Vehicle[k].end_time, Vehicle[k].end_time - Vehicle[k].start_time, Vehicle[k].Cost.travel_cost);

		//	printf("Request	Load	ArrivalTime	StartofService	EndofService	DepartureTime	WaitingTime	RideTime	TravelCost	SATC	Battery_level	ExcessRideTime	ETW	LTW	TwV	RTV\n");
			printf("req	d	A	B	C	D	W	RT	T.T.	SATC	Btry	ERT	ETW	LTW	TwV	RTV\n");

			printf("%d	0	0.00	0.00	0.00	0.00	0.00	0.00	0.00	%0.02f	0.00	0.00	0.00	0.00	0.00\n",
				Vehicle[k].StartingPoint + 1, Vehicle[k].start_time);

			//printf("TravelTime[%d][%d] = %f\n", Vehicle[k].StartingPoint, Vehicle[k].path[0]-1, Problem.TravelTime[Vehicle[k].StartingPoint][Vehicle[k].path[0] - 1]);
			float TravelCost = 0;
			float CurrentTravelCost = TravelCost;

			int depart_location = Vehicle[k].StartingPoint; // 2*n + k;
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				CurrentTravelCost = TravelCost;

				//if (Vehicle[k].path[j] <= n / 2)
				{
					int reqID = Vehicle[k].path[j];
					int node = reqID - 1;

					if (j == 0)
					{
						TravelCost += Problem.TravelTime[depart_location][node];
					}
					else if (j > 0)
					{
						TravelCost += Problem.TravelTime[depart_location][node];
					}


					Request[reqID - 1].D = Request[reqID - 1].C;
					if ((reqID - 1) >= n)
					{
						printf("%d	%d	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f\n",
							reqID, Request[reqID - 1].p, Request[reqID - 1].A, Request[reqID - 1].B, Request[reqID - 1].C, Request[reqID - 1].D, Request[reqID - 1].W, Request[reqID - 1].L, TravelCost, TravelCost - CurrentTravelCost, Request[reqID - 1].L - Problem.TravelTime[reqID - 1 - n][reqID - 1], Problem.Request[reqID - 1 - n].dropoff.EarliestTime, Problem.Request[reqID - 1 - n].dropoff.LatestTime, Request[reqID - 1].time_window_violation, Request[reqID - 1].ride_time_violation);
					}
					else
					{
						printf("%d	%d	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	0.00	%0.02f	%0.02f	%0.02f	%0.02f\n",
							reqID, Request[reqID - 1].p, Request[reqID - 1].A, Request[reqID - 1].B, Request[reqID - 1].C, Request[reqID - 1].D, Request[reqID - 1].W, Request[reqID - 1].L, TravelCost, TravelCost - CurrentTravelCost, Problem.Request[reqID - 1].pickup.EarliestTime, Problem.Request[reqID - 1].pickup.LatestTime, Request[reqID - 1].time_window_violation, Request[reqID - 1].ride_time_violation);
					}


					//printf("Problem.TravelTime[%d][%d] = %f\n", depart_location, node, Problem.TravelTime[depart_location][node]);

					depart_location = node;
				}
			}
			CurrentTravelCost = TravelCost;

			int arrival_location = Vehicle[k].EndingPoint; // 2 * n + m + k;
			TravelCost += Problem.TravelTime[depart_location][arrival_location];
			printf("%d	%d	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f	%0.02f\n",
				Vehicle[k].EndingPoint + 1, 0, Vehicle[k].end_time, 0, 0, 0, 0, 0, TravelCost, TravelCost - CurrentTravelCost, 0, 0, 0, 0, 0);

			printf("Travelcost:	%f\n", TravelCost);
			printf("Actual Travelcost:	%f\n", Vehicle[k].Cost.travel_cost);
			printf("Actual BatteryViol:	%f\n", Vehicle[k].Cost.batteryViolation);
			printf("\n\n");

			if (Vehicle[k].size == 0)
			{
				TotalTravelCost += 0;
			}
			else
				TotalTravelCost += TravelCost;// Vehicle[k].Cost.travel_cost;
		}

		for (int k = 0; k < TotalVehicles; k++)
			printf("Vehicle[%d].Cost.travel_cost = %0.02f\n", k, Vehicle[k].Cost.travel_cost);
		printf("Object.Cost.travel_cost = %0.02f\n", Cost.travel_cost);
		printf("\nTotal Travel Cost: %f\n", TotalTravelCost);
		for (int k = 0; k < TotalVehicles; k++)
			printf("Vehicle[%d].Cost.batteryViol = %0.02f\n", k, Vehicle[k].Cost.batteryViolation);
		printf("Object.Cost.batteryViolation = %0.02f\n", Cost.batteryViolation);
		for (int k = 0; k < TotalVehicles; k++)
			printf("Vehicle[%d].Cost.max_battery_Violation = %0.02f\n", k, Vehicle[k].Cost.max_battery_Violation);
		printf("Object.Cost.max_battery_Violation = %0.02f\n", Cost.max_battery_Violation);

		// Find Waiting Time
		float TotalWaitingTime = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			for (int j = 0; j < Vehicle[k].size; j++)
				TotalWaitingTime += Request[Vehicle[k].path[j] - 1].W;
		}

		printf("Total Vehicle Waiting Time: %f\n", TotalWaitingTime);

		// Find Ride Time
		float TotalRideTime = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			for (int j = 0; j < Vehicle[k].size; j++)
				TotalRideTime += Request[Vehicle[k].path[j] - 1].L;
		}

		printf("Total Passenger Ride Time: %f\n", TotalRideTime);
		printf("Total Feasible Req_Count: %d\n", Feasible_Req_Count);

		int total_requests = 0;
		for (int k = 0; k < m; k++)
			total_requests += (Vehicle[k].size / 2);
		printf("Total Req_Count: %d\n", total_requests);

		float final_excessrideTime = 0;
		for (int k = 0; k < m; k++)
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				int node = Vehicle[k].path[j] - 1;
				if (node >= n)
					final_excessrideTime += Request[node].L - Problem.TravelTime[node - n][node];

			}
		printf("Total ExcessRideTime (Actual):	%f\n", final_excessrideTime);

		printf("==> ultimate cost: %0.02f", w1*Cost.travel_cost + w2*final_excessrideTime);
		printf("\n\n");


		/*for (int i = 0; i < 2 * n + 2*m; i++)
		for (int j = i; j < 2 * n + 2 * m; j++)
		{
		printf("TravelTime[%d][%d]=	%f\n", i, j, Problem.TravelTime[i][j]);
		}*/

		/*printf("Nodes	Earliest	Latest\n");
		for (int i = 0; i < 2 * n; i++)
		{
		if (i>=n)
		printf("%d	%f	%f\n", i, Problem.Request[i - n].dropoff.EarliestTime, Problem.Request[i - n].dropoff.LatestTime);
		else
		printf("%d	%f	%f\n", i, Problem.Request[i].pickup.EarliestTime, Problem.Request[i].pickup.LatestTime);
		}*/

	}

	__host__ __device__
		void printcost()
	{
		printf("isFeasible %s\n", Cost.isFeasible == 0 ? "false" : "true");
		printf("travelcost %f (%s)\n", Cost.travel_cost, Cost.travel_cost == 0 ? "true" : "false");
		printf("timeWindow %f (%s) \n", Cost.time_window, Cost.time_window == 0 ? "true" : "false");
		printf("rideTimeViol %f (%s) \n", Cost.ride_time, Cost.ride_time == 0 ? "true" : "false");
		printf("loadViol %d (%s) \n", Cost.load, Cost.load == 0 ? "true" : "false");
		printf("durationViol %f (%s) \n", Cost.duration, Cost.duration == 0 ? "true" : "false");
	}

	__host__ __device__
		void Copy(solution &inSolution)
	{
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].Copy(inSolution.Vehicle[k]);

		for (int i = 0; i < 2 * TotalRequests; i++)
			Request[i].Copy(inSolution.Request[i]);

		Cost.Copy(inSolution.Cost);

		Problem.Copy(inSolution.Problem);

		allinsertion = inSolution.allinsertion;

	}

	__host__ __device__
		void fetchinput()
	{

		//Dummy Variables
		int totalvehicle, totreq;
		printf("Enter the totalvehicles:\n");
		cin >> totalvehicle >> totreq;

		printf("Enter the vehiclepaths:\n");

		//Set VehiclePaths
		for (int k = 0; k < TotalVehicles; k++)
		{
			Vehicle[k].id = k;
			 printf("sz:");
			int sz; cin >> sz;
			Vehicle[k].initialize(sz);
			cout << " path: ";
			if (sz > 0)
			{

				for (int j = 0; j < sz; j++)
				{
					int val; cin >> val;
					Vehicle[k].path[j] = val;
				}
			}
			cout << "\n";
		}

		//Set RequestIds
		for (int i = 0; i<TotalRequests; i++)
			Request[i].id = i + 1;
	}

	__host__ __device__
		void print_inputstyle()
	{

		printf("\n-------------------\n");

		//Dummy Variables
		printf("%d\n%d\n\n", n, m);


		//Set VehiclePaths
		for (int k = 0; k < m; k++)
		{
			printf("%d\n\n", Vehicle[k].size);
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				printf("%d\n", Vehicle[k].path[j]);
			}
			printf("\n");
		}

		printf("\n-------------------\n");


	}

	__host__ __device__
		void generateInput()
	{
		int a[3] = { 12, 12, 24 };

		int val = 1;
		//Set VehiclePaths
		for (int k = 0; k < TotalVehicles; k++)
		{

			Vehicle[k].id = k;

			int sz = a[k];// 2 * TotalRequests / TotalVehicles;
			Vehicle[k].initialize(sz);

			for (int j = 0; j < sz; j += 2, val++)
			{
				Vehicle[k].path[j] = val;
				Vehicle[k].path[j + 1] = val + TotalRequests;
			}

			if (val == TotalRequests)
				break;
		}

		//Set RequestIds
		for (int i = 0; i<TotalRequests; i++)
			Request[i].id = i + 1;
	}

	__host__ __device__
		void randomizeInput()
	{


		int a[3] = { 12, 12, 24 };

		int val = 1;
		//Set VehiclePaths
		for (int k = 0; k < TotalVehicles; k++)
		{

			Vehicle[k].id = k;

			int sz = a[k];// 2 * TotalRequests / TotalVehicles;
			Vehicle[k].initialize(sz);

			for (int j = 0; j < sz; j += 2, val++)
			{
				Vehicle[k].path[j] = val;
				Vehicle[k].path[j + 1] = val + TotalRequests;
			}

			if (val == TotalRequests)
				break;
		}

		//Set RequestIds
		for (int i = 0; i<TotalRequests; i++)
			Request[i].id = i + 1;
	}

	__host__ __device__
		void Boot()
	{
		for (int k = 0; k < TotalVehicles; k++)
		{
			Vehicle[k].id = k + 1;
			Vehicle[k].initialize(0);
			Vehicle[k].path[0] = 0;
		}

		for (int i = 0; i < TotalRequests; i++)
		{
			Request[i].id = i + 1;
			allinsertion = 1;
		}
	}

	__host__ __device__
		void ReBoot()
	{

		//1) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		Vehicle[0].cumulativeInsertions = 0;
		int cumulativeInsertions = 0;
		for (int k = 1; k < TotalVehicles; k++)
		{
			Vehicle[k].cumulativeInsertions = cumulativeInsertions + Vehicle[k-1].insertways;
			cumulativeInsertions += Vehicle[k-1].insertways;
		}

		//2) Update allinsertion
		int sum = 0;
		for (int k = 0; k < TotalVehicles; k++)
			sum += Vehicle[k].insertways;

		allinsertion = sum;

	}

/*	__host__ __device__
		bool isComplete()
	{
		bool Bool = false;

		int count = 0;
		for (int k = 0; k < TotalVehicles; k++)
			count += Vehicle[k].size;

		if (count == 2 * TotalRequests)
			Bool = true;

		return Bool;
	}*/

	__host__ __device__
		void RandomShuffle()
	{

		//Pick Vehicle length randomly
		int a[TotalVehicles];
		float divident = 0;
		for (int i = 0; i<TotalVehicles; i++)
		{
		//	srand(time_now++);
			a[i] = rand() % 100;
			//printf("%d\n", a[i]);
			divident += a[i];
		}

		float factor = (2 * n / divident);
		//printf("%f = %d / %f \n", factor, 2 * n, divident);

		int sum = 0;
		for (int i = 0; i<TotalVehicles - 1; i++)
		{
			a[i] *= factor;
			if ((a[i] % 2 != 0))
				a[i] += 1;

			sum += a[i];
		}

		a[TotalVehicles - 1] = 2 * n - sum;



		//////////////////

		//Initialization
		for (int k = 0; k < TotalVehicles; k++)
		{
			Vehicle[k].id = k;
			Vehicle[k].initialize(a[k]);
		}

		//Set RequestIds
		for (int i = 0; i<TotalRequests; i++)
			Request[i].id = i + 1;

		/////////////////

		//Copy
		int sz[TotalVehicles];
		for (int i = 0; i < TotalVehicles; i++)
			sz[i] = a[i];


		int isSelected[2 * n];
		for (int i = 0; i < 2 * n; i++)
			isSelected[i] = 0;

		//Initialize
		for (int i = 0; i < TotalVehicles; i++)
			for (int j = 0; j < sz[i]; j++)
			{
				Vehicle[i].path[j] = 0;
				Vehicle[i].size = sz[i];
			}

		//Filling loop
		for (int i = 0; i < TotalVehicles; i++)
		{


			//cout << "Vehicle " << i << "\n";

			int fill_count = 0;
			int j = 0;

			while (fill_count < Vehicle[i].size)
			{
				//Select Request
				int req;
				int found = 0;
				while (!found)
				{
					//srand(time_now++);
					req = RandomNumber(0, TotalRequests - 1);

					if (isSelected[req])
						continue;
					else
					{
						found = 1;
						isSelected[req] = 1;
					}

				}


				//Select Positions to insert that req.
				//srand(time_now++);=
				int ed1 = Vehicle[i].size - 2;
				int start = RandomNumber(0, ed1);
				int found_start = 0;

				int roomExist = 0;
				for (int l = start + 1; l < sz[i]; l++)
					if (Vehicle[i].path[l] == 0)
						roomExist = 1;

				while (!found_start)
				{

					if (!roomExist || Vehicle[i].path[start])
					{
						roomExist = 0;

						//srand(time_now++);
						start = RandomNumber(0, ed1);

						for (int l = start + 1; l < sz[i]; l++)
							if (Vehicle[i].path[l] == 0)
								roomExist = 1;
					}

					if (Vehicle[i].path[start] == 0 && roomExist == 1)
					{
						Vehicle[i].path[start] = req + 1;
						found_start = 1;
					}
				}

				//srand(time_now++);
				int st = start + 1;
				int ed = sz[i] - 1;
				int gap = RandomNumber(st, ed);//RandomNumber(start + 1, Vehicle[i].size - 1);
				int found_gap = 0;

				while (!found_gap)
				{
					if (Vehicle[i].path[gap])
					{
						//srand(time_now++);
						int st = start + 1;
						int ed = sz[i] - 1;
						gap = RandomNumber(st, ed);//RandomNumber(start + 1, Vehicle[i].size - 1);
					}
					else
					{
						Vehicle[i].path[gap] = req + 1 + n;
						found_gap = 1;
					}
				}

				fill_count += 2;
				j++;

			}
		}



		IdentifyCurrvehicle();

		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].id = k + 1;

		for (int i = 0; i < TotalRequests; i++)
			Request[i].id = i + 1;

	}

	__host__ __device__
		void IdentifyCurrvehicle()
	{

		//IdentifyCurrVehicle
		int count = 0;
		for (int k = 0; k < TotalVehicles; k++)
		{
			for (int j = 0; j < Vehicle[k].size; j++)
				if (Vehicle[k].path[j] <= TotalRequests)
				{
					Request[Vehicle[k].path[j] - 1].currvehicle = k; //Set currvehicle
					Request[Vehicle[k].path[j] - 1].totinsertion = 0; //Initialize totinsertion
					count++;
				}

			if (count > TotalRequests)
				break;
		}


		//ReqInsertionways
		for (int i = 0; i < TotalRequests; i++)
		{
			Request[i].totinsertion = 0;

			for (int k = 0; k < TotalVehicles; k++)
				if (k != Request[i].currvehicle)
					Request[i].totinsertion += Vehicle[k].insertways;
		}

		//Intra_ReqInsertionWays
		for (int i = 0; i < TotalRequests; i++)
		{
			int current = Request[i].currvehicle;
			Request[i].intra_totinsertion = (Vehicle[current].size - 2 + 1) * (Vehicle[current].size - 2 + 2) / 2;
		}


	}

	__host__ __device__
		void UpdateInformation(insertion *Insertion)
	{
		//Align
		int i = Insertion[0].request;
		int k = Insertion[0].vehid;

		//1) CurrVehicle
		Request[i].currvehicle = k;

		//2) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		//3) ReqInsertionways
		for (int i = 0; i < TotalRequests; i++)
		{
			Request[i].totinsertion = 0;

			for (int k = 0; k < TotalVehicles; k++)
				if (k != Request[i].currvehicle)
					Request[i].totinsertion += Vehicle[k].insertways;
		}

		//4) Intra_ReqInsertionWays
		for (int i = 0; i < TotalRequests; i++)
		{
			int k = Request[i].currvehicle;
			Request[i].intra_totinsertion = ((Vehicle[k].size - 2) + 1) * (((Vehicle[k].size - 2)) + 2) / 2;
		}
	}

	__host__ __device__
		void printinput()
	{

		printf("\nTotalRequests: %d TotalVehicles: %d\n", TotalRequests, TotalVehicles);
		for (int k = 0; k < TotalVehicles; k++)
		{
			Vehicle[k].id = k;
			printf("sz: %d Path: ", Vehicle[k].size);
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				printf("%d ", Vehicle[k].path[j]);
			}
			printf("\n");
		}
		printf("\n");
	}

	__host__ __device__
		void printinput_with_details(problem *Problem)
	{
		int Totalsize = 0;
		for (int k = 0; k < TotalVehicles; k++)
			Totalsize += Vehicle[k].size;

		printf("\nTotalRequests: %d TotalVehicles: %d\n", TotalRequests, TotalVehicles);
		printf("Totalsize: %d\n"), Totalsize;
		for (int k = 0; k < TotalVehicles; k++)
		{
			Vehicle[k].id = k;
			printf("sz: %d \n", Vehicle[k].size);
			printf("node	Serving_TW	Original_TW\n");
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				int node = Vehicle[k].path[j] - 1;
				if (node < n)
					printf("%d	[%0.0f	%0.0f]	[%0.0f	%0.0f]\n", node+1, Request[node].B, Request[node].C, Problem[0].Request[node].pickup.EarliestTime, Problem[0].Request[node].pickup.LatestTime);
				else
					printf("%d	[%0.0f	%0.0f]	[%0.0f	%0.0f]\n", node+1, Request[node].B, Request[node].C, Problem[0].Request[node-n].dropoff.EarliestTime, Problem[0].Request[node-n].dropoff.LatestTime);
			}
			printf("\n");
		}
		printf("\n");
	}

	__host__ __device__
		void printinsertdetails()
	{

		//Print
		for (int i = 0; i < TotalRequests; i++)
			printf("ReqId %d CurrVehicle %d\n", Request[i].id, Request[i].currvehicle);

		int ALLinsertion = 0;
		for (int i = 0; i < TotalRequests; i++)
		{
			printf("Insertion of %d onto other vehicles: Ways: %d\n", Request[i].id, Request[i].totinsertion);
			ALLinsertion += Request[i].totinsertion;
		}
		printf("\nALL Insertion Possibilities: %d\n", ALLinsertion);

		if (ALLinsertion / n > Expectation)
			printf("WARNING !!!!! Increase the EXPECTATION !!!!!!\n");
	}

	__host__ __device__
		int fetchALLinsertion()
	{
		int ALLinsertion = 0;

		for (int i = 0; i < TotalRequests; i++)
			ALLinsertion += Request[i].totinsertion;

		return ALLinsertion;
	}

	__host__ __device__
		void EV_based_evaluation()
	{
		for (int k = 0; k < m; k++)
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				//
			}

	}

	__host__ __device__
		void update_overall_cost()
	{
		Temp.reset();
		for (int k = 0; k < TotalVehicles; k++)
		{

			Temp.travel_cost += Vehicle[k].Cost.travel_cost;
			Temp.load += Vehicle[k].Cost.load;
			Temp.time_window += Vehicle[k].Cost.time_window;
			Temp.ride_time += Vehicle[k].Cost.ride_time;
			Temp.duration += Vehicle[k].Cost.duration;
			Temp.excessrideTime += Vehicle[k].Cost.excessrideTime;

			Temp.batteryViolation += Vehicle[k].Cost.batteryViolation;

			Cost.Copy(Temp);
			Cost.getFeasibility();
		}
	}

	__host__ __device__
		void Offline_eight_step_evaluation()
	{
		Feasible_Req_Count = 0;
		total_req_served = 0;

		Temp.reset();
		for (int k = 0; k < TotalVehicles; k++)
		{
			if (Vehicle[k].size == 0)
			{
				Vehicle[k].Cost.reset();
				Vehicle[k].Cost.isFeasible = true;
			}
			else
			{
				//OfflineEvaluateByRoute(k + 1);
				Offline_route_evaluation(k + 1);
				//Offline_new_route_evaluation(k);
				//Offline_new_route_evaluation_for_DARP(k);
				Vehicle[k].Cost.getFeasibility();
			}

			Temp.on_board_user_wait_time += Vehicle[k].Cost.on_board_user_wait_time;
			Temp.route_duration += Vehicle[k].Cost.route_duration;
			Temp.early_arrival_time += Vehicle[k].Cost.early_arrival_time;

			Temp.travel_cost += Vehicle[k].Cost.travel_cost;
			Temp.load += Vehicle[k].Cost.load;
			Temp.time_window += Vehicle[k].Cost.time_window;
			Temp.ride_time += Vehicle[k].Cost.ride_time;
			Temp.duration += Vehicle[k].Cost.duration;
			Temp.excessrideTime += Vehicle[k].Cost.excessrideTime;

			Temp.batteryViolation += Vehicle[k].Cost.batteryViolation;
			Temp.max_battery_Violation = max(Temp.max_battery_Violation, Vehicle[k].Cost.batteryViolation);

			Temp.overall_ideal_cost += Vehicle[k].Cost.overall_ideal_cost;
			total_req_served += (Vehicle[k].size / 2);
		}

		Cost = Temp;
		Cost.getFeasibility();

		if (Cost.batteryViolation < 0.001)
			Cost.isOverall_batteryfeasible = true;
		else
			Cost.isOverall_batteryfeasible = false;

		// update max_battery_violation
		/*Cost.max_battery_Violation = Vehicle[0].Cost.max_battery_Violation;
		for (int k = 1; k < TotalVehicles; k++)
		if (Vehicle[k].Cost.max_battery_Violation > Cost.max_battery_Violation)
		Cost.max_battery_Violation = Vehicle[k].Cost.max_battery_Violation;*/


		for (int k = 0; k < TotalVehicles; k++)
		{
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				// Here we are not checking load violation, considering no violation due to previously feasible route provided by the solver before . THat's makes sense !!!!

				int node = Vehicle[k].path[j] - 1;
				if (node < n)
				{
					if (Request[node].time_window_violation < 0.01)
						if (Request[node].ride_time_violation  < 0.01)
							if (Request[node + n].time_window_violation < 0.01)
								if (Request[node + n].ride_time_violation < 0.01)
									Feasible_Req_Count++;
				}
			}
		}
	}

	__host__ __device__
		void Offline_as_early_as_possible_evaluation()
	{
		Feasible_Req_Count = 0;
		total_req_served = 0;

		Temp.reset();
		for (int k = 0; k < TotalVehicles; k++)
		{

			if (Vehicle[k].size == 0)
			{
				Vehicle[k].Cost.reset();
				Vehicle[k].Cost.isFeasible = true;
			}
			else
			{
				Offline_route_as_early_as_possible_evaluation(k);
				Vehicle[k].Cost.getFeasibility();
			}


			Temp.travel_cost += Vehicle[k].Cost.travel_cost;
			Temp.load += Vehicle[k].Cost.load;
			Temp.time_window += Vehicle[k].Cost.time_window;
			Temp.ride_time += Vehicle[k].Cost.ride_time;
			Temp.duration += Vehicle[k].Cost.duration;
			Temp.excessrideTime += Vehicle[k].Cost.excessrideTime;

			Temp.batteryViolation += Vehicle[k].Cost.batteryViolation;
			Temp.max_battery_Violation = max(Temp.max_battery_Violation, Vehicle[k].Cost.max_battery_Violation);

			total_req_served += (Vehicle[k].size / 2);
		}

		Cost = Temp;
		Cost.getFeasibility();

		if (Cost.batteryViolation < 0.001)
			Cost.isOverall_batteryfeasible = true;
		else
			Cost.isOverall_batteryfeasible = false;

		// update max_battery_violation
		/*Cost.max_battery_Violation = Vehicle[0].Cost.max_battery_Violation;
		for (int k = 1; k < TotalVehicles; k++)
		if (Vehicle[k].Cost.max_battery_Violation > Cost.max_battery_Violation)
		Cost.max_battery_Violation = Vehicle[k].Cost.max_battery_Violation;*/


		for (int k = 0; k < TotalVehicles; k++)
		{
			for (int j = 0; j < Vehicle[k].size; j++)
			{

				int node = Vehicle[k].path[j] - 1;
				if (node < n)
				{
					if (Request[node].time_window_violation < 0.01)
						if (Request[node].ride_time_violation  < 0.01)
							if (Request[node + n].time_window_violation < 0.01)
								if (Request[node + n].ride_time_violation < 0.01)
									Feasible_Req_Count++;
				}
			}
		}
	}

	__host__ __device__
		cost Offline_fetch_route_evaluation(int k)
	{
		//here k starts from zero.

		if (Vehicle[k].size == 0)
		{
			Vehicle[k].Cost.reset();
			Vehicle[k].Cost.isFeasible = true;
		}
		else
		{
			Offline_new_route_evaluation(k + 1);
			Vehicle[k].Cost.getFeasibility();
		}
		return Vehicle[k].Cost;
	}

	__host__ __device__
		void update_A_B_C_W(int k, int current_index, float depart_time, int depart_location, int arrival_location)
	{
		// passed k value starts from zero
		// update from this current_index (this index excluded, dont update it.)

		vertex VertexNow;

		float current_time = depart_time;
		int current_location = depart_location;

		for (int i = current_index + 1; i < Vehicle[k].size; i++)
		{
			int j = Vehicle[k].path[i] - 1;
			VertexNow = find_vertex(Problem, j);

			Request[j].travel_cost = Problem.TravelTime[current_location][j];
			Request[j].A = current_time + Problem.TravelTime[current_location][j];
			Request[j].B = max(Request[j].A, VertexNow.EarliestTime);
			Request[j].C = Request[j].B + VertexNow.ServiceTime;
			Request[j].W = Request[j].B - Request[j].A;

			current_time = Request[j].C;
			current_location = j;
		}

		Vehicle[k].end_time = current_time + Problem.TravelTime[current_location][arrival_location];
	}

	__host__ __device__
		void update_load(int k)
	{
		int current_load = 0;
		for (int i = 0; i < Vehicle[k].size; i++)
		{
			int j = Vehicle[k].path[i] - 1;

			if (j < n)
				current_load += Problem.Request[j].load;
			else
				current_load -= Problem.Request[j - n].load;

			Request[j].p = current_load;
		}
	}

	__host__ __device__
		void update_L_p(int k)
	{
		for (int i = 0; i < Vehicle[k].size; i++)
		{
			int j = Vehicle[k].path[i] - 1;

			if (j < n)
				Request[j].L = 0;
			else
				Request[j].L = Request[j].B - Request[j - n].C;
		}
	}

	__host__ __device__
		void compute_slack_for_vertex(int k, int start_vertex_index, float &min_slack, float &cumulative_waiting_time)
	{
		float *factor_1 = new float[ExpectedPath]; // summation Wp
		float *factor_2 = new float[ExpectedPath]; // l{j} - B{j}
		float *factor_3 = new float[ExpectedPath]; // RTC - P{j}; if j=>dropoff, else 0.

		float *slack = new float[ExpectedPath];

		// 1) computing factor_1 i.e., running Wp

		cumulative_waiting_time = 0;
		for (int i = start_vertex_index; i < Vehicle[k].size; i++)
		{
			int node = Vehicle[k].path[i] - 1;
			if (i == 0)
				cumulative_waiting_time += Request[node].W;
			else if (i == start_vertex_index)
				cumulative_waiting_time = 0;
			else
				cumulative_waiting_time += Request[node].W;



			//printf("cumulative_waiting_time	=	%f,	Request[%d].W	=	%f\n", cumulative_waiting_time, node, Request[node].W);

			factor_1[i] = cumulative_waiting_time;
		}

		// 2) compute factor_2 and factor_3

		for (int i = start_vertex_index; i < Vehicle[k].size; i++)
		{
			int node = Vehicle[k].path[i] - 1;
			if (node < n) //pickup
			{
				factor_2[i] = max(0, Problem.Request[node].pickup.LatestTime - Request[node].B);
				factor_3[i] = Problem.Request[node].RideTimeConstraint;
			}
			else //dropoff
			{
				factor_2[i] = max(0, Problem.Request[node - n].dropoff.LatestTime - Request[node].B);
				factor_3[i] = max(0, Problem.Request[node - n].RideTimeConstraint - Request[node].L);
			}
		}

		// 3) compute minimum_slack_time

		for (int i = start_vertex_index; i < Vehicle[k].size; i++)
			slack[i] = factor_1[i] + min(factor_2[i], factor_3[i]);

		int min_i = start_vertex_index;
		float minimum_slack = slack[start_vertex_index];
		for (int i = start_vertex_index + 1; i < Vehicle[k].size; i++)
		{
			if (slack[i] < minimum_slack)
			{
				// if (Vehicle[k].path[i] - 1 < n) // only pickups are considered!!!
				{
					minimum_slack = slack[i];
					min_i = i;
				}
			}
		}

		min_slack = minimum_slack;


		/*printf("For vehicle k: %d\n", k);
		printf("i	node	factor_1[i]	factor_2[i]	factor_3[i]	slack[i]\n");
		for (int i = start_vertex_index; i < Vehicle[k].size; i++)
		{
		if (Vehicle[k].path[i] - 1 < n)
		printf("%d	%d	%f	%f	%f	%f	pickup\n", i, Vehicle[k].path[i], factor_1[i], factor_2[i], factor_3[i], slack[i]);
		else
		printf("%d	%d	%f	%f	%f	%f\n", i, Vehicle[k].path[i], factor_1[i], factor_2[i], factor_3[i], slack[i]);
		}
		printf("\nmin_slack = %f, at index %d\n", min_slack, min_i);*/

		delete[] factor_1;
		delete[] factor_2;
		delete[] factor_3;
		delete[] slack;
	}

	__host__ __device__
		void compute_violations(int k, int arrival_location)
	{
		Vehicle[k].Cost.reset();
		for (int i = 0; i < Vehicle[k].size; i++)
		{
			int j = Vehicle[k].path[i] - 1;

			if (j < n) //pickup
			{
				Request[j].ride_time_violation = 0;
				Request[j].time_window_violation = max((Request[j].B - Problem.Request[j].pickup.LatestTime), 0.0);
			}
			else //dropoff
			{
				Request[j].ride_time_violation = max((Request[j].L - Problem.Request[j - n].RideTimeConstraint), 0.0);
				Request[j].time_window_violation = max((Request[j].B - Problem.Request[j - n].dropoff.LatestTime), 0.0);
			}

			if (Request[j].p > Problem.Vehicle[k].capacity)
				Vehicle[k].Cost.load = max(Vehicle[k].Cost.load, Request[j].p - Problem.Vehicle[k].capacity);

			Vehicle[k].Cost.ride_time += Request[j].ride_time_violation;
			Vehicle[k].Cost.time_window += Request[j].time_window_violation;
			Vehicle[k].Cost.travel_cost += Request[j].travel_cost;
		}

		int last_visited_node = Vehicle[k].path[Vehicle[k].size - 1] - 1;
		Vehicle[k].Cost.travel_cost += Problem.TravelTime[last_visited_node][arrival_location];

		Vehicle[k].Cost.duration = max((Vehicle[k].end_time - Vehicle[k].start_time) - Problem.Vehicle[k].duration_constraint, 0.0);
		Vehicle[k].Cost.duration += max(Vehicle[k].end_time - Problem.Vehicle[k].end_time, 0.0);


	}

	__host__ __device__
		void Offline_route_evaluation(int k)
	{
		//here, passed k value starts from 1.

		int vehid = k - 1;
		k--;

		if (Vehicle[k].size < 1)
		{
			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = Problem.Vehicle[vehid].start_time;
			float depart_time = Vehicle[k].start_time;
			int depart_location = 2 * n + vehid;
			Vehicle[k].end_time = Vehicle[k].start_time + Problem.TravelTime[depart_location][2 * n + m + vehid];
			Vehicle[k].Cost.duration = max(Vehicle[k].end_time - Vehicle[k].start_time - Problem.Vehicle[vehid].duration_constraint, 0.0);
			Vehicle[k].Cost.duration += max(Vehicle[k].end_time - Problem.Vehicle[vehid].end_time, 0.0);
		2	Vehicle[k].Cost.travel_cost += Problem.TravelTime[depart_location][2 * n + m + vehid];
		}
		else
		{

			//  1. Set D0 : ¼? e0.

			int depart_location = 2 * n + k;
			int arrival_location = 2 * n + m + k;

			float depart_time = Problem.Vehicle[vehid].start_time;

			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = depart_time;

			//  2. Compute Ai, Wi, Bi and Di and for each vertex vi in the route.

			for (int i = 0; i < Vehicle[k].size; i++)
			{
				int node = Vehicle[k].path[i] - 1;

				//printf("node: %d\n", node);

				Request[node].A = 0;
				Request[node].B = 0;
				Request[node].C = 0;
				Request[node].L = 0;
				Request[node].W = 0;
			}

			//print_with_battery_level();
			//printf("ABOVE IS BEFORE FIRST UPDATE for k = %d!!!!!!!!\n", k);
			// system("PAUSE");

			update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);

			update_load(k);
			update_L_p(k);

			//print_with_battery_level();
			//printf("ABOVE IS AFTER FIRST UPDATE for k = %d!!!!!!!!\n", k);
			// system("PAUSE");

			//	3. Compute F0.

			float cumulative_waiting_time;
			float F_slack[ExpectedPath];
			compute_slack_for_vertex(k, 0, F_slack[0], cumulative_waiting_time);

			//printf("k = %d, F_slack[0] = %f\n cumulative_waiting_time = %f\n", k, F_slack[0], cumulative_waiting_time);

			//	4. Set D0 : ¼? e0 þ? minfF0;0<p<q Wpg.

			depart_time = Problem.Vehicle[vehid].start_time + min(F_slack[0], cumulative_waiting_time);
			Vehicle[k].start_time = depart_time;

			//	5. Update Ai, Wi, Bi and Di for each vertex vi in the route.

			//print_with_battery_level();
			//printf("ABOVE IS BEFORE SLACK AT DEPOT - UPDATE for k = %d!!!!!!!!\n", k);
			//// system("PAUSE");

			update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);

			//	6. Compute Li for each request assigned to the route.

			update_load(k);
			update_L_p(k);

			//print_with_battery_level();
			//printf("ABOVE IS AFTER SLACK AT DEPOT - UPDATE for k = %d!!!!!!!!\n", k);
			//// system("PAUSE");

			compute_violations(k, arrival_location);
			//Vehicle[k].Cost.print();
			//printf("ABOVE IS FIRST VIOLATION COMPUTATION bef eight-step for k = %d\n", k);
			// system("PAUSE");

			//	7. For everyvertex vj that corresponds to the origin of a request j

			for (int index = 0; index < Vehicle[k].size; index++)
			{
				vertex Curr_Vertex;

				int j = Vehicle[k].path[index] - 1;

				if (j < n) //pickup
				{
					if (0)// (Request[j].p <= 1) && (Request[j + n].ride_time_violation == 0) && (j + n == Vehicle[k].path[index + 1] - 1) )
					{
						printf("Before slacking: REQ %d has RTD violation of %f\n", j, Request[j].ride_time_violation);
						//// system("PAUSE");

						Curr_Vertex = find_vertex(Problem, j);

						//	(a) Compute Fj.

						compute_slack_for_vertex(k, index, F_slack[index], cumulative_waiting_time);

						//	(b) Set Bj : ¼? Bj þ? minfFj; j<p<q Wpg; Dj:¼? Bj þ? dj.

						Request[j].B = Request[j].B + min(F_slack[index], cumulative_waiting_time);
						Request[j].C = Request[j].B + Curr_Vertex.ServiceTime;
						Request[j].W = Request[j].B - Request[j].A;

						printf("UPDATED VALUES:\n");
						printf("ETW	LTW	Request[j].B	Request[j].C	Request[j].W\n");
						printf("%0.02f	%0.02f	%0.02f	%0.02f	%0.02f\n", Curr_Vertex.EarliestTime, Curr_Vertex.LatestTime, Request[j].B, Request[j].C, Request[j].W);

						//	(c) Update Ai, Wi, Bi and Di, for each vertex vi that comes after vj in the route.

						update_A_B_C_W(k, index, Request[j].C, j, arrival_location);

						//	(d) Update the ride time Li for each request i whose destination vertex is after vertex vj.

						update_L_p(k);


						print_with_battery_level();
						printf("ABOVE IS AFTER SLACKING OFF req = %d at vehicle %d!!!!!!!!\n", j, k);
						compute_violations(k, arrival_location);
						Vehicle[k].Cost.print();
						printf("After slacking: REQ %d has RTD violation of %f\n", j, Request[j + n].ride_time_violation);
						printf("SPECIAL!!!!!!!!!!!!!!!!!!!\n");
						// system("PAUSE");
					}
					else if (0)//Request[j + n].ride_time_violation == 0)
					{
						printf("PROCEDURES IS SKIPPED FOR REQ %d since zero RTD violation of %f\n", j, Request[j].ride_time_violation);
						// system("PAUSE");
					}
					/*else if (Request[j].p > 1)
					{
					printf("PROCEDURES IS SKIPPED FOR REQ %d since A CUSTOMER IS INSIDE ALREADY (p = %d), btw RTD of REQ is %f\n", j, Request[j].p, Request[j].ride_time_violation);
					// system("PAUSE");
					}*/
					else
					{
						//printf("Before slacking: REQ %d has RTD violation of %f\n", j, Request[j].ride_time_violation);
						//// system("PAUSE");

						Curr_Vertex = find_vertex(Problem, j);

						//	(a) Compute Fj.

						compute_slack_for_vertex(k, index, F_slack[index], cumulative_waiting_time);

						//	(b) Set Bj : ¼? Bj þ? minfFj; j<p<q Wpg; Dj:¼? Bj þ? dj.

						Request[j].B = Request[j].B + min(F_slack[index], cumulative_waiting_time);
						Request[j].C = Request[j].B + Curr_Vertex.ServiceTime;
						Request[j].W = Request[j].B - Request[j].A;

						//printf("UPDATED VALUES:\n");
						//printf("ETW	LTW	Request[j].B	Request[j].C	Request[j].W\n");
						//printf("%0.02f	%0.02f	%0.02f	%0.02f	%0.02f\n", Curr_Vertex.EarliestTime, Curr_Vertex.LatestTime, Request[j].B, Request[j].C, Request[j].W);

						//	(c) Update Ai, Wi, Bi and Di, for each vertex vi that comes after vj in the route.

						update_A_B_C_W(k, index, Request[j].C, j, arrival_location);

						//	(d) Update the ride time Li for each request i whose destination vertex is after vertex vj.

						update_L_p(k);


						//print_with_battery_level();
						//printf("ABOVE IS AFTER SLACKING OFF req = %d at vehicle %d!!!!!!!!\n", j, k);
						compute_violations(k, arrival_location);
						//Vehicle[k].Cost.print();
						//printf("After slacking: REQ %d has RTD violation of %f\n", j, Request[j + n].ride_time_violation);
						// system("PAUSE");
					}
				}


			}

			//	8. Compute changes in violations of vehicle load, route duration, time window and ride time constraints.

			compute_violations(k, arrival_location);
			//Vehicle[k].Cost.print();

			// calculate excessrideTime
			Vehicle[k].Cost.excessrideTime = 0;
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				int node = Vehicle[k].path[j] - 1;
				if (node >= n)
					Vehicle[k].Cost.excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);
			}

			//printf("End of evaluation for vehicle k = %d, with excessrideTime = %f:!!\n", k, Vehicle[k].Cost.excessrideTime);
			// system("PAUSE");

			// Modified objective function
			Vehicle[k].Cost.on_board_user_wait_time = 0;
			Vehicle[k].Cost.early_arrival_time = 0;
			Vehicle[k].Cost.excessrideTime = 0;
			Vehicle[k].Cost.route_duration = 0;
			for (int index = 0; index < Vehicle[k].size; index++)
			{
				int node = Vehicle[k].path[index] - 1;

				Vehicle[k].Cost.on_board_user_wait_time += (Request[node].W * (Request[node].p-1));

				if (node < n)
				{
					if (Problem.Request[node].pickup.EarliestTime > Request[node].A)
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node].pickup.EarliestTime - Request[node].A);
				}
				else
				{
					Vehicle[k].Cost.excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);

					if (Problem.Request[node - n].dropoff.EarliestTime > Request[node].A)
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node - n].dropoff.EarliestTime - Request[node].A);
				}

			}
			Vehicle[k].Cost.route_duration = Vehicle[k].end_time - Vehicle[k].start_time;

		}

		Vehicle[k].Cost.overall_ideal_cost = (w1 * Vehicle[k].Cost.travel_cost + w2 * Vehicle[k].Cost.excessrideTime)
							+ 1* Vehicle[k].Cost.load
							+ 1* Vehicle[k].Cost.time_window
							+ 1* Vehicle[k].Cost.ride_time
							+ 1* Vehicle[k].Cost.duration;

	}


	__host__ __device__
		void Offline_mod_route_evaluation(int k)
	{
		//here, passed k value starts from 1.

		int vehid = k - 1;
		k--;

		if (Vehicle[k].size < 1)
		{
			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = Problem.Vehicle[vehid].start_time;
			float depart_time = Vehicle[k].start_time;
			int depart_location = 2 * n + vehid;
			Vehicle[k].end_time = Vehicle[k].start_time + Problem.TravelTime[depart_location][2 * n + m + vehid];
			Vehicle[k].Cost.duration = max(Vehicle[k].end_time - Vehicle[k].start_time - Problem.Vehicle[vehid].duration_constraint, 0.0);
			Vehicle[k].Cost.duration += max(Vehicle[k].end_time - Problem.Vehicle[vehid].end_time, 0.0);
			Vehicle[k].Cost.travel_cost += Problem.TravelTime[depart_location][2 * n + m + vehid];
		}
		else
		{

			//  1. Set D0 : ¼? e0.
			int depart_location = 2 * n + k;
			int arrival_location = 2 * n + m + k;
			float depart_time = Problem.Vehicle[vehid].start_time;

			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = depart_time;

			//  2. Compute Ai, Wi, Bi and Di and for each vertex vi in the route.
			for (int i = 0; i < Vehicle[k].size; i++)
			{
				int node = Vehicle[k].path[i] - 1;
				Request[node].A = 0;
				Request[node].B = 0;
				Request[node].C = 0;
				Request[node].L = 0;
				Request[node].W = 0;
			}

			update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);
			update_load(k);
			update_L_p(k);

			//	3. Compute F0.

			float cumulative_waiting_time;
			float F_slack[ExpectedPath];
			compute_slack_for_vertex(k, 0, F_slack[0], cumulative_waiting_time);

			//	4. Set D0 : ¼? e0 þ? minfF0;0<p<q Wpg.

			depart_time = Problem.Vehicle[vehid].start_time + min(F_slack[0], cumulative_waiting_time);
			Vehicle[k].start_time = depart_time;

			update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);

			//	6. Compute Li for each request assigned to the route.
			update_load(k);
			update_L_p(k);
			compute_violations(k, arrival_location);

			//	7. For everyvertex vj that corresponds to the origin of a request j

			for (int index = 0; index < Vehicle[k].size; index++)
			{
				vertex Curr_Vertex;

				int j = Vehicle[k].path[index] - 1;

				if (j < n) //pickup
				{

					Curr_Vertex = find_vertex(Problem, j);

					//	(a) Compute Fj.
					compute_slack_for_vertex(k, index, F_slack[index], cumulative_waiting_time);

					//	(b) Set Bj : ¼? Bj þ? minfFj; j<p<q Wpg; Dj:¼? Bj þ? dj.
					Request[j].B = Request[j].B + min(F_slack[index], cumulative_waiting_time);
					Request[j].C = Request[j].B + Curr_Vertex.ServiceTime;
					Request[j].W = Request[j].B - Request[j].A;

					//	(c) Update Ai, Wi, Bi and Di, for each vertex vi that comes after vj in the route.
					update_A_B_C_W(k, index, Request[j].C, j, arrival_location);
					update_L_p(k);
					compute_violations(k, arrival_location);
				}
			}

			//	8. Compute changes in violations of vehicle load, route duration, time window and ride time constraints.

			// delay departure at the depot

			int first_node = Vehicle[k].path[0] - 1;

			Vehicle[k].start_time += Request[first_node].W;
			Request[first_node].A = Request[first_node].B;
			Request[first_node].W = Request[first_node].B - Request[first_node].A;

			compute_violations(k, arrival_location);

			// Modified objective function
			Vehicle[k].Cost.on_board_user_wait_time = 0;
			Vehicle[k].Cost.early_arrival_time = 0;
			Vehicle[k].Cost.excessrideTime = 0;
			Vehicle[k].Cost.route_duration = 0;
			for (int index = 0; index < Vehicle[k].size; index++)
			{
				int node = Vehicle[k].path[index] - 1;

				Vehicle[k].Cost.on_board_user_wait_time += (Request[node].W * (Request[node].p-1));

				if (node < n)
				{
					if (Problem.Request[node].pickup.EarliestTime > Request[node].A)
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node].pickup.EarliestTime - Request[node].A);
				}
				else
				{
					Vehicle[k].Cost.excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);

					if (Problem.Request[node - n].dropoff.EarliestTime > Request[node].A)
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node - n].dropoff.EarliestTime - Request[node].A);
				}

			}
			Vehicle[k].Cost.route_duration = Vehicle[k].end_time - Vehicle[k].start_time;

		}

		Vehicle[k].Cost.overall_ideal_cost = (w1 * Vehicle[k].Cost.travel_cost + w2 * Vehicle[k].Cost.excessrideTime)
							+ 1* Vehicle[k].Cost.load
							+ 1* Vehicle[k].Cost.time_window
							+ 1* Vehicle[k].Cost.ride_time
							+ 1* Vehicle[k].Cost.duration;

	}



	__host__ __device__
		void Offline_modified_route_evaluation(int k)
	{
		k--;

		int vehid = k;

		if (Vehicle[k].size < 1)
		{
			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = Problem.Vehicle[vehid].start_time;
			float depart_time = Vehicle[k].start_time;
			int depart_location = 2 * n + vehid;
			Vehicle[k].end_time = Vehicle[k].start_time + Problem.TravelTime[depart_location][2 * n + m + vehid];
			Vehicle[k].Cost.duration = max(Vehicle[k].end_time - Vehicle[k].start_time - Problem.Vehicle[vehid].duration_constraint, 0.0);
			Vehicle[k].Cost.duration += max(Vehicle[k].end_time - Problem.Vehicle[vehid].end_time, 0.0);
			Vehicle[k].Cost.travel_cost += Problem.TravelTime[depart_location][2 * n + m + vehid];
		}
		else
		{

			//  1. Set D0 : ¼ e0.

			int depart_location = 2 * n + vehid;;
			int arrival_location = 2 * n + m + vehid;

			float depart_time = Problem.Vehicle[vehid].start_time;

			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = depart_time;

			//  2. Compute Ai, Wi, Bi and Di and for each vertex vi in the route.

			float A[2 * n],
				B[2 * n],
				C[2 * n],
				W[2 * n],
				L[2 * n],
				travel_cost[2 * n];
			int	p[2 * n];

			for (int i = 0; i < Vehicle[k].size; i++)
			{
				int j = Vehicle[k].path[i] - 1;

				A[j] = 0;
				B[j] = 0;
				C[j] = 0;
				L[j] = 0;
				W[j] = 0;
				p[j] = 0;
			}

			update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);
			update_L_p(k);

			//	3. Compute F0.

			float cumulative_waiting_time;
			float F_slack;
			compute_slack_for_vertex(k, 0, F_slack, cumulative_waiting_time);

			//	4. Set D0 : ¼ e0 þ minfF0;0<p<q Wpg.

			depart_time = Problem.Vehicle[vehid].start_time + min(F_slack, cumulative_waiting_time);
			Vehicle[k].start_time = depart_time;

			//	5, 6. Update Ai, Wi, Bi, Di, Li, p for each vertex vi in the route.

			update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);
			update_L_p(k);

			compute_violations(k, arrival_location);

			//	7. For everyvertex vj that corresponds to the origin of a request j

			for (int index = 0; index < Vehicle[k].size; index++)
			{

				int j = Vehicle[k].path[index] - 1;

				if (j < n) //pickup
				{
					float Serv_time = Problem.Request[j].pickup.ServiceTime;

					//	(a) Compute Fj.

					compute_slack_for_vertex(k, index, F_slack, cumulative_waiting_time);

					//	(b) Set Bj : ¼ Bj þ minfFj; j<p<q Wpg; Dj:¼ Bj þ dj.

					B[j] = B[j] + min(F_slack, cumulative_waiting_time);
					C[j] = B[j] + Serv_time;
					W[j] = B[j] - A[j];

					//	(c) Update Ai, Wi, Bi and Di, for each vertex vi that comes after vj in the route.

					update_A_B_C_W(k, index, C[j], j, arrival_location);

					//	(d) Update the ride time Li for each request i whose destination vertex is after vertex vj.

					update_L_p(k);

					compute_violations(k, arrival_location);
				}


			}


			// delay departure at the depot

			int first_node = Vehicle[k].path[0] - 1;

			Vehicle[k].start_time += W[first_node];
			A[first_node] = B[first_node];
			W[first_node] = B[first_node] - A[first_node];


			//	8. Compute changes in violations of Vehicle[k] load, route duration, time window and ride time constraints.

			compute_violations(k, arrival_location);


			// Modified objective function
			Vehicle[k].Cost.on_board_user_wait_time = 0;
			Vehicle[k].Cost.early_arrival_time = 0;
			Vehicle[k].Cost.excessrideTime = 0;
			Vehicle[k].Cost.route_duration = 0;
			for (int index = 0; index < Vehicle[k].size; index++)
			{
				int node = Vehicle[k].path[index] - 1;

				Vehicle[k].Cost.on_board_user_wait_time += (W[node] * (p[node]-1));

				if (node < n)
				{
					if (Problem.Request[node].pickup.EarliestTime > A[node])
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node].pickup.EarliestTime - A[node]);
				}
				else
				{
					Vehicle[k].Cost.excessrideTime += (L[node] - Problem.TravelTime[node - n][node]);

					if (Problem.Request[node - n].dropoff.EarliestTime > A[node])
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node - n].dropoff.EarliestTime - A[node]);
				}

			}
			Vehicle[k].Cost.route_duration = Vehicle[k].end_time - Vehicle[k].start_time;

			Vehicle[k].Cost.overall_ideal_cost = (w1 * Vehicle[k].Cost.travel_cost + w2 * Vehicle[k].Cost.excessrideTime)
								+ 1* Vehicle[k].Cost.load
								+ 1* Vehicle[k].Cost.time_window
								+ 1* Vehicle[k].Cost.ride_time
								+ 1* Vehicle[k].Cost.duration;

			// copy back stuffs

			for (int i = 0; i < Vehicle[k].size; i++)
			{
				int j = Vehicle[k].path[i] - 1;

				Request[j].A = A[j];
				Request[j].B = B[j];
				Request[j].C = C[j];
				Request[j].L = L[j];
				Request[j].W = W[j];
				Request[j].p = p[j];
			}

		}

	}
	__host__ __device__
		void sequential_slacking(int first_pickup_index, int zero_split_index, int k)
	{
		// j => current node :: (i < j <= zero_split_index)

		// Go through every node j from i to zero_split_index
		// if request[j].W, then shift as much waiting time as possible to the node i (first_pickup_node)

		//printf("Enter sequential slack\n");

		// feasibly shift all waiting times of the nodes between (i+1)th node and zero_split_node into the ith node.

		int i = first_pickup_index;

		for (int j = i + 1; j <= zero_split_index; j++)
		{
			int node = Vehicle[k].path[j] - 1;
			if (Request[node].W > 0)
			{
				//printf("GOING TO SLACK FOR: Request[%d].W ==> %f\n", node, Request[node].W);
				//	// system("PAUSE");

				perform_slack(i, j, k, zero_split_index);
			}
			else
			{
				//printf("No slack performed: Request[%d].W ==> %f\n", node, Request[node].W);
				//	// system("PAUSE");
			}
		}
	}

	__host__ __device__
		void perform_slack(int i, int j, int k, int zero_split_index)
	{
		// i => firstpickup node, j => current_node

		//printf("BELOW: Slacking for i = %d, j = %d, k = %d\n", i, j, k);
		//// system("PAUSE");

		// calculate slack difference: LTW - B from node i to node j-1

		float slack[ExpectedPath];

		vertex Curr_Vertex;
		for (int l = i; l <= j - 1; l++)
		{
			int node = Vehicle[k].path[l] - 1;

			Curr_Vertex = find_vertex(Problem, node);
			slack[l] = max(0, Curr_Vertex.LatestTime - Request[node].B);
		}

		/*printf("index	LT	B	slack\n");
		for (int l = i; l <= j - 1; l++)
		{
		int node = Vehicle[k].path[l] - 1;

		Curr_Vertex = find_vertex(Problem, node);
		printf("%d	%f	%f	%f\n", l, Curr_Vertex.LatestTime, Request[node].B, slack[l]);
		}*/

		// find min slack difference

		float min_slack = slack[i];
		for (int l = i + 1; l <= j - 1; l++)
			if (slack[l] < min_slack)
				min_slack = slack[l];

		//printf("min_slack = %f\n", min_slack);

		// find min slack difference = min(min_slack, request[j].W)

		bool is_j_rep_ending_depot = false;
		if (j - Vehicle[k].size == 0)
			is_j_rep_ending_depot = true;

		float variable_2;
		if (is_j_rep_ending_depot) // if (j ==> index implicating the arrival depot i.e. the node after the final customer dropoff)
		{
			variable_2 = max(0, Problem.Vehicle[k].duration_constraint - Vehicle[k].end_time);

			min_slack = min(min_slack, variable_2);
		}
		else
		{
			int node = Vehicle[k].path[j] - 1;
			variable_2 = Request[node].W;

			min_slack = min(min_slack, variable_2);
		}


		//printf("min_slack = %f ==> min(min_slack, %f)\n", min_slack, variable_2);
		//// system("PAUSE");

		// adjust the times at node i (first_pickup_index)

		int node = Vehicle[k].path[i] - 1;

		Request[node].B += min_slack;
		Request[node].C += min_slack;
		Request[node].W += min_slack;

		// adjust the times from node i+1 to j-1 (in-between-nodes)

		for (int l = i + 1; l <= j - 1; l++)
		{
			int node = Vehicle[k].path[l] - 1;

			Request[node].A += min_slack;
			Request[node].B += min_slack;
			Request[node].C += min_slack;
		}

		// adjust the times at node j (zero_split_index)

		if (is_j_rep_ending_depot) // if (j ==> index implicating the arrival depot i.e. the node after the final cutomer dropoff)
		{
			Vehicle[k].end_time += min_slack;
		}
		else
		{
			int node = Vehicle[k].path[j] - 1;

			Request[node].A += min_slack;
			Request[node].W -= min_slack;
		}


		update_L_p(k);

		//print_with_battery_level();
		/*if (!is_j_rep_ending_depot)
		printf("ABOVE is after internal-slack at currentnode %d, between firstpickupnode = %d, and zeroSplitnode = %d at k = %d\n", Vehicle[k].path[j], Vehicle[k].path[i], Vehicle[k].path[zero_split_index], k);
		else
		printf("ABOVE is after internal-slack at currentnode %d, between firstpickupnode = %d, and Endingdepot at k = %d\n", Vehicle[k].path[j], Vehicle[k].path[i], k);*/

		//	// system("PAUSE");

	}

	__host__ __device__
		void Offline_new_route_evaluation(int vehicle_id)
	{
		//here, passed k value starts from 1.

		int vehid = vehicle_id;// -1;
		int k = vehicle_id;// -1;

		//int vehid = k - 1;
		//k--;
		//printf("k=%d\n", k);


		Vehicle[k].StartingPoint = 2 * n + k;
		Vehicle[k].EndingPoint = 2 * n + m + k;


		if (Vehicle[k].size < 1)
		{
			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = Problem.Vehicle[k].start_time;
			float depart_time = Vehicle[k].start_time;
			int depart_location = Vehicle[k].StartingPoint;// 2 * n + k;//
			Vehicle[k].end_time = Vehicle[k].start_time + Problem.TravelTime[depart_location][Vehicle[k].EndingPoint];
			Vehicle[k].Cost.duration = max(Vehicle[k].end_time - Vehicle[k].start_time - Problem.Vehicle[k].duration_constraint, 0.0);
			Vehicle[k].Cost.duration += max(Vehicle[k].end_time - Problem.Vehicle[k].end_time, 0.0);
			Vehicle[k].Cost.travel_cost += Problem.TravelTime[depart_location][Vehicle[k].EndingPoint];
		}
		else
		{

			//  1. Set D0 : ¼? e0.

			int depart_location = Vehicle[k].StartingPoint;
			int arrival_location = Vehicle[k].EndingPoint;


			float depart_time = Problem.Vehicle[vehid].start_time;

			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = depart_time;


			//  2. Compute Ai, Wi, Bi and Di and for each vertex vi in the route.

			for (int i = 0; i < Vehicle[k].size; i++)
			{
				int node = Vehicle[k].path[i] - 1;

				//printf("node: %d\n", node);

				Request[node].A = 0;
				Request[node].B = 0;
				Request[node].C = 0;
				Request[node].L = 0;
				Request[node].p = 0;
				Request[node].W = 0;
			}

			//print_with_battery_level();
			//printf("ABOVE IS BEFORE FIRST UPDATE for k = %d!!!!!!!!\n", k);
			//// system("PAUSE");

			// scheduling - start service at every node ASAP

			update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);

			update_load(k);
			update_L_p(k);

			//print_with_battery_level();
			//printf("ABOVE IS AFTER FIRST UPDATE for k = %d!!!!!!!!\n", k);
			//// system("PAUSE");

			compute_violations(k, arrival_location);
			//Vehicle[k].Cost.print();
			//printf("ABOVE IS FIRST VIOLATION COMPUTATION bef sequential-slacking for k = %d\n", k);
			//// system("PAUSE");


			//	3. Compute F0.

			//	4. Set D0 : ¼? e0 þ? minfF0;0<p<q Wpg.

			//	5. Update Ai, Wi, Bi and Di for each vertex vi in the route.

			//	6. Compute Li for each request assigned to the route.

			//	7. For everyvertex vj that corresponds to the origin of a request j

			// From every first pickup node to next zero-split node: perform sequential-slacking.

			for (int index = 0; index < Vehicle[k].size; index++)
			{
				if (index != Vehicle[k].size - 1)
				{
					int first_pickup_index = index;

					int zero_split_index;
					for (int slide = index + 1; slide < Vehicle[k].size; slide++)
						if (Request[Vehicle[k].path[slide] - 1].p == 0)
						{
							zero_split_index = slide;
							break;
						}

					//printf("\n\nFIRST PICKUP NODE:	%d, NEXT_ZERO_SPLIT_NODE:	%d\n", Vehicle[k].path[first_pickup_index], Vehicle[k].path[zero_split_index]);
					//// system("PAUSE");

					sequential_slacking(first_pickup_index, zero_split_index, k);
					//printf("Exit sequential slack\n");

					index = zero_split_index;
				}
			}

			// delay departure from the starting depot

			int first_node = Vehicle[k].path[0] - 1;

			Vehicle[k].start_time += Request[first_node].W;
			Request[first_node].A = Request[first_node].B;
			Request[first_node].W = 0;

			//	8. Compute changes in violations of vehicle load, route duration, time window and ride time constraints.

			compute_violations(k, arrival_location);
			//Vehicle[k].Cost.print();

			// calculate excessrideTime
			Vehicle[k].Cost.excessrideTime = 0;
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				int node = Vehicle[k].path[j] - 1;
				if (node >= n)
					Vehicle[k].Cost.excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);
			}

			//printf("End of evaluation for vehicle k = %d, with excessrideTime = %f:!!\n", k, Vehicle[k].Cost.excessrideTime);
			//// system("PAUSE");


			// Modified objective function
			Vehicle[k].Cost.on_board_user_wait_time = 0;
			Vehicle[k].Cost.early_arrival_time = 0;
			Vehicle[k].Cost.excessrideTime = 0;
			Vehicle[k].Cost.route_duration = 0;
			for (int index = 0; index < Vehicle[k].size; index++)
			{
				int node = Vehicle[k].path[index] - 1;

				Vehicle[k].Cost.on_board_user_wait_time += (Request[node].W * (Request[node].p-1));

				if (node < n)
				{
					if (Problem.Request[node].pickup.EarliestTime > Request[node].A)
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node].pickup.EarliestTime - Request[node].A);
				}
				else
				{
					Vehicle[k].Cost.excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);

					if (Problem.Request[node - n].dropoff.EarliestTime > Request[node].A)
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node - n].dropoff.EarliestTime - Request[node].A);
				}

			}
			Vehicle[k].Cost.route_duration = Vehicle[k].end_time - Vehicle[k].start_time;

		}

		Vehicle[k].Cost.overall_ideal_cost = (w1 * Vehicle[k].Cost.travel_cost + w2 * Vehicle[k].Cost.excessrideTime)
							+ 1* Vehicle[k].Cost.load
							+ 1* Vehicle[k].Cost.time_window
							+ 1* Vehicle[k].Cost.ride_time
							+ 1* Vehicle[k].Cost.duration;


	}

	__host__ __device__
		void Offline_route_as_early_as_possible_evaluation(int vehicle_id)
	{
		//here, passed k value starts from 1.

		int vehid = vehicle_id;// -1;
		int k = vehicle_id;// -1;

		Vehicle[k].StartingPoint = 2 * n + k;
		Vehicle[k].EndingPoint = 2 * n + m + k;


		if (Vehicle[k].size < 1)
		{
			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = Problem.Vehicle[k].start_time;
			float depart_time = Vehicle[k].start_time;
			int depart_location = Vehicle[k].StartingPoint;// 2 * n + k;//
			Vehicle[k].end_time = Vehicle[k].start_time + Problem.TravelTime[depart_location][Vehicle[k].EndingPoint];
			Vehicle[k].Cost.duration = max(Vehicle[k].end_time - Vehicle[k].start_time - Problem.Vehicle[k].duration_constraint, 0.0);
			Vehicle[k].Cost.duration += max(Vehicle[k].end_time - Problem.Vehicle[k].end_time, 0.0);
			Vehicle[k].Cost.travel_cost += Problem.TravelTime[depart_location][Vehicle[k].EndingPoint];
		}
		else
		{

			//  1. Set D0 : ¼? e0.

			int depart_location = Vehicle[k].StartingPoint;
			int arrival_location = Vehicle[k].EndingPoint;


			float depart_time = Problem.Vehicle[vehid].start_time;

			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = depart_time;


			//  2. Compute Ai, Wi, Bi and Di and for each vertex vi in the route.

			for (int i = 0; i < Vehicle[k].size; i++)
			{
				int node = Vehicle[k].path[i] - 1;

				//printf("node: %d\n", node);

				Request[node].A = 0;
				Request[node].B = 0;
				Request[node].C = 0;
				Request[node].L = 0;
				Request[node].p = 0;
				Request[node].W = 0;
			}

			//print_with_battery_level();
			//printf("ABOVE IS BEFORE FIRST UPDATE for k = %d!!!!!!!!\n", k);
			//// system("PAUSE");

			// scheduling - start service at every node ASAP

			update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);

			update_load(k);
			update_L_p(k);

			//print_with_battery_level();
			//printf("ABOVE IS AFTER FIRST UPDATE for k = %d!!!!!!!!\n", k);
			//// system("PAUSE");

			compute_violations(k, arrival_location);
			//Vehicle[k].Cost.print();
			//printf("ABOVE IS FIRST VIOLATION COMPUTATION bef sequential-slacking for k = %d\n", k);
			//// system("PAUSE");


			//	3. Compute F0.

			//	4. Set D0 : ¼? e0 þ? minfF0;0<p<q Wpg.

			//	5. Update Ai, Wi, Bi and Di for each vertex vi in the route.

			//	6. Compute Li for each request assigned to the route.

			//	7. For everyvertex vj that corresponds to the origin of a request j

			// From every first pickup node to next zero-split node: perform sequential-slacking.

/*
			for (int index = 0; index < Vehicle[k].size; index++)
			{
				if (index != Vehicle[k].size - 1)
				{
					int first_pickup_index = index;

					int zero_split_index;
					for (int slide = index + 1; slide < Vehicle[k].size; slide++)
						if (Request[Vehicle[k].path[slide] - 1].p == 0)
						{
							zero_split_index = slide;
							break;
						}

					//printf("\n\nFIRST PICKUP NODE:	%d, NEXT_ZERO_SPLIT_NODE:	%d\n", Vehicle[k].path[first_pickup_index], Vehicle[k].path[zero_split_index]);
					//// system("PAUSE");

					sequential_slacking(first_pickup_index, zero_split_index, k);
					//printf("Exit sequential slack\n");

					index = zero_split_index;
				}
			}
			*/
			// delay departure from the starting depot

			int first_node = Vehicle[k].path[0] - 1;

			Vehicle[k].start_time += Request[first_node].W;
			Request[first_node].A = Request[first_node].B;
			Request[first_node].W = 0;

			//	8. Compute changes in violations of vehicle load, route duration, time window and ride time constraints.

			compute_violations(k, arrival_location);
			//Vehicle[k].Cost.print();


			// calculate excessrideTime
			Vehicle[k].Cost.excessrideTime = 0;
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				int node = Vehicle[k].path[j] - 1;
				if (node >= n)
					Vehicle[k].Cost.excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);
			}

			//printf("End of evaluation for vehicle k = %d, with excessrideTime = %f:!!\n", k, Vehicle[k].Cost.excessrideTime);
			//// system("PAUSE");

		}


	}

	__host__ __device__
		void Offline_new_route_evaluation_for_DARP(int vehicle_id)
	{
		//here, passed k value starts from 1.

		int vehid = vehicle_id;// -1;
		int k = vehicle_id;// -1;

		Vehicle[k].StartingPoint = 2 * n + k;
		Vehicle[k].EndingPoint = 2 * n + m + k;


		//  1. Set D0 : ¼? e0.

		int depart_location = Vehicle[k].StartingPoint;
		int arrival_location = Vehicle[k].EndingPoint;


		float depart_time = Problem.Vehicle[vehid].start_time;

		Vehicle[k].Cost.reset();
		Vehicle[k].start_time = depart_time;


		//  2. Compute Ai, Wi, Bi and Di and for each vertex vi in the route.

		for (int i = 0; i < Vehicle[k].size; i++)
		{
			int node = Vehicle[k].path[i] - 1;

			//printf("node: %d\n", node);

			Request[node].A = 0;
			Request[node].B = 0;
			Request[node].C = 0;
			Request[node].L = 0;
			Request[node].p = 0;
			Request[node].W = 0;
		}

		update_A_B_C_W(k, -1, depart_time, depart_location, arrival_location);
		update_load(k);
		update_L_p(k);
		compute_violations(k, arrival_location);

		//	3. Compute F0.
		//	4. Set D0 : ¼? e0 þ? minfF0;0<p<q Wpg.
		//	5. Update Ai, Wi, Bi and Di for each vertex vi in the route.
		//	6. Compute Li for each request assigned to the route.
		//	7. For everyvertex vj that corresponds to the origin of a request j
		// From every first pickup node to next zero-split node: perform sequential-slacking.

		for (int index = 0; index < Vehicle[k].size; index++)
		{
			if (index != Vehicle[k].size - 1)
			{
				int first_pickup_index = index;

				int zero_split_index;
				for (int slide = index + 1; slide < Vehicle[k].size; slide++)
					if (Request[Vehicle[k].path[slide] - 1].p == 0)
					{
						zero_split_index = slide;
						break;
					}

				sequential_slacking(first_pickup_index, zero_split_index, k);

				index = zero_split_index;
			}
		}

		// delay departure from the starting depot

		int first_node = Vehicle[k].path[0] - 1;

		Vehicle[k].start_time += Request[first_node].W;
		Request[first_node].A = Request[first_node].B;
		Request[first_node].W = 0;

		//	8. Compute changes in violations of vehicle load, route duration, time window and ride time constraints.

		compute_violations(k, arrival_location);

		// calculate excessrideTime
		Vehicle[k].Cost.excessrideTime = 0;
		for (int j = 0; j < Vehicle[k].size; j++)
		{
			int node = Vehicle[k].path[j] - 1;
			if (node >= n)
				Vehicle[k].Cost.excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);
		}



		// perform battery evaluation

		//Offline_battery_evaluation(k);

	}

	__host__ __device__
		void find_zero_split_node_indices(int k, int *zero_split_node_index, int &zsp_size)
	{
		zsp_size = 0;

		for (int index = 0; index < Vehicle[k].size; index++)
		{
			int j = Vehicle[k].path[index] - 1;

			if (Request[j].p == 0)
			{
				zero_split_node_index[zsp_size] = index;
				zsp_size++;
			}
		}
	}

	__host__ __device__
		void perform_delay_departure_from_at_pickup_node(int loop, int *zero_split_node_index, int zsp_size, int k)
	{
		// loop => index of first_zero_split_index present in zero_split_node_index[] array.

		//printf("ENTER INTO delay_departure_from_at_pickup_node with loop = %d\n", loop);
		//// system("PAUSE");

		for (int l = zsp_size - 2; l >= loop; l--)
		{
			int j = zero_split_node_index[l + 1] + 1; // second_pick_up_node_index
			int i = zero_split_node_index[l] + 1; // first_pick_up_node_index

												  /*for (int pp = 0; pp < zsp_size; pp++)
												  printf("zero_split_node_index[%d] = %d ==> corr.node. = %d\n", pp, zero_split_node_index[pp], Vehicle[k].path[zero_split_node_index[pp]]);*/

												  //printf("loop = %d, l = %d, ZSP_node_bef_i = %d, ZSP_node_bef_j = %d\n", loop, l, Vehicle[k].path[zero_split_node_index[l]], Vehicle[k].path[zero_split_node_index[l + 1]]);
												  //// system("PAUSE");

												  //printf("first_pick_up_node_index = %d, second_pick_up_node_index = %d\n", i, j);
												  //// system("PAUSE");

			shift_waiting_times_from_j_to_i(i, j, k);

		}

		/*	for (int pp = 0; pp < zsp_size; pp++)
		printf("zero_split_node_index[%d] = %d ==> corr.node. = %d\n", pp, zero_split_node_index[pp], Vehicle[k].path[zero_split_node_index[pp]]);
		printf("EXITED from delay_departure_from_at_pickup_node with loop = %d\n", loop);
		printf("Backward shifted all the waiting times of pick up nodes back to the pickup node %d\n", Vehicle[k].path[zero_split_node_index[loop] + 1]);*/
		//// system("PAUSE");


	}

	__host__ __device__
		void shift_waiting_times_from_j_to_i(int i, int j, int k)
	{
		// i => first_pick_up_node_index node, j => second_pick_up_node_index

		//printf("BELOW: Shifting waiting time from node = %d + 1 to node = %d + 1, k = %d\n", Vehicle[k].path[j -1], Vehicle[k].path[i - 1], k);
		//// system("PAUSE");

		// calculate slack difference: LTW - B from node i to node j-1

		float slack[ExpectedPath];

		vertex Curr_Vertex;
		for (int l = i; l <= j - 1; l++)
		{
			int node = Vehicle[k].path[l] - 1;

			Curr_Vertex = find_vertex(Problem, node);
			slack[l] = max(0, Curr_Vertex.LatestTime - Request[node].B);
		}

		/*printf("nodes	LT	B	slack\n");
		for (int l = i; l <= j - 1; l++)
		{
		int node = Vehicle[k].path[l] - 1;

		Curr_Vertex = find_vertex(Problem, node);
		printf("%d	%f	%f	%f\n", Vehicle[k].path[l], Curr_Vertex.LatestTime, Request[node].B, slack[l]);
		}*/

		// find min slack difference

		float min_slack = slack[i];
		for (int l = i + 1; l <= j - 1; l++)
			if (slack[l] < min_slack)
				min_slack = slack[l];

		//printf("min_slack = %f\n", min_slack);

		// find min slack difference = min(min_slack, request[j].W)

		bool is_j_rep_ending_depot = false;
		if (j - Vehicle[k].size == 0)
			is_j_rep_ending_depot = true;

		float variable_2;
		if (is_j_rep_ending_depot) // if (j ==> index implicating the arrival depot i.e. the node after the final cutomer dropoff)
		{
			variable_2 = max(0, Problem.Vehicle[k].duration_constraint - Vehicle[k].end_time);

			min_slack = min(min_slack, variable_2);
		}
		else
		{
			int node = Vehicle[k].path[j] - 1;
			variable_2 = Request[node].W;

			min_slack = min(min_slack, variable_2);
		}


		//printf("min_slack = %f ==> min(min_slack, %f)\n", min_slack, variable_2);
		//// system("PAUSE");

		// adjust the times at node i (first_pickup_index)

		int node = Vehicle[k].path[i] - 1;

		Request[node].B += min_slack;
		Request[node].C += min_slack;
		Request[node].W += min_slack;

		// adjust the times from node i+1 to j-1 (in-between-nodes)

		for (int l = i + 1; l <= j - 1; l++)
		{
			int node = Vehicle[k].path[l] - 1;

			Request[node].A += min_slack;
			Request[node].B += min_slack;
			Request[node].C += min_slack;
		}

		// adjust the times at node j (zero_split_index)

		if (is_j_rep_ending_depot) // if (j ==> index implicating the arrival depot i.e. the node after the final cutomer dropoff)
		{
			Vehicle[k].end_time += min_slack;
		}
		else
		{
			int node = Vehicle[k].path[j] - 1;

			Request[node].A += min_slack;
			Request[node].W -= min_slack;
		}


		update_L_p(k);

		//print_with_battery_level();
		/*if (!is_j_rep_ending_depot)
		printf("ABOVE is after internal-slack at currentnode %d, between firstpickupnode = %d, and zeroSplitnode at k = %d\n", Vehicle[k].path[j], Vehicle[k].path[i], k);*/

		//// system("PAUSE");

	}

	__host__ __device__
		void OfflineEvaluateCost()
	{
		Feasible_Req_Count = 0;
		total_req_served = 0;

		Temp.reset();
		for (int k = 0; k < TotalVehicles; k++)
		{

			if (Vehicle[k].size == 0)
			{
				Vehicle[k].Cost.reset();
				Vehicle[k].Cost.isFeasible = true;
			}
			else
			{
				//OfflineEvaluateByRoute(k + 1);
				//Offline_route_evaluation(k + 1);
				Offline_mod_route_evaluation(k + 1);
				//Offline_new_route_evaluation(k + 1);
				//Offline_modified_route_evaluation(k + 1);
				Vehicle[k].Cost.getFeasibility();
			}


			Temp.on_board_user_wait_time += Vehicle[k].Cost.on_board_user_wait_time;
			Temp.route_duration += Vehicle[k].Cost.route_duration;
			Temp.early_arrival_time += Vehicle[k].Cost.early_arrival_time;

			Temp.travel_cost += Vehicle[k].Cost.travel_cost;
			Temp.load += Vehicle[k].Cost.load;
			Temp.time_window += Vehicle[k].Cost.time_window;
			Temp.ride_time += Vehicle[k].Cost.ride_time;
			Temp.duration += Vehicle[k].Cost.duration;
			Temp.excessrideTime += Vehicle[k].Cost.excessrideTime;

			Temp.batteryViolation += Vehicle[k].Cost.batteryViolation;
			Temp.max_battery_Violation = max(Temp.max_battery_Violation, Vehicle[k].Cost.batteryViolation);

			Temp.overall_ideal_cost += Vehicle[k].Cost.overall_ideal_cost;
			total_req_served += (Vehicle[k].size / 2);

		}

		Cost = Temp;
		Cost.isFeasible = Cost.getFeasibility();

		// update battery feasibility
		if (Cost.batteryViolation < 0.001)
			Cost.isOverall_batteryfeasible = true;
		else
			Cost.isOverall_batteryfeasible = false;

		for (int k = 0; k < TotalVehicles; k++)
		{
			for (int j = 0; j < Vehicle[k].size; j++)
			{
				// Here we are not checking load violation, considering no violation due to previously feasible route provided by the solver before . THat's makes sense !!!!

				int node = Vehicle[k].path[j] - 1;
				if (node < n)
				{
					if (Request[node].time_window_violation < 0.01)
						if (Request[node].ride_time_violation  < 0.01)
							if (Request[node + n].time_window_violation < 0.01)
								if (Request[node + n].ride_time_violation < 0.01)
									Feasible_Req_Count++;
				}
			}
		}

	}

	__host__ __device__
		cost OfflineFetchRouteCost(int k)
	{
		//here k starts from zero.

		if (Vehicle[k].size == 0)
		{
			Vehicle[k].Cost.reset();
			Vehicle[k].Cost.isFeasible = true;
		}
		else
		{
			OfflineEvaluateByRoute(k + 1);
			Vehicle[k].Cost.getFeasibility();
		}


		return Vehicle[k].Cost;
	}

	__host__ __device__
		void OfflineEvaluateByRoute(int k)
	{
		//here, passed k value starts from 1.

		int vehid = k - 1;
		k--;
		int routesize = Vehicle[k].size;
		int *route = &Vehicle[k].path[0];

		//if (routesize == 0)
		//{
		//	Vehicle[k].Cost.reset();
		//	goto looper;
		//}


		//leverage
		for (int i = 0; i < routesize; i++)
			route[i] -= 1;

		// k = vehicle number = {integer: 1<=k<=m}, where {m = number of vehicle}
		if (routesize<1)
		{
			Vehicle[k].Cost.reset();
			Vehicle[k].start_time = Problem.Vehicle[vehid].start_time;
			float depart_time = Vehicle[k].start_time;
			int depart_location = 2 * n + vehid;
			Vehicle[k].end_time = Vehicle[k].start_time + Problem.TravelTime[depart_location][2 * n + m + vehid];
			Vehicle[k].Cost.duration = max(Vehicle[k].end_time - Vehicle[k].start_time - Problem.Vehicle[vehid].duration_constraint, 0.0);
			Vehicle[k].Cost.duration += max(Vehicle[k].end_time - Problem.Vehicle[vehid].end_time, 0.0);
			Vehicle[k].Cost.travel_cost += Problem.TravelTime[depart_location][2 * n + m + vehid];

			Vehicle[k].Cost.overall_ideal_cost = (w1 * Vehicle[k].Cost.travel_cost + w2 * Vehicle[k].Cost.excessrideTime)
								+ 1* Vehicle[k].Cost.load
								+ 1* Vehicle[k].Cost.time_window
								+ 1* Vehicle[k].Cost.ride_time
								+ 1* Vehicle[k].Cost.duration;
		}
		else
		{

			// 1 Find out who is in the vehicle now (NOT needed for static case)
			int ppl_in_vehicle[n];

			int p = 0;
			for (int i = 0; i<n; i++)
				ppl_in_vehicle[i] = 0;


			// 2 Departure from depot / get depot&start time

			int depart_location = 2 * n + vehid;
			float depart_time = Problem.Vehicle[vehid].start_time;
			Vehicle[k].start_time = depart_time;

			Vehicle[k].Cost.reset();

			int j;
			vertex VertexNow;

			for (int i = 0; i<routesize; i++)
			{
				//cout << "part 3.1\n";
				j = route[i];
				VertexNow = find_vertex(Problem, j);


				Request[j].travel_cost = Problem.TravelTime[depart_location][j];
				Request[j].A = depart_time + Problem.TravelTime[depart_location][j];
				Request[j].B = max(Request[j].A, VertexNow.EarliestTime);
				depart_time = Request[j].B + VertexNow.ServiceTime;
				depart_location = j;
				Request[j].C = depart_time;
				Request[j].W = Request[j].B - Request[j].A;
				//p = p + find_request_load(Problem, j);

				if (j < n)
					p += Problem.Request[j].load;
				else
					p -= Problem.Request[j - n].load;

				Request[j].p = p;

				//cout << "part 3.2\n";

				if (j >= n) //dropoff
				{
					/*if (Problem.Request[j - n].isPickUpDone)
					Request[j].L = Request[j].B - Problem.Request[j - n].departureTime;
					else if (ppl_in_vehicle[j - n])
					Request[j].L = Request[j].B - Request[j - n].C;
					else
					Request[j].L = 0;*/

					Request[j].L = Request[j].B - Request[j - n].C;
					Request[j].ride_time_violation = max((Request[j].L - Problem.Request[j - n].RideTimeConstraint), 0.0);
				}
				else //pickup
				{
					ppl_in_vehicle[j] = 1;
					Request[j].L = 0;
					Request[j].ride_time_violation = 0;
				}

				Request[j].time_window_violation = max((Request[j].B - VertexNow.LatestTime), 0.0);

				if (p > Problem.Vehicle[vehid].capacity)
					Vehicle[k].Cost.load = max(Vehicle[k].Cost.load, p - Problem.Vehicle[vehid].capacity);

				Vehicle[k].Cost.ride_time += Request[j].ride_time_violation;
				Vehicle[k].Cost.time_window += Request[j].time_window_violation;
				Vehicle[k].Cost.travel_cost += Request[j].travel_cost;
			}

			/*
			if (mode == 0)
			rideTimeCompress(k, route[0], routesize, mode);
			else*/

			k++;
			OfflineDurationCompress(k); // Parallel DurationCompression for all possible route_seq on Veh 'k'
			k--;

			if (Request[route[0]].W > 0) //PUT this condition inside expression itself.
			{
				Vehicle[k].start_time += Request[route[0]].W;
				Request[route[0]].W = 0.0;
				Request[route[0]].A = Request[route[0]].B;
			}

			Vehicle[k].end_time = depart_time + Problem.TravelTime[depart_location][2 * n + m + vehid];
			Vehicle[k].Cost.duration = max(Vehicle[k].end_time - Vehicle[k].start_time - Problem.Vehicle[vehid].duration_constraint, 0.0);
			Vehicle[k].Cost.duration += max(Vehicle[k].end_time - Problem.Vehicle[vehid].end_time, 0.0);
			Vehicle[k].Cost.travel_cost += Problem.TravelTime[depart_location][2 * n + m + vehid];


			///////////////////////////////////
			/// Modified objective function ///
			Vehicle[k].Cost.excessrideTime = 0;
			Vehicle[k].Cost.early_arrival_time = 0;
			Vehicle[k].Cost.on_board_user_wait_time = 0;
			Vehicle[k].Cost.route_duration = Vehicle[k].end_time - Vehicle[k].start_time;
			for (int index = 0; index < Vehicle[k].size; index++)
			{
				int node = Vehicle[k].path[index];

				Vehicle[k].Cost.on_board_user_wait_time += Request[node].W * (Request[node].p-1);

				if (node < n)
				{
					if (Problem.Request[node].pickup.EarliestTime > Request[node].A)
						Vehicle[k].Cost.early_arrival_time += Problem.Request[node].pickup.EarliestTime - Request[node].A;
				}
				else
				{
					Vehicle[k].Cost.excessrideTime += Request[node].L - Problem.TravelTime[node][node - n];

					if (Problem.Request[node - n].dropoff.EarliestTime > Request[node].A)
						Vehicle[k].Cost.early_arrival_time += (Problem.Request[node - n].dropoff.EarliestTime - Request[node].A);
				}
			}
			/////////////////////////////////
		}

		//recovery
		for (int i = 0; i < routesize; i++)
			route[i] += 1;

		// calculate excessrideTime
		/*Vehicle[k].Cost.excessrideTime = 0;
		for (int j = 0; j < Vehicle[k].size; j++)
		{
			int node = Vehicle[k].path[j] - 1;
			if (node >= n)
				Vehicle[k].Cost.excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);
		}*/



		Vehicle[k].Cost.overall_ideal_cost = (w1 * Vehicle[k].Cost.travel_cost + w2 * Vehicle[k].Cost.excessrideTime)
							+ 1* Vehicle[k].Cost.load
							+ 1* Vehicle[k].Cost.time_window
							+ 1* Vehicle[k].Cost.ride_time
							+ 1* Vehicle[k].Cost.duration;

	}

	__host__ __device__
		void OfflineDurationCompress(int k)
	{
		k--;
		int routesize_2 = Vehicle[k].size;
		int *route_2 = &Vehicle[k].path[0];

		if (routesize_2>0)
		{
			//int* route = &route_sequence;

			// check who gets in the vehicle
			int ppl_in_vehicle_2[n];
			float all_delay[ExpectedPath];


			for (int i = 0; i<n; i++)
				ppl_in_vehicle_2[i] = 0; // INCLUDED IT LATER...

			for (int i = 0; i<routesize_2; i++)
			{
				if (route_2[i] < n)
					ppl_in_vehicle_2[route_2[i]] = 1;
				//                else
				//                    ppl_in_vehicle[route[i]] = 0;
			}

			for (int i = 0; i<routesize_2; i++)
				all_delay[i] = 0;

			{

				for (int i = routesize_2 - 2; i >= 0; i--)
				{
					float W_prev = Request[route_2[i + 1]].W;
					float E = max(find_vertex(Problem, route_2[i]).LatestTime - Request[route_2[i]].B, 0.0);
					all_delay[i] = min(W_prev + all_delay[i + 1], E);
				}

				bool change = 1;



				while (change)
				{

					//cout << "boom\n";
					change = 0;

					for (int i = 0; i<routesize_2 - 1; i++)
					{
						if (route_2[i] >= n)
						{
							if (Problem.Request[route_2[i] - n].isPickUpDone)
							{
								all_delay[i] = min(all_delay[i], Problem.Request[route_2[i] - n].RideTimeConstraint - Request[route_2[i]].L);
								all_delay[i] = max(all_delay[i], 0.0);
							}
							else
								if (ppl_in_vehicle_2[route_2[i] - n])
								{
									int request_num = route_2[i] - n;
									int j = -1;

									for (int h = i - 1; h >= 0; h--)
										if (route_2[h] == request_num)
										{
											j = h;
											break;
										}

									if (j >= 0)
									{
										all_delay[i] = min(all_delay[i], all_delay[j] + Problem.Request[route_2[i] - n].RideTimeConstraint - Request[route_2[i]].L);
										all_delay[i] = max(all_delay[i], 0.0);
									}
								}
						}
					}

					for (int i = routesize_2 - 2; i >= 0; i--)
					{
						float W_prev = Request[route_2[i + 1]].W;
						float E = max(find_vertex(Problem, route_2[i]).LatestTime - Request[route_2[i]].B, 0.0);
						float delay = min(W_prev + all_delay[i + 1], E);

						if (delay < all_delay[i])
						{
							all_delay[i] = delay;
							change = 1;
						}
					}
				}
			}

			for (int h = 0; h<routesize_2 - 1; h++)
			{
				float possible_delay;
				if (route_2[h]<n)
					possible_delay = all_delay[h];
				else if (Request[route_2[h]].A > Request[route_2[h]].B)
					possible_delay = Request[route_2[h]].A - Request[route_2[h]].B;
				else
					continue;

				Request[route_2[h]].W += possible_delay;

				if (Request[route_2[h]].W < 0)
				{
					possible_delay -= Request[route_2[h]].W;
					Request[route_2[h]].W = 0;
				}

				Request[route_2[h]].B += possible_delay;
				Request[route_2[h]].C += possible_delay;
				Request[route_2[h + 1]].A += possible_delay;
				Request[route_2[h + 1]].W -= possible_delay;


				if (route_2[h] >= n)
				{
					if ((Problem.Request[route_2[h] - n].isPickUpDone) || (ppl_in_vehicle_2[route_2[h] - n]))
					{
						Request[route_2[h]].L += possible_delay;
						float ride_time_violation = max((Request[route_2[h]].L - Problem.Request[route_2[h] - n].RideTimeConstraint), 0.0);

						if (ride_time_violation < 1e-10)
							ride_time_violation = 0.0;

						Vehicle[k].Cost.ride_time += ride_time_violation;
						Vehicle[k].Cost.ride_time -= Request[route_2[h]].ride_time_violation;
						Request[route_2[h]].ride_time_violation = ride_time_violation;
					}
				}
				else
				{
					Request[route_2[h] + n].L -= possible_delay;
					float ride_time_violation = max((Request[route_2[h] + n].L - Problem.Request[route_2[h]].RideTimeConstraint), 0.0);

					if (ride_time_violation < 1e-10)
						ride_time_violation = 0.0;

					Vehicle[k].Cost.ride_time += ride_time_violation;
					Vehicle[k].Cost.ride_time -= Request[route_2[h] + n].ride_time_violation;
					Request[route_2[h] + n].ride_time_violation = ride_time_violation;
				}
			}

			if (Vehicle[k].Cost.ride_time < 1e-10)
				Vehicle[k].Cost.ride_time = 0.0;
		}
	}

	__host__ __device__
		float Calculate_excess_rideTime(int k)
	{
		//here, passed k value starts from ZERO.

		float excessrideTime = 0;

		for (int j = 0; j < Vehicle[k].size; j++)
		{
			int node = Vehicle[k].path[j] - 1;
			if (node >= n)
				excessrideTime += (Request[node].L - Problem.TravelTime[node - n][node]);
		}

		return excessrideTime;
	}

	__host__ __device__
		void RemoveRequest(int req, int k)
	{
		//req = request , k = currvehicle
		//req,k starts from zero.
		//Task: Remove the request & alter the respective vehicle size.


		//printf("remove: %d from %d\n",req,k);
		//spl case-1
		if (Vehicle[k].size == 2)
		{

			//clearing for o/p formatting.
			for (int index = 0; index < Vehicle[k].size; index++)
				Vehicle[k].path[index] = 0;

			Vehicle[k].size -= 2;
			//goto looper;
		}

		//spl case-2
		else if (Vehicle[k].size == 0)
		{
			//goto looper;
		}

		else
		{
			//req i, currvehicle & k starts frm zero
			int *temp = new int[Vehicle[k].size - 2];


			//remove i & store path in temp.
			//for (int index = 0; index < Vehicle[k].size - 2; index++)


			int index = 0;
			int entry = 0;
			while (index < Vehicle[k].size - 2)
			{
				for (int j = 0; j < Vehicle[k].size; j++)
					if (Vehicle[k].path[j] != req + 1)
						if (Vehicle[k].path[j] != req + 1 + n)
							if (Vehicle[k].path[j] != 0)
							{
								temp[index] = Vehicle[k].path[j];
								index++;
							}
				//printf("	loop index: %d\n", index);
				entry++;

				if (entry > 300)
				{
					printf("infinite loop at void RemoveRequest():\n system EXITED!!\n");
					exit(0);
					// system("PAUSE");
				}
			}


			//clearing for o/p formatting.
			for (int index = 0; index < Vehicle[k].size; index++)
				Vehicle[k].path[index] = 0;

			//copy altered path to vehicle.
			for (int index = 0; index < Vehicle[k].size - 2; index++)
			{
				Vehicle[k].path[index] = temp[index];
			}

			Vehicle[k].size -= 2;

			delete[] temp;

		}

	looper:


		////////// Criticial Updates /////////////

		//1) Update New CurrVehicle (Never Mind for here).

		//2) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		//3) Update allinsertion
		allinsertion = 0;
		for (int k = 0; k < TotalVehicles; k++)
			allinsertion += Vehicle[k].insertways;

		///////////////////////////////////////////

		//printf("Do nothing\n");
		//printf("-----\n");
	}

	__host__ __device__
		void RemoveRequest_Spl_Purpose(int req, int k)
	{
		if (Vehicle[k].size == 2)
		{
			//clearing for o/p formatting.
			for (int index = 0; index < 2; index++)
				Vehicle[k].path[index] = 0;

			Vehicle[k].size = 0;
			//goto looper;
		}
		else
		{

			bool BoolOne = 0;
			bool BoolTwo = 0;

			for (int j = 0; j < Vehicle[k].size - 2; j++)
			{

				if (Vehicle[k].path[j] == req + 1)
					BoolOne = 1;
				else if (Vehicle[k].path[j] == req + 1 + n)
					BoolTwo = 1;

				Vehicle[k].path[j] = Vehicle[k].path[j + BoolOne + BoolTwo];
			}
		}


	}

	__host__ __device__
		void LocateRequest(position *Position, int req)
	{
		//req starts from zero.

		Position[0].request = req;
		Position[0].request = Request[req].currvehicle;

		int i = req;
		int k = Request[req].currvehicle;

		for (int j = 0; j < Vehicle[k].size; j++)
		{
			int node = Vehicle[k].path[j] - 1;
			if (node == req)
				Position[0].start = j;

			if (node == req + n)
			{
				Position[0].gap = j - Position[0].start;
				break;
			}
		}
	}

	__host__ __device__
		void ReInsertRequest(int REQ, int VEH, int starting, int ending)
	{
		//printf("reinsertion: %d into %d at pos(%d,%d)\n", REQ, VEH, starting, ending);
		//Align

					//printf("Here 0\n");

		int req = REQ;
		int toVeh = VEH;
		int start = starting;
		int end = ending;

		//here, req and toVeh starts from zero.
		int k = toVeh;

					//printf("Here 1\n");

		//backup the path
		int *temp = new int[Vehicle[k].size];
		for (int j = 0; j < Vehicle[k].size; j++)
			temp[j] = Vehicle[k].path[j];

					/*printf("size: %d Path: ", Vehicle[k].size);
					for (int j = 0; j < Vehicle[k].size; j++)
						printf("%d ",temp[j]);
					printf("\n");

					printf("Here 2\n");*/

		//increment size
		Vehicle[k].size += 2;

		//clear up
		for (int i = 0; i < Vehicle[k].size; i++)
			Vehicle[k].path[i] = 0;

					/*printf("size: %d Path: ", Vehicle[k].size);
					for (int j = 0; j < Vehicle[k].size; j++)
						printf("%d ",Vehicle[k].path[j]);
					printf("\n");

					printf("Here 3\n");*/

		//copy
		Vehicle[k].path[start] = req + 1;
		Vehicle[k].path[end] = req + 1 + n;

					/*printf("size: %d Path: ", Vehicle[k].size);
					for (int j = 0; j < Vehicle[k].size; j++)
						printf("%d ",Vehicle[k].path[j]);
					printf("\n");

					printf("Here 4\n");*/

		//copyback temp
		int j = 0;
		for (int i = 0; i < Vehicle[k].size; i++)
			if (Vehicle[k].path[i] == 0)
			{
				Vehicle[k].path[i] = temp[j];
				j++;
			}

					/*printf("size: %d Path: ", Vehicle[k].size);
					for (int j = 0; j < Vehicle[k].size; j++)
						printf("%d ",Vehicle[k].path[j]);
					printf("\n");

					printf("Here 5\n");*/

		delete[] temp;

					//printf("Here 6\n");

		////////// Criticial Updates /////////////

		//1) Update New Current Vehicle
		Request[req].currvehicle = toVeh;

		//2) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		//3) Update allinsertion
		allinsertion = 0;
		for (int k = 0; k < TotalVehicles; k++)
			allinsertion += Vehicle[k].insertways;

					//printf("Here 7\n");

		///////////////////////////////////////////
		//printf("-----\n");
	}

	__host__ __device__
		void ReInsertRequest_PickUp_Point(int REQ, int VEH, int starting)
	{
		//printf("reinsertion: %d into %d at pos(%d,%d)\n", REQ, VEH, starting, ending);
		//Align
		int req = REQ;
		int toVeh = VEH;
		int start = starting;

		//here, req and toVeh starts from zero.
		int k = toVeh;


		//backup the path
		int *temp = new int[Vehicle[k].size];
		for (int j = 0; j < Vehicle[k].size; j++)
			temp[j] = Vehicle[k].path[j];

		//increment size
		Vehicle[k].size += 1;

		//clear up
		for (int i = 0; i < Vehicle[k].size; i++)
			Vehicle[k].path[i] = 0;

		//copy
		Vehicle[k].path[start] = req + 1;

		//copyback temp
		int j = 0;
		for (int i = 0; i < Vehicle[k].size; i++)
			if (Vehicle[k].path[i] == 0)
			{
				Vehicle[k].path[i] = temp[j];
				j++;
			}

		delete[] temp;



		////////// Criticial Updates /////////////

		//1) Update New Current Vehicle
		Request[req].currvehicle = toVeh;

		//2) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		//3) Update allinsertion
		allinsertion = 0;
		for (int k = 0; k < TotalVehicles; k++)
			allinsertion += Vehicle[k].insertways;

		///////////////////////////////////////////
		//printf("-----\n");
	}

	__host__ __device__
		void ReInsertRequest_DropOff_Point(int REQ, int VEH, int ending)
	{
		//printf("reinsertion: %d into %d at pos(%d,%d)\n", REQ, VEH, starting, ending);
		//Align
		int req = REQ;
		int toVeh = VEH;
		int end = ending;

		//here, req and toVeh starts from zero.
		int k = toVeh;


		//backup the path
		int *temp = new int[Vehicle[k].size];
		for (int j = 0; j < Vehicle[k].size; j++)
			temp[j] = Vehicle[k].path[j];

		//increment size
		Vehicle[k].size += 1;

		//clear up
		for (int i = 0; i < Vehicle[k].size; i++)
			Vehicle[k].path[i] = 0;

		//copy
		Vehicle[k].path[end] = req + 1 + n;

		//copyback temp
		int j = 0;
		for (int i = 0; i < Vehicle[k].size; i++)
			if (Vehicle[k].path[i] == 0)
			{
				Vehicle[k].path[i] = temp[j];
				j++;
			}

		delete[] temp;



		////////// Criticial Updates /////////////

		//1) Update New Current Vehicle
		Request[req].currvehicle = toVeh;

		//2) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		//3) Update allinsertion
		allinsertion = 0;
		for (int k = 0; k < TotalVehicles; k++)
			allinsertion += Vehicle[k].insertways;

		///////////////////////////////////////////
		//printf("-----\n");
	}

	__host__ __device__
		void ReInsertRequest_PickUp_Point_Reverse(int REQ, int VEH, int index_location)
	{
		//printf("reinsertion: %d into %d at pos(%d,%d)\n", REQ, VEH, starting, ending);
		//Align
		int req = REQ;
		int toVeh = VEH;

		//here, req and toVeh starts from zero.
		int k = toVeh;


		//backup the path
		int *temp = new int[Vehicle[k].size];
		for (int j = 0; j < Vehicle[k].size; j++)
			temp[j] = Vehicle[k].path[j];

		//increment size
		Vehicle[k].size += 1;

		//clear up
		for (int i = 0; i < Vehicle[k].size; i++)
			Vehicle[k].path[i] = 0;

		//copy /////////////////////
		Vehicle[k].path[index_location] = req + 1 + n; // insertion of dropOff first (after which, pickup will be inserted!)

		//copyback temp
		int j = 0;
		for (int i = 0; i < Vehicle[k].size; i++)
			if (Vehicle[k].path[i] == 0)
			{
				Vehicle[k].path[i] = temp[j];
				j++;
			}

		delete[] temp;



		////////// Criticial Updates /////////////

		//1) Update New Current Vehicle
		Request[req].currvehicle = toVeh;

		//2) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		//3) Update allinsertion
		allinsertion = 0;
		for (int k = 0; k < TotalVehicles; k++)
			allinsertion += Vehicle[k].insertways;

		///////////////////////////////////////////
		//printf("-----\n");
	}

	__host__ __device__
		void ReInsertRequest_DropOff_Point_Reverse(int REQ, int VEH, int index_location)
	{
		//printf("reinsertion: %d into %d at pos(%d,%d)\n", REQ, VEH, starting, ending);
		//Align
		int req = REQ;
		int toVeh = VEH;

		//here, req and toVeh starts from zero.
		int k = toVeh;


		//backup the path
		int *temp = new int[Vehicle[k].size];
		for (int j = 0; j < Vehicle[k].size; j++)
			temp[j] = Vehicle[k].path[j];

		//increment size
		Vehicle[k].size += 1;

		//clear up
		for (int i = 0; i < Vehicle[k].size; i++)
			Vehicle[k].path[i] = 0;

		//copy //////////////////////////
		Vehicle[k].path[index_location] = req + 1; // insertion of pickup (after dropoff)

		//copyback temp
		int j = 0;
		for (int i = 0; i < Vehicle[k].size; i++)
			if (Vehicle[k].path[i] == 0)
			{
				Vehicle[k].path[i] = temp[j];
				j++;
			}

		delete[] temp;



		////////// Criticial Updates /////////////

		//1) Update New Current Vehicle
		Request[req].currvehicle = toVeh;

		//2) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		//3) Update allinsertion
		allinsertion = 0;
		for (int k = 0; k < TotalVehicles; k++)
			allinsertion += Vehicle[k].insertways;

		///////////////////////////////////////////
		//printf("-----\n");
	}

	__host__ __device__
		void InsertRequest(insertion *Insertion)
	{

		//Align
		int req = Insertion[0].request;
		int fromVeh = Request[Insertion[0].request].currvehicle;
		int toVeh = Insertion[0].vehid;
		int start = Insertion[0].start;
		int gap = Insertion[0].gap;

		//here, req, fromVeh & toVeh starts from zero.
		int k = toVeh;

		RemoveRequest(req, fromVeh);

		//backup the path
		int *temp = new int[Vehicle[k].size];
		for (int j = 0; j < Vehicle[k].size; j++)
			temp[j] = Vehicle[k].path[j];

		//increment size
		Vehicle[k].size += 2;

		//clear up
		for (int i = 0; i < Vehicle[k].size; i++)
			Vehicle[k].path[i] = 0;

		//copy
		Vehicle[k].path[start] = req + 1;
		Vehicle[k].path[start + gap + 1] = req + 1 + n;

		//copyback temp
		int j = 0;
		for (int i = 0; i < Vehicle[k].size; i++)
			if (Vehicle[k].path[i] == 0)
			{
				Vehicle[k].path[i] = temp[j];
				j++;
			}

		delete[] temp;


		////////// Criticial Updates /////////////

		//1) Update New Current Vehicle
		Request[req].currvehicle = toVeh;

		//2) Vehicle - Insertways
		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].insertways = (Vehicle[k].size + 1) * (Vehicle[k].size + 2) / 2;

		//3) Update allinsertion
		allinsertion = 0;
		for (int k = 0; k < TotalVehicles; k++)
			allinsertion += Vehicle[k].insertways;

		///////////////////////////////////////////

	}

	__host__ __device__
		~solution() {}
};
