
__host__ __device__
location find_location(problem& Problem, int i)
{
	location Location;

	if (i < n)
	{
		Location = Problem.Request[i].pickup.Location;
	}
	else if (i < 2 * n)
	{
		Location = Problem.Request[i - n].dropoff.Location;
	}
	else if (i < 2 * n + m)
	{
		Location = Problem.Vehicle[i - 2 * n].start;
	}
	else
	{
		Location = Problem.Vehicle[i - 2 * n - m].end;
	}

	return Location;
}

__host__ __device__
vertex find_vertex(problem& Problem, int i)
{
	// the node value i represents req => starts from 0

	vertex Vertex;

	if (i < n)
		Vertex = Problem.Request[i].pickup;
	else
		Vertex = Problem.Request[i - n].dropoff;


/*
	bool coeff = i < n;
	Vertex.EarliestTime = (coeff * Problem.Request[i].pickup.EarliestTime) + ((1-coeff) * Problem.Request[i - n].dropoff.EarliestTime);
	Vertex.LatestTime = (coeff * Problem.Request[i].pickup.LatestTime) + ((1-coeff) * Problem.Request[i - n].dropoff.LatestTime);
	Vertex.ServiceTime = (coeff * Problem.Request[i].pickup.ServiceTime) + ((1-coeff) * Problem.Request[i - n].dropoff.ServiceTime);
	Vertex.Location.ID = (coeff * Problem.Request[i].pickup.Location.ID) + ((1-coeff) * Problem.Request[i - n].dropoff.Location.ID);
	Vertex.Location.x = (coeff * Problem.Request[i].pickup.Location.x) + ((1-coeff) * Problem.Request[i - n].dropoff.Location.x);
	Vertex.Location.y = (coeff * Problem.Request[i].pickup.Location.y) + ((1-coeff) * Problem.Request[i - n].dropoff.Location.y);
*/

	return Vertex;
}

__host__ __device__
int find_request_load(problem& Problem, int i)
{
	int Load;

	if (i < n)
		Load = Problem.Request[i].load;
	else
		Load = -Problem.Request[i - n].load;

	return Load;
}

struct container
{
	int start, gap, vehid;

	__host__ __device__
		void reset()
	{
		start = 0;
		gap = 0;
		vehid = 0;
	}
};

struct holder
{
	int ID;
	int total_insertions;

	int start[Expectation];
	int gap[Expectation];

	int vehid; //vehicleid, starts from 0.
	int vehid_size; //vehicle size

	holder()
	{
		reset();
	}

	~holder() {}

	__host__ __device__
	void reset()
	{
		ID = 0; total_insertions = 0; vehid = 0; vehid_size = 0;
	}

	__host__ __device__
	void initialize(int vehSize, int totInsertions)
	{
		if (vehSize == 0)
		{
			start[0] = 0;
			gap[0] = 0;
			total_insertions = 1;
			vehid = 0;
		}
		else if (vehSize > 0)
		{
			ID = vehSize;
			vehid_size = vehSize;
			total_insertions = totInsertions;
			vehid = 0; //can be changed anytime.


			/*int indices = 0;
			for (int limiter = size; limiter >= 0; limiter--)
				for (int i_start = 0; i_start <= limiter; i_start++)
				{
					start[indices] = i_start;
					gap[indices] = size - limiter;
					indices++;
				}*/

			int indices = 0; int i_start = 0;
			for (int bound = vehid_size+1; bound>0; bound--, i_start++)
			{
				for (int i_gap = 0; i_gap < bound; i_gap++, indices++)
				{
					start[indices] = i_start;
					gap[indices] = i_gap;

				}
			}


			// print utility
			if (0)
			{
				printf("---------------\n");
				if (total_insertions == indices)
					printf("All is well\n");

				printf("holder ID: %d, vehid_size: %d, total_ins: %d\n", ID, vehid_size, total_insertions);

				printf("str: ");
				indices = 0, i_start = 0;
				for (int bound = vehid_size+1; bound>0; bound--, i_start++)
				{
					for (int i_gap = 0; i_gap < bound; i_gap++, indices++)
						printf("%d ", start[indices]);
					printf("|");
				}

				printf("\ngap: ");
				indices = 0, i_start = 0;
				for (int bound = vehid_size+1; bound>0; bound--, i_start++)
				{
					for (int i_gap = 0; i_gap < bound; i_gap++, indices++)
						printf("%d ", gap[indices]);
					printf("|");
				}


				printf("\n---------------\n");
			}


		}

	}
};

struct skeleton
{
	int req;
	int vehid;
	int start;
	int gap;

	skeleton() {}

	~skeleton() {}
};

struct sol
{
	int solID;
	vehicle_detail Vehicle;
	//request_detail Request[2 * TotalRequests];
	float current_cost;

	//int entrycount;

	bool isFeasible;

	bool ppl_in_vehicle[n];
	//float all_delay[ExpectedPath];

	int route[ExpectedPath];
	int routesize;

	int request; //starts from zero here.
	int vehid; //starts from zero here.

	int start; //best_i
	int gap; //best_j


	__host__ __device__
		void clean()
	{
		solID = 0;
		Vehicle.Cost.reset();
		//for (int i = 0; i < 2 * TotalRequests; i++)
		//	Request[i].clean();
		current_cost = 0;

		isFeasible = 0;

		for (int i = 0; i < n; i++)
			ppl_in_vehicle[i] = 0;

		/*for (int i = 0; i < ExpectedPath; i++)
		all_delay[i] = 0;*/

		for (int i = 0; i < ExpectedPath; i++)
			route[i] = 0;
		routesize = 0;

		request = 0;
		vehid = 0;

		start = 0;
		gap = 0;
	}

	__host__ __device__ sol() {}

	__host__ __device__
		void fixation(int size)
	{
		current_cost = 0;
		routesize = size + 2;
		for (int i = 0; i < routesize; i++)
			route[i] = 0;
	}

	__host__ __device__
		~sol() {}

	__host__ __device__
		void Copy(sol &inSol)
	{
		solID = inSol.solID;
		Vehicle.Copy(inSol.Vehicle);
		current_cost = inSol.current_cost;

		isFeasible = inSol.isFeasible;

		for (int i = 0; i < n; i++)
			ppl_in_vehicle[i] = inSol.ppl_in_vehicle[i];


		routesize = inSol.routesize;
		for (int i = 0; i < routesize; i++)
			route[i] = inSol.route[i];

		request = inSol.request;
		vehid = inSol.vehid;

		start = inSol.start;
		gap = inSol.gap;
	}

	__host__ __device__
		void evaluateByRoute_SM(problem *Problem)
	{

		int i, j;

		int routesequence[ExpectedPath];
		int size = routesize;

		//leverage
		for (i = 0; i < size; i++)
			routesequence[i] = route[i] - 1;

		// k = vehicle number = {integer: 1<=k<=m}, where {m = number of vehicle}
/*		if (size<1)
		{
			Vehicle.Cost.reset();
			Vehicle.start_time = Problem[0].Vehicle[vehid].start_time;
			float depart_time = Vehicle.start_time;
			int depart_location = 2 * n + vehid;
			Vehicle.end_time = Vehicle.start_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];
		}
		else*/
		{
			float RTV[2 * TotalRequests];
			float TwV[2 * TotalRequests];
			float tc[2 * TotalRequests];
			int Req_p[2 * TotalRequests];

			float A[2 * TotalRequests],
				B[2 * TotalRequests],
				C[2 * TotalRequests],
				D[2 * TotalRequests],
				W[2 * TotalRequests],
				L[2 * TotalRequests];

			float V_RT = 0;
			float V_TW = 0;
			float V_tc = 0;


			//vertex VertexNow;
			float EarliestTime;
			float ServiceTime;
			float LatestTime;
			int load;

			// 1 Find out who is in the vehicle now (NOT needed for static case)

			int p = 0;
			// 2 Departure from depot / get depot&start time

			int depart_location = 2 * n + vehid;
			float depart_time = Problem[0].Vehicle[vehid].start_time;
			Vehicle.start_time = depart_time;

			Vehicle.Cost.reset();

			bool coeff;

			for (i = 0; i<size; i++)
			{
				//cout << "part 3.1\n";
				j = routesequence[i];
				//VertexNow = find_vertex(Problem[0], j);

				coeff = j < n;

				EarliestTime = (coeff)*  Problem[0].Request[j].pickup.EarliestTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.EarliestTime;
				ServiceTime = (coeff)*  Problem[0].Request[j].pickup.ServiceTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.ServiceTime;
				LatestTime = (coeff)*  Problem[0].Request[j].pickup.LatestTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.LatestTime;
				load = (coeff)*  Problem[0].Request[j].load; + (1 - coeff) * -Problem[0].Request[j - n].load;



				/*if (j < n)
				{
					EarliestTime = Problem[0].Request[j].pickup.EarliestTime;
					ServiceTime = Problem[0].Request[j].pickup.ServiceTime;
					LatestTime = Problem[0].Request[j].pickup.LatestTime;
					load = Problem[0].Request[j].load;
				}
				else
				{
					EarliestTime = Problem[0].Request[j - n].dropoff.EarliestTime;
					ServiceTime = Problem[0].Request[j - n].dropoff.ServiceTime;
					LatestTime = Problem[0].Request[j - n].dropoff.LatestTime;
					load = -Problem[0].Request[j - n].load;
				}*/


				tc[j] = Problem[0].TravelTime[depart_location][j];
				A[j] = depart_time + Problem[0].TravelTime[depart_location][j];
				B[j] = max(A[j], EarliestTime);
				depart_time = B[j] + ServiceTime;
				depart_location = j;
				C[j] = depart_time;
				W[j] = B[j] - A[j];
				p = p + load;
				Req_p[j] = p;

				//cout << "part 3.2\n";

				coeff = j >= n;

				L[j] = (coeff) * (B[j] - C[j - n]);
				RTV[j] = (coeff)*  max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
				ppl_in_vehicle[j] = (coeff) * 0 + (1 - coeff) * 1;


				/*if (j >= n) //dropoff
				{
					L[j] = B[j] - C[j - n];
					RTV[j] = max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
				}
				else //pickup
				{
					ppl_in_vehicle[j] = 1;
					L[j] = 0;
					RTV[j] = 0;
				}*/

				TwV[j] = max((B[j] - LatestTime), 0.0);

				if (p > Problem[0].Vehicle[vehid].capacity)
					Vehicle.Cost.load = max(Vehicle.Cost.load, p - Problem[0].Vehicle[vehid].capacity);

				/*RTV[j] = Request[j].ride_time_violation;
				TwV[j] = Request[j].time_window_violation;
				tc[j] = Request[j].travel_cost;*/

				/*Vehicle.Cost.ride_time += Request[j].ride_time_violation;
				Vehicle.Cost.time_window += Request[j].time_window_violation;
				Vehicle.Cost.travel_cost += Request[j].travel_cost;*/
			}

			/*for (i = 0; i < size; i++)
			{
			j = routesequence[i];
			Request[j].A = A[j];
			Request[j].B = B[j];
			Request[j].C = C[j];
			Request[j].W = W[j];
			Request[j].L = L[j];
			}*/



			for (i = 0; i < size; i++)
			{
				j = routesequence[i];
				V_RT += RTV[j];
				V_TW += TwV[j];
				V_tc += tc[j];
			}

			Vehicle.Cost.ride_time = V_RT;
			Vehicle.Cost.time_window = V_TW;
			Vehicle.Cost.travel_cost = V_tc;


			//durationCompress_SM(Problem, routesequence, size, A, B, C, W, L, V_RT, RTV); // Parallel DurationCompression for all possible route_seq on Veh 'k'

			if (W[routesequence[0]] > 0) //PUT this condition inside expression itself.
			{
				Vehicle.start_time += W[routesequence[0]];
				W[routesequence[0]] = 0.0;
				A[routesequence[0]] = B[routesequence[0]];
			}

			Vehicle.end_time = depart_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];

			// Modified objective function
			/*
			float excess_ride_time = 0, passenger_waiting_time = 0, route_duration = 0, early_arrival_time = 0;

			for (i = 0; i < size; i++)
			{
			int node = routesequence[i];
			if (node < n)
			{
			passenger_waiting_time += W[node] * (Req_p[node] - Problem[0].Request[node].load);
			early_arrival_time += (Problem[0].Request[node].pickup.EarliestTime - A[node]);
			}
			else
			{
			excess_ride_time += B[node] - D[node - n] - Problem[0].TravelTime[node][node - n];

			passenger_waiting_time += W[node] * (Req_p[node] - (-Problem[0].Request[node - n].load));
			early_arrival_time += (Problem[0].Request[node].dropoff.EarliestTime - A[node]);
			}

			}
			route_duration = Vehicle.end_time - Vehicle.start_time;
			*/

			//recovery
			//for (int i = 0; i < routesize; i++)
				//routesequence[i] += 1;

			// calculate excessrideTime
			Vehicle.Cost.excessrideTime = 0;
			for (int j = 0; j < routesize; j++)
			{
				int node = route[j] - 1;
				if (node >= n)
					Vehicle.Cost.excessrideTime += (L[node] - Problem[0].TravelTime[node - n][node]);
			}
		}

		Vehicle.Cost.getFeasibility();
	}

	__host__ __device__
		void durationCompress_SM(problem *Problem, int *routesequence, int size, float *A, float *B, float *C, float *W, float *L, float &V_RT, float *RTV)
	{

			float all_delay[ExpectedPath];

			float W_prev, E, delay;
			int i, j, h, request_num;
			float possible_delay, ride_time_violation;

			for (i = 0; i <size; i++)
				all_delay[i] = 0;


			for (i = size - 2; i >= 0; i--)
			{
				W_prev = W[routesequence[i + 1]];
				E = max(find_vertex(Problem[0], routesequence[i]).LatestTime - B[routesequence[i]], 0.0);
				all_delay[i] = min(W_prev + all_delay[i + 1], E);
			}

			for (i = 0; i<size - 1; i++)
			{
				if (routesequence[i] >= n)
				{
					// i = dropoff, j = pickup

					request_num = routesequence[i] - n;
					j = -1;

					for (h = i - 1; h >= 0; h--)
						if (routesequence[h] == request_num)
						{
							j = h;
							break;
						}

					if (j >= 0)
					{
						all_delay[i] = min(all_delay[i], all_delay[j] + Problem[0].Request[routesequence[i] - n].RideTimeConstraint - L[routesequence[i]]);
						all_delay[i] = max(all_delay[i], 0.0);
					}

				}
			}

			for (i = size - 2; i >= 0; i--)
			{
				W_prev = W[routesequence[i + 1]];
				E = max(find_vertex(Problem[0], routesequence[i]).LatestTime - B[routesequence[i]], 0.0);
				delay = min(W_prev + all_delay[i + 1], E);

				bool coeff = delay < all_delay[i];

				all_delay[i] = coeff * delay;

			}



			for (h = 0; h<size - 1; h++)
			{

				if (routesequence[h]<n)
					possible_delay = all_delay[h];
				else if (A[routesequence[h]] > B[routesequence[h]])
					possible_delay = A[routesequence[h]] - B[routesequence[h]];
				else
					continue;

				W[routesequence[h]] += possible_delay;

				if (W[routesequence[h]] < 0)
				{
					possible_delay -= W[routesequence[h]];
					W[routesequence[h]] = 0;
				}

				B[routesequence[h]] += possible_delay;
				C[routesequence[h]] += possible_delay;
				A[routesequence[h + 1]] += possible_delay;
				W[routesequence[h + 1]] -= possible_delay;

				if (routesequence[h] >= n)
				{
					L[routesequence[h]] += possible_delay;
					ride_time_violation = max((L[routesequence[h]] - Problem[0].Request[routesequence[h] - n].RideTimeConstraint), 0.0);

					if (ride_time_violation < 1e-10)
						ride_time_violation = 0.0;

					Vehicle.Cost.ride_time += ride_time_violation;
					Vehicle.Cost.ride_time -= RTV[routesequence[h]];
					RTV[routesequence[h]] = ride_time_violation;
				}
				else
				{
					L[routesequence[h] + n] -= possible_delay;
					ride_time_violation = max((L[routesequence[h] + n] - Problem[0].Request[routesequence[h]].RideTimeConstraint), 0.0);

					if (ride_time_violation < 1e-10)
						ride_time_violation = 0.0;

					Vehicle.Cost.ride_time += ride_time_violation;
					Vehicle.Cost.ride_time -= RTV[routesequence[h] + n];
					RTV[routesequence[h] + n] = ride_time_violation;
				}
			}

			if (Vehicle.Cost.ride_time < 1e-10)
				Vehicle.Cost.ride_time = 0.0;

	}

	__host__ __device__
		void evaluateByRoute_standardDARP(problem *Problem)
	{

		int i, j;

		int routesequence[ExpectedPath];
		int size = routesize;

		//leverage
		for (i = 0; i < size; i++)
			routesequence[i] = route[i] - 1;

		// k = vehicle number = {integer: 1<=k<=m}, where {m = number of vehicle}
		if (size<1)
		{
			Vehicle.Cost.reset();
			Vehicle.start_time = Problem[0].Vehicle[vehid].start_time;
			float depart_time = Vehicle.start_time;
			int depart_location = 2 * n + vehid;
			Vehicle.end_time = Vehicle.start_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];
		}
		else
		{
			float RTV[2 * TotalRequests];
			float TwV[2 * TotalRequests];
			float tc[2 * TotalRequests];
			int Req_p[2 * TotalRequests];

			float A[2 * TotalRequests],
				B[2 * TotalRequests],
				C[2 * TotalRequests],
				D[2 * TotalRequests],
				W[2 * TotalRequests],
				L[2 * TotalRequests];

			float V_RT = 0;
			float V_TW = 0;
			float V_tc = 0;


			//vertex VertexNow;
			float EarliestTime;
			float ServiceTime;
			float LatestTime;
			int load;

			// 1 Find out who is in the vehicle now (NOT needed for static case)

			int p = 0;
			// 2 Departure from depot / get depot&start time

			int depart_location = 2 * n + vehid;
			float depart_time = Problem[0].Vehicle[vehid].start_time;
			Vehicle.start_time = depart_time;

			Vehicle.Cost.reset();

			bool coeff;

			for (i = 0; i<size; i++)
			{
				//cout << "part 3.1\n";
				j = routesequence[i];
				//VertexNow = find_vertex(Problem[0], j);

				coeff = j < n;

				/*EarliestTime = (coeff)*  Pblm[0].Request[j].pickup.EarliestTime + (1 - coeff) * Pblm[0].Request[j - n].dropoff.EarliestTime;
				ServiceTime = (coeff)*  Pblm[0].Request[j].pickup.ServiceTime + (1 - coeff) * Pblm[0].Request[j - n].dropoff.ServiceTime;
				LatestTime = (coeff)*  Pblm[0].Request[j].pickup.LatestTime + (1 - coeff) * Pblm[0].Request[j - n].dropoff.LatestTime;
				load = (coeff)*  Pblm[0].Request[j].load; + (1 - coeff) * -Pblm[0].Request[j - n].load;*/


				if (j < n)
				{
					EarliestTime = Problem[0].Request[j].pickup.EarliestTime;
					ServiceTime = Problem[0].Request[j].pickup.ServiceTime;
					LatestTime = Problem[0].Request[j].pickup.LatestTime;
					load = Problem[0].Request[j].load;
				}
				else
				{
					EarliestTime = Problem[0].Request[j - n].dropoff.EarliestTime;
					ServiceTime = Problem[0].Request[j - n].dropoff.ServiceTime;
					LatestTime = Problem[0].Request[j - n].dropoff.LatestTime;
					load = -Problem[0].Request[j - n].load;
				}


				tc[j] = Problem[0].TravelTime[depart_location][j];
				A[j] = depart_time + Problem[0].TravelTime[depart_location][j];
				B[j] = max(A[j], EarliestTime);
				depart_time = B[j] + ServiceTime;
				depart_location = j;
				C[j] = depart_time;
				W[j] = B[j] - A[j];
				p = p + load;
				Req_p[j] = p;

				//cout << "part 3.2\n";

				/*coeff = j >= n;

				L[j] = (coeff) * (B[j] - C[j - n]);
				RTV[j] = (coeff)*  max((L[j] - Pblm[0].Request[j - n].RTC), 0.0);
				ppl_in_vehicle[j] = (coeff) * 0 + (1 - coeff) * 1;*/


				if (j >= n) //dropoff
				{
					L[j] = B[j] - C[j - n];
					RTV[j] = max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
				}
				else //pickup
				{
					ppl_in_vehicle[j] = 1;
					L[j] = 0;
					RTV[j] = 0;
				}

				TwV[j] = max((B[j] - LatestTime), 0.0);

				if (p > Problem[0].Vehicle[vehid].capacity)
					Vehicle.Cost.load = max(Vehicle.Cost.load, p - Problem[0].Vehicle[vehid].capacity);

				/*RTV[j] = Request[j].ride_time_violation;
				TwV[j] = Request[j].time_window_violation;
				tc[j] = Request[j].travel_cost;*/

				/*Vehicle.Cost.ride_time += Request[j].ride_time_violation;
				Vehicle.Cost.time_window += Request[j].time_window_violation;
				Vehicle.Cost.travel_cost += Request[j].travel_cost;*/
			}

			/*for (i = 0; i < size; i++)
			{
			j = routesequence[i];
			Request[j].A = A[j];
			Request[j].B = B[j];
			Request[j].C = C[j];
			Request[j].W = W[j];
			Request[j].L = L[j];
			}*/



			for (i = 0; i < size; i++)
			{
				j = routesequence[i];
				V_RT += RTV[j];
				V_TW += TwV[j];
				V_tc += tc[j];
			}


			Vehicle.Cost.ride_time = V_RT;
			Vehicle.Cost.time_window = V_TW;
			Vehicle.Cost.travel_cost = V_tc;





			durationCompress(/*vehid - 1, route[0], routesize,*/ Problem, routesequence, size, A, B, C, W, L, V_RT, RTV); // Parallel DurationCompression for all possible route_seq on Veh 'k'

			if (W[routesequence[0]] > 0) //PUT this condition inside expression itself.
			{
				Vehicle.start_time += W[routesequence[0]];
				W[routesequence[0]] = 0.0;
				A[routesequence[0]] = B[routesequence[0]];
			}

			Vehicle.end_time = depart_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];

			// Modified objective function

			for (int index = 0; index < size; index++)
			{
				int node = routesequence[index];

				Vehicle.Cost.on_board_user_wait_time += (W[node] * (Req_p[node] - 1));

				if (node < n)
				{
					if (Problem[0].Request[node].pickup.EarliestTime > A[node])
						Vehicle.Cost.early_arrival_time += (Problem[0].Request[node].pickup.EarliestTime - A[node]);

				}
				else
				{
					Vehicle.Cost.excessrideTime += (L[node] - Problem[0].TravelTime[node][node - n]);

					if (Problem[0].Request[node - n].dropoff.EarliestTime > A[node])
						Vehicle.Cost.early_arrival_time += (Problem[0].Request[node - n].dropoff.EarliestTime - A[node]);
				}
			}
			Vehicle.Cost.route_duration = Vehicle.end_time - Vehicle.start_time;
		}
		Vehicle.Cost.getFeasibility();
		Vehicle.Cost.overall_ideal_cost = (w1 * Vehicle.Cost.travel_cost + w2 * Vehicle.Cost.excessrideTime)
				+ 1* Vehicle.Cost.load
				+ 1* Vehicle.Cost.time_window
				+ 1* Vehicle.Cost.ride_time
				+ 1* Vehicle.Cost.duration;
	}

	__host__ __device__
		void evaluateByRoute(problem *Problem)
	{

		int i, j;

		int routesequence[ExpectedPath];
		int size = routesize;

		//leverage
		for (i = 0; i < size; i++)
			routesequence[i] = route[i] - 1;

		// k = vehicle number = {integer: 1<=k<=m}, where {m = number of vehicle}
		if (size<1)
		{
			Vehicle.Cost.reset();
			Vehicle.start_time = Problem[0].Vehicle[vehid].start_time;
			float depart_time = Vehicle.start_time;
			int depart_location = 2 * n + vehid;
			Vehicle.end_time = Vehicle.start_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];
		}
		else
		{
			float RTV[2 * TotalRequests];
			float TwV[2 * TotalRequests];
			float tc[2 * TotalRequests];
			int Req_p[2 * TotalRequests];

			float A[2 * TotalRequests],
				B[2 * TotalRequests],
				C[2 * TotalRequests],
				D[2 * TotalRequests],
				W[2 * TotalRequests],
				L[2 * TotalRequests];

			float V_RT = 0;
			float V_TW = 0;
			float V_tc = 0;


			//vertex VertexNow;
			float EarliestTime;
			float ServiceTime;
			float LatestTime;
			int load;

			// 1 Find out who is in the vehicle now (NOT needed for static case)

			int p = 0;
			// 2 Departure from depot / get depot&start time

			int depart_location = 2 * n + vehid;
			float depart_time = Problem[0].Vehicle[vehid].start_time;
			Vehicle.start_time = depart_time;

			Vehicle.Cost.reset();

			bool coeff;

			for (i = 0; i<size; i++)
			{
				//cout << "part 3.1\n";
				j = routesequence[i];
				//VertexNow = find_vertex(Problem[0], j);

				coeff = j < n;

				/*EarliestTime = (coeff)*  Problem[0].Request[j].pickup.EarliestTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.EarliestTime;
				ServiceTime = (coeff)*  Problem[0].Request[j].pickup.ServiceTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.ServiceTime;
				LatestTime = (coeff)*  Problem[0].Request[j].pickup.LatestTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.LatestTime;
				load = (coeff)*  Problem[0].Request[j].load; + (1 - coeff) * -Problem[0].Request[j - n].load;*/



				if (j < n)
				{
					EarliestTime = Problem[0].Request[j].pickup.EarliestTime;
					ServiceTime = Problem[0].Request[j].pickup.ServiceTime;
					LatestTime = Problem[0].Request[j].pickup.LatestTime;
					load = Problem[0].Request[j].load;
				}
				else
				{
					EarliestTime = Problem[0].Request[j - n].dropoff.EarliestTime;
					ServiceTime = Problem[0].Request[j - n].dropoff.ServiceTime;
					LatestTime = Problem[0].Request[j - n].dropoff.LatestTime;
					load = -Problem[0].Request[j - n].load;
				}


				tc[j] = Problem[0].TravelTime[depart_location][j];
				A[j] = depart_time + Problem[0].TravelTime[depart_location][j];
				B[j] = max(A[j], EarliestTime);
				depart_time = B[j] + ServiceTime;
				depart_location = j;
				C[j] = depart_time;
				W[j] = B[j] - A[j];
				p = p + load;
				Req_p[j] = p;

				//cout << "part 3.2\n";

				/*coeff = j >= n;

				L[j] = (coeff) * (B[j] - C[j - n]);
				RTV[j] = (coeff)*  max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
				ppl_in_vehicle[j] = (coeff) * 0 + (1 - coeff) * 1;*/


				if (j >= n) //dropoff
				{
					L[j] = B[j] - C[j - n];
					RTV[j] = max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
				}
				else //pickup
				{
					ppl_in_vehicle[j] = 1;
					L[j] = 0;
					RTV[j] = 0;
				}

				TwV[j] = max((B[j] - LatestTime), 0.0);

				if (p > Problem[0].Vehicle[vehid].capacity)
					Vehicle.Cost.load = max(Vehicle.Cost.load, p - Problem[0].Vehicle[vehid].capacity);

				/*RTV[j] = Request[j].ride_time_violation;
				TwV[j] = Request[j].time_window_violation;
				tc[j] = Request[j].travel_cost;*/

				/*Vehicle.Cost.ride_time += Request[j].ride_time_violation;
				Vehicle.Cost.time_window += Request[j].time_window_violation;
				Vehicle.Cost.travel_cost += Request[j].travel_cost;*/
			}

			/*for (i = 0; i < size; i++)
			{
			j = routesequence[i];
			Request[j].A = A[j];
			Request[j].B = B[j];
			Request[j].C = C[j];
			Request[j].W = W[j];
			Request[j].L = L[j];
			}*/



			for (i = 0; i < size; i++)
			{
				j = routesequence[i];
				V_RT += RTV[j];
				V_TW += TwV[j];
				V_tc += tc[j];
			}

			Vehicle.Cost.ride_time = V_RT;
			Vehicle.Cost.time_window = V_TW;
			Vehicle.Cost.travel_cost = V_tc;


			durationCompress(Problem, routesequence, size, A, B, C, W, L, V_RT, RTV); // Parallel DurationCompression for all possible route_seq on Veh 'k'

			if (W[routesequence[0]] > 0) //PUT this condition inside expression itself.
			{
				Vehicle.start_time += W[routesequence[0]];
				W[routesequence[0]] = 0.0;
				A[routesequence[0]] = B[routesequence[0]];
			}

			Vehicle.end_time = depart_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];

			// Modified objective function
			/*
			float excess_ride_time = 0, passenger_waiting_time = 0, route_duration = 0, early_arrival_time = 0;

			for (i = 0; i < size; i++)
			{
			int node = routesequence[i];
			if (node < n)
			{
			passenger_waiting_time += W[node] * (Req_p[node] - Problem[0].Request[node].load);
			early_arrival_time += (Problem[0].Request[node].pickup.EarliestTime - A[node]);
			}
			else
			{
			excess_ride_time += B[node] - D[node - n] - Problem[0].TravelTime[node][node - n];

			passenger_waiting_time += W[node] * (Req_p[node] - (-Problem[0].Request[node - n].load));
			early_arrival_time += (Problem[0].Request[node].dropoff.EarliestTime - A[node]);
			}

			}
			route_duration = Vehicle.end_time - Vehicle.start_time;
			*/

			//recovery
			for (int i = 0; i < routesize; i++)
				routesequence[i] += 1;

			// calculate excessrideTime
			Vehicle.Cost.excessrideTime = 0;
			for (int j = 0; j < routesize; j++)
			{
				int node = route[j] - 1;
				if (node >= n)
					Vehicle.Cost.excessrideTime += (L[node] - Problem[0].TravelTime[node - n][node]);
			}


			////////////////////////////////////
			/// Modified objective function ////
			Vehicle.Cost.excessrideTime = 0;
			Vehicle.Cost.early_arrival_time = 0;
			Vehicle.Cost.on_board_user_wait_time = 0;
			Vehicle.Cost.route_duration = Vehicle.end_time - Vehicle.start_time;
			for (int index = 0; index < size; index++)
			{
				int node = routesequence[index];

				Vehicle.Cost.on_board_user_wait_time += (W[node] * (Req_p[node] - 1));

				if (node < n)
				{
					if (Problem[0].Request[node].pickup.EarliestTime > A[node])
						Vehicle.Cost.early_arrival_time += (Problem[0].Request[node].pickup.EarliestTime - A[node]);

				}
				else
				{
					Vehicle.Cost.excessrideTime += (L[node] - Problem[0].TravelTime[node][node - n]);

					if (Problem[0].Request[node - n].dropoff.EarliestTime > A[node])
						Vehicle.Cost.early_arrival_time += (Problem[0].Request[node - n].dropoff.EarliestTime - A[node]);
				}
			}
			//////////////////////////////////////

		}

		Vehicle.Cost.getFeasibility();

		Vehicle.Cost.overall_ideal_cost = (w1 * Vehicle.Cost.travel_cost + w2 * Vehicle.Cost.excessrideTime)
				+ 1* Vehicle.Cost.load
				+ 1* Vehicle.Cost.time_window
				+ 1* Vehicle.Cost.ride_time
				+ 1* Vehicle.Cost.duration;


	}

	__host__ __device__
		void durationCompress(problem *Problem, int *routesequence, int size, float *A, float *B, float *C, float *W, float *L, float &V_RT, float *RTV)
	{
		if (size>0)
		{

			float all_delay[ExpectedPath];



			float W_prev, E, delay;
			int i, j, h, request_num;
			float possible_delay, ride_time_violation;

			for (i = 0; i <size; i++)
				all_delay[i] = 0;

			{
				for (i = size - 2; i >= 0; i--)
				{
					W_prev = W[routesequence[i + 1]];
					E = max(find_vertex(Problem[0], routesequence[i]).LatestTime - B[routesequence[i]], 0.0);
					all_delay[i] = min(W_prev + all_delay[i + 1], E);
				}

				bool change = 1;
				int entrycount = 0;

				while (change)
				{


					entrycount++;
					if (entrycount > 1)
						break;

					//cout << "boom\n";
					change = 0;

					for (i = 0; i<size - 1; i++)
					{
						if (routesequence[i] >= n)
						{
							// i = dropoff, j = pickup

							request_num = routesequence[i] - n;
							j = -1;

							for (h = i - 1; h >= 0; h--)
								if (routesequence[h] == request_num)
								{
									j = h;
									break;
								}

							if (j >= 0)
							{
								all_delay[i] = min(all_delay[i], all_delay[j] + Problem[0].Request[routesequence[i] - n].RideTimeConstraint - L[routesequence[i]]);
								all_delay[i] = max(all_delay[i], 0.0);
							}

						}
					}

					for (i = size - 2; i >= 0; i--)
					{
						W_prev = W[routesequence[i + 1]];
						E = max(find_vertex(Problem[0], routesequence[i]).LatestTime - B[routesequence[i]], 0.0);
						delay = min(W_prev + all_delay[i + 1], E);

						if (delay < all_delay[i])
						{
							all_delay[i] = delay;
							change = 1;
						}
					}
				}
			}

			for (h = 0; h<size - 1; h++)
			{

				if (routesequence[h]<n)
					possible_delay = all_delay[h];
				else if (A[routesequence[h]] > B[routesequence[h]])
					possible_delay = A[routesequence[h]] - B[routesequence[h]];
				else
					continue;

				W[routesequence[h]] += possible_delay;

				if (W[routesequence[h]] < 0)
				{
					possible_delay -= W[routesequence[h]];
					W[routesequence[h]] = 0;
				}

				B[routesequence[h]] += possible_delay;
				C[routesequence[h]] += possible_delay;
				A[routesequence[h + 1]] += possible_delay;
				W[routesequence[h + 1]] -= possible_delay;

				if (routesequence[h] >= n)
				{
					L[routesequence[h]] += possible_delay;
					ride_time_violation = max((L[routesequence[h]] - Problem[0].Request[routesequence[h] - n].RideTimeConstraint), 0.0);

					if (ride_time_violation < 1e-10)
						ride_time_violation = 0.0;

					Vehicle.Cost.ride_time += ride_time_violation;
					Vehicle.Cost.ride_time -= RTV[routesequence[h]];
					RTV[routesequence[h]] = ride_time_violation;
				}
				else
				{
					L[routesequence[h] + n] -= possible_delay;
					ride_time_violation = max((L[routesequence[h] + n] - Problem[0].Request[routesequence[h]].RideTimeConstraint), 0.0);

					if (ride_time_violation < 1e-10)
						ride_time_violation = 0.0;

					Vehicle.Cost.ride_time += ride_time_violation;
					Vehicle.Cost.ride_time -= RTV[routesequence[h] + n];
					RTV[routesequence[h] + n] = ride_time_violation;
				}
			}

			if (Vehicle.Cost.ride_time < 1e-10)
				Vehicle.Cost.ride_time = 0.0;
		}
	}

	__host__ __device__
		void CPU_evaluateByRoute(problem *Problem)
	{


		int routesequence[ExpectedPath];
		int size = routesize;

		//leverage
		for (int i = 0; i < size; i++)
			routesequence[i] = route[i] - 1;

		// k = vehicle number = {integer: 1<=k<=m}, where {m = number of vehicle}
		if (size<1)
		{
			Vehicle.Cost.reset();
			Vehicle.start_time = Problem[0].Vehicle[vehid].start_time;
			float depart_time = Vehicle.start_time;
			int depart_location = 2 * n + vehid;
			Vehicle.end_time = Vehicle.start_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];
		}
		else
		{
			float RTV[2 * TotalRequests];
			float TwV[2 * TotalRequests];
			float tc[2 * TotalRequests];
			int Req_p[2 * TotalRequests];

			float A[2 * TotalRequests],
				B[2 * TotalRequests],
				C[2 * TotalRequests],
				D[2 * TotalRequests],
				W[2 * TotalRequests],
				L[2 * TotalRequests];

			float V_RT = 0;
			float V_TW = 0;
			float V_tc = 0;


			//vertex VertexNow;
			float EarliestTime;
			float ServiceTime;
			float LatestTime;
			int load;

			// 1 Find out who is in the vehicle now (NOT needed for static case)

			int p = 0;
			// 2 Departure from depot / get depot&start time

			int depart_location = 2 * n + vehid;
			float depart_time = Problem[0].Vehicle[vehid].start_time;
			Vehicle.start_time = depart_time;

			Vehicle.Cost.reset();

			bool coeff;

			for (int i = 0; i<size; i++)
			{
				//cout << "part 3.1\n";
				int j = routesequence[i];
				//VertexNow = find_vertex(Problem[0], j);

				coeff = j < n;

				/*EarliestTime = (coeff)*  Problem[0].Request[j].pickup.EarliestTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.EarliestTime;
				ServiceTime = (coeff)*  Problem[0].Request[j].pickup.ServiceTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.ServiceTime;
				LatestTime = (coeff)*  Problem[0].Request[j].pickup.LatestTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.LatestTime;
				load = (coeff)*  Problem[0].Request[j].load; + (1 - coeff) * -Problem[0].Request[j - n].load;*/



				if (j < n)
				{
					EarliestTime = Problem[0].Request[j].pickup.EarliestTime;
					ServiceTime = Problem[0].Request[j].pickup.ServiceTime;
					LatestTime = Problem[0].Request[j].pickup.LatestTime;
					load = Problem[0].Request[j].load;
				}
				else
				{
					EarliestTime = Problem[0].Request[j - n].dropoff.EarliestTime;
					ServiceTime = Problem[0].Request[j - n].dropoff.ServiceTime;
					LatestTime = Problem[0].Request[j - n].dropoff.LatestTime;
					load = -Problem[0].Request[j - n].load;
				}


				tc[j] = Problem[0].TravelTime[depart_location][j];
				A[j] = depart_time + Problem[0].TravelTime[depart_location][j];
				B[j] = max(A[j], EarliestTime);
				depart_time = B[j] + ServiceTime;
				depart_location = j;
				C[j] = depart_time;
				W[j] = B[j] - A[j];
				p = p + load;
				Req_p[j] = p;

				//cout << "part 3.2\n";

				/*coeff = j >= n;

				L[j] = (coeff) * (B[j] - C[j - n]);
				RTV[j] = (coeff)*  max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
				ppl_in_vehicle[j] = (coeff) * 0 + (1 - coeff) * 1;*/


				if (j >= n) //dropoff
				{
					L[j] = B[j] - C[j - n];
					RTV[j] = max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
				}
				else //pickup
				{
					ppl_in_vehicle[j] = 1;
					L[j] = 0;
					RTV[j] = 0;
				}

				TwV[j] = max((B[j] - LatestTime), 0.0);

				if (p > Problem[0].Vehicle[vehid].capacity)
					Vehicle.Cost.load = max(Vehicle.Cost.load, p - Problem[0].Vehicle[vehid].capacity);

				/*RTV[j] = Request[j].ride_time_violation;
				TwV[j] = Request[j].time_window_violation;
				tc[j] = Request[j].travel_cost;*/

				/*Vehicle.Cost.ride_time += Request[j].ride_time_violation;
				Vehicle.Cost.time_window += Request[j].time_window_violation;
				Vehicle.Cost.travel_cost += Request[j].travel_cost;*/
			}

			/*for (i = 0; i < size; i++)
			{
			j = routesequence[i];
			Request[j].A = A[j];
			Request[j].B = B[j];
			Request[j].C = C[j];
			Request[j].W = W[j];
			Request[j].L = L[j];
			}*/



			for (int i = 0; i < size; i++)
			{
				int j = routesequence[i];
				V_RT += RTV[j];
				V_TW += TwV[j];
				V_tc += tc[j];
			}


			Vehicle.Cost.ride_time = V_RT;
			Vehicle.Cost.time_window = V_TW;
			Vehicle.Cost.travel_cost = V_tc;


			CPU_durationCompress(/*vehid - 1, route[0], routesize,*/ Problem, routesequence, size, A, B, C, W, L, V_RT, RTV); // Parallel DurationCompression for all possible route_seq on Veh 'k'

			if (W[routesequence[0]] > 0) //PUT this condition inside expression itself.
			{
				Vehicle.start_time += W[routesequence[0]];
				W[routesequence[0]] = 0.0;
				A[routesequence[0]] = B[routesequence[0]];
			}

			Vehicle.end_time = depart_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];

			//recovery
			for (int j = 0; j < routesize; j++)
				routesequence[j] += 1;

			// calculate excessrideTime
			Vehicle.Cost.excessrideTime = 0;
			for (int j = 0; j < routesize; j++)
			{
				int node = route[j] - 1;
				if (node >= n)
					Vehicle.Cost.excessrideTime += (L[node] - Problem[0].TravelTime[node - n][node]);
			}
		}


		Vehicle.Cost.getFeasibility();
	}

	__host__ __device__
		void CPU_durationCompress(/*int k, int& route_sequence, int routesize,*/ problem *Problem, int *routesequence, int size, float *A, float *B, float *C, float *W, float *L, float &V_RT, float *RTV)
	{
		if (size>0)
		{

			float all_delay[ExpectedPath];

			for (int i = 0; i < size; i++)
				all_delay[i] = 0;

			{
				for (int i = size - 2; i >= 0; i--)
				{
					float W_prev = W[routesequence[i + 1]];
					float E = max(find_vertex(Problem[0], routesequence[i]).LatestTime - B[routesequence[i]], 0.0);
					all_delay[i] = min(W_prev + all_delay[i + 1], E);
				}

				bool change = 1;
				int entrycount = 0;

				while (change)
				{


					entrycount++;
					if (entrycount > 1)
						break;

					//cout << "boom\n";
					change = 0;

					for (int i = 0; i<size - 1; i++)
					{
						if (routesequence[i] >= n)
						{
							// i = dropoff, j = pickup

							int request_num = routesequence[i] - n;
							int j = -1;

							for (int h = i - 1; h >= 0; h--)
								if (routesequence[h] == request_num)
								{
									j = h;
									break;
								}

							if (j >= 0)
							{
								all_delay[i] = min(all_delay[i], all_delay[j] + Problem[0].Request[routesequence[i] - n].RideTimeConstraint - L[routesequence[i]]);
								all_delay[i] = max(all_delay[i], 0.0);
							}

						}
					}

					for (int i = size - 2; i >= 0; i--)
					{
						float W_prev = W[routesequence[i + 1]];
						float E = max(find_vertex(Problem[0], routesequence[i]).LatestTime - B[routesequence[i]], 0.0);
						float delay = min(W_prev + all_delay[i + 1], E);

						if (delay < all_delay[i])
						{
							all_delay[i] = delay;
							change = 1;
						}
					}
				}
			}

			for (int h = 0; h<size - 1; h++)
			{
				float possible_delay;

				if (routesequence[h]<n)
					possible_delay = all_delay[h];
				else if (A[routesequence[h]] > B[routesequence[h]])
					possible_delay = A[routesequence[h]] - B[routesequence[h]];
				else
					continue;

				W[routesequence[h]] += possible_delay;

				if (W[routesequence[h]] < 0)
				{
					possible_delay -= W[routesequence[h]];
					W[routesequence[h]] = 0;
				}

				B[routesequence[h]] += possible_delay;
				C[routesequence[h]] += possible_delay;
				A[routesequence[h + 1]] += possible_delay;
				W[routesequence[h + 1]] -= possible_delay;

				if (routesequence[h] >= n)
				{
					L[routesequence[h]] += possible_delay;
					float ride_time_violation = max((L[routesequence[h]] - Problem[0].Request[routesequence[h] - n].RideTimeConstraint), 0.0);

					if (ride_time_violation < 1e-10)
						ride_time_violation = 0.0;

					Vehicle.Cost.ride_time += ride_time_violation;
					Vehicle.Cost.ride_time -= RTV[routesequence[h]];
					RTV[routesequence[h]] = ride_time_violation;
				}
				else
				{
					L[routesequence[h] + n] -= possible_delay;
					float ride_time_violation = max((L[routesequence[h] + n] - Problem[0].Request[routesequence[h]].RideTimeConstraint), 0.0);

					if (ride_time_violation < 1e-10)
						ride_time_violation = 0.0;

					Vehicle.Cost.ride_time += ride_time_violation;
					Vehicle.Cost.ride_time -= RTV[routesequence[h] + n];
					RTV[routesequence[h] + n] = ride_time_violation;
				}
			}

			if (Vehicle.Cost.ride_time < 1e-10)
				Vehicle.Cost.ride_time = 0.0;
		}
	}

	__host__ __device__
		void print()
	{
		printf("Printing Demo \n");
		printf("Demo[0].isFeasible %s\n", Vehicle.Cost.isFeasible == 0 ? "false" : "true");
		printf("Demo[0].travelcost %f (%s)\n", Vehicle.Cost.travel_cost, Vehicle.Cost.travel_cost == 0 ? "true" : "false");
		printf("Demo[0].timeWindow %f (%s) \n", Vehicle.Cost.time_window, Vehicle.Cost.time_window == 0 ? "true" : "false");
		printf("Demo[0].rideTimeViol %f (%s) \n", Vehicle.Cost.ride_time, Vehicle.Cost.ride_time == 0 ? "true" : "false");
		printf("Demo[0].loadViol %d (%s) \n", Vehicle.Cost.load, Vehicle.Cost.load == 0 ? "true" : "false");
		printf("Demo[0].durationViol %f (%s) \n", Vehicle.Cost.duration, Vehicle.Cost.duration == 0 ? "true" : "false");
		printf("Demo[0].currentcost %0.02f \n", current_cost);
	}

	__host__ __device__
		void PrintInformation()
	{
		printf("----\n");
		printf("Here's the Printed Information:\n");
		printf("Req %d Veh %d Start %d Gap %d solID %d currentCost %f\n", request, vehid, start, gap, solID, current_cost);
		printf("isFeasible %s\n", Vehicle.Cost.isFeasible == 0 ? "false" : "true");
		printf("isWholeRouteFeasible? %s\n", isFeasible == 0 ? "false" : "true");
		printf("travelcost %f (%s)\n", Vehicle.Cost.travel_cost, Vehicle.Cost.travel_cost == 0 ? "true" : "false");
		printf("timeWindow %f (%s) \n", Vehicle.Cost.time_window, Vehicle.Cost.time_window == 0 ? "true" : "false");
		printf("rideTimeViol %f (%s) \n", Vehicle.Cost.ride_time, Vehicle.Cost.ride_time == 0 ? "true" : "false");
		printf("loadViol %d (%s) \n", Vehicle.Cost.load, Vehicle.Cost.load == 0 ? "true" : "false");
		printf("durationViol %f (%s) \n", Vehicle.Cost.duration, Vehicle.Cost.duration == 0 ? "true" : "false");
		printf("----\n");
	}


	__host__ __device__
		void update_A_B_C_W(problem *Problem, int current_index, float depart_time, int depart_location, int arrival_location, float *travel_cost, float *A, float *B, float *C, float *W, int *routesequence, int size)
	{
		// passed k value starts from zero
		// update from this current_index (this index excluded, dont update it.)

		int index, j;

		float current_time = depart_time;
		int current_location = depart_location;

		float EarliestTime, ServiceTime;

		bool coeff;

		for (index = current_index + 1; index < size; index++)
		{
			j = routesequence[index];

			coeff = j < n;

			EarliestTime = coeff * Problem[0].Request[j].pickup.EarliestTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.EarliestTime;
			ServiceTime = coeff * Problem[0].Request[j].pickup.ServiceTime + (1 - coeff) * Problem[0].Request[j - n].dropoff.ServiceTime;

			/*if (j < n)
			{
			EarliestTime = Problem[0].Request[j].pickup.EarliestTime;
			ServiceTime = Problem[0].Request[j].pickup.ServiceTime;
			}
			else
			{
			EarliestTime = Problem[0].Request[j - n].dropoff.EarliestTime;
			ServiceTime = Problem[0].Request[j - n].dropoff.ServiceTime;
			}*/

			travel_cost[j] = Problem[0].TravelTime[current_location][j];
			A[j] = current_time + Problem[0].TravelTime[current_location][j];

			coeff = A[j] > EarliestTime;
			//B[j] = coeff * A[j] + (1 - coeff) * EarliestTime;
			B[j] = max(A[j], EarliestTime);
			C[j] = B[j] + ServiceTime;
			W[j] = B[j] - A[j];



			current_time = C[j];
			current_location = j;
		}

		Vehicle.end_time = current_time + Problem[0].TravelTime[current_location][arrival_location];
	}

	__host__ __device__
		void update_L_p(problem *Problem, float *B, float *C, float *L, int *p, int *routesequence, int size)
	{
		int current_load = 0;

		int j;
		bool coeff;

		for (int i = 0; i < size; i++)
		{
			j = routesequence[i];

			/*coeff = j < n;

			current_load = coeff * (current_load + 1) + (1 - coeff) * (current_load - 1);
			L[j] = coeff * 0 + (1 - coeff) * (B[j] - C[j - n]);*/

			if (j < n)
			{
				current_load += 1;// Problem[0].Request[j].load;
				L[j] = 0;

			}
			else
			{
				current_load -= 1;// Problem[0].Request[j - n].load;
				L[j] = B[j] - C[j - n];
			}

			p[j] = current_load;

		}
	}

	__host__ __device__
		void compute_slack_for_vertex(problem *Problem, int start_vertex_index, float *min_slack, float *cumulative_waiting_time, float *B, float *L, float *W, int *routesequence, int size)
	{
		float factor_1[3 * n]; // summation Wp
		float factor_2[3 * n]; // l{j} - B{j}
		float factor_3[3 * n]; // RTC - P{j}; if j=>dropoff, else 0.

		float slack[3 * n];

		int node;

		// 1) computing factor_1 i.e., running Wp

		cumulative_waiting_time[0] = 0;
		for (int i = start_vertex_index; i < size; i++)
		{
			node = routesequence[i];
			if (i == 0)
				cumulative_waiting_time[0] += W[node];
			else if (i == start_vertex_index)
				cumulative_waiting_time[0] = 0;
			else
				cumulative_waiting_time[0] += W[node];

			factor_1[i] = cumulative_waiting_time[0];
		}

		// 2) compute factor_2 and factor_3

		for (int i = start_vertex_index; i < size; i++)
		{
			node = routesequence[i];
			if (node < n) //pickup
			{
				factor_2[i] = max(0, Problem[0].Request[node].pickup.LatestTime - B[node]);
				factor_3[i] = Problem[0].Request[node].RideTimeConstraint;
			}
			else //dropoff
			{
				factor_2[i] = max(0, Problem[0].Request[node - n].dropoff.LatestTime - B[node]);
				factor_3[i] = max(0, Problem[0].Request[node - n].RideTimeConstraint - L[node]);
			}
		}

		// 3) compute minimum_slack_time

		for (int i = start_vertex_index; i < size; i++)
			slack[i] = factor_1[i] + min(factor_2[i], factor_3[i]);

		int min_i = start_vertex_index;
		float minimum_slack = slack[start_vertex_index];
		for (int i = start_vertex_index + 1; i < size; i++)
		{
			if (slack[i] < minimum_slack)
			{
				// if (Vehicle.path[i] - 1 < n) // only pickups are considered!!!
				{
					minimum_slack = slack[i];
					min_i = i;
				}
			}
		}

		min_slack[0] = minimum_slack;

	}

	__host__ __device__
		void compute_violations(problem *Problem, int arrival_location, float *travel_cost, float *B, float *L, int *p, int *routesequence, int size)
	{
		float ride_time_violation, time_window_violation;

		int index, j;

		Vehicle.Cost.reset();

		float total_ride_time_violation = 0;
		float total_time_window_violation = 0;
		float total_travel_cost = 0;

		bool coeff;

		for (index = 0; index < size; index++)
		{
			j = routesequence[index];

			coeff = j < n;

			total_ride_time_violation += (1 - coeff) *  max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
			total_time_window_violation += coeff * max((B[j] - Problem[0].Request[j].pickup.LatestTime), 0.0) + (1 - coeff) * max((B[j] - Problem[0].Request[j - n].dropoff.LatestTime), 0.0);

			//if (j < n) //pickup
			//{
			//	ride_time_violation = 0;
			//	time_window_violation = max((B[j] - Problem[0].Request[j].pickup.LatestTime), 0.0);
			//}
			//else //dropoff
			//{
			//	ride_time_violation = max((L[j] - Problem[0].Request[j - n].RideTimeConstraint), 0.0);
			//	time_window_violation = max((B[j] - Problem[0].Request[j - n].dropoff.LatestTime), 0.0);
			//}

			if (p[j] > Problem[0].Vehicle[vehid].capacity)
				Vehicle.Cost.load = max(Vehicle.Cost.load, p[j] - Problem[0].Vehicle[vehid].capacity);

			//total_ride_time_violation += ride_time_violation;
			//total_time_window_violation += time_window_violation;
			total_travel_cost += travel_cost[j];
		}

		Vehicle.Cost.ride_time = total_ride_time_violation;
		Vehicle.Cost.time_window = total_time_window_violation;
		Vehicle.Cost.travel_cost = total_travel_cost;

		int last_visited_node = routesequence[size - 1];
		Vehicle.Cost.travel_cost += Problem[0].TravelTime[last_visited_node][arrival_location];

		Vehicle.Cost.duration = max((Vehicle.end_time - Vehicle.start_time) - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
		Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);


	}

	__host__ __device__
		void Offline_route_evaluation(problem *Problem)
	{
		int i, j;

		int routesequence[ExpectedPath];
		int size = routesize;

		//leverage
		for (i = 0; i < size; i++)
			routesequence[i] = route[i] - 1;

		if (size < 1)
		{
			Vehicle.Cost.reset();
			Vehicle.start_time = Problem[0].Vehicle[vehid].start_time;
			float depart_time = Vehicle.start_time;
			int depart_location = 2 * n + vehid;
			Vehicle.end_time = Vehicle.start_time + Problem[0].TravelTime[depart_location][2 * n + m + vehid];
			Vehicle.Cost.duration = max(Vehicle.end_time - Vehicle.start_time - Problem[0].Vehicle[vehid].duration_constraint, 0.0);
			Vehicle.Cost.duration += max(Vehicle.end_time - Problem[0].Vehicle[vehid].end_time, 0.0);
			Vehicle.Cost.travel_cost += Problem[0].TravelTime[depart_location][2 * n + m + vehid];
		}
		else
		{

			//  1. Set D0 : ¼ e0.

			int depart_location = 2 * n + vehid;;
			int arrival_location = 2 * n + m + vehid;

			float depart_time = Problem[0].Vehicle[vehid].start_time;

			Vehicle.Cost.reset();
			Vehicle.start_time = depart_time;

			//  2. Compute Ai, Wi, Bi and Di and for each vertex vi in the route.

			float A[2 * n],
				B[2 * n],
				C[2 * n],
				W[2 * n],
				L[2 * n],
				travel_cost[2 * n];
			int	p[2 * n];

			/*for (int i = 0; i < size; i++)
			{
			int j = routesequence[i];

			A[j] = 0;
			B[j] = 0;
			C[j] = 0;
			L[j] = 0;
			W[j] = 0;
			p[j] = 0;
			}*/

			update_A_B_C_W(Problem, -1, depart_time, depart_location, arrival_location, travel_cost, A, B, C, W, routesequence, size);
			update_L_p(Problem, B, C, L, p, routesequence, size);

			//	3. Compute F0.

			float cumulative_waiting_time[1];
			float F_slack[1];
			compute_slack_for_vertex(Problem, 0, F_slack, cumulative_waiting_time, B, L, W, routesequence, size);

			//	4. Set D0 : ¼ e0 þ minfF0;0<p<q Wpg.

			depart_time = Problem[0].Vehicle[vehid].start_time + min(F_slack[0], cumulative_waiting_time[0]);
			Vehicle.start_time = depart_time;

			//	5, 6. Update Ai, Wi, Bi, Di, Li, p for each vertex vi in the route.

			update_A_B_C_W(Problem, -1, depart_time, depart_location, arrival_location, travel_cost, A, B, C, W, routesequence, size);
			update_L_p(Problem, B, C, L, p, routesequence, size);

			compute_violations(Problem, arrival_location, travel_cost, B, L, p, routesequence, size);

			//	7. For everyvertex vj that corresponds to the origin of a request j

			for (i = 0; i < size; i++)
			{

				j = routesequence[i];

				if (j < n) //pickup
				{
					float Serv_time = Problem[0].Request[j].pickup.ServiceTime;

					//	(a) Compute Fj.

					compute_slack_for_vertex(Problem, i, F_slack, cumulative_waiting_time, B, L, W, routesequence, size);

					//	(b) Set Bj : ¼ Bj þ minfFj; j<p<q Wpg; Dj:¼ Bj þ dj.

					B[j] = B[j] + min(F_slack[0], cumulative_waiting_time[0]);
					C[j] = B[j] + Serv_time;
					W[j] = B[j] - A[j];

					//	(c) Update Ai, Wi, Bi and Di, for each vertex vi that comes after vj in the route.

					update_A_B_C_W(Problem, i, C[j], j, arrival_location, travel_cost, A, B, C, W, routesequence, size);

					//	(d) Update the ride time Li for each request i whose destination vertex is after vertex vj.

					update_L_p(Problem, B, C, L, p, routesequence, size);

					compute_violations(Problem, arrival_location, travel_cost, B, L, p, routesequence, size);
				}


			}


			// delay departure at the depot

			j = routesequence[0];

			Vehicle.start_time += W[j];
			A[j] = B[j];
			W[j] = B[j] - A[j];


			//	8. Compute changes in violations of vehicle load, route duration, time window and ride time constraints.

			compute_violations(Problem, arrival_location, travel_cost, B, L, p, routesequence, size);


			// Modified objective function
			Vehicle.Cost.excessrideTime = 0;
			Vehicle.Cost.early_arrival_time = 0;
			Vehicle.Cost.on_board_user_wait_time = 0;
			for (i = 0; i < size; i++)
			{
				j = routesequence[i];
				Vehicle.Cost.on_board_user_wait_time += (W[j] * (p[j]-1));
				if (j < n)
				{
					if (Problem[0].Request[j].pickup.EarliestTime > A[j])
						Vehicle.Cost.early_arrival_time += (Problem[0].Request[j].pickup.EarliestTime - A[j]);
				}
				else
				{
					Vehicle.Cost.excessrideTime += (L[j] - Problem[0].TravelTime[j - n][j]);

					if (Problem[0].Request[j - n].dropoff.EarliestTime > A[j])
						Vehicle.Cost.early_arrival_time += (Problem[0].Request[j - n].dropoff.EarliestTime - A[j]);
				}
			}
			Vehicle.Cost.route_duration = Vehicle.end_time - Vehicle.start_time;


		}

	}


};
