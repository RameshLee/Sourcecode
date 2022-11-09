
struct position
{
	int start;
	int gap;

	int request;
	int vehicle;
};

struct vehicle_detail
{
	//evaluation
	float start_time;
	float end_time;
	cost Cost;

	//insertion
	int path[ExpectedPath];
	int id;
	int size;
	int insertways;
	int cumulativeInsertions;

	int StartingPoint;
	int EndingPoint;

	__host__ __device__ vehicle_detail()
	{
		StartingPoint = 0;
		EndingPoint = 0;
	}

	__host__ __device__
		void initialize(int size)
	{

		if (size + 2 > ExpectedPath)
		{
			printf("ExpectedPath Reached. Overflow!!!");
			exit(0);
		}
		this->size = size;
		//path = new int[size];
		//CHECK(cudaMallocManaged((void**)&path, size * sizeof(int)));
		insertways = (size + 1) * (size + 2) / 2;
		cumulativeInsertions = 0;
	}

	__host__ __device__
		void Copy(vehicle_detail &inVehicle)
	{
			StartingPoint = inVehicle.StartingPoint;
			EndingPoint = inVehicle.EndingPoint;

			id = inVehicle.id;
			size = inVehicle.size;

			for (int j = 0; j < size; j++)
				path[j] = inVehicle.path[j];

			insertways = inVehicle.insertways;
			cumulativeInsertions = inVehicle.cumulativeInsertions;
	}

	__host__ __device__
		~vehicle_detail() {}
};

struct request_detail
{
	position Position;

	//evaluation
	float A, B, C, D, W, L;
	int p;
	float time_window_violation; //twV
	float ride_time_violation; //rtV
	float travel_cost;

	//insertion
	int id;
	int currvehicle; //starts from zero
	int totinsertion; //inter
	int intra_totinsertion; //intra
	int allinsertion; //intra + inter

	__host__ __device__ request_detail() {}

	__host__ __device__
		void clean()
	{
		A, B, C, D, W, L = 0;
		p = 0;
		time_window_violation = 0;
		ride_time_violation = 0;
		travel_cost = 0;
	}

	__host__ __device__
		void initialize()
	{
		totinsertion = 0;
	}

	__host__ __device__
		void Copy(request_detail &inRequest)
	{
		id = inRequest.id;
		currvehicle = inRequest.currvehicle;
		totinsertion = inRequest.totinsertion;
		intra_totinsertion = inRequest.intra_totinsertion;
		allinsertion = inRequest.allinsertion;
	}

	__host__ __device__ ~request_detail() {}
};

struct location
{
	int ID;
	float x, y;

	__host__ __device__
		void Copy(location &inLocation)
	{
		ID = inLocation.ID;
		x = inLocation.x;
		y = inLocation.y;
	}
};

struct vertex
{
	location Location;
	float ServiceTime,
		EarliestTime,
		LatestTime;

	__host__ __device__
		void Copy(vertex &inVertex)
	{
		Location.Copy(inVertex.Location);
		ServiceTime = inVertex.ServiceTime;
		EarliestTime = inVertex.EarliestTime;
		LatestTime = inVertex.LatestTime;

	}

	__host__ __device__
		void print_vertex()
	{
		printf("%d	%0.03f	%0.03f	%0.0f	%d	%0.0f	%0.0f\n", Location.ID, Location.x, Location.y, ServiceTime, 0, EarliestTime, LatestTime);
	}
};

struct request
{
	int ID;
	bool isDropoffCritical;
	float RideTimeConstraint;
	int load;

	vertex pickup;
	vertex dropoff;

	bool isPickUpDone;
	int vehicle;
	float departureTime;

	__host__ __device__
		request() : isPickUpDone(false) {}

	__host__ __device__
		void Copy(request &inRequest)
	{
		ID = inRequest.ID;
		isDropoffCritical = inRequest.isDropoffCritical;
		RideTimeConstraint = inRequest.RideTimeConstraint;
		load = inRequest.load;

		pickup.Copy(inRequest.pickup);
		dropoff.Copy(inRequest.dropoff);
	}

	__host__ __device__
		void initalize()
	{
		isPickUpDone = false;
	}

	__host__ __device__
		~request() {}
};

struct vehicle
{
	int ID;
	int capacity;

	location start;
	location end;

	float start_time;
	float end_time;
	float duration_constraint;

	__host__ __device__
		vehicle()
	{
		reset();
	}

	__host__ __device__
		void Copy(vehicle &inVehicle)
	{
		ID = inVehicle.ID;
		capacity = inVehicle.capacity;

		start.Copy(inVehicle.start);
		end.Copy(inVehicle.end);

		start_time = inVehicle.start_time;
		end_time = inVehicle.end_time;
		duration_constraint = inVehicle.duration_constraint;
	}


	__host__ __device__
		void initialize()
	{
		reset();
	}

	__host__ __device__
		void reset()
	{
		//ID = 0;
		capacity = 0;
		start.ID = 0;
		start.x = 0;
		start.y = 0;
		end = start;
		start_time = 0;
		end_time = 1440;
		duration_constraint = 0;
	}

	__host__ __device__
		~vehicle() {}

};

struct request_info
{
	int load;
	vertex pickup;
	vertex dropoff;
	float RTC;
};

struct vehicle_info
{
	int capacity;
	float start_time;
	float end_time;
	float duration_constraint;
};


struct problem
{
	//EV-aspect

	request Request[TotalRequests];
	vehicle Vehicle[TotalVehicles];
	float TravelTime[Var][Var];r
	//float TravelCost[Var][Var];

	__host__ __device__
		problem()
	{
		reset();
	}

	__host__ __device__
		void Copy(problem &inProblem)
	{

		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].Copy(inProblem.Vehicle[k]);

		for (int i = 0; i < TotalRequests; i++)
			Request[i].Copy(inProblem.Request[i]);

		for (int i = 0; i < Var; i++)
			for (int j = 0; j < Var; j++)
			{
				TravelTime[i][j] = inProblem.TravelTime[i][j];
				//TravelCost[i][j] = inProblem.TravelCost[i][j];
			}
	}

		void Copy_Except_SampledRequests(problem &inProblem, int s, myList *SampledRequests)
	{

		for (int k = 0; k < TotalVehicles; k++)
			Vehicle[k].Copy(inProblem.Vehicle[k]);

		for (int i = 0; i < TotalRequests; i++)
		{
			if (SampledRequests[s].isFound(i) == false)
				Request[i].Copy(inProblem.Request[i]);
		}

		for (int i = 0; i < Var; i++)
			for (int j = 0; j < Var; j++)
			{
				TravelTime[i][j] = inProblem.TravelTime[i][j];
				//TravelCost[i][j] = inProblem.TravelCost[i][j];
			}
	}

	__host__ __device__
		void initialize()
	{
		reset();
		for (int i = 0; i < TotalRequests; i++)
		{
			Request[i].initalize();
			Vehicle[i].initialize();
		}

	}

	__host__ __device__
		void reset()
	{
		for (int i = 0; i < Var; i++)
			for (int j = 0; j < Var; j++)
			{
				TravelTime[i][j] = 0;
				TravelTime[i][j] = 0;
			}
	}

	__host__ __device__
		~problem() {}

		void print()
	{
		for (int i = 0; i < TotalRequests; i++)
		{
			Request[i].pickup.Location.ID = i+1;
			printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Request[i].pickup.Location.ID, Request[i].pickup.Location.x, Request[i].pickup.Location.y, Request[i].pickup.ServiceTime,
				Request[i].load, Request[i].pickup.EarliestTime, Request[i].pickup.LatestTime);
		}

		for (int i = 0; i < TotalRequests; i++)
		{
			Request[i].dropoff.Location.ID = n+i+1;
			printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Request[i].dropoff.Location.ID, Request[i].dropoff.Location.x, Request[i].dropoff.Location.y, Request[i].dropoff.ServiceTime,
				0, Request[i].dropoff.EarliestTime, Request[i].dropoff.LatestTime);
		}
	}
};
