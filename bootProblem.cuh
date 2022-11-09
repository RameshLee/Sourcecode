void SortRequestsBasedOnTW(problem &Original)
{
	vertex PickUp;
	vertex DropOff;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			if (i!=j)
				if (Original.Request[i].pickup.EarliestTime < Original.Request[j].pickup.EarliestTime)
				{
					PickUp.Copy(Original.Request[i].pickup);
					Original.Request[i].pickup.Copy(Original.Request[j].pickup);
					Original.Request[j].pickup.Copy(PickUp);

					DropOff.Copy(Original.Request[i].dropoff);
					Original.Request[i].dropoff.Copy(Original.Request[j].dropoff);
					Original.Request[j].dropoff.Copy(DropOff);
				}
	}


	for (int i = 0; i < n; i++)
	{
		Original.Request[i].dropoff.EarliestTime = 0;
		Original.Request[i].dropoff.LatestTime = 1440;
	}

}


//ExternalFunction
__host__ __device__
void update_relationship(problem& Problem)
{
	//Set Req ID
	for (int i = 0; i < TotalRequests; i++)
	{
		Problem.Request[i].ID = i + 1;
	}

	//Based on Eucledian Distance
	for (int i = 0; i<Var; i++)
		for (int j = 0; j<i; j++)
		{
			float x_diff = find_location(Problem, i).x - find_location(Problem, j).x;
			float y_diff = find_location(Problem, i).y - find_location(Problem, j).y;
			Problem.TravelTime[i][j] = sqrt(pow(x_diff, 2) + pow(y_diff, 2));
			Problem.TravelTime[j][i] = sqrt(pow(x_diff, 2) + pow(y_diff, 2));
		}



	//TravelCost = TravelTime
/*	for (int i = 0; i<Var; i++)
		for (int j = 0; j < Var; j++)
		{
			if (i != j)
				Problem.TravelCost[i][j] = Problem.TravelTime[i][j];
			else
				Problem.TravelCost[i][j] = 0;
		}*/

}

__host__ __device__
void time_window_adjustment(problem& Problem)
{
	for (int i = 0; i<TotalRequests; i++)
	{
		if ((Problem.Request[i].dropoff.LatestTime - Problem.Request[i].dropoff.EarliestTime)<(Problem.Request[i].pickup.LatestTime - Problem.Request[i].pickup.EarliestTime))
			Problem.Request[i].isDropoffCritical = true;
		else
			Problem.Request[i].isDropoffCritical = false;

		float distance = Problem.TravelTime[i][i + n];

		Problem.Request[i].pickup.EarliestTime = max(Problem.Request[i].pickup.EarliestTime, Problem.Request[i].dropoff.EarliestTime
			- distance - Problem.Request[i].RideTimeConstraint - Problem.Request[i].pickup.ServiceTime);
		Problem.Request[i].pickup.LatestTime = min(Problem.Request[i].pickup.LatestTime, Problem.Request[i].dropoff.LatestTime - distance - Problem.Request[i].pickup.ServiceTime);

		Problem.Request[i].dropoff.EarliestTime = max(Problem.Request[i].dropoff.EarliestTime, Problem.Request[i].pickup.EarliestTime + distance + Problem.Request[i].pickup.ServiceTime);
		Problem.Request[i].dropoff.LatestTime = min(Problem.Request[i].dropoff.LatestTime, Problem.Request[i].pickup.LatestTime + distance + Problem.Request[i].RideTimeConstraint + Problem.Request[i].pickup.ServiceTime);
	}
}

void BootProblem(problem &Original)
{
	string filename = DeclaredProblem;


	ifstream myFile(filename.c_str());

	if (!myFile)
	{
		cerr << "File could not be opened" << endl;
	}

	int* temp = new int[5];
	int* dummy = new int;
	int* mm = new int;
	int* nn = new int;
	location* Depot = new location;

	myFile >> *mm >> *nn >> temp[0] >> temp[1] >> temp[2];

	*nn = *nn / 2;

	myFile >> Depot[0].ID >> Depot[0].x >> Depot[0].y >> *dummy >> *dummy >> temp[3] >> temp[4];

	for (int i = 0; i < TotalRequests; i++)
	{
		myFile >> Original.Request[i].pickup.Location.ID >> Original.Request[i].pickup.Location.x >> Original.Request[i].pickup.Location.y >> Original.Request[i].pickup.ServiceTime
			>> Original.Request[i].load >> Original.Request[i].pickup.EarliestTime >> Original.Request[i].pickup.LatestTime;
	}

	for (int i = 0; i < TotalRequests; i++)
	{
		myFile >> Original.Request[i].dropoff.Location.ID >> Original.Request[i].dropoff.Location.x >> Original.Request[i].dropoff.Location.y >> Original.Request[i].dropoff.ServiceTime
			>> *dummy >> Original.Request[i].dropoff.EarliestTime >> Original.Request[i].dropoff.LatestTime;
	}

	for (int k = 0; k < TotalVehicles; k++)
	{
		Original.Vehicle[k].capacity = temp[1];
		Original.Vehicle[k].duration_constraint = temp[0];
		Original.Vehicle[k].start_time = temp[3];
		Original.Vehicle[k].end_time = temp[4];
		Original.Vehicle[k].start = Depot[0];
		Original.Vehicle[k].end = Depot[0];
	}

	for (int i = 0; i < TotalRequests; i++)
	{
		Original.Request[i].RideTimeConstraint = temp[2];
	}

	// PART II
	time_window_adjustment(Original);

	// print utility
/*	if (1)
	{
		printf("%d %d %d %d %d\n", *mm, 2*(*nn), temp[0], temp[1], temp[2]);

		printf("0	%0.02f	%0.02f	0	0	%d	%d\n", Depot[0].x, Depot[0].y, temp[3], temp[4]);

		for (int i = 0; i < TotalRequests; i++)
		{
			printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Original.Request[i].pickup.Location.ID, Original.Request[i].pickup.Location.x, Original.Request[i].pickup.Location.y, Original.Request[i].pickup.ServiceTime,
				Original.Request[i].load, Original.Request[i].pickup.EarliestTime, Original.Request[i].pickup.LatestTime);
		}

		for (int i = 0; i < TotalRequests; i++)
		{
			printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Original.Request[i].dropoff.Location.ID, Original.Request[i].dropoff.Location.x, Original.Request[i].dropoff.Location.y, Original.Request[i].dropoff.ServiceTime,
				0, Original.Request[i].dropoff.EarliestTime, Original.Request[i].dropoff.LatestTime);
		}
	}*/

	/////////////////////////////////////////////////////
	/////// OSCO: Sort requests and print problem ///////
	////////////////////////////////////////////////////
	SortRequestsBasedOnTW(Original);
	if (0)
	{
		printf("%d %d %d %d %d\n", *mm, 2*(*nn), temp[0], temp[1], temp[2]);

		printf("0	%0.02f	%0.02f	0	0	%d	%d\n", Depot[0].x, Depot[0].y, temp[3], temp[4]);

		for (int i = 0; i < TotalRequests; i++)
		{
			Original.Request[i].pickup.Location.ID = i+1;
			printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Original.Request[i].pickup.Location.ID, Original.Request[i].pickup.Location.x, Original.Request[i].pickup.Location.y, Original.Request[i].pickup.ServiceTime,
				Original.Request[i].load, Original.Request[i].pickup.EarliestTime, Original.Request[i].pickup.LatestTime);
		}

		for (int i = 0; i < TotalRequests; i++)
		{
			Original.Request[i].dropoff.Location.ID = n+i+1;
			printf("%d	%0.02f	%0.02f	%0.0f	%d	%0.0f	%0.0f\n", Original.Request[i].dropoff.Location.ID, Original.Request[i].dropoff.Location.x, Original.Request[i].dropoff.Location.y, Original.Request[i].dropoff.ServiceTime,
				0, Original.Request[i].dropoff.EarliestTime, Original.Request[i].dropoff.LatestTime);
		}
	}
	/////////////////////////////////////////////////////


	// part-III
	update_relationship(Original);

	delete[] temp;
	delete dummy;
	delete mm;
	delete nn;
	delete Depot;
	myFile.close();


}





/*void Update_ConstantMemory(problem *Problem)
{
	pblm h_Pblm[1];

	for (int i = 0; i < TotalRequests; i++)
	{
		h_Pblm[0].Request[i].load = Problem[0].Request[i].load;
		h_Pblm[0].Request[i].pickup.Copy(Problem[0].Request[i].pickup);
		h_Pblm[0].Request[i].dropoff.Copy(Problem[0].Request[i].dropoff);
		h_Pblm[0].Request[i].RTC = Problem[0].Request[i].RideTimeConstraint;
	}

	for (int k = 0; k < TotalVehicles; k++)
	{
		h_Pblm[0].Vehicle[k].capacity = Problem[0].Vehicle[k].capacity;
		h_Pblm[0].Vehicle[k].start_time = Problem[0].Vehicle[k].start_time;
		h_Pblm[0].Vehicle[k].end_time = Problem[0].Vehicle[k].end_time;
		h_Pblm[0].Vehicle[k].duration_constraint = Problem[0].Vehicle[k].duration_constraint;
	}

	CHECK(cudaMemcpyToSymbol(Pblm, h_Pblm, sizeof(pblm)));
}*/

