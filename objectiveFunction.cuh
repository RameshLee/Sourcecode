
struct cost
{

	float overall_current_cost;
	float overall_ideal_cost;

	//////////////////////////
	bool isMED_feasible;
	bool isOverall_batteryfeasible;
	float total_feasible_charging_time;

	//EV-aspect
	float excessrideTime;


	int start;
	int vehicle;
	float waitingTime;
	float PassengerRideTime;
	float early_arrival_time;
	float late_arrival_time;
	float passenger_idle_time;
	float passenger_waiting_time;


	//just for fleet minimsation
	int vehid;
	//////////////////////////

	float travel_cost;


	// float excess_ride_time, passenger_waiting_time, route_duration, early_arrival_time;


	//to check violations
	float max_battery_Violation;
	float batteryViolation; // bV storage
	float time_window; //twV storage
	float ride_time; //rtV storage
	int load; //CaV
	float duration; // duV
	bool isFeasible;

	float on_board_user_wait_time;
	float route_duration;

	float violation;

	__host__ __device__
		cost()
	{
		reset();
	}

	__host__ __device__
		void Copy(cost &inCost)
	{

		on_board_user_wait_time = inCost.on_board_user_wait_time;
		route_duration = inCost.route_duration;
		early_arrival_time = inCost.early_arrival_time;

		overall_current_cost = inCost.overall_current_cost;
		overall_ideal_cost = inCost.overall_ideal_cost;

		total_feasible_charging_time = inCost.total_feasible_charging_time;
		isMED_feasible = inCost.isMED_feasible;

		isOverall_batteryfeasible = inCost.isOverall_batteryfeasible;

		max_battery_Violation = inCost.max_battery_Violation;
		batteryViolation = inCost.batteryViolation;

		excessrideTime = inCost.excessrideTime;
		travel_cost = inCost.travel_cost;

		time_window = inCost.time_window;
		ride_time = inCost.ride_time;
		load = inCost.load;
		duration = inCost.duration;

		isFeasible = inCost.isFeasible;


		waitingTime = inCost.waitingTime;
		PassengerRideTime = inCost.PassengerRideTime;
		early_arrival_time = inCost.early_arrival_time;
		late_arrival_time = inCost.late_arrival_time;
		passenger_idle_time = inCost.passenger_idle_time;
		passenger_waiting_time = inCost.passenger_waiting_time;

	}

	__host__ __device__
		void reset()
	{

		on_board_user_wait_time = 0;
		route_duration = 0;
		early_arrival_time = 0;

		overall_current_cost = 0;
		overall_ideal_cost = 0;

		total_feasible_charging_time = 0;

		// excess_ride_time = 0, passenger_waiting_time = 0, route_duration = 0, early_arrival_time = 0;
		max_battery_Violation = 0;
		batteryViolation = 0;

		excessrideTime = 0;
		travel_cost = 0;
		time_window = 0;
		ride_time = 0;
		load = 0;
		duration = 0;
		isFeasible = false;

		isMED_feasible = false;
		isOverall_batteryfeasible = false;
	}

	__host__ __device__
		bool getFeasibility()
	{
		//if all zero means, no violation i.e. feasible.
		isFeasible = false;

		if (time_window == 0)
			if (ride_time == 0)
				if (load == 0)
					if (duration == 0)
						isFeasible = true;

		if (time_window < 0.001)
			if (ride_time < 0.001)
				if (load < 0.001)
					if (duration < 0.001)
						isFeasible = true;

		/*if (travel_cost < 0.01)
		isFeasible = false;*/

		return isFeasible;
	}

	__host__ __device__
		bool get_Overall_battery_feasibility()
	{
		//if all zero means, no violation i.e. feasible.

		if (batteryViolation < 0.001)
			isOverall_batteryfeasible = true;
		else
			isOverall_batteryfeasible = false;

		return isOverall_batteryfeasible;
	}

	__host__ __device__
		bool getPossibility()
	{
		bool Bool = false;
		if (time_window == 0)
			if (load == 0)
				Bool = true;

		return Bool;
	}

	__host__ __device__
		void print()
	{
		printf("Cost Information:\n");
		printf("isFeasible %s\n", isFeasible == 0 ? "false" : "true");
		printf("travelcost %f (%s)\n", travel_cost, travel_cost == 0 ? "true" : "false");
		printf("timeWindow %f (%s) \n", time_window, time_window == 0 ? "true" : "false");
		printf("rideTimeViol %f (%s) \n", ride_time, ride_time == 0 ? "true" : "false");
		printf("loadViol %d (%s) \n", load, load == 0 ? "true" : "false");
		printf("durationViol %f (%s) \n", duration, duration == 0 ? "true" : "false");
		printf("batteryViol %f (%s) \n", batteryViolation, batteryViolation == 0 ? "true" : "false");
		printf("\n");
	}

	__host__ __device__
		~cost() {}
};

struct coefficient
{
	float delta_max;
	float lambda_max;

	float delta;
	float lambda;
	float alpha;
	float beta;
	float tau;
	float gamma;

	float spl_parameter;

	__host__ __device__
		coefficient() {}

	__host__ __device__
		~coefficient() {}

	__host__ __device__
		void initialize()
	{
		delta_max = max(0.1, sqrt((float)((n - 1)*(m - 1))) / 10.0);
		lambda_max = max(0.01, sqrt((float)((n - 1)*(m - 1))) / 70.0);

		delta_max = 10;// 0.5;
		lambda_max = 0.015*sqrt(n*m);

		reset();
	}

	__host__ __device__
		void reset()
	{
		delta = delta_max;
		lambda = lambda_max;

		alpha = 1;
		beta = 1;
		tau = 1;
		gamma = 1;

		spl_parameter = 1;
	}

	__host__ __device__
		void randomize()
	{
		srand(time_now++);
		delta = rand() % 21 * delta_max / 20;

		srand(time_now++);
		lambda = rand() % 21 * lambda_max / 20;

	}

/*	__host__ __device__
		float cost_function(cost &Cost)
	{
		float Double;

		Double = (w1 * Cost.travel_cost + w2 * Cost.excessrideTime)
			+ alpha* Cost.load
			+ beta* Cost.time_window
			+ tau* Cost.ride_time
			+ gamma* Cost.duration;

		return Double;
	}*/

	__host__ __device__
		float cost_function(cost &Cost)
	{
		float Double;

		Double = (w1 * Cost.travel_cost + w2 * Cost.excessrideTime + w3 *Cost.on_board_user_wait_time + w4 * Cost.route_duration + w5 * Cost.early_arrival_time)
			+ alpha* Cost.load
			+ beta* Cost.time_window
			+ tau* Cost.ride_time
			+ gamma* Cost.duration;

		return Double;
	}


	__host__ __device__
		float cost_function_with_batteryViolation(cost &Cost)
	{
		float Double;

		//float a = RandomNumber(0, 1);

		/*Double = (w1 * Cost.travel_cost + w2 * Cost.excessrideTime) + (1 * Cost.max_battery_Violation)
		+ alpha* Cost.load
		+ beta* Cost.time_window
		+ tau* Cost.ride_time
		+ gamma* Cost.duration;*/

		Double = (w1 * Cost.travel_cost + w2 * Cost.excessrideTime) + (spl_parameter * Cost.batteryViolation)
			+ alpha* Cost.load
			+ beta* Cost.time_window
			+ tau* Cost.ride_time
			+ gamma* Cost.duration;

		return Double;
	}

	__host__ __device__
		float cost_function_with_batteryViolation_2(cost &Cost)
	{
		float Double;

		//float a = RandomNumber(0, 1);

		Double = (w1 * Cost.travel_cost + w2 * Cost.excessrideTime) + (spl_parameter * Cost.max_battery_Violation)
			+ alpha* Cost.load
			+ beta* Cost.time_window
			+ tau* Cost.ride_time
			+ gamma* Cost.duration;

		return Double;
	}

	__host__ __device__
		float cost_function_with_batteryViolation_and_MED_travel_cost(cost &Cost, cost &MED_Solution_Cost)
	{
		float Double;

		//float a = RandomNumber(0, 1);

		Double = (500 * Cost.batteryViolation) + MED_Solution_Cost.travel_cost;

		/*Double = (w1 * Cost.travel_cost + w2 * Cost.excessrideTime) + (1 * Cost.max_battery_Violation)
		+ alpha* Cost.load
		+ beta* Cost.time_window
		+ tau* Cost.ride_time
		+ gamma* Cost.duration
		+ MED_Solution_Cost.travel_cost + alpha * MED_Solution_Cost.time_window;*/

		return Double;
	}

	__host__ __device__
		float cost_function_only_for_MED_Solution(cost &MED_Solution_Cost)
	{
		float Double;

		Double = MED_Solution_Cost.travel_cost;

		return Double;
	}

	__host__ __device__
		float get_violation(cost &Cost)
	{
		float Double;

		Double =
			+alpha* Cost.load
			+ beta* Cost.time_window
			+ tau* Cost.ride_time
			+ gamma* Cost.duration;

		return Double;
	}

	__host__ __device__
		void update(cost Cost)
	{
		if (Cost.load)
			alpha *= (1 + delta);
		else
			alpha /= (1 + delta);

		if (Cost.time_window)
			beta *= (1 + delta);
		else
			beta /= (1 + delta);

		if (Cost.ride_time)
			tau *= (1 + delta);
		else
			tau /= (1 + delta);

		if (Cost.duration)
			gamma *= (1 + delta);
		else
			gamma /= (1 + delta);

		if (Cost.batteryViolation)
			spl_parameter *= (1 + delta);
		else
			spl_parameter /= (1 + delta);

		alpha = max(alpha, 1e-100);
		beta = max(beta, 1e-100);
		gamma = max(gamma, 1e-100);
		tau = max(tau, 1e-100);
		spl_parameter = max(spl_parameter, 1e-100);
	}


};
