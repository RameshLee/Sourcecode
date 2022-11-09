struct myList
{
	int list[n];
	int size;

	myList()
	{
		size = 0;
		for (int i = 0; i < n; i++)
			list[i] = 0;
	}

	void clear()
	{
		size = 0;
		for (int i = 0; i < n; i++)
			list[i] = 0;
	}

	~myList() {}

	void print()
	{
		printf("size: %d\n", size);
		for (int i = 0; i < size; i++)
			printf("list[%d]=%d\n", i, list[i]);
		printf("---\n");
	}

	void simple_print()
	{
		printf("size: %d =>[", size);
		for (int i = 0; i < size; i++)
			printf("%d ", list[i]);
		printf("]\n");
	}

	void push(int element)
	{
        if (size >= n)
        {
            printf("ERROR on utility: push() size goes greater than n: causing heap corruption!:\n");
        }

        list[size] = element;
        size++;

        //sort();
	}

	void popout(int element)
	{
        if (size == 0)
        {
            printf("ERROR on utility: No items to pop out from the list:\n");
        }

       if (size == 1)
        {
            printf("ERROR on utility: pop() size goes zero: causing heap corruption!:\n");
        }


		int *req_list = new int[size-1];

		int index = 0;
		for (int i = 0; i < size; i++)
			if (list[i] != element)
			{
				req_list[index] = list[i];
				index++;
			}

		size--;
		for (int i = 0; i < size; i++)
			list[i] = req_list[i];

		delete[] req_list;

		//sort();

		/*int Activate = 0;
		for (int i = 0; i < size; i++)
		{
			if (list[i] == element)
			{
				Activate = 1;
				size--;
			}

			if (Activate == 1)
				list[i] = list[i + 1];
		}*/
	}

    void sort()
    {
        for (int i=0; i<size; i++)
            for (int j=i+1; j<size; j++)
                if (i != j && list[i] > list[j])
                {
                    int temp = list[i];
                    list[i] = list[j];
                    list[j] = temp;
                }

    }

	void randomize()
	{
		if (size > 0)
		{
			random_shuffle(&list[0], &list[size]);
		}
		else
		{
			printf("ERROR ON UTILITY: CANNOT SORT A ZERO ARRAY!!!!\n");
		}
	}

	bool isFound(int element)
	{
		for (int i = 0; i < size; i++)
			if (list[i] == element)
				return 1;

		return 0;
	}

};

size: 1 =>[96 ]
Control enters the scenario 15
-- Before removing request 96 from vehicle 4

TotalRequests: 120 TotalVehicles: 11
sz: 14 Path: 10 23 11 6 130 143 131 126 37 157 61 181 45 165
sz: 10 Path: 2 122 22 31 142 151 18 138 54 174
sz: 16 Path: 7 8 1 127 4 128 121 46 124 166 52 172 17 137 55 175
sz: 10 Path: 21 24 41 141 144 161 49 169 5 125
sz: 18 Path: 3 28 123 148 43 40 19 15 39 32 139 58 160 152 163 159 178 135
sz: 8 Path: 14 13 35 134 133 155 20 140
sz: 14 Path: 44 164 48 168 38 158 33 153 50 170 9 129 60 180
sz: 4 Path: 27 147 62 182
sz: 6 Path: 26 146 59 179 34 154
sz: 8 Path: 30 150 12 132 42 162 29 149
sz: 16 Path: 16 25 136 145 47 167 57 177 53 173 51 171 36 156 56 176

Available: size: 36 =>[14 47 48 33 38 19 44 56 51 37 32 8 11 57 18 58 17 49 53 55 61 41 28 39 60 59 46 35 31 16 36 42 54 50 52 4 ]
UnAvailable: size: 58 =>[63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 62 ]
Sampled: size: 0 =>[]
Rejected: size: 1 =>[62 ]