This document has the code base for GPU-based algorithm for online stochastic dial-a-ride problem.

Degree of dynamism and expected scenarios are set as constants.
The problem is defined under the constant named Str1.
Total requests and total vehicles are set as constants.

To run the program:

Import the cuda module in ComputeCanada using the command: "module load cuda/11.0".
Compile and run the code using the command: "nvcc sourcecode.cu -o TestRun && ./TestRun".
