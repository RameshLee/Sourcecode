This document has the code base for GPU-based algorithm for online stochastic dial-a-ride problem.

Degree of dynamism and expected scenarios are set as constants.
The problem is defined under the constant named Str1.
Total requests and total vehicles are set as constants.

To run the program:

Login to the compute canada cluster.
Start an interactive session using the command: salloc --account=<account_name> --gres=gpu:<gpu_name>:1 --cpus-per-task=4 --mem=define --time=0-00:define
Import the cuda module in ComputeCanada using the command: "module load cuda/versionNumber".
Compile and run the code using the command: "nvcc sourcecode.cu -o Run && ./Run".
