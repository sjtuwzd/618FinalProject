# 618FinalProject

Usage

Squential Version:

make
./port_analyzer 50 10000 EXP best
50 is the number of trading days to run
10000 is the number of runs of Mento Carlo Simulation
EXP is the suffix for the datafiles
best is one of the two mode this program runs, if it's something other than best, it will simply run the performance evaluation
if it is exactly "best" then the program would solve the single-transaction problem

Parallel Version:
mpirun -np mpirun -np 12 ./port_analyzer_parallel 50 10000 EXP best

Special thanks to 
https://pbpython.com/monte-carlo.html
