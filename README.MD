The implementation for MPI is very similar to how omp and normal parallel version were implemented. You just add a barrier after each move/check/update.
The formula for dividing the nr of people across the processes has also been rewritten for clarity, but it behaves functionally the same(first ids will
take+1 person untill the rest is 0. Ex: 6 processes,10 people: proc 0: 0,1 ; proc 1: 2,3 ; proc 2: 4,5 ; proc 3: 6,7 ; proc 4: 8 ; proc 5: 9 where format is: 
proc //processNumber// : personID,personID where is the case.). 

Each process beside the main one will first receive its respective number of people like above. The main process(rank 0) will be the one to send the respective
people to each process, and also will send his respective people to himself. Each process will go through the algorithm mentioned above, moving/checking/updating
its chunk. At each iteration of sim time, there will be an "Allreduce" for the infection counter. This will get the max values for each index in each process, 
effectively computing the value for the infection counter. 

EX: Say proc 0 has infection counter to values:(0,0,1) ; proc 1: (1,1,1) ; proc 2: (0,0,0) ; for these values, the infection counter will become (1,1,1) : max_infec1=max(proc_0_infec_1,proc_1_infec_1,proc_2_infec_1) = max(0,1,0) = 1 ; the same is done for all the other values in the infection counter.

The logic for measuring the time in MPI version is done using the course steps: 
    1. Barrier to have all processes start together 
    2. Barrier to know when allprocesses have finished
    3. Each process measureslocaltime between barriers
    4. Send localtimes to Master
    5. Master finds longestlocaltime

Time measurments/Values:

Sending a large nr of processes on one computer will lead to the program crashing because it uses up too much memory. In my case, the limit for when I run it from the PC is 100 processes for 10K file, and 25 processes for the laptop version for the 10K file. A small number of processes generally gives a better speedup compared to a larger number of processes. What's pretty interesting is that for the file where there are 10 people, the best version is the one with only 1 process. 

All the other versions are inferior to the serial version. For epidemics10K, the best calculated value was the one with 5 processes. All the results and tries with different processes can be seen in the files written like this: epidemicsNrPeople_NrProcesses_proc_results.txt (ex: epidemics10_2_proc_results contains the information of 10 people, 2 processes).

