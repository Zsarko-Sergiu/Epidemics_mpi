#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define FILE_BUFFER_SIZE 256


/*

    To measure runtime of a
    program composed by several
    processes:
    1. Barrier to have all processes
    start together
    2. Barrier to know when all
    processes have finished
    3. Each process measures
    localtime between barriers
    4. Send localtimes to Master
    5. Master finds longest
    localtime
*/

typedef struct person
{
    int id,x,y,init_status,movement,amp;
}person;

typedef struct person_thread
{
    int id,first,last;
}person_thread;

//main values used in the program.
int rows,cols,N,sim_time,nr_threads;
person* people;
int* infection_counter;
int** matr;
int size;

//files for output and debugging
FILE* f2;
FILE* f3;
FILE* f4;
int choice;

//time structures
struct timespec start,finish;
struct timespec start_output,finish_output;
double elapsed,elapsed_output,Tserial;
double start_time,stop_time,elapsed_time,total_time;

#define INFECTED_DURATION 2
#define IMMUNE_DURATION 11
#define CHUNK_SIZE 1000

void initialize(char* total_sim_time,char* file_name,char* nr_threads);
void start_simulation_serial(int sim_time,int n,int* infection_counter,person* people,int** matr);
void move_person(person* p,int** matr);
void print_people(person* people,int n,int* infection_counter);
void copy_vector(int* infection_counter_serial,int* infection_counter);
void copy_people(person* people_serial,person* people);
void copy_matrix(int** matr_serial,int** matr);
void write_vars_to_file(person* people,FILE* out);
void start_mpi(int sim_time,int n);
void print_matrix(int** matr);


int main(int argc, char* argv[]) {
    int rank,size;
    MPI_Init(NULL, NULL);

    

    // Initialize rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Start debugging
    printf("Rank %d: MPI initialized. Total ranks: %d\n", rank, size);
    fflush(stdout);  // Ensure all output is flushed

    if (rank == 0) {
        // printf("in rank 0.");
        if (argc < 4) {
            fprintf(stderr, "Error: Not enough arguments. Expected <sim_time> <file_name> <nr_threads>\n");
            MPI_Abort(MPI_COMM_WORLD, 1);  // Terminate all MPI processes
        }

        // Initialize only on rank 0
        printf("Rank 0: Initializing simulation...\n");
        fflush(stdout);
        initialize(argv[1], argv[2], argv[3]);
    }

    // MPI_Finalize();
    // Synchronize after initialization
    // printf("Rank %d: Reached pre-Broadcast barrier.\n", rank);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // Broadcast necessary simulation parameters
    MPI_Bcast(&sim_time, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Broadcast the people array

    // printf("Rank %d: Broadcast before important data structures complete: Sim time: %d, N: %d, rows: %d, cols: %d\n", rank, sim_time, N, rows, cols);

    if(rank!=0)
        people=malloc(N * sizeof(person));    
    // printf("Broadcast before people\n");\
    fflush(stdout);
    MPI_Bcast(people, N * sizeof(person), MPI_BYTE, 0, MPI_COMM_WORLD);

    // if (rank != 0) {
    //     // printf("Rank %d: Received the people array.\n", rank);
    //     // fflush(stdout);
    //     // printf("%d\n",N);
    //     // for (int i = 0; i < N; i++) {
    //     //     printf("Person %d: id=%d, x=%d, y=%d, status=%d, movement=%d, amp=%d\n",
    //     //            i, people[i].id, people[i].x, people[i].y,
    //     //            people[i].init_status, people[i].movement, people[i].amp);
    //     // }
    // }
    fflush(stdout);

    // Broadcast the matrix (matr)
    // printf("Broadcast before matrix\n");
    fflush(stdout);
    if(rank!=0)
    {
        // printf("entered rank\n");
        fflush(stdout);
        matr=(int**)malloc(rows* sizeof(int*));
        for(int i=0;i<rows;i++)
        {
            matr[i]=(int*)malloc(sizeof(int)*cols);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("right before bc\n");
    fflush(stdout);
    for (int i = 0; i < rows; i++) {
        // printf("Trying broadcast ; at i:%d ; rank:%d\n",i,rank);
        // fflush(stdout);
        MPI_Bcast(matr[i], cols, MPI_INT, 0, MPI_COMM_WORLD);
    }
    // printf("finished bc\n");
    fflush(stdout);
    if(rank!=0)
    {
        // printf("Printing matrix from rank %d",rank);
        // fflush(stdout);
        // print_matrix(matr);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank!=0)
    {
        infection_counter=(int*)malloc(sizeof(int) * N);
    }
    // Broadcast the infection counter
    // printf("Broadcast before infection counter\n");
    fflush(stdout);
    MPI_Bcast(infection_counter, N, MPI_INT, 0, MPI_COMM_WORLD);
    // if(rank!=0)
    // {
    //     print_people(people,N,infection_counter);
    // }

    // printf("Rank %d: Broadcast complete. Sim time: %d, N: %d, rows: %d, cols: %d\n", rank, sim_time, N, rows, cols);
    fflush(stdout);
    
    // MPI_Finalize();
    // Call the parallel simulation
    start_mpi(sim_time, N);

    printf("Rank %d: Simulation completed.\n", rank);
    fflush(stdout);

    MPI_Finalize();
    return 0;
}

void print_matrix(int** matr)
{
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            printf("%d ",matr[i][j]);
        }
        printf("\n");
    }
}


void init_matrix(int** matr)
{
    for (int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            matr[i][j]=0;
        }
    }
}

void move_person(person* p,int** matr)
{
    int direction=p->movement;
    //left
    if(direction==3)
    {
        //y decrease, x stays same
        int new_y=p->y-p->amp;
        if(new_y<=0)
        {
            new_y=0;
            p->movement=2;
        }
        matr[p->x][p->y]=0;
        matr[p->x][new_y]=p->id;
        p->y=new_y;
    }
    //right
    if(direction==2)
    {
        //y decrease, x stays same
        int new_y=p->y+p->amp;
        if(new_y>=cols)
        {
            new_y=cols-1;
            p->movement=3;
        }
        matr[p->x][p->y]=0;
        matr[p->x][new_y]=p->id;
        p->y=new_y;
    }
    //down
    if(direction==1)
    {
        //x decrease, y stays same
        int new_x=p->x-p->amp;
        if(new_x<=0)
        {
            new_x=0;
            p->movement=0;
        }
        matr[p->x][p->y]=0;
        // printf("%d %d\n",new_x,p->id);
        matr[new_x][p->y]=p->id;
        p->x=new_x;
    }
    //up
    if(direction==0)
    {
        //x decrease, y stays same
        int new_x=p->x+p->amp;
        if(new_x>=rows)
        {
            new_x=rows-1;
            p->movement=1;
        }
        matr[p->x][p->y]=0;
        matr[new_x][p->y]=p->id;
        p->x=new_x;
    }
}

void check_for_infections(person* people,int n,person* p)
{
    for(int i=0;i<n;i++)
    {
        if(people[i].x==p->x && people[i].y==p->y && people[i].id!=p->id)
        {
            //got one person who is on the same square as the person
            person* p2=&people[i];
            
            //check for immunity ; if either of the people are immune, we don't need to check if they can be infected.
            if(p2->init_status>=2 || p->init_status>=2)
                continue;
            //check for infections for p2 or p
            //if p2 is infected(0)
            if(p2->init_status<=0)
                {
                    //if the other person on the same space is susceptible, theyu can be infected (set to 0). Only if they are susceptible will the infection counter increase and the other person
                    //become infected. If both people are already infected, won't increase infection counter or set p->init_status to 0 since that will effectively reset the infected_duration for 
                    //that other person.
                    if(p->init_status==1)
                        {
                            p->init_status=0;
                            
                        }
                }
            //if p is infected
            else if(p->init_status<=0)
                {
                    if(p2->init_status==1)
                    {
                        p2->init_status=0;
                        
                    }
                }
        } 
    }   
}

void update_infections(person* people, int n,int* infection_counter,person* p) {
        if(p->init_status==0)
        {
            infection_counter[p->id-1]++;
        }
        //infected. check if its value is also bigger than the infected duration
        if(p->init_status<=0 && p->init_status>-INFECTED_DURATION)
        {
            p->init_status--;
        }
            
        else if(p->init_status==-INFECTED_DURATION)
        {
            //finished being infected; make person immune now.
            //1 is hardcoded here because we have two base cases. Either "0" as infected or "1" as susceptible. When we reach "1" we know that the person has finished being immune and is now susceptible again
            p->init_status=1+IMMUNE_DURATION;
        }
        else if(p->init_status>1)
        {
            p->init_status--;
        }
        //if people[i] didn't go through any of the ifs, that means that the person is susceptible (1) value ; nothing needs to be changed in this case.
}


void start_simulation_serial(int sim_time,int n,int* infection_counter,person* people,int** matr)
{
    // printf("%d\n",choice);
    while(sim_time>0)
    {
        // printf("%d\n",sim_time);
        // printf("Simulation nr %d\n",sim_time);
        for(int i=0;i<n;i++)
        {
            move_person(&people[i],matr);
            
        }
        for(int i=0;i<n;i++)
        {
            check_for_infections(people,n,&people[i]);
        }
        for(int i=0;i<n;i++)
            update_infections(people,n,infection_counter,&people[i]);
        // print_matrix(matr);
        // print_people(people,n,infection_counter);
        if(choice==1)
        {
            write_vars_to_file(people,f2);
            fprintf(f2,"\n");
        }
        
        sim_time--;
    }
    clock_gettime(CLOCK_MONOTONIC,&start_output);
    if(choice==0)
    {
        write_vars_to_file(people,f2);
    }
    clock_gettime(CLOCK_MONOTONIC,&finish_output);
    elapsed_output=(finish_output.tv_sec-start_output.tv_sec);
    elapsed_output+=(finish_output.tv_nsec-start_output.tv_nsec)/pow(10,9);
    printf("Finished simulation for serial\n");
    // print_matrix(matr);
    // print_people(people,n,infection_counter);
}


void simulate_mpi(person_thread *my_thread) 
{
    int sim_time_threads=sim_time; // Each process has its own simulation time
    
    MPI_Barrier(MPI_COMM_WORLD);
    start_time=MPI_Wtime();
    // printf("Starting simulation from thread nr:%d\n",my_thread->id);
    // printf("first %d last %d id %d\n",my_thread->first,my_thread->last,my_thread->id);
    fflush(stdout);
    while (sim_time_threads>0) 
    {
        // 1. Move people
        // printf("Sim time from id %d:%d\n",my_thread->id,sim_time_threads);
        for (int i=my_thread->first;i<=my_thread->last;i++) 
        {
            move_person(&people[i],matr);
        }

        // Synchronize all processes after moving people
        // printf("Rank %d reached barrier after moving people.\n", my_thread->id);
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        // 2. Check for infections
        // printf("Rank %d reached barrier after checking people.\n", my_thread->id);
        fflush(stdout);
        for (int i=my_thread->first;i<=my_thread->last;i++) 
        {
            check_for_infections(people,N,&people[i]);
        }

        // Synchronize all processes after checking for infections
        MPI_Barrier(MPI_COMM_WORLD);

        // 3. Update infections
        // printf("Rank %d reached barrier after updating people.\n", my_thread->id);
        fflush(stdout);
        for (int i=my_thread->first;i<=my_thread->last;i++) {
            update_infections(people,N,infection_counter,&people[i]);
        }

        // Synchronize all processes after updating infections
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Allreduce(MPI_IN_PLACE,infection_counter,N,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

        // rank 0 process writes to output file like in serial vers
        if (choice==1 && my_thread->id==0) {
            write_vars_to_file(people,f3);
            fprintf(f3,"\n");
        }

        // Decrease simulation time for the process
        sim_time_threads--;
    }

    MPI_Allreduce(MPI_IN_PLACE,infection_counter,N,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

    stop_time=MPI_Wtime();
    elapsed_time=stop_time-start_time;
    //reduce all values into total time, total time will be the maximum value of time from the processes; variable total time will be sent to main process (rank 0)
    MPI_Reduce(&elapsed_time,&total_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(my_thread->id==0)
    {
        fprintf(f4,"MPI Time:%lf\n",total_time);
    }
}

/*
void start_mpi(int sim_time, int N) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate load distribution
    int base_people = N / size;
    int remainder = N % size;
    int local_count = base_people + (rank < remainder ? 1 : 0);
    int start_index = rank * base_people + (rank < remainder ? rank : remainder);
    int end_index = start_index + local_count;

    if(rank == 0)
        fprintf(f4, "Nr processes:%d\n", size);

    int* local_infection = malloc(sizeof(int) * local_count);
    if (!local_infection) {
        printf("Failed mem alloc for local_infection at rank %d", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Scatter using MPI_Send/Recv
    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            int count = base_people + (i < remainder ? 1 : 0);
            int offset = i * base_people + (i < remainder ? i : remainder);
            
            if (i == 0) {
                // Root process directly uses the global people array
                // No need for local_people, root keeps the data in the global `people`
            } else {
                MPI_Send(&people[offset], count * sizeof(person), MPI_BYTE, i, 0, MPI_COMM_WORLD);
            }
        }
    } else {
        // Non-root processes receive their portion of the people array
        MPI_Recv(&people[start_index], local_count * sizeof(person), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Simulate
    person_thread thread_info;
    thread_info.id = rank;
    thread_info.first = start_index;
    thread_info.last = end_index - 1;

    simulate_mpi(&thread_info);

    MPI_Barrier(MPI_COMM_WORLD);

    // Gather results
    if (rank == 0) {
        for (int source_rank = 1; source_rank < size; source_rank++) {
            int start_index = source_rank * base_people + (source_rank < remainder ? source_rank : remainder);
            int end_index = start_index + base_people + (source_rank < remainder ? 1 : 0) - 1;
            MPI_Recv(&people[start_index], (end_index - start_index + 1) * sizeof(person), MPI_BYTE, source_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        MPI_Send(&people[start_index], (end_index - start_index) * sizeof(person), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
    }

    free(local_infection);

    if (rank == 0) {
        printf("Simulation completed; printing to file.\n");
        write_vars_to_file(people, f3);
        double efficiency = Tserial / (size * total_time);
        double speedup = efficiency * size;
        fprintf(f4, "SPEEDUP:%lf\nEFFICIENCY:%lf", efficiency, speedup);
    }
}
*/

void start_mpi(int sim_time, int N) {
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate load distribution
    int base_people=N/size;
    int remainder=N%size;
    int local_count=base_people+(rank < remainder ? 1 : 0);
    int start_index=rank*base_people+(rank < remainder ? rank : remainder);
    int end_index=start_index+local_count;

    if(rank==0)
        fprintf(f4,"Nr processes:%d\n",size);

    // int* local_infection=malloc(sizeof(int)*local_count);
    // if (!local_infection) 
    // {
    //     printf("Failed mem alloc for local_people/local infec at rank %d",rank);
    //     MPI_Abort(MPI_COMM_WORLD, 1);
    // }


    // printf("Before scatter\n");
    // Scatter using MPI_Send/Recv
    if(rank==0) 
    {
        // Root process sends the data to other processes
        for (int i=1;i<size;i++) {
            int count=base_people+(i < remainder ? 1 : 0);
            int offset=i*base_people+(i < remainder ? i : remainder);
            
            MPI_Send(&people[offset],count*sizeof(person),MPI_BYTE,i,0,MPI_COMM_WORLD);
            // MPI_Send(&infection_counter[offset],count*sizeof(int),MPI_INT,i,0,MPI_COMM_WORLD);
        }
    } 
    else 
    {
        MPI_Recv(&people[start_index],local_count*sizeof(person),MPI_BYTE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        // MPI_Recv(local_infection,local_count*sizeof(int),MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    // Simulate
    person_thread thread_info;
    thread_info.id=rank;
    thread_info.first=start_index;
    thread_info.last=end_index-1;

    // printf("Before simulate\n");
    simulate_mpi(&thread_info);
    // printf("After simulate\n");
    // printf("For rank:%d\n",rank);
    // fflush(stdout);
    // for(int i=start_index;i<end_index;i++)
    // {
    //     printf("I value: %d ; infection counter value:%d\n",i,infection_counter[i]);
    //     fflush(stdout);
    // }

    MPI_Barrier(MPI_COMM_WORLD);
    // Gather results
        // printf("Final people for rank %d:\n",rank);
        // fflush(stdout);
        // print_people(people,N,infection_counter);
        // fflush(stdout);
        // MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank==0) 
    {
        // Root process receives data from other processes
        for(int source_rank=1;source_rank<size;source_rank++) 
        {
            int start_index=source_rank*base_people+(source_rank < remainder ? source_rank : remainder);
            int end_index=start_index+base_people+(source_rank < remainder ? 1 : 0)-1;

            MPI_Recv(&people[start_index],(end_index-start_index+1)*sizeof(person),MPI_BYTE,source_rank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // MPI_Recv(&infection_counter[start_index],(end_index-start_index+1)*sizeof(int),MPI_INT,source_rank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    } 
    else 
    {
        // Non-root processes send their portion of people
        MPI_Send(&people[start_index],(end_index-start_index)*sizeof(person),MPI_BYTE,0,1,MPI_COMM_WORLD);
        // MPI_Send(&infection_counter[start_index],(end_index-start_index)*sizeof(int),MPI_INT,0,1,MPI_COMM_WORLD);
    }

    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // free(local_infection);
    

    if(rank==0) 
    {
        printf("Simulation completed; printing to file.\n");
        write_vars_to_file(people,f3);
        fflush(stdout);
        double efficiency=Tserial/(size*total_time);
        double speedup=efficiency*size;
        fprintf(f4,"SPEEDUP:%lf\nEFFICIENCY:%lf",efficiency,speedup);
    }
}

void copy_vector(int* infection_counter_serial,int* infection_counter)
{
    // printf("%d\n",N);
    
    for(int i=0;i<N;i++)
    {
        // printf("%d\n",i);
        infection_counter_serial[i]=infection_counter[i];
    }
    // printf("test\n");
}

void copy_people(person* people_serial,person* people)
{
    // printf("testing\n");
    for(int i=0;i<N;i++)
    {
        people_serial[i]=people[i];
    }
}

void copy_matrix(int** matr_serial,int** matr)
{
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            matr_serial[i][j]=matr[i][j];
        }
    }
}

void write_vars_to_file(person* people,FILE* out)
{
    // printf("Trying to write vars to file\n");
    for(int i=0;i<N;i++)
    {
        int tmp=people[i].init_status < 0 ? 0 : people[i].init_status;
        char* status;
        if(tmp==0)
        {
            status="Infected (0)";
        }
        if(tmp==1)
        {
            status="Susceptible (1)";
        }
        if(tmp>1)
        {
            status="Immune (>1)";
        }
        // printf(" ID:%d CoordX:%d CoordY:%d Status:%s Infection Counter:%d\n",
        //     people[i].id,
        //     people[i].x,
        //     people[i].y,
        //     status,
        //     infection_counter[people[i].id-1]);
        fprintf(out," ID:%d CoordX:%d CoordY:%d Status:%s Infection Counter:%d\n",
            people[i].id,
            people[i].x,
            people[i].y,
            status,
            infection_counter[people[i].id-1]);

    }   
}



int check_results_people(person* people_serial)
{
    int result=1;
    for(int i=0;i<N;i++)
    {
        if(
            people_serial[i].amp!=people[i].amp || 
            people_serial[i].id!=people[i].id || 
            people_serial[i].init_status!=people[i].init_status || 
            people_serial[i].x!=people[i].x || 
            people_serial[i].y!=people[i].y || 
            people_serial[i].movement!=people[i].movement
        )
            {
                printf("Error ; people differ: %d\n",i);
                result=0;
            }
            
    }
    return result;
}

int check_results_infection(int* infection_counter_serial)
{
    int result=1;
    for(int i=0;i<N;i++)
    {
        if(infection_counter_serial[i]!=infection_counter[i])
        {
            printf("Error ; infections differ : %d\n",i);
            result=0;;
        }
    }
    return result;
}

void check_results(person* people_serial,int** matr_serial,int* infection_counter_serial,int val)
{
    //NOTE: Here you only check for the people and the infections. You don't check for the matrix also since for higher number of threads, the results will most likely differ.
    //At the end, you will have some people that will be at the same position lets say a person with id nr 1, nr 6 and nr 7. In the serial version, the person there will be for example the person
    //with id nr 1 or 6 or 7 while in the parallel version, there can also be those values. The values in the matrix do not account for overlap. Because of this, only the x and y coordinates of the
    //people in the people array will only be accounted for in terms of coordinates.
    if(check_results_people(people_serial) && check_results_infection(infection_counter_serial))
    {
        printf("Values in serial and parallel are the same.\n");
    }
}

void initialize(char* total_sim_time,char* file_name,char* nr_threads_str)
{
    printf("In initialize..\n");
    fflush(stdout);

    FILE* f=fopen(file_name,"r");
    if(!f)
    {
        printf("Error opening file\n");
        exit(1);
    }
    fscanf(f,"%d %d",&rows,&cols);
    
    char line[FILE_BUFFER_SIZE];
    fgets(line,FILE_BUFFER_SIZE,f);
    fscanf(f,"%d",&N);

    printf("Welcome. Select your preffered mode:\n 0 -> Normal mode (evolution of each person not printed)\n 1 -> Debug mode: (evolution of each person printed.)\n");
    fflush(stdout);
    scanf("%d",&choice);

    // printf("%d\n",choice);

    //allocate dynamically both people array and matrix for bigger values.
    people=(person*)malloc(sizeof(person) * N);
    if(!people)
    {
        printf("Error at mem alloc\n");
        return;
    }
    
    int cnt=0;
    //N by N matrix ; values in the matrix will represent the id's of the people
    matr=(int**)malloc(sizeof(int*)*rows);
    if(!matr)
    {
        printf("error allocating memory for matrix 2\n");
        return;
    }
    for(int i=0;i<rows;i++)
    {
        matr[i]=(int*)malloc(sizeof(int) * cols);
        if(!matr[i])
        {
            printf("Error; ran out of memory allocating matrix 1\n");
            return;
        }
    }

    // printf("Before init matr\n");
    // fflush(stdout);
    init_matrix(matr);
    
    fgets(line,FILE_BUFFER_SIZE,f);
    while(fgets(line,FILE_BUFFER_SIZE,f))
    {
        sscanf(line,"%d %d %d %d %d %d",
            &people[cnt].id,
            &people[cnt].x,
            &people[cnt].y,
            &people[cnt].init_status,
            &people[cnt].movement,
            &people[cnt].amp);
        
        matr[people[cnt].x][people[cnt].y]=people[cnt].id;
        cnt++;
    }

    int** matr_reset=(int**)malloc(sizeof(int*)*rows);
    if(!matr_reset)
    {
        printf("error allocating memory for matrix 2\n");
        return;
    }
    for(int i=0;i<rows;i++)
    {
        matr_reset[i]=(int*)malloc(sizeof(int) * cols);
        if(!matr_reset[i])
        {
            printf("Error; ran out of memory allocating matrix 1\n");
            return;
        }
    }

    person* people_reset=(person*)malloc(sizeof(person) * N);
    if(!people_reset)
    {
        printf("Error at mem alloc\n");
        return;
    }
    // printf("Before copy people\n");
    // fflush(stdout);
    copy_people(people_reset,people);
    copy_matrix(matr_reset,matr);

    infection_counter=(int*)malloc(sizeof(int)*N);
    if(!infection_counter)
    {
        printf("Error allocating memory for infection counter\n");
        return;
    }
    for(int i=0;i<N;i++)
    {
        
        infection_counter[i]=0;
    }

    // print_people(people,N,infection_counter);
    // printf("\n");
    
    sim_time=atoi(total_sim_time);
    nr_threads=atoi(nr_threads_str);

    int* infection_counter_serial=(int*)malloc(sizeof(int)*N);
    if(!infection_counter_serial)
    {
        printf("Error allocating memory\n");
    }
    int** matr_serial=(int**)malloc(sizeof(int*)*rows);
    for(int i=0;i<rows;i++)
    {
        matr_serial[i]=(int*)malloc(sizeof(int) * cols);
        if(!matr[i])
        {
            printf("Error; ran out of memory allocating matrix 1\n");
            return;
        }
    }
    person* people_serial=(person*)malloc(sizeof(person)*N);
    // printf("Before creating paths\n");
    // fflush(stdout);
    //_serial_out.txt -> 12 characters (leave 1 for '/0')
    //_parallel_out.txt -> 14 characters (leave 1 for '/0')
    char* path_serial=(char*)malloc(sizeof(char) * (strlen(file_name)+12));
    char path_mpi[50];
    // printf("After creating mpi file\n");
    // fflush(stdout);

    strcpy(path_serial,file_name);
    path_serial=strtok(path_serial,".txt");
    snprintf(path_mpi,50,"%s_mpi_out.txt",path_serial);
    strcat(path_serial,"_serial_out.txt");
    // printf("Path mpi:%s\n",path_mpi);
    fflush(stdout);

    char* results_file=(char*)malloc(sizeof(char)*(strlen(file_name)+9));
    strcpy(results_file,file_name);
    results_file=strtok(results_file,".txt");
    strcat(results_file,"_results.txt");

    // printf("Before opening files\n");
    // fflush(stdout);
    // printf("MPI path:%s",path_mpi);
    // fflush(stdout);
    f2 = fopen(path_serial,"w");
    f4=fopen(results_file,"w");
    f3=fopen(path_mpi,"w");

    fprintf(f4,"nr simulations:%d\n",sim_time,nr_threads);
    //serial
    //start:
    // printf("\nBefore starting simulation\n");
    // fflush(stdout);
    clock_gettime(CLOCK_MONOTONIC,&start);
    start_simulation_serial(sim_time,N,infection_counter,people,matr);
    clock_gettime(CLOCK_MONOTONIC,&finish);

    elapsed=(finish.tv_sec-start.tv_sec);
    elapsed+=(finish.tv_nsec-start.tv_nsec)/pow(10,9);
    Tserial=elapsed-elapsed_output;
    fprintf(f4,"SERIAL:%lf\n",Tserial);

    //store the values for matr,infection counter and people somewhere and reset the values to work on the mpi version.
    
        copy_vector(infection_counter_serial,infection_counter);
        
        copy_matrix(matr_serial,matr);
        
        copy_people(people_serial,people);
        
    //reset the values in matr,people and infection_counter

    copy_people(people,people_reset);
    copy_matrix(matr,matr_reset);

    for(int i=0;i<N;i++)
    {
        infection_counter[i]=0;
    }

}

void print_people(person* people,int N,int* infection_counter)
{
    for(int i=0;i<N;i++)
    {
        int tmp=people[i].init_status < 0 ? 0 : people[i].init_status;
        // if(tmp>1)
        // {
        //     tmp=1;
        // }
        printf("Person nr %d: ID:%d CoordX:%d CoordY:%d Status:%d Movement:%d Amp:%d Infection Counter:%d\n",
            i+1,
            people[i].id,
            people[i].x,
            people[i].y,
            tmp,
            people[i].movement,
            people[i].amp,
            infection_counter[people[i].id-1]);

    }   
}