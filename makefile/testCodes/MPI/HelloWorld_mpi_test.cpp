#include <mpi.h>
#include <stdio.h>
#include <vector>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);

 
   int * values = new int[16];
 
   for(int i = 0; i<4; i++)
          values[4*world_rank + i] =  world_rank;           
   
  
   {
      printf("Print in rank %d \n", world_rank);
      for(int j=0; j<4; j++)
	printf("v[ %d ]= %d \n", 4*world_rank+ j, values[4*world_rank+j]); 
   }
  
 
   MPI_Barrier(MPI_COMM_WORLD);
   if(world_rank==1)
   {
      printf("Print from master rank \n");
      for(int i=0; i<16; i++)
	printf("v[%d]=%d \n", i, values[i]);
  }


    // Finalize the MPI environment.
    MPI_Finalize();
}

