#include <mpi.h>
#include <stdio.h>
int main(int argc, char** argv)
{ 
    MPI_Init(&argc, &argv); 
    int num;
    MPI_Comm_rank(MPI_COMM_WORLD,&num);

    if(num==0) { // "master"
    std::cout<<"Hello from rank "<< num<<" recieveing messages " << std::endl;
    MPI_Status status;
    char txt[100];
    int ierror = MPI_Recv(txt, 100, MPI_CHAR,1, 42, MPI_COMM_WORLD, &status);
    std::cout << txt << "\n";
    }
    else {
    std::cout<<"Hello from rank "<<num<<" sending message " << std::endl;
    std::string text="Hello world!";
    MPI_Request request;
    int ierror = MPI_Isend(const_cast<char*>(text.c_str()), text.size()+1, MPI_CHAR,0, 42, MPI_COMM_WORLD, &request);
    }

    MPI_Finalize();

   return 0;
}

