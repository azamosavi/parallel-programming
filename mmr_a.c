//compile: mpicc mmr_a  mmr_a.c
//Ru : mpirun -np 32 -npernode 8 -hostfile ~/nodes  ./mmr_a 240
//
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#define  prob_size ((argc > 1) ? atoi(argv[1]) : 64)//the size of matrix is the input argument from user 
#define amat(I,J) A[I+N*J]
#define bmat(I,J) B[I+N*J]
#define cmat(I,J) C[I+N*J]
#define Atemp(I,J) Atemp[I+N*J]


void multp (double *A, double *B, double *C, int N, int Nc, int index);

int main(int argc, char *argv[])

{
  int myid, numprocs;
  int  N,Nc, sum;
  int i,j,k,index,p,t,l;
  double *A, *B, *C,*Atemp;
  int to_send, recv_from;
  double start_time, end_time,tp;
  

  MPI_Status status;
  MPI_Request request, request2;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

 
  // Calculate number of rows each processor holds
   N=prob_size;
   Nc=N/numprocs;
  

  // Allocate Space for matrces A, B and C
   A=(double *)malloc (N*Nc* sizeof(double));
   assert (A!=NULL);
   B=(double *)malloc (N*Nc* sizeof(double));
   assert (B!=NULL);
   C=(double *)malloc (N*Nc* sizeof(double));
   assert (C!=NULL);
   Atemp=(double *)malloc (N*Nc* sizeof(double));
   assert (Atemp!=NULL);
   

  //Inititalize in order to test, using all ones will lead to the result all Ns.
 // set_all amat and bmat to zero and all cmat to one

  for (i=0; i<N*Nc; i++){
       A[i]=(int) rand()%10;
       B[i]=(int) rand()%10 ;
       C[i]=0;
       Atemp[i]=A[i];
  //     printf("the values are %f,%f, %f .\n", amat(j,myid+i),bmat(j,myid+i), cmat(j,myid+i));
    }
   
  // Finding neighbours for each processor
      to_send=(myid-1+numprocs)% numprocs;
      recv_from= (myid+1) % numprocs;
 
    

  // Start the timer
  MPI_Barrier (MPI_COMM_WORLD);
  start_time=MPI_Wtime();

  // Perform  matrix multiplication here
  //First step : Ck=Ck+Ak*Bkk
     multp ( A, B, C,  N,  Nc,myid);

  //Next step is to do multiplication in a ring among processors
  //So that each processor gives N*nc of A to others to complete 
  //the multiplication

       index=myid;
     for (p=1; p<numprocs;p++){
        index=(index+1) % numprocs;

      //Synchronous message passing
      MPI_Sendrecv(Atemp, N*Nc, MPI_DOUBLE, to_send, to_send,Atemp,N*Nc, MPI_DOUBLE,recv_from, myid, MPI_COMM_WORLD, NULL);
       multp (Atemp, B, C,  N,  Nc, index);
    
   }


  // Stop timer
  MPI_Barrier(MPI_COMM_WORLD);
  end_time=MPI_Wtime();
 //printing the results to test
 /*
 for(j=0; j<Nc ; j++)
      for(i=0;i<N;i++){
          printf("process %d ,a_ element(%d,%d)=%f\n",myid,i,j,amat(i,j) );
          printf("process %d ,b_ element(%d,%d)=%f\n",myid,i,j,bmat(i,j) );
          printf("process %d ,c_ element(%d,%d)=%f\n",myid,i,j,cmat(i,j) );
      }
*/
 //Printing the time it takse to do the multiplication
 if (myid==0) {
   tp=end_time-start_time;
   printf("Total time of multiplication is %7.3e\n",tp);
   }




  free(A);
  free(B);
  free(C);
  free(Atemp);
  MPI_Finalize();
return;
}

void multp (double *A, double *B, double *C, int N, int Nc, int index){
 int sum;
 int i,j,k;
 for (i=0; i<N; i++){
    for(j=0; j<Nc; j++){
       sum=0;
       for (k=0; k<Nc; k++){
           sum=sum+amat(i,k)*bmat(k+index*Nc,j);
          }
         cmat(i,j)=cmat(i,j)+sum;
        }
      }       
}
