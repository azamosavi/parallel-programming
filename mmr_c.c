//compile: mpicc mmr_c  mmr_c.c
//Run: mpirun -np 32 -npernode 8 -hostfile ~/nodes  ./mmr_c 240
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define  prob_size ((argc > 1) ? atoi(argv[1]) : 64)//getting the size of problem from user
#define amat(I,J) A[I+N*J]
#define bmat(I,J) B[I+N*J]
#define cmat(I,J) C[I+N*J]
#define Atemp(I,J) Atemp[I+N*J]
#define Atemp_bf(I,J) Atemp_bf[I+N*J]

void multp (double *A, double *B, double *C, int N, int Nc, int index,double * Atemp_bf);

int main(int argc, char *argv[])

{
  int myid, numprocs;
  int N, Nc, sum;
  int i,j,k,index,p,t,l;
  double *A, *B, *C,*Atemp, *Atemp_bf;
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
    Atemp_bf=(double *)malloc (N*Nc* sizeof(double)); //Different buffer for sending in asynchronous message passing
   assert (Atemp_bf!=NULL);


  //Inititalize in order to test, using all ones will lead to the result all Ns.
 // set_all amat and bmat to zero and all cmat to one

  for (i=0; i<N*Nc; i++){
       A[i]=(double) rand();
       B[i]=(double) rand() ;
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
     multp ( A, B, C,  N,  Nc,myid, Atemp_bf);

  //Next step is to do multiplication in a ring among processors
  //So that each processor gives N*nc of A to others to complete 
  //the multiplication

       index=myid;
     for (p=1; p<numprocs;p++){
        index=(index+1) % numprocs;
      
      // Asunchronous message passing
      // Using two different buffers for send and recieve
     // MPI_Irecv(Atemp,N*Nc, MPI_DOUBLE, recv_from,recv_from, MPI_COMM_WORLD,&request);
       MPI_Isend(Atemp_bf,N*Nc, MPI_DOUBLE, to_send, myid, MPI_COMM_WORLD,&request2);
       MPI_Irecv(Atemp,N*Nc, MPI_DOUBLE, recv_from,recv_from, MPI_COMM_WORLD,&request);
       MPI_Wait (&request, &status);
       MPI_Wait (&request2, &status);


       multp (Atemp, B, C,  N,  Nc, index, Atemp_bf);
    
   }


  // Stop timer
  MPI_Barrier(MPI_COMM_WORLD);
  end_time=MPI_Wtime();

// Printing the result to test
/*

 for(j=0; j<Nc ; j++)
      for(i=0;i<N;i++){
          printf("process %d ,a_ element(%d,%d)=%f\n",myid,i,j,amat(i,j) );
          printf("process %d ,b_ element(%d,%d)=%f\n",myid,i,j,bmat(i,j) );
          printf("process %d ,c_ element(%d,%d)=%f\n",myid,i,j,cmat(i,j) );
      }
*/
    
if (myid==0) {
   tp=end_time-start_time;
   printf("Total time of multiplication is %7.3e\n",tp);
   }
   

  free(A);
  free(B);
  free(C);
  free(Atemp_bf);
  free(Atemp);
  MPI_Finalize();
return;
}

void multp (double *A, double *B, double *C, int N, int Nc, int index, double *Atemp_bf){
 int sum;
 int i,j,k;
 for (i=0; i<N; i++){
    for(j=0; j<Nc; j++){
       sum=0;
       Atemp_bf(i,j)=amat(i,j);
       for (k=0; k<Nc; k++){
           sum=sum+amat(i,k)*bmat(k+index*Nc,j);
         }
         cmat(i,j)=cmat(i,j)+sum;
        }
      }       
}

