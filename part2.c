#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define NS 64


int main(int argc,char* argv[])
{
    int i,j,k,kp;
    int from,to;
    double matrix_a[NS][NS];
    double matrix_b[NS][NS];
    double matrix_c[NS][NS];
    int myrank, p, index;
    int to_send, recv_from;
    double starttime, endtime;

    
    MPI_Request request, request2;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if ((NS % p)!=0)
    {
    printf ("n should be devisible by number of processors \n");
    return -1;
    }
    MPI_Status status;
    // Define new MPI data type for column decomposition
    MPI_Datatype acol;
    MPI_Datatype acoltype;
    MPI_Type_vector(NS,NS/p,NS,MPI_DOUBLE,&acol);       
    MPI_Type_commit(&acol);
    MPI_Type_create_resized(acol, 0, (NS/p)*sizeof(double), &acoltype);
    MPI_Type_commit(&acoltype); 

    if(myrank == 0)
    { 
        //Initialize global matrices
          for(i=0; i<NS; ++i){
            for(j=0; j<NS; ++j){
                matrix_a[i][j] = rand() % 10;
                matrix_b[i][j] = rand() % 10;
                matrix_c[i][j]=0;
           }
         }

/*
        printf("\n matrix A \n");
        printf("----------------\n");
        for(i=0; i<NS; ++i)
        {
            for(j=0; j<NS; ++j)
                printf("%.1f ", matrix_a[i][j]);
            printf("\n");
        }
        
        printf("\n matrixB \n");
        printf("----------------\n");
         for(i=0; i<NS; ++i)
        {
            for(j=0; j<NS; ++j)
                printf("%.1f ", matrix_b[i][j]);
            printf("\n");
        }
*/ 
  

      }
      // defining local matrices for each processor 
      double my_c[NS/p][NS];
      double my_a[NS/p][NS];
      double my_b[NS][NS/p];
      double temp_c[NS/p][NS/p];
   
      starttime= MPI_Wtime();
      // scatter matrix A and B
      MPI_Scatter (matrix_a, NS*NS/p, MPI_DOUBLE, my_a, NS*NS/p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Scatter (matrix_b, 1, acoltype, my_b, NS*(NS/p), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
      MPI_Barrier (MPI_COMM_WORLD); 
     // Initialize local matrices of each processor
      for (i=0; i<NS/p; i++) {
        for (j=0; j<NS; j++) {
         my_c[i][j]=0;
        }
      }
    
      // Finding neighbours for each processor
      to_send=(myrank-1+p)% p;
      recv_from= (myrank+1) % p;
      from = myrank * NS/p;
      to = (myrank+1) * NS/p;
     
     
//multiplication starts here
      for (kp=0; kp<p; kp++){
        index=(myrank+kp) % p;     
        for (i=0; i<NS/p; i++) {
          for (j=0; j<NS/p; j++) {
            temp_c[i][j]=0;
            for (k=0; k<NS; k++){
               temp_c[i][j] += my_a[i][k]*my_b[k][j];
               my_c[i][j+index*NS/p]=temp_c[i][j];
       }
    // printf ("value of C from %d i: %d,j: %d is %d\n",myrank,i,j+index*n/p, my_c[i][j+index*n/p]);
    }
  }   
     
      MPI_Isend(my_b,NS*NS/p, MPI_DOUBLE, to_send, myrank, MPI_COMM_WORLD,&request2);
      MPI_Irecv(my_b,NS*NS/p, MPI_DOUBLE, recv_from,recv_from, MPI_COMM_WORLD,&request);
      MPI_Wait (&request, &status);
      MPI_Wait (&request2, &status);      
} 
   
//end of multiplication

      MPI_Barrier (MPI_COMM_WORLD);   
     //gathering the results of all proceesors
      MPI_Gather(my_c,NS*(NS/p),MPI_DOUBLE, matrix_c[from], NS*(NS/p), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);
      endtime   = MPI_Wtime();

    
    //Display area
   if(myrank == 0)
    {
 /*
       printf("\nRESULT: matrix c \n");
        printf("----------------\n");
         for (i=0; i<NS; i++)
         {
          printf("\n");

          for (j=0; j<NS; j++)
           printf("%.1f  ", matrix_c[i][j]);
         }
*/
        printf("\n\nParellel Time %f seconds\n",endtime-starttime);
    }
  
    
    MPI_Finalize();
    return 0;
}



