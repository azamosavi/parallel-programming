#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define n 64



 
int main(int argc,char* argv[])
{
    int i,j,k;
    double matrix_a[n][n];
    double matrix_b[n][n];
    double matrix_c[n][n];
    MPI_Datatype acol;
    MPI_Datatype acoltype;
    int myrank, p;
    double starttime, endtime;
 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if ((n % p)!=0)
    {
    printf ("n should be devisible by number of processors \n");
    return -1;
    }
    //MPI_Status status;
    // Define new MPI data type for column decomposition
    MPI_Type_vector(n,n/p,n,MPI_DOUBLE,&acol);       
    MPI_Type_commit(&acol);
    MPI_Type_create_resized(acol, 0, (n/p)*sizeof(double), &acoltype);
    MPI_Type_commit(&acoltype); 
    
    
    starttime   = MPI_Wtime();
       if(myrank == 0)
    { 
        // Initialize global matrices
        for(i=0; i<n; ++i){
            for(j=0; j<n; ++j){
                matrix_a[i][j] = rand() % 10;
                matrix_b[i][j] = rand() % 10;
                matrix_c[i][j]=0;
           }
         }
/*
        printf("\n matrix A \n");
        printf("----------------\n");
        for(i=0; i<n; ++i)
        {
            for(j=0; j<n; ++j)
                printf("%.1f  ", matrix_a[i][j]);
            printf("\n");
        }

        printf("\n matrixB \n");
        printf("----------------\n");
         for(i=0; i<n; ++i)
        {
            for(j=0; j<n; ++j)
                printf("%.1f  ", matrix_b[i][j]);
            printf("\n");
        }
*/ 
    }
     
      // defining local matrices for each processor
      double my_b[n][n/p];
      double my_c[n][n/p];

     //Broadcast A to all processors and scatter B to all processors with their chunks of columns
      MPI_Bcast (&matrix_a, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
      MPI_Scatter (matrix_b, 1, acoltype, my_b, n*(n/p), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
      MPI_Barrier (MPI_COMM_WORLD);
   
     // start of multiplication
      for (i=0; i<n; i++) {
        for (j=0; j<n/p; j++) {
          my_c[i][j]=0;
          for (k=0; k<n; k++){
	         my_c[i][j] += matrix_a[i][k]*my_b[k][j];
       }
    // printf ("value of C from %d is th %d and %d is %.1f\n",myrank,i,j, my_c[i][j]);
    }
  }
      MPI_Barrier (MPI_COMM_WORLD);
     // Gather the results from all processors
      MPI_Gather(my_c,n*n/p, MPI_DOUBLE,matrix_c[0]+myrank*n/p, 1, acoltype, 0, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);
      endtime   = MPI_Wtime();
    //end of multiplication
    
    //Display area
    if(myrank == 0)
    {
/*     
        printf("\nRESULT: matrix c \n");
        printf("----------------\n");
         for (i=0; i<n; i++)
         {
          printf("\n");

          for (j=0; j<n; j++)
           printf("%.1f ", matrix_c[i][j]);
         }
*/
         printf("\n\nParellel Time %f seconds\n",endtime-starttime);
    }

    MPI_Finalize();
    return 0;
}


