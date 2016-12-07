#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#define n 64


int main(int argc,char* argv[])
{
    int i,j,k,kp;
    int from,to;
    double matrix_a[n][n];
    double matrix_b[n][n];
    double matrix_c[n][n];
    

    MPI_Datatype ablock;
    MPI_Datatype ablocktype;
    int myrank, p, index;
    int to_send, recv_from;
    double starttime, endtime;
 
    
    MPI_Request request, request2;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int q=sqrt(p);
    if ((q % 2)!=0)
    {
    printf ("The number of processors should be power of four\n");
    return -1;
    }
    if ((n % q)!=0)
    {
    printf ("n should be devisible by square root of number of processors \n");
    return -1;
    }
     //Define new MPI data type for block decomposition
    MPI_Type_vector(n/q,n/q,n, MPI_DOUBLE, &ablock);       
    MPI_Type_commit(&ablock);
    MPI_Type_create_resized(ablock, 0,sizeof(double), &ablocktype);
    MPI_Type_commit(&ablocktype); 
    
     //Define 2D communication of processors
    MPI_Comm comm_2d;
    int dims[2], free_coords[2], periods[2]; 
    int my2drank, mycoords[2]; 
    dims[0] = dims[1] = q;
    periods[0] = periods[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d); 
    MPI_Comm_rank(comm_2d, &my2drank); 
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords); 
    MPI_Comm row_comm;  
    MPI_Comm col_comm;  
    // Set up row communicators 
    free_coords[0] = 0; 
    free_coords[1] = 1;
    MPI_Cart_sub(comm_2d, free_coords, &row_comm);
    // Set up column communicators 
    free_coords[0] = 1; 
    free_coords[1] = 0;
    MPI_Cart_sub(comm_2d, free_coords, &col_comm);
    
    if(myrank == 0)
    { 
        //Initialize global matrices
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
                printf("%0.1f", matrix_a[i][j]);
            printf("\n");
        }
        printf("\n matrix B \n");
        printf("----------------\n");
        for(i=0; i<n; ++i)
        {
            for(j=0; j<n; ++j)
                printf("%0.1f", matrix_b[i][j]);
            printf("\n");
        }
*/ 
    }
    
      
      
      // defining local matrices for each processor
      double my_c[n/q][n/q];
      double my_a[n/q][n/q];
      double my_b[n/q][n/q];
      double temp_a[n/q][n/q];
       
      starttime   = MPI_Wtime();      
     // define displacement and counts for the scatterv of blocks
      int disps[q*q];
      int counts[q*q];
      for (i=0; i<q; i++) {
        for (j=0; j<q; j++) {
            disps[i*q+j] = (i*n*(n/q))+(j*(n/q));
            counts [i*q+j] = 1;
        }
    }
     
      MPI_Scatterv (matrix_a, counts, disps, ablocktype, my_a, (n/q)*(n/q),MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Scatterv (matrix_b, counts, disps, ablocktype, my_b, (n/q)*(n/q),MPI_DOUBLE, 0, MPI_COMM_WORLD);
      // Initialize local matrices of each processor
      for (i=0; i<n/q; i++){
            for(j=0; j<n/q; j++){
              temp_a[i][j]=my_a[i][j];
              my_c[i][j]=0;
           }
         }
      MPI_Barrier (MPI_COMM_WORLD);

           
     int source = (mycoords[0] + 1) % q;
     int dest = (mycoords[0]+ q - 1) % q;
     int stage;
     for (stage=0; stage < q; stage++) {
        kp = (mycoords[0] + stage) % q;
        if (kp == mycoords[1]) {
              MPI_Bcast(my_a,(n/q)*(n/q), MPI_DOUBLE, kp, row_comm);
              //local multiplication
              for (i=0; i<n/q; i++){
                for(j=0; j<n/q; j++){
                 for(k=0; k<n/q; k++){
                   my_c[i][j]+=my_a[i][k]*my_b[k][j];
                  }
               }
             }


        } else {
              MPI_Bcast(temp_a,(n/q)*(n/q), MPI_DOUBLE, kp, row_comm);
             //local multiplication
              for (i=0; i<n/q; i++){
                for(j=0; j<n/q; j++){
                 for(k=0; k<n/q; k++){
                   my_c[i][j]+=temp_a[i][k]*my_b[k][j];
                  }
               }
             }

        }
      
         MPI_Sendrecv_replace(my_b, (n/q)*(n/q), MPI_DOUBLE, dest, 0, source, 0, col_comm, &status);
 }
      
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_Gatherv(my_c,(n/q)*(n/q), MPI_DOUBLE, matrix_c, counts, disps, ablocktype, 0, MPI_COMM_WORLD);
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
           printf("%0.1f ", matrix_c[i][j]);
         }
*/
         printf("\n\nParellel Time %f seconds\n",endtime-starttime);
    }
    
    MPI_Finalize();
    return 0;
}




