#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int **allocate_board(int Rows, int Cols)
{
	// allocate Rows rows, each row is a pointer to int
	int **board = (int **)malloc(Rows * sizeof(int *));
	int row;
	// for each row allocate Cols ints
	for (row = 0; row < Rows; row++) {
		board[row] = (int *)malloc(Cols * sizeof(int));
	}
	return board;
}

void print_matrix(int **A,int n){
	int i,j;
	for( i=0;i < n;i++)
	{
		for( j=0;j < n;j++)
			printf("value of matrix at i %d , j %d is %d  \n ",i,j, A[i][j]);
	}

}

int** sumofMatrix(int n, int **A, int **B){
	int i,j;
	int **C=allocate_board(n,n);
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			C[i][j]=A[i][j]+B[i][j];
		}
	}
	return C;
}

int** subofMatrix(int n, int **A, int **B){
	int i,j;
	int **C=allocate_board(n,n);
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			C[i][j]=A[i][j]-B[i][j];
		}
	}
	return C;
}

int** mulofMatrix(int n, int **A, int **B){
	int i,j,k;
	int **C=allocate_board(n,n);
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			C[i][j]=0;
			for (k=0; k<n; k++){
				C[i][j]=C[i][j]+A[i][k]*B[k][j];
			}
		}
	}
	return C;
}


int** comp_strassen(int n, int **A, int **B){
	int i,j;
	int **C=allocate_board(n,n);
	if (n<=16)
		C=mulofMatrix(n, A,B);
	else
	{

		int N=n;
		n=n/2;
		int **A11=allocate_board(n,n);
		int **A12=allocate_board(n,n);
		int **A21=allocate_board(n,n);
		int **A22=allocate_board(n,n);
		int **B11=allocate_board(n,n);
		int **B12=allocate_board(n,n);
		int **B21=allocate_board(n,n);
		int **B22=allocate_board(n,n);
		int **C11=allocate_board(n,n);
		int **C12=allocate_board(n,n);
		int **C21=allocate_board(n,n);
		int **C22=allocate_board(n,n);
		int **P1=allocate_board(n,n);
		int **P2=allocate_board(n,n);
		int **P3=allocate_board(n,n);
		int **P4=allocate_board(n,n);
		int **P5=allocate_board(n,n);
		int **P6=allocate_board(n,n);
		int **P7=allocate_board(n,n);
		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
				A11[i][j]=A[i][j];
				A12[i][j]=A[i][j+N/2];
				A21[i][j]=A[i+N/2][j];
				A22[i][j]=A[i+N/2][j+N/2];
				B11[i][j] = B[i    ][j];
				B12[i][j] = B[i    ][j+N/2];
				B21[i][j] = B[i+N/2][j];
				B22[i][j] = B[i+N/2][j+N/2];
			}
		}




#pragma opm task
		{
			C11=subofMatrix(n,B12,B22);
			P1=comp_strassen(n,A11,C11);
		}
#pragma opm task
		{
			C11=sumofMatrix(n,A11,A12);
			P2=comp_strassen(n,C11,B22);
		}

#pragma opm task
		{
			C11=sumofMatrix(n,A21,A22);
			P3=comp_strassen(n,C11,B11);
		}

#pragma opm task
		{
			C11=subofMatrix(n,B21,B11);
			P4=comp_strassen(n,A22,C11);
		}


#pragma opm task
		{
			C11=sumofMatrix(n,A11,A22);
			C12=sumofMatrix(n,B11,B22);
			P5=comp_strassen(n ,C11,C12);
		}



#pragma opm task
		{
			C11=subofMatrix(n,A12,A22);
			C12=sumofMatrix(n,B21,B22);
			P6=comp_strassen(n,C11,C12);
		}
#pragma opm task
		{
			C11=subofMatrix(n,A21,A11);
			C12=sumofMatrix(n,B11,B12);
			P7=comp_strassen(n,C11,C12);
		}

#pragma omp taskwait

		C11=sumofMatrix(n,P5,P4);
		C11=subofMatrix(n,C11,P2);
		C11=sumofMatrix(n,C11,P6);

		C12=sumofMatrix(n,P1,P2);

		C21=sumofMatrix(n,P3,P4);

		C22=sumofMatrix(n,P5,P1);
		C22=subofMatrix(n,C22,P3);
		C22=sumofMatrix(n,C22,P7);

		for(i=0;i < n ;i++)
			for(j=0;j < n ;j++)
			{
				C[i    ][j] = C11[i][j];
				C[i    ][j+N/2] = C12[i][j];
				C[i+N/2][j] = C21[i][j];
				C[i+N/2][j+N/2] = C22[i][j];
			}


		free(A11);free(A12);free(A21);free(A22);
		free(B11);free(B12);free(B21);free(B22);
		free(C11);free(C12);free(C21);free(C22);
		free(P1);free(P2);free(P3);free(P4);
		free(P5);free(P6);free(P7);

	}   
	return C;
}


int isPowerOfTwo (unsigned int x)
{
	unsigned int powerOfTwo = 1;

	while (powerOfTwo < x && powerOfTwo < 2147483648)
		powerOfTwo *= 2;
	return (x == powerOfTwo);
}

int main(int argc, char *argv[])
{
	int n,threads;
	if (argc==3){
		n=strtol(argv[1],NULL, 10);
		threads=strtol(argv[2],NULL, 10);   

	}
	else{
		printf("Please enter the dimensions of matrices and the number of threads  \n");
		scanf("%ld", &n);
		scanf("%ld", &threads);
	}
	while (isPowerOfTwo(n)!=1  && isPowerOfTwo(threads)!=1 ){

		printf("they should be power of two please enter the dimensions of matrices and the number of threads  \n");
		scanf("%ld", &n);
		scanf("%ld", &threads);


	}

	//printf ("result of power %.8f\n", log10(n)/log10(2));
	//printf ("the dimension is %d \n", n);
	int **A=allocate_board(n,n);
	int **B=allocate_board(n,n);
	int **C=allocate_board(n,n);   
	double ostart, oend;
	int i,j;
	double begV, endV;

	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			A[i][j]=rand() % 10 + 1;
			B[i][j]=rand() % 10 + 1;

		}
	}
	printf("Strasen Algorithm\n");
	ostart=omp_get_wtime();

	omp_set_num_threads(threads);
#pragma omp parallel 
	//#pragma omp parallel num_threads(threads)
	{
#pragma omp single
		{
			C=comp_strassen(n,A,B);

		}
	}
	oend=omp_get_wtime();
	//print_matrix(C,n);
	printf("time of omp for strassen is %.8f \n",((double) oend-ostart));
	printf("normal multiplication algorithm \n");

	int **C1=allocate_board(n,n);
	ostart=omp_get_wtime();
	C1=mulofMatrix(n,A,B);
	oend=omp_get_wtime();
	//print_matrix(C1,n);
	printf("time of normal is %.8f \n",((double) oend-ostart));
	return 0;
}

