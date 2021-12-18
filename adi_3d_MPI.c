#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define  Max(a,b) ((a)>(b)?(a):(b))

void 
init(int n, float *A)
{ 
	for (int i = 0; i < n; i++) {
	    for (int j = 0; j < n; j++) {
	        for (int k = 0; k < n; k++) {
	            if (i == 0 || i == n - 1 || j == 0 || j == n - 1 || k == 0 || k == n - 1) {
	                A[i * n * n + n * j + k] = 0.;
                } else {
                    A[i * n * n + n * j + k] = (4. + i + j + k);
                }
            }
        }
    }
} 

float 
relax_mpi(int n, float *A, int myrank, float *myA1, float *myA2, int *cnt1, int *disp1, int *cnt2, int *disp2, MPI_Datatype COLRES)
{
    int myn = cnt1[myrank];
    int arr_size = myn * n;
    float eps = 0.;
    float eps2 = 0.;
    MPI_Scatterv(A, cnt1, disp1, COLRES, myA1, arr_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    int cnt = 0;
    for (int j = 0; j < myn; j++) {
        for (int i = 0; i < n; i++) {
            if (i != 0 && i != n - 1) {
                myA1[cnt] = (myA1[cnt - 1] + myA1[cnt + 1]) / 2.; 
            }
            cnt++;
        }
    }
    MPI_Gatherv(myA1, arr_size, MPI_FLOAT, A, cnt1, disp1, COLRES, 0, MPI_COMM_WORLD);
    arr_size = cnt2[myrank];
    myn = arr_size / (n * n);
    MPI_Scatterv(A, cnt2, disp2, MPI_FLOAT, myA2, arr_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    cnt = 0;
    for (int i = 0; i < myn; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (j != 0 && j != n - 1) {
                    myA2[cnt] = (myA2[cnt - n] + myA2[cnt + n]) / 2.;
                }
                cnt++;
            }
        }
    }
    cnt = 0;
    for (int i = 0; i < myn; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (k != 0 && k != n - 1) {
                    float e = myA2[cnt];
                    myA2[cnt] = (myA2[cnt - 1] + myA2[cnt + 1]) / 2.;
                    eps = Max(eps, fabs(e - myA2[cnt]));
                }
                cnt++;
            }
        }
    }
    MPI_Gatherv(myA2, arr_size, MPI_FLOAT, A, cnt2, disp2, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Reduce(&eps, &eps2, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eps2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    return eps2;
}

void 
wrap(int n, float * A, int myrank, int sz, int itmax, float mineps) 
{
    int *cnt1, * disp1;
    cnt1 = calloc(sz, sizeof(*cnt1));
    disp1 = calloc(sz, sizeof(*disp1));
    cnt1[0] = (n * n / sz + (0 < (n * n) % sz));
    disp1[0] = 0;
    for (int i = 1; i < sz; i++) {
        cnt1[i] = cnt1[i - 1];
        if (i == (n * n) % sz) {
            cnt1[i]--;
        }
        disp1[i] = disp1[i - 1] + cnt1[i - 1];
    }
    int * cnt2, * disp2;
    cnt2 = calloc(sz, sizeof(*cnt2));
    disp2 = calloc(sz, sizeof(*disp2));
    cnt2[0] = (n / sz + (0 < n % sz)) * n * n;
    disp2[0] = 0;
    for (int i = 1; i < sz; i++) {
        cnt2[i] = cnt2[i - 1];
        if (i == n % sz) {
            cnt2[i] -= n * n;
        }
        disp2[i] = disp2[i - 1] + cnt2[i - 1];
    }
    int arr_size = cnt1[myrank] * n;
    float *myA1 = calloc(arr_size, sizeof(*myA1));
    arr_size = cnt2[myrank];
    float *myA2 = calloc(arr_size, sizeof(*myA2));
    MPI_Datatype COL, COLRES;
    MPI_Type_vector(n, 1, n * n, MPI_FLOAT, &COL);
    MPI_Type_commit(&COL);
    MPI_Type_create_resized(COL, 0, sizeof(*myA1), &COLRES);
    MPI_Type_commit(&COLRES);
    for(int it = 0; it < itmax; it++)
    {
	    float eps = 0.;
	    eps = relax_mpi(n, A, myrank, myA1, myA2, cnt1, disp1, cnt2, disp2, COLRES);
	    if (eps < mineps) {
            break;
        }
    }
    MPI_Type_free(&COL);
    MPI_Type_free(&COLRES);
    free(myA1);
    free(myA2);
    free(cnt1);
    free(disp1);
    free(cnt2);
    free(disp2);
}

void 
verify(int n, float *A)
{
	float s = 0.;
	for (int i = 1; i < n - 1; i++) {
	    for (int j = 1; j < n - 1; j++) {
	        for (int k = 1; k < n - 1; k++) {
		        s += A[i * n * n + n * j + k] * (i + 1) * (j + 1) * (k + 1) / (n * n * n);
	        }
        }
    }
	printf("S = %f\n", s);
}

int main (int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    const float mineps = 0.1e-7;
    const int itmax = 100;
    int numproc = 64;
    FILE * f = fopen("output_proc64", "w");
    for (int n = 16; n <= 256; n += 16) {
        double start, end;
        int myrank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        float *A;
        if (myrank == 0) {
            A = calloc(n * n * n, sizeof(*A));
            init(n, A);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == 0) {
            start = MPI_Wtime();
        }
        wrap(n, A, myrank, size, itmax, mineps);
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == 0) {
            end = MPI_Wtime();
            fprintf(f, "%d\t%d\t%lf\n", n, numproc, end - start);
            // verify(n, A); we comment so as not to waste time on checking as it has already been checked
            free(A);
        }
    }
    fclose(f);
    MPI_Finalize();
	return 0;
}
