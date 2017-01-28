#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

void random_numbers (int n, double a[n]){
    srand(time(NULL)); 
    for(int i=0; i<n; i++) 
        a[i] = ((double) rand())/RAND_MAX;
}

void static_numbers (int n, double a[n]){
    srand(time(NULL)); 
    for(int i=0; i<n; i++) 
        a[i] = (double) n - i;
}

void print_numbers (int n, double a[n]){
    for(int i=0; i<n; i++)
        printf("%.5e\n", a[i]);
}

void print_inumbers (int n, int a[n]){
    for(int i=0; i<n; i++)
        printf("%d\n", a[i]);
}

int compare (const void *e1, const void *e2){
    double *i1 = (double*)e1;
    double *i2 = (double*)e2;
    return ((*i1 < *i2) ? -1 : (*i1 > *i2) ? +1 : 0);
}

void distribute(const int num_proc, const int block_size, const int last_block_size, int displ[num_proc], int send_counts[num_proc]){
	displ[0] = 0;
	send_counts[0] = block_size;
	
	int i;
	for (i = 1; i < num_proc; i++){
		displ[i] = displ[i-1] + block_size;
		
		if (i != num_proc - 1) {
			send_counts[i] = block_size;			
		} else {
			send_counts[i] = last_block_size;
		}
	}
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s N\n", argv[0]);
        return 1;
    }

    int asize = atoi(argv[1]);
    int rank, procs, N_block, N_block_last;
    int root = 0;
	
	int *send_counts, *displ, *recv_counts;
	
	int j, k; // buckets
	
	double *a, *local_arr, *splitter, *all_splitters, *buckets, *bucketBuffer, *localBucket;
	double *OutputBuffer, *Output;
	double start, stop; // time
	
	a = (double*)calloc(asize, sizeof(double));
	// static_numbers(asize, a);
	random_numbers(asize, a);
	
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

	// start time measurement
	if (rank == root){
		start = omp_get_wtime();
	}
	
	N_block = asize / procs;
	N_block_last = asize / procs + asize % procs;
    if (rank == procs - 1) {
		N_block = N_block_last;
	}
	
    local_arr = (double *) malloc (N_block * sizeof (double));
	displ = (int*)calloc(procs, sizeof(int));
	send_counts = (int*)calloc(procs, sizeof(int));
	distribute(procs, asize / procs, N_block_last, displ, send_counts);
	MPI_Scatterv(a, send_counts, displ, MPI_DOUBLE, local_arr, N_block, MPI_DOUBLE, root, MPI_COMM_WORLD);
	
	// printf("[Scatterv]\n");
	// for (int i = 0; i < procs; i++) {
		// MPI_Barrier(MPI_COMM_WORLD);
		// if (i == rank) {
			// printf("Rank: %d\n", rank);
			// print_numbers(N_block, local_arr);
		// }
	// }
	// return 0;
	
	qsort ((double *) local_arr, N_block, sizeof(double), compare);

	splitter = (double *) malloc (sizeof(double) * (procs - 1));
	for (int i = 0; i < (procs - 1); i++){
		splitter[i] = local_arr[asize/(procs*procs) * (i+1)];
    }
	if (rank == root) {
		all_splitters = (double *) malloc (sizeof (double) * procs * (procs-1));
	}
	MPI_Gather(splitter, procs-1, MPI_DOUBLE, all_splitters, procs-1, MPI_DOUBLE, root, MPI_COMM_WORLD);	
	
    if (rank == root){
        qsort ((double *) all_splitters, procs*(procs-1), sizeof(double), compare);
        for (int i = 0; i < (procs-1); i++)
			splitter[i] = all_splitters[(procs-1)*(i+1)];
    }
	
    MPI_Bcast(splitter, procs-1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	
    // local buckets
    buckets = (double *) malloc (sizeof (double) * (asize + procs));
	
	// distributing local array among local buckets
    j = 0;
    k = 1;
	for (int i = 0; i < N_block; i++){
		if (j < procs - 1){
			if (local_arr[i] < splitter[j]){
				buckets[((N_block + 1) * j) + k++] = local_arr[i]; 
			} else {
				buckets[(N_block + 1) * j] = k-1;
				k=1;
				j++;
				i--;
			}
		} else 
			buckets[((N_block + 1) * j) + k++] = local_arr[i]; // if we are behind of the last splitter we can place items there without additional calculations
	}
	buckets[(N_block + 1) * j] = k-1;
	
	// printf("After local buckets sorting\n");
	// for (int i = 0; i < procs; i++) {
		// MPI_Barrier(MPI_COMM_WORLD);
		// if (i == rank) {
			// printf("Rank: %d\n", rank);
			// print_numbers(asize + procs, buckets);
		// }
	// }
	
	
	bucketBuffer = (double*)calloc(2 * asize + procs, sizeof(double));
	distribute(procs, N_block + 1, N_block_last + 1, displ, send_counts);
	recv_counts = (int*)calloc(procs, sizeof(int));
	for(int i = 0; i < procs; i++)
		recv_counts[i] = N_block_last + 1;
	MPI_Alltoallv(buckets, send_counts, displ, MPI_DOUBLE, bucketBuffer, recv_counts, displ, MPI_DOUBLE, MPI_COMM_WORLD);

	// printf("After Alltoall\n");
	// for (int i = 0; i < procs; i++) {
		// MPI_Barrier(MPI_COMM_WORLD);
		// if (i == rank) {
			// printf("Rank: %d\n", rank);
			// print_numbers(asize + procs, bucketBuffer);
		// }
	// }
	
	/**** Rearranging BucketBuffer ****/
	localBucket = (double *) malloc (sizeof (double) * 2 * N_block);

	int count = 1;

	for (j=0; j<procs; j++) {
	k = 1;
		for (int i=0; i< (int)bucketBuffer[(asize / procs + 1) * j]; i++)
			localBucket[count++] = bucketBuffer[(asize / procs + 1) * j + k++];
	}
	localBucket[0] = count-1;
	
	// printf("After rearraging\n");
	// for (int i = 0; i < procs; i++) {
		// MPI_Barrier(MPI_COMM_WORLD);
		// if (i == rank) {
			// printf("Rank: %d\n", rank);
			// print_numbers(asize + 1, localBucket);
		// }
	// }
    
	/**** Sorting Local Buckets using qsort ****/
	int NoElementsToSort = (int)localBucket[0];
	qsort ((double *) &localBucket[1], NoElementsToSort, sizeof(double), compare); 

	// printf("After local sorting\n");
	// for (int i = 0; i < procs; i++) {
		// MPI_Barrier(MPI_COMM_WORLD);
		// if (i == rank) {
			// printf("Rank: %d\n", rank);
			// print_numbers(asize + 1, localBucket);
		// }
	// }

	
	/**** Gathering sorted sub blocks at root ****/
	if(rank == root) {
		OutputBuffer = (double *) malloc (sizeof(double) * 2 * asize);
		Output = (double *) malloc (sizeof (double) * asize);
		recv_counts = (int*)calloc(procs, sizeof(int));
		displ = (int*)calloc(procs, sizeof(int));
		displ[0] = 0;
		for(int i = 0; i < procs; i++){
			recv_counts[i] = 2 * N_block + 1;
			if (i > 0) displ[i] = displ[i-1] + 2 * N_block + 1;
		}
	}

	MPI_Gatherv(localBucket, count, MPI_DOUBLE, OutputBuffer, recv_counts, displ, MPI_DOUBLE, root, MPI_COMM_WORLD);

	// CHECKED
	
	if (rank == root){
		count = 0;
		for(int j = 0 ; j < procs; j++){
			k = 1;
			for(int i = 0; i<OutputBuffer[(2 * N_block + 1) * j]; i++) 
				 Output[count++] = OutputBuffer[(2 * N_block + 1) * j + k++];
    	}
		// print_numbers(asize, Output);
		
		stop = omp_get_wtime();
		printf("Parallel: %.5f seconds\n", stop - start);
		
		start = omp_get_wtime();
		qsort ((double *) a, asize, sizeof(double), compare);
		stop = omp_get_wtime();
		printf("Serial: %.5f seconds\n", stop - start);
		
		for(int i = 0; i < asize; i++)
			if(Output[i] != a[i]){
				printf("Error! Sorted data isn't correct at %d (%.5f != %.5f)\n", i, Output[i], a[i]);
				break;
			}
		
		free(OutputBuffer);
		free(Output);
		free(all_splitters);
	}
	
	free(localBucket);
	free(bucketBuffer);
	free(buckets);
	free(splitter);
	free(local_arr);
	
    MPI_Finalize();
	
	return 0;
}