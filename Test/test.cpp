#include <iostream>
#include <stdlib.h>
#include <omp.h>

int main() {
	#ifdef _OPENMP
		std::cout << "OpenMP is on" << std::endl;
	#endif
	int a[10];
	int i =0;
	int id, proc_nums, thread_num;

	omp_set_num_threads(5);

	proc_nums = omp_get_num_procs();
	printf("Num of procs: %d\n", proc_nums);
	thread_num = omp_get_num_threads();

	for (i = 0; i < 10; i++) {
		a[i] = i;
	}
	id = omp_get_thread_num();
	printf("Thread num in consequtive part 1: %d\n", id);
	
	#pragma omp parallel shared(a) private(id, i)
	{
		id = omp_get_thread_num();
		printf("Thread num in parallel part: %d\n", id);
	}

	id = omp_get_thread_num();
	printf("Thread num in consequtive part 2: %d\n", id);
	return 0;
}