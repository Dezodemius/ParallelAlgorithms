//
// Created by dezodemius on 19.04.2020.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <limits.h>

float dot(int n, float* partA, float* partB);
float* generateArray(int n);
float dot_optimized(int n, float * partA, float * partB);

// Точка входа в программу.
void main(int argc, char** argv)
{
  int np, IDp, N = 1e3;
  double sum = 0.0;
  
  float* a = generateArray(N);
  float* b = generateArray(N);
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &IDp);
  int offset = (N / np) * IDp;
  
  double t0 = MPI_Wtime();
  double tempSum = dot(N / np, a + offset, b + offset);
  MPI_Allreduce(&tempSum, &sum, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
  double t1 = MPI_Wtime();
  
  if (IDp == 0){
    double t_1 = 0.000003;
    double t_p = t1 - t0;
    double s = t_1 / t_p;
    double e = s / np * 100;
    printf("%d (\\verb\"Allreduce\") & %.8f & %.8f & %.4f \\%% & & & \\%%\\\\ \\hline\n", np, t_p, s, e);
  }
  
  MPI_Finalize();
}

// Скалярное произведение
float dot(int n, float * partA, float * partB)
{
  int i;
  float sum = 0.0;
  
  for(i = 0; i < n; i++)
    sum += partA[i] * partB[i];
  
  return sum;
}

// Скалярное произведение
float dot_optimized(int n, float * partA, float * partB)
{
  int i;
  float sum = 0.0;
  
  for(i = 0; i < n; i+=4)
    sum += partA[i] * partB[i] + partA[i + 1] * partB[i + 1] +
           partA[i + 2] * partB[i + 2] + partA[i + 3] * partB[i + 3];
  
  return sum;
}

// Сгенерировать массив.
float * generateArray(int n)
{
  int i;
  float *array = (float*)malloc(n*sizeof(float));
  
  for(i = 0; i < n; i++)
    array[i] = i / 100.0;
  
  return array;
}
