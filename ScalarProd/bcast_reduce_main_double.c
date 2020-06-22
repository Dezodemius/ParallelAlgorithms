//
// Created by dezodemius on 19.04.2020.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <limits.h>

double dot(int n, double* partA, double* partB);
double* generateArray(int n);

// Точка входа в программу.
void main(int argc, char** argv)
{
  int np, IDp, N = 1e3;
  double sum = 0.0;
  
  double *a = (double*)malloc(N*sizeof(double));
  double *b = (double*)malloc(N*sizeof(double));

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &IDp);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  a = generateArray(N);
  b = generateArray(N);
  MPI_Bcast(a, N, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(b, N, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  int offset = (N / np) * IDp;
  
  double t0 = MPI_Wtime();
  double tempSum = dot(N / np, a + offset, b + offset);
  MPI_Reduce(&tempSum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (IDp == 0){
    double t_1 = 0.000004;
    double t_p = MPI_Wtime() - t0;
    double s = t_1 / t_p;
    double e = s / np * 100;
    printf("& %.8f & %.4f & %.4f \\\\ \\hline\n", t_p, s, e);
  }
  
  MPI_Finalize();
  free(a);
  free(b);
}

// Скалярное произведение
double dot(int n, double * partA, double * partB)
{
  int i;
  double sum = 0.0;
  
  for(i = 0; i < n; i++)
    sum += partA[i] * partB[i];
  
  return sum;
}

// Скалярное произведение
double dot_optimized(int n, double * partA, double * partB)
{
  int i;
  double sum = 0.0;
  
  for(i = 0; i < n; i+=4)
    sum += partA[i] * partB[i] + partA[i + 1] * partB[i + 1] +
           partA[i + 2] * partB[i + 2] + partA[i + 3] * partB[i + 3];
  
  return sum;
}

// Сгенерировать массив.
double * generateArray(int n)
{
  int i;
  double *array = (double*)malloc(n*sizeof(double));
  
  for(i = 0; i < n; i++)
    array[i] = 1.0;
  
  return array;
}
