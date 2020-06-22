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
double dot_optimized(int n, double * partA, double * partB);

// Точка входа в программу.
void main(int argc, char** argv)
{
  int np, IDp, N = 1e8;
  double sum = 0.0;
  
  double* a = generateArray(N);
  double* b = generateArray(N);
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &IDp);
  int offset = (N / np) * IDp;
  
  double t0 = MPI_Wtime();
  double tempSum = dot_optimized(N / np, a + offset, b + offset);
  MPI_Allreduce(&tempSum, &sum, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
  double t1 = MPI_Wtime();
  
  if (IDp == 0){
    double t_1 = 0.176161;
    double t_p = t1 - t0;
    double s = t_1 / t_p;
    double e = s / np * 100;
    printf("& %.4f & %.4f & %.2f \\%%\\\\ \\hline\n", t_p, s, e);
  }
  
  MPI_Finalize();
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
    array[i] = i / 100.0;
  
  return array;
}
