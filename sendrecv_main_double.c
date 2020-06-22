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
double dot_optimized(int n, double* partA, double* partB);
double* generateArray(int n);

// Точка входа в программу.
void main(int argc, char** argv)
{
  int processesNumber, processId, N = 1e8;
  double sum = 0.0, tempSum = 0.0;
  double* a = generateArray(N);
  double* b = generateArray(N);
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &processesNumber);
  MPI_Comm_rank(MPI_COMM_WORLD, &processId);
  
  int n = N / processesNumber;
  int offset = n * processId;
  double *buff_a = (double*)malloc(n*sizeof(double));
  double *buff_b = (double*)malloc(n*sizeof(double));
      
  MPI_Sendrecv(a + offset, n, MPI_DOUBLE, (processId + 1) % processesNumber, 1,
    buff_a, n, MPI_DOUBLE, (processId + processesNumber - 1) % processesNumber, 1, MPI_COMM_WORLD, &status);

  MPI_Sendrecv(b + offset, n, MPI_DOUBLE, (processId + 1) % processesNumber, 1,
    buff_b, n, MPI_DOUBLE, (processId + processesNumber - 1) % processesNumber, 1, MPI_COMM_WORLD, &status);
  double t0 = MPI_Wtime();
  if (processId > 0){
    tempSum = dot_optimized(n, buff_a, buff_b);
    MPI_Send(&tempSum, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
  } else {
    sum = dot_optimized(n, buff_a, buff_b);
    int i;
    for (i = 1; i < processesNumber; i++){
      MPI_Recv(&tempSum, n, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
      sum += tempSum;  
    }
    
    double t_1 = 0.1761;
    double t_p = MPI_Wtime() - t0;
    double s = t_1 / t_p;
    double e = s / processesNumber * 100;
    printf("& %.8f & %.4f & %.4f \\\\ \\hline\n", t_p, s, e);
  }  

  MPI_Finalize();
}

// Скалярное произведение
double dot(int n, double partA[], double partB[])
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
