#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <sstream> 

#include "mpi.h"

typedef double real;

real** getSymmetricMatrix(int n) {
  real** matrix = (real**)malloc(n * sizeof(real*));
  for (int i = 0; i < n; i++)
    matrix[i] = (real*)malloc(n * sizeof(real));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      matrix[i][j] = matrix[j][i] = i / 100.0;
    }
  }
  return matrix;
}

real* getVector(int n) {
  real* array = (real*)malloc(n * sizeof(real));

  for (int i = 0; i < n; ++i)
    array[i] = i / 100.0;

  return array;
}

void showVector(real* vector, int n) {
  for (int j = 0; j < n; j++) {
    std::cout << vector[j] << "\t";
  }
}

void showMatrix(real** matrix, int n) {
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      std::cout << matrix[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}

real* multiply(real** M, int n, int m, real* V, int k) {
  real* resultVector = (real*)malloc(n * sizeof(real));
  real sum;
  for (int i = 0; i < m; i++) {
    sum = 0.0;
    for (int j = 0; j < n; j++) {
      sum += M[i][j] * V[j];
    }
    resultVector[k + i] = sum;
  }
  return resultVector;
}

// MPI_Allgather
int main(int argc, char** argv)
{
  std::stringstream convert(argv[1]);

  int n;
  if (!(convert >> n))
    n = 0;
  double start, stop;

  real** A, * b, ** A_;
  real* result = (real*)malloc(n * sizeof(real));
  int IDp, np;
  int root = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &IDp);

  double m = n / np;

  A = getSymmetricMatrix(n);
  b = getVector(n);
  A_ = getSymmetricMatrix(n * m);

  MPI_Scatter(A, n * m, MPI_DOUBLE, A_, n * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  start = MPI_Wtime();
  real* x = multiply(A, n, m, b, IDp * m);
  stop = MPI_Wtime();
  MPI_Allgather(x, m, MPI_DOUBLE, result, m, MPI_DOUBLE, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  free(x);

  if (IDp == 0)
  {
    double t = stop - start;
    double t0 = 0.0025865;
    double Sp = t0 / t;
    double Ep = Sp * 100.0 / np;
    std::cout << t << " & " << Sp << " & " << Ep << std::endl;
  }

  free(result);
  free(b);

  free(A);
  free(A_);
  MPI_Finalize();
  return 0;
}

//// MPI_Allreduce
//int main(int argc, char** argv)
//{
//  std::stringstream convert(argv[1]);
//
//  int n;
//  if (!(convert >> n))
//    n = 0;
//  double start, stop;
//
//  real** A, * b, ** A_;
//  real* result = (real*)malloc(n * sizeof(real));
//  int IDp, np;
//  int root = 0;
//
//  MPI_Init(&argc, &argv);
//  MPI_Comm_size(MPI_COMM_WORLD, &np);
//  MPI_Comm_rank(MPI_COMM_WORLD, &IDp);
//
//  double m = n / np;
//
//  A = getSymmetricMatrix(n);
//  b = getVector(n);
//  A_ = getSymmetricMatrix(n * m);
//
//  MPI_Scatter(A, m, MPI_FLOAT, A_, m, MPI_FLOAT, 0, MPI_COMM_WORLD);
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  start = MPI_Wtime();
//  real* x = multiply(A, n, m, b, IDp * m);
//  stop = MPI_Wtime();
//  MPI_Allreduce(x, result, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  free(x);
//
//  if (IDp == 0)
//  {
//    double t = stop - start;
//    double t0 = 6.55e-05;
//    double Sp = t0 / t;
//    double Ep = Sp * 100.0 / np;
//    std::cout << t << " & " << Sp << " & " << Ep << std::endl;
//  }
//
//  free(result);
//  free(b);
//
//  free(A);
//  free(A_);
//  MPI_Finalize();
//  return 0;
//}

// MPI_Reduce
//int main(int argc, char** argv)
//{
//  std::stringstream convert(argv[1]);
//
//  int n;
//  if (!(convert >> n))
//    n = 0;
//  double start, stop;
//
//  real** A, * b, ** A_;
//  real* result = (real*)malloc(n * sizeof(real));
//  int IDp, np;
//  int root = 0;
//
//  MPI_Init(&argc, &argv);
//  MPI_Comm_size(MPI_COMM_WORLD, &np);
//  MPI_Comm_rank(MPI_COMM_WORLD, &IDp);
//
//  double m = n / np;
//
//  A = getSymmetricMatrix(n);
//  b = getVector(n);
//  A_ = getSymmetricMatrix(n * m);
//
//  MPI_Scatter(A, m, MPI_FLOAT, A_, m, MPI_FLOAT, 0, MPI_COMM_WORLD);
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  start = MPI_Wtime();
//  real* x = multiply(A, n, m, b, IDp * m);
//  stop = MPI_Wtime();
//  MPI_Reduce(x, result, n, MPI_FLOAT, MPI_SUM, root, MPI_COMM_WORLD);
//
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  free(x);
//
//  if (IDp == 0)
//  {
//    double t = stop - start;
//    double t0 = 6.55e-05;
//    double Sp = t0 / t;
//    double Ep = Sp * 100.0 / np;
//    std::cout << t << " & " << Sp << " & "  << Ep << std::endl;
//  }
//
//  free(result);
//  free(b);
//
//  free(A);
//  free(A_);
//  MPI_Finalize();
//  return 0;
//}

// Последовательно.
//int main(int argc, char** argv)
//{
//  std::stringstream convert(argv[1]);
//
//  int n;
//  if (!(convert >> n))
//    n = 0;
//  double start, stop;
//
//  real **A = getSymmetricMatrix(n);
//  real *b = getVector(n);
//  
//  MPI_Init(&argc, &argv);
//
//  start = MPI_Wtime();
//  real* x = multiply(A, n, n, b, 0);
//  stop = MPI_Wtime();
//
//  std::cout << (stop - start) << std::endl;
//
//  MPI_Finalize();
//
//  delete[] A[0];
//  delete[] A;
//  delete b;
//  delete x;
//
//  return 0;
//}
