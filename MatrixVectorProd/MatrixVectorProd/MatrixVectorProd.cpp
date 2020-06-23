#include <iostream>
#include <ctime>
#include <cstdlib>
#include <sstream> 

#include "mpi.h"

typedef float real;

real** getSymmetricMatrix(int n) {
  real** matrix = new real*[n];
  for (int i = 0; i < n; i++)
    matrix[i] = new real[n];

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      matrix[i][j] = matrix[j][i] = i / 100.0;
    }
  }
  return matrix;
}

real* getVector(int n) {
  real* array = new real[n];

  for (int i = 0; i < n; ++i)
    array[i] = i / 100.0;

  return array;
}

void showVector(real* vector, int n) {
  for (int j = 0; j < n; j++) {
    std::cout << vector[j] << "\t";
  }
}

real* multiply(real** M, int n, int m, real* V, int k) {
  real* resultVector = new real[n];
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

//// MPI_Reduce
//int main(int argc, char** argv)
//{
//  int n = 1000;
// /* std::stringstream convert(argv[1]);
//
//  int n;
//  if (!(convert >> n))
//    n = 0;*/
//  double start, stop;
//
//  real* A, * b, * A_;
//  real* result = (real*)malloc(n * sizeof(real));
//  int IDp, np;
//  double root = 0.0;
//
//  MPI_Init(&argc, &argv);
//  MPI_Comm_size(MPI_COMM_WORLD, &np);
//  MPI_Comm_rank(MPI_COMM_WORLD, &IDp);
//
//  A = getSymmetricMatrix(n);
//  b = getVector(n);
//  std::cout << "awdawd" << std::endl;
//
//  A_ = (real*)malloc((n * n / np) * sizeof(real));
//
//  std::cout << "awdawd" << std::endl;
//  MPI_Scatter(A, n * n / np, MPI_FLOAT, A_, n * n / np, MPI_FLOAT, 0, MPI_COMM_WORLD);
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  std::cout << "awdawd2" << std::endl;
//  start = MPI_Wtime();
//  real* x = multiply(A, n, n / np, b, IDp * n / np);
//  stop = MPI_Wtime();
//  MPI_Reduce(x, result, n, MPI_FLOAT, MPI_SUM, root, MPI_COMM_WORLD);
//
//  std::cout << "awdawd3" << std::endl;
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  if (IDp == 0)
//  {
//    double t = stop - start;
//    std::cout << t << std::endl;
//    double t0 = 0.0025865;
//    double Sp = t0 / t;
//    double Ep = Sp * 100.0 / np;
//    std::cout << Sp << Ep << std::endl;
//  }
//
//  std::cout << "awdawd4" << std::endl;
//  free(A);
//  free(A_);
//  free(x);
//  free(b);
//  free(result);
//
//  MPI_Finalize();
//  return 0;
//}

// Последовательно.
int main(int argc, char** argv)
{
  std::stringstream convert(argv[1]);

  int n;
  if (!(convert >> n))
    n = 0;
  double start, stop;

  real **A = getSymmetricMatrix(n);
  real *b = getVector(n);
  
  MPI_Init(&argc, &argv);

  start = MPI_Wtime();
  real* x = multiply(A, n, n, b, 0);
  stop = MPI_Wtime();

  std::cout << (stop - start) << std::endl;

  MPI_Finalize();

  //delete[] A[0];
  //delete[] A;
  //delete b;
  //delete x;

  return 0;
}
