//
// Created by dezodemius on 15.04.2020.
//
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

/**
 * Сгенерировать случайный массив float.
 * @param N - Размер массива.
 * @param min - Минимальное число в массиве.
 * @param max - Максимальное число в массиве.
 * @return Массив, заполненный случайными данными.
 */
float* gen_float_array(int N, const int min, const int max){
  int i;
  srand(time(NULL));
  float *array = (float*)malloc(N*sizeof(float));
  
  for(i = 0; i < N; i++)
    array[i] = rand() / (float) RAND_MAX   * (max - min) + min;
  
  return array;
}

/**
 * Сгенерировать случайный массив double.
 * @param N - Размер массива.
 * @param min - Минимальное число в массиве.
 * @param max - Максимальное число в массиве.
 * @return Массив, заполненный случайными данными.
 */
double * gen_double_array(int N, const int min, const int max){
  srand(time(NULL));
  double *array = (double*)malloc(N*sizeof(double));
  
  int i;
  for(i = 0; i < N; i++)
    array[i] = rand() / (double) RAND_MAX   * (max - min) + min;
  
  return array;
}

/**
 * Обычное скалярное произведение векторов для float.
 * @param a - Первый вектор.
 * @param b - Второй вектор.
 * @param N - Количество элементов в векторах.
 * @return Результат скалярного произведения.
 */
float dot_float(float *a, float *b, int N){
  int i;
  float sum = 0.0;
  for (i = 0; i < N; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

/**
 * Обычное скалярное произведение векторов для double.
 * @param a - Первый вектор.
 * @param b - Второй вектор.
 * @param N - Количество элементов в векторах.
 * @return Результат скалярного произведения.
 */
double dot_double(double *a, double *b, int N){
  int i;
  double sum = 0.0;
  for (i = 0; i < N; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

/**
 * Скалярное произведение векторов для float c развёртыванием цикла с шагом 2.
 * @param a - Первый вектор.
 * @param b - Второй вектор.
 * @param N - Количество элементов в векторах.
 * @return Результат скалярного произведения.
 */
float dot_float_deploy2(float *a, float *b, int N){
  int i;
  float sum = 0.0;
  for (i = 0; i < N; i+=2) {
    sum += (a[i] * b[i]) + (a[i + 1] * b[i + 1]);
  }
  return sum;
}

/**
 * Скалярное произведение векторов для float c развёртыванием цикла с шагом 4.
 * @param a - Первый вектор.
 * @param b - Второй вектор.
 * @param N - Количество элементов в векторах.
 * @return Результат скалярного произведения.
 */
float dot_float_deploy4(float *a, float *b, int N){
  int i;
  float sum = 0.0;
  for (i = 0; i < N; i+=4) {
    sum += (a[i] * b[i]) + (a[i + 1] * b[i + 1]) +
        (a[i + 2] * b[i + 2]) + (a[i + 3] * b[i + 3]);
  }
  return sum;
}

/**
 * Скалярное произведение векторов для double c развёртыванием цикла с шагом 2.
 * @param a - Первый вектор.
 * @param b - Второй вектор.
 * @param N - Количество элементов в векторах.
 * @return Результат скалярного произведения.
 */
double dot_double_deploy2(double *a, double *b, int N){
  int i;
  double sum = 0.0;
  for (i = 0; i < N; i+=2) {
    sum += (a[i] * b[i]) + (a[i + 1] * b[i + 1]);
  }
  return sum;
}

/**
 * Скалярное произведение векторов для double c развёртыванием цикла с шагом 4.
 * @param a - Первый вектор.
 * @param b - Второй вектор.
 * @param N - Количество элементов в векторах.
 * @return Результат скалярного произведения.
 */
double dot_double_deploy4(double *a, double *b, int N){
  int i;
  double sum = 0.0;
  for (i = 0; i < N; i+=4) {
    sum += (a[i] * b[i]) + (a[i + 1] * b[i + 1]) +
           (a[i + 2] * b[i + 2]) + (a[i + 3] * b[i + 3]);
  }
  
  return sum;
}

/**
 * Засечь время выполнения обычного скалярного произведения float.
 * @param N - Объём тестируемых данных.
 * @return Время выполнения.
 */
double float_scalar_prod_pinpoint_runtime(float *a, float *b, int N){
  double t0 = clock();
  float result = dot_float(a, b, N);
  double t1 = clock();
  
  return (t1 - t0) / (double)CLOCKS_PER_SEC;
}

/**
 * Засечь время выполнения обычного скалярного произведения double.
 * @param N - Объём тестируемых данных.
 * @return Время выполнения.
 */
double double_scalar_prod_pinpoint_runtime(double *a, double *b, int N){
  double t0 = clock();
  double result = dot_double(a, b, N);
  double t1 = clock();
  
  return (t1 - t0) / (double)CLOCKS_PER_SEC;
}

/**
 * Засечь время выполнения скалярного произведения с развёртыванием цикла с шагом 2 double.
 * @param N - Объём тестируемых данных.
 * @return Время выполнения.
 */
double deploy2_double_scalar_prod_pinpoint_runtime(double *a, double *b, int N){
  double t0 = clock();
  double result = dot_double_deploy2(a, b, N);
  double t1 = clock();
  
  return (t1 - t0) / (double)CLOCKS_PER_SEC;
}

/**
 * Засечь время выполнения скалярного произведения с развёртыванием цикла с шагом 4 double.
 * @param N - Объём тестируемых данных.
 * @return Время выполнения.
 */
double deploy4_double_scalar_prod_pinpoint_runtime(double *a, double *b, int N){
  double t0 = clock();
  double result = dot_double_deploy4(a, b, N);
  double t1 = clock();
  
  return (t1 - t0) / (double)CLOCKS_PER_SEC;
}

/**
 * Засечь время выполнения скалярного произведения с развёртыванием цикла с шагом 2 float.
 * @param N - Объём тестируемых данных.
 * @return Время выполнения.
 */
double deploy2_float_scalar_prod_pinpoint_runtime(float *a, float *b, int N){
  double t0 = clock();
  float result = dot_float_deploy2(a, b, N);
  double t1 = clock();
  
  return (t1 - t0) / (double)CLOCKS_PER_SEC;
}

/**
 * Засечь время выполнения скалярного произведения с развёртыванием цикла с шагом 4 float.
 * @param N - Объём тестируемых данных.
 * @return Время выполнения.
 */
double deploy4_float_scalar_prod_pinpoint_runtime(float *a, float *b, int N){
  double t0 = clock();
  float result = dot_float_deploy4(a, b, N);
  double t1 = clock();
  
  return (t1 - t0) / (double)CLOCKS_PER_SEC;
}

/**
 * Получить имя оптимизации
 * @param argv - Аргументы командной строки.
 * @return Имя текущей оптимизации.
 */
char *get_optimization_name(char **argv){
  char* optimization;
  
  if (strcmp(argv[1], "-O0") == 0)
    optimization = "Без оптимизации";
  if (strcmp(argv[1], "-O1") == 0)
    optimization = "Оптимизация -O1";
  if (strcmp(argv[1], "-O2") == 0)
    optimization = "Оптимизация -O2";
  if (strcmp(argv[1], "-O3") == 0)
    optimization = "Оптимизация -O3";
  
  return optimization;
}

void main(int argc, char** argv){
  if (argc <= 1)
    return;
  char *optimization = get_optimization_name(argv);
  
  int i, N1 = 1e3, N2 = 1e6, N3 = 1e8, min = -10, max = 10;
  double f_n1 = 0.0, f_n2 = 0.0, f_n3 = 0.0,
    d_n1 = 0.0, d_n2 = 0.0, d_n3 = 0.0,
    f_deploy2_n1 = 0.0, f_deploy2_n2 = 0.0, f_deploy2_n3 = 0.0,
    f_deploy4_n1 = 0.0, f_deploy4_n2 = 0.0, f_deploy4_n3 = 0.0,
    d_deploy2_n1 = 0.0, d_deploy2_n2 = 0.0, d_deploy2_n3 = 0.0,
    d_deploy4_n1 = 0.0, d_deploy4_n2 = 0.0, d_deploy4_n3 = 0.0,
    N = 6.0;
  
  double *a_n1 = gen_double_array(N1, min, max);
  double *a_n2 = gen_double_array(N2, min, max);
  double *a_n3 = gen_double_array(N3, min, max);
  
  float *b_n1 = gen_float_array(N1, min, max);
  float *b_n2 = gen_float_array(N2, min, max);
  float *b_n3 = gen_float_array(N3, min, max);
  
  for (i = 0; i < N; ++i) {
    f_n1 += float_scalar_prod_pinpoint_runtime(b_n1, b_n1, N1);
    f_n2 += float_scalar_prod_pinpoint_runtime(b_n2, b_n2, N2);
    f_n3 += float_scalar_prod_pinpoint_runtime(b_n3, b_n3, N3);
  
    d_n1 += double_scalar_prod_pinpoint_runtime(a_n1, a_n1, N1);
    d_n2 += double_scalar_prod_pinpoint_runtime(a_n2, a_n2, N2);
    d_n3 += double_scalar_prod_pinpoint_runtime(a_n3, a_n3, N3);
  
    f_deploy2_n1 += deploy2_float_scalar_prod_pinpoint_runtime(b_n1, b_n1, N1);
    f_deploy2_n2 += deploy2_float_scalar_prod_pinpoint_runtime(b_n2, b_n2, N2);
    f_deploy2_n3 += deploy2_float_scalar_prod_pinpoint_runtime(b_n3, b_n3, N3);
  
    f_deploy4_n1 += deploy4_float_scalar_prod_pinpoint_runtime(b_n1, b_n1, N1);
    f_deploy4_n2 += deploy4_float_scalar_prod_pinpoint_runtime(b_n2, b_n2, N2);
    f_deploy4_n3 += deploy4_float_scalar_prod_pinpoint_runtime(b_n3, b_n3, N3);
  
    d_deploy2_n1 += deploy2_double_scalar_prod_pinpoint_runtime(a_n1, a_n1, N1);
    d_deploy2_n2 += deploy2_double_scalar_prod_pinpoint_runtime(a_n2, a_n2, N2);
    d_deploy2_n3 += deploy2_double_scalar_prod_pinpoint_runtime(a_n3, a_n3, N3);
  
    d_deploy4_n1 += deploy4_double_scalar_prod_pinpoint_runtime(a_n1, a_n1, N1);
    d_deploy4_n2 += deploy4_double_scalar_prod_pinpoint_runtime(a_n2, a_n2, N2);
    d_deploy4_n3 += deploy4_double_scalar_prod_pinpoint_runtime(a_n3, a_n3, N3);
  }
  
  FILE *file;
  if((file= fopen("tables.txt", "a"))==NULL)
  {
    perror("Error occured while opening file");
    return;
  }
  
  fprintf(file, "\\begin{table}[H]\n"
         "  \\caption{Скалярное произведение. %s}\n"
         "  \\centering\n"
         "  \\begin{tabular}{|c|*7{c|}}\n"
         "  \\hline\n"
         "  \\textnumero & Вариант &\\multicolumn{6}{|c|}{ Время, c}\\\\ \\cline{3-8}\n"
         "  & оптимизации &\\multicolumn{3}{|c|}{float}&\\multicolumn{3}{|c|}{double}\n"
         "  \\\\ \\cline{3-8}\n"
         "  &  & $N=10^3$ & $N=10^6$& $N=10^8$& $N=10^3$ &$N=10^6$ & $N=10^8$\\\\\n"
         "  \\hline 1. & Без оптимизации & %lf & %lf  & %lf  & %lf  & %lf & %lf\\\\ \\hline\n"
         "2. & Цикл с шагом 2 & %lf & %lf  & %lf  & %lf  & %lf & %lf\\\\ \\hline\n"
         "3. & Цикл с шагом 4 & %lf & %lf  & %lf  & %lf  & %lf & %lf\\\\ \\hline\n"
         "\t\\end{tabular}\n"
         "\\end{table}\n\n",
      optimization,
      f_n1 / N, f_n2 / N, f_n3 / N, d_n1 / N, d_n2 / N, d_n3 / N,
      f_deploy2_n1 / N, f_deploy2_n2 / N, f_deploy2_n3 / N, d_deploy2_n1 / N, d_deploy2_n2 / N, d_deploy2_n3 / N,
      f_deploy4_n1 / N, f_deploy4_n2 / N, f_deploy4_n3 / N, d_deploy4_n1 / N, d_deploy4_n2 / N, d_deploy4_n3 / N);

  fclose(file);
  free(a_n1);
  free(a_n2);
  free(a_n3);
  free(b_n1);
  free(b_n2);
  free(b_n3);
  return;
}
