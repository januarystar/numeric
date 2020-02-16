#include <stdlib.h>
#include "simple.h"

void Diagonal(double *A[], double *b, double *x, int N)
{
  int k;

  for (k = 0; k < N; k++)
    x[k] = b[k] / A[k][k];

}

void Left(double **A, double *b, double *x, int N)
{
  int k;

  x[0] = b[0] / A[0][0];

  for (k = 1; k < N; k++)
  {
    int j;
    double p = 0;

    for (j = 0; j < k; j++)
      p += A[k][j] * x[j];

    x[k] = (b[k] - p) / A[k][k];
  }
}

void Right(double *A[], double *b, double *x, int N)
{
  int k;

  x[N - 1] = b[N - 1] / A[N - 1][N - 1];

  for (k = N - 2; k >= 0; k--)
  {
    int j;
    double p = 0;

    for (j = k + 1; j < N; j++)
      p += A[k][j] * x[j];

    x[k] = (b[k] - p) / A[k][k];
  }
}