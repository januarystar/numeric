#include <stdlib.h>

#include "ldlt.h"
#include "simple.h"
#include "matrix.h"


void LDLTDecomposition(double **A, double ***L, double ***D, double ***LT, int N)
{
  int i, k;

  /* Variables for LDLT */
  double
    *f = malloc(sizeof(double) * N),
    *x = malloc(sizeof(double) * N);
  double beta;

  *L = malloc(sizeof(double *) * N);
  for (i = 0; i < N; i++)
    (*L)[i] = malloc(sizeof(double) * N);

  *D = malloc(sizeof(double *) * N);
  for (i = 0; i < N; i++)
    (*D)[i] = malloc(sizeof(double) * N);

  *LT = malloc(sizeof(double *) * N);
  for (i = 0; i < N; i++)
    (*LT)[i] = malloc(sizeof(double) * N);

  /* A1 = L1 * D1 * R1 */

  (*L)[0][0] = 1;
  (*D)[0][0] = A[0][0];
  (*LT)[0][0] = 1;

  for (k = 2; k <= N; k++)
  {
    /* Fill f_k-1 */
    for (i = 0; i < k - 1; i++)
      f[i] = A[i][k - 1];

    /*** Set L and LT ***/
    /* Solve X */
    Left(MatrMulMatr(Transposing(*LT, k - 1, k - 1), \
      Transposing(*D, k - 1, k - 1), \
      k - 1, k - 1, k - 1),
      f, x, k - 1);

    for (i = 0; i < k - 1; i++)
    {
      (*L)[k - 1][i] = x[i];
      (*L)[i][k - 1] = 0;
      (*LT)[i][k - 1] = (*L)[k - 1][i];
      (*LT)[k - 1][i] = 0;
    }
    (*L)[k - 1][k - 1] = 1;
    (*LT)[k - 1][k - 1] = 1;


    /*** Set Dk ***/
    /* Solve beta */
    beta = A[k - 1][k - 1] - DotProduct(VecMulMatr(x, *D, k - 1, k - 1), x, k - 1);
    for (i = 0; i < k - 1; i++)
    {
      (*D)[i][k - 1] = 0;
      (*D)[k - 1][i] = 0;
    }
    (*D)[k - 1][k - 1] = beta;
  }

  free(f);
  free(x);
} 

void LDLTSolve(double **A, double *b, double *x, int N)
{
  double **L, **D, **LT;
  double
    *y = malloc(sizeof(double) * N),
    *z = malloc(sizeof(double) * N);

  LDLTDecomposition(A, &L, &D, &LT, N);

  /* Ly = b */
  Left(L, b, y, N);
  /* Dz = y */
  Diagonal(D, y, z, N);
  /* LTx = z */
  Right(LT, z, x, N);

  free(y);
  free(z);
} /* End of 'LDRSolve' function */