#define _CRT_SECURE_NO_WARNINGS
#define _CRTDBG_MAP_ALLOC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <crtdbg.h>
#include "ldlt.h"
#include "matrix.h"


int LoadFromFile(FILE *F, double ***A, double **b, int *N)
{
  int i, j;

  if (fscanf(F, "%i", N) < 1)
    return 0;

  *A = malloc(*N * sizeof(double *));

  for (i = 0; i < *N; i++)
    (*A)[i] = malloc(sizeof(double) * *N);

  for (i = 0; i < *N; i++)
    for (j = 0; j < *N; j++)
      fscanf(F, "%lf", &((*A)[i][j]));

  *b = malloc(sizeof(double) * *N);

  for (i = 0; i < *N; i++)
    fscanf(F, "%lf", &((*b)[i]));

  return 1;
}

int main(void)
{
  _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
  double **A = NULL, *b = NULL;
  double *b0 = NULL;
  double *x0 = NULL;
  double *deltax = NULL;
  double *deltab = NULL;
  double *x = NULL;
  int N, i;

  N = 13;

  b0 = malloc(N * sizeof(double));
  x0 = malloc(N * sizeof(double));
  deltax = malloc(N * sizeof(double));
  deltab = malloc(N * sizeof(double));
  x = malloc(N * sizeof(double));

  FILE *F = fopen("solve1.out", "w");
  FILE *Fin = fopen("m1.in", "r");
  FILE* Fdeltax = fopen("deltax1.txt", "w");
  FILE* Fepsx = fopen("epsx.txt", "w");
  FILE *Fdeltab = fopen("deltab1.txt", "w");

  for (i = 0; i < N; i++)
    b0[i] = 1;

  //b0[0] = 2;
  //b0[1] = 2.0002;

  if (F == NULL || Fin == NULL || Fdeltax == NULL || Fdeltab == NULL || Fepsx == NULL)
    return;

  while (LoadFromFile(Fin, &A, &b, &N))
  {
    //x = malloc(sizeof(double) * N); 

    LDLTSolve(A, b, x, N);
    LDLTSolve(A, b0, x0, N);

	for (i = 0; i < N; i++)
	  deltab[i] = b[i] - b0[i];

    if (IsEqual(MatrMulVec(A, x, N, N), b, N) == 0)
      printf("OK\n");
    else
      printf("BAD: %i\n", IsEqual(MatrMulVec(A, x, N, N), b, N));



    for (i = 0; i < N; i++)
    {
      deltax[i] = (x[i] - x0[i]);
      fprintf(F, "%.16f     %.16f     %.16f\n", x[i], x0[i], deltax[i]);
      
    }
    fprintf(Fdeltax, "%.15f\n", SqrNorm(deltax, N));
    fprintf(Fepsx, "%.15f\n", SqrNorm(deltax, N) / SqrNorm(x0, N));
	//fprintf(Fdeltax, "%.15f AAA %.15f AAA %.15f\n", SqrNorm(deltab, N), SqrNorm(b, N), SqrNorm(deltab, N) / SqrNorm(b, N));
    fprintf(F, "\n\n");
  }

  fclose(F);
  fclose(Fin);


  fprintf(Fdeltab, "%.10f %.10f %.10f", SqrNorm(deltab, N), SqrNorm(b0, N), (SqrNorm(deltab, N)) / SqrNorm(b0, N));

  return 0;
} /* End of 'main' function */