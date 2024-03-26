/* libtwocop.c Main library file.
 * This file contains the functions to calculate the test, we separate the
 * functions into two libraries, the first is libtwocopU.h which handles all
 * the computations where only one of the U or V matrix is needed.
 * and the second is libtwocopUV.h where the computations that use U and V
 * are needed, we note here that some functions are of the form f(U,V) while
 * other expect f(V,U), since it is completely symmetrics, we only need to 
 * implement the f(U,V) case and call the function with U = V and V = U.
 */
#define  Pi 3.14159265358979323846264338
#include <R.h>
#include "libtwocop.h"
#include <Rmath.h>
#define DEBUG1 0


/* Used from utils.c */ 
double mean(double *x, int n)
{
  int i; 
  double sum = 0.0;
  
  for(i=0;i<n;i++)
    sum += x[i];
  
  return sum/((double) n);
}

/* Used from rang.c */

void rang2(double *x, int n, double *y)
  
{
  int i, j;
  double  rr;
  
  
  
  
  for(i=0;i<n;i++)
  {
    rr = 1.0;
    for(j=0;j<n;j++)
    {
      if(x[j] < x[i])
        rr += 1.0;
    }
    y[i] = rr/((double) (1.0+n));
  }
}


/* Original litwocop.c */

double geth(int n)
{
  /* return 3.5 * sqrt(1. / 12) * exp(-.33333333333333333 * log(n)); */
  return 1. / sqrt(n);
}



void pvalue(double **X, double **Y, int n, int m, int d, int Nsim, int paired, double *S, double *Shat, double *pval)
{
  
  double **U = NULL;
  double **V = NULL;
  double **U1 = NULL;
  double **V1 = NULL;
  double **xi = NULL;
  double **eta = NULL;
  double xibar, etabar;
  int i, j;
  
  /* Initialize U and V matrices */
  /* U is d x n */
  U = malloc(d * sizeof(double *));
  for (i = 0; i < d; i++)
    U[i] = malloc(n * sizeof(double));
  
  /* V is d x m */
  V = malloc(d * sizeof(double *));
  for (i = 0; i < d; i++)
    V[i] = malloc(m * sizeof(double));
  
  /* U1 is d x m */		
  U1 = malloc(d * sizeof(double *));
  for (i = 0; i < d; i++)
    U1[i] = malloc(n * sizeof(double));
  
  /* V1 is d x m */
  V1 = malloc(d * sizeof(double *));
  for (i = 0; i < d; i++)
    V1[i] = malloc(m * sizeof(double));
  
  
  
  
  for (i = 0; i < d; i++)
  {
    for (j = 0; j < n; j++)
      U1[i][j] = X[j][i];
    
    rang2(U1[i], n, U[i]);
    
    for (j = 0; j < m; j++)
      V1[i][j] = Y[j][i];
    
    rang2(V1[i], m, V[i]);
  }
  
  
  
  /* Generate xi and eta vectors */
  
  
  
  /* xi is Nsim x n */
  xi = malloc(Nsim * sizeof(double *));
  for (i = 0; i < Nsim; i++)
    xi[i] = malloc(n * sizeof(double));
  
  /* eta is Nsim x m */
  eta = malloc(Nsim * sizeof(double *));
  for (i = 0; i < Nsim; i++)
    eta[i] = malloc(m * sizeof(double));
  
  /*---------------------------------------*/
  
  GetRNGstate(); 
  
  if(paired==0)
  {
    
    
    for (i = 0; i < Nsim; i++) {
      for (j = 0; j < n; j++)
        xi[i][j] = rnorm(0.,1.); 
      
      xibar = mean(xi[i],n);    
      
      
      for (j = 0; j < m; j++)
        eta[i][j] = rnorm(0.,1.); 
      
      
      etabar = mean(eta[i],m);  
      
      /* Center the vectors as we always use the centered version */
      for (j = 0; j < n; j++)
        xi[i][j] -= xibar;
      
      for (j = 0; j < m; j++)
        eta[i][j] -= etabar;
    }
  }
  else
  {
    for (i = 0; i < Nsim; i++) {
      for (j = 0; j < n; j++)
        xi[i][j] = rnorm(0.,1.); 
      
      xibar = mean(xi[i],n);    
      
      
      /* Center the vectors as we always use the centered version */
      for (j = 0; j < n; j++)
        xi[i][j] -= xibar;
      
      for (j = 0; j < m; j++)
        eta[i][j] = xi[i][j];
    }
  }
  /*---------------------------------------*/
  
  PutRNGstate();
  
  twoCopTest(U, V, n, m, d, xi, eta, Nsim, S, Shat, pval);
  
  
  
  
  for (i = 0; i < d; i++)
  { 
    free(V1[i]); free(U1[i]);  free(V[i]);  free(U[i]); 
  }
  free(V1); free(U1); free(V); free(U);
  
  for (i = 0; i < Nsim; i++) 
  {
    free(xi[i]); free(eta[i]);
  }
  free(xi); free(eta);
  
  
}









void twoCopTest(double **U, double **V, int n1, int n2, int d,
                double **xi_c, double **eta_c, int N, double *S, double *Shat, double *pval)
{
  int i, l,m;
  int counter;
  double h1, h2;
  double C, D, CD;
  double **Muu, **Mvv, **Muv, **Mvu;
  double ***H1uu, ***H1vv, ***H1uv, ***H1vu;
  double ***H2uu, ***H2vv, ***H2uv;
  double ***H3uu, ***H3vv;
  double ****H3uv;
  h1 = geth(n1);
  h2 = geth(n2);
  m = 0;
  counter = 0;
  
  /*========================================================================*
   *                     PREPARATION OF THE MATRICES                        *
   *========================================================================*/
  Muu = malloc(n1 * sizeof(double *));
  Muv = malloc(n1 * sizeof(double *));
  for (i = 0; i < n1; ++i) {
    Muu[i] = malloc(n1 * sizeof(double));
    Muv[i] = malloc(n2 * sizeof(double));
  }
  Mvv = malloc(n2 * sizeof(double *));
  Mvu = malloc(n2 * sizeof(double *));
  for (i = 0; i < n2; ++i) {
    Mvv[i] = malloc(n2 * sizeof(double));
    Mvu[i] = malloc(n1 * sizeof(double));
  }
  
  /* Initialization of the H matrices */
  H1uu = malloc(d * sizeof(double **));
  H1vv = malloc(d * sizeof(double **));
  H1uv = malloc(d * sizeof(double **));
  H1vu = malloc(d * sizeof(double **));
  
  H2uu = malloc(d * sizeof(double **));
  H2vv = malloc(d * sizeof(double **));
  H2uv = malloc(d * sizeof(double **));
  
  H3uu = malloc((d * (d - 1) / 2) * sizeof(double **));
  H3vv = malloc((d * (d - 1) / 2) * sizeof(double **));
  
  H3uv = malloc(d*sizeof(double***));
  for(l = 0;l<d;++l){
    H3uv[l] = malloc(d*sizeof(double**));
    for(m = 0;m<d;++m){
      if(l  == m) 
        continue;
      H3uv[l][m] = malloc(n1*sizeof(double*));
      for(i = 0;i<n1;i++)
        H3uv[l][m][i] = malloc(n2*sizeof(double));
    }
  }
  
  
  
  
  for (l = 0; l < d; ++l) {
    H1uu[l] = malloc(n1 * sizeof(double *));
    H1vv[l] = malloc(n2 * sizeof(double *));
    H1uv[l] = malloc(n1 * sizeof(double *));
    H1vu[l] = malloc(n2 * sizeof(double *));
    H2uu[l] = malloc(n1 * sizeof(double *));
    H2uv[l] = malloc(n1 * sizeof(double *));
    H2vv[l] = malloc(n2 * sizeof(double *));
    for (i = 0; i < n1; ++i) {
      H1uu[l][i] = malloc(n1 * sizeof(double));
      H2uu[l][i] = malloc(n1 * sizeof(double));
      H1uv[l][i] = malloc(n2 * sizeof(double));
      H2uv[l][i] = malloc(n2 * sizeof(double));
    }
    for (i = 0; i < n2; ++i) {
      H1vv[l][i] = malloc(n2 * sizeof(double));
      H2vv[l][i] = malloc(n2 * sizeof(double));
      H1vu[l][i] = malloc(n1 * sizeof(double));
    }
  }
  for (l = 0; l < (d * (d - 1) / 2); ++l) {
    H3uu[l] = malloc(n1 * sizeof(double *));
    H3vv[l] = malloc(n2 * sizeof(double *));
    for (i = 0; i < n1; ++i) 
      H3uu[l][i] = malloc(n1 * sizeof(double));
    
    for (i = 0; i < n2; ++i)
      H3vv[l][i] = malloc(n2 * sizeof(double));
  }
  
  /*========================================================================*
   *                  MATRICES AVAILABLE BEYOND THIS POINT:                 *
   *========================================================================*/
  makeMmats(Muu, Mvv, Muv, Mvu, U, V, n1, n2, d);
  
  
  S[0] = calcS(Muu, Mvv, Muv, n1, n2);
  
  
  makeHuu(H1uu, H2uu, H3uu, Muu, U, n1, d, h1);
  
  
  makeHuu(H1vv, H2vv, H3vv, Mvv, V, n2, d, h2);
  
  
  makeHuv(H1uv, H1vu, H2uv, H3uv, Muv, U, V, n1, n2, d, h1, h2);
  
  
  /* Begin the computation of the P-value indicator function. */
  for (i = 0; i < N; ++i) {
    C = calcC(H1uu, H2uu, H3uu, Muu, xi_c[i], n1, d, h1);
    D = calcC(H1vv, H2vv, H3vv, Mvv, eta_c[i], n2, d, h2);
    CD = calcCD(H1uv, H1vu, H2uv, H3uv, Muv,
                xi_c[i], eta_c[i], n1, n2, d, h1, h2);
    
    Shat[i] = (n2 * C - CD + n1 * D) / (double) (n1 + n2);
    
    if (Shat[i] > S[0])
      counter++;
  }
  
  /*========================================================================*
   *                  Garbage Collection, we free all the matrices          *
   *========================================================================*/
  for (i = 0; i < n1; ++i) {
    free(Muu[i]);
    free(Muv[i]);
  }
  free(Muu);
  free(Muv);
  for (i = 0; i < n2; ++i) {
    free(Mvv[i]);
    free(Mvu[i]);
  }
  free(Mvv);
  free(Mvu);
  
  /* Free all 3D matrices */
  for (l = 0; l < d; ++l) {
    for (i = 0; i < n1; ++i) {
      free(H1uu[l][i]);
      free(H2uu[l][i]);
      free(H1uv[l][i]);
      free(H2uv[l][i]);
    }
    for (i = 0; i < n2; ++i) {
      free(H1vv[l][i]);
      free(H2vv[l][i]);
      free(H1vu[l][i]);
    }
    free(H1uu[l]);
    free(H1uv[l]);
    
    free(H2uu[l]);
    free(H2uv[l]);
    
    free(H1vv[l]);
    free(H1vu[l]);
    free(H2vv[l]);
  }
  
  for (l = 0; l < d * (d - 1) / 2; ++l) {
    for (i = 0; i < n1; ++i) {
      free(H3uu[l][i]);
    }
    for (i = 0; i < n2; ++i)
      free(H3vv[l][i]);
    
    free(H3uu[l]);
    free(H3vv[l]);
  }
  
  
  for (l = 0;l< d;++l){ 
    for (m = 0;m < d;++m){
      if(l == m) 
        continue;
      for(i = 0;i < n1;++i) 
        free(H3uv[l][m][i]); 
      free(H3uv[l][m]); 
    } 
    free(H3uv[l]); 
  } 
  
  free(H1uu);
  free(H2uu);
  free(H3uu);
  free(H1vv);
  free(H2vv);
  free(H3vv);
  free(H1uv);
  free(H1vu);
  free(H2uv);
  free(H3uv);
  /*========================================================================*
   *                  We may now return the statistic:                      *
   *========================================================================*/
  pval[0] = (double) counter / (double) N;
}


void makeHuu(double ***H1uu, double ***H2uu, double ***H3uu,
             double **Muu, double **U, int n, int d, double h)
{
  int i, j, k, p, l, m;
  double ***D;
  double ***G1;
  double sumH1, sumH2, sumH3;
  double a, b;
  D = malloc(d * sizeof(double **));
  for (i = 0; i < d; i++) {
    D[i] = malloc(n * sizeof(double *));
    for (j = 0; j < n; j++)
      D[i][j] = malloc(n * sizeof(double));
  }
  makeDuu(D, U, n, d);
  G1 = malloc(d * sizeof(double **));
  for (i = 0; i < d; i++) {
    G1[i] = malloc(n * sizeof(double *));
    for (j = 0; j < n; j++)
      G1[i][j] = malloc(n * sizeof(double));
  }
  makeG1(G1, U, n, d, h);
  for (l = 0; l < d; ++l) {
    
    /* Make H1uu and H2uu */
    /*----------------------------------------------------------------------*/
    for (i = 0; i < n; ++i) {
      for (j = 0; j < n; ++j) {
        sumH1 = 0.0;
        sumH2 = 0.0;
        for (k = 0; k < n; ++k) {
          /* Make H1 */
          sumH1 +=
          Muu[i][k] * D[l][i][k] * G1[l][((U[l][i] > U[l][j]) ? i : j)][k];
          for (p = 0; p < n; ++p) {
            /* Make H2 */
            a = min(1.0, h + min(U[l][k], U[l][p]));
            b = max(max(U[l][i], U[l][j]), max(U[l][k], U[l][p]) - h);
            if (a > b)
              sumH2 += Muu[k][p] * D[l][k][p] * (a - b);
            /* Make H3 */
          }
        }
        H1uu[l][i][j] = sumH1;
        H2uu[l][i][j] = sumH2;
      }
    }
    /*----------------------------------------------------------------------*/
    
    /* Make the H3uu matrix */
    /*----------------------------------------------------------------------*/
    for (m = 0; m < l; ++m) {
      for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
          sumH3 = 0.0;
          for (k = 0; k < n; ++k) {
            for (p = 0; p < n; ++p) {
              sumH3 +=
                Muu[k][p] * D[l][k][p] * D[m][k][p] *
                G1[l][(U[l][i] > U[l][p] ? i : p)][k] *
                G1[m][(U[m][j] > U[m][k] ? j : k)][p];
            }
          }
          H3uu[XTOIJ(l, m)][i][j] = sumH3;
          
        }
      }
    }
  }
  /* We need to Free, G1 and D */
  for (l = 0; l < d; l++) {
    for (i = 0; i < n; i++)
      free(G1[l][i]);
    free(G1[l]);
  }
  free(G1);
  /* We free D */
  for (l = 0; l < d; l++) {
    for (i = 0; i < n; i++)
      free(D[l][i]);
    free(D[l]);
  }
  free(D);
}

void makeHuv(double ***H1uv, double ***H1vu, double ***H2uv, double ****H3uv,
             double **Muv, double **U, double **V, int n1, int n2, int d,
             double h1, double h2)
{
  int i, j, k, p, l, m;
  double ***D;
  double sumH1, sumH2, sumH3;
  double a, b, c, e;
  double u;
  D = malloc(d * sizeof(double **));
  for (l = 0; l < d; ++l) {
    D[l] = malloc(n1 * sizeof(double *));
    for (i = 0; i < n1; ++i) {
      D[l][i] = malloc(n2 * sizeof(double));
    }
  }
  
  makeDuv(D, U, V, n1, n2, d);
  
  for (l = 0; l < d; ++l) {
    
    for (i = 0; i < n1; ++i) {
      for (j = 0; j < n2; ++j) {
        sumH1 = 0;
        u = max(U[l][i], V[l][j]);
        /* Do H1 uv */
        for (k = 0; k < n2; ++k) {
          /* Compute G1 */
          a = min(1.0, V[l][k] + h2);
          b = max(u, V[l][k] - h2);
          if (a > b)
            sumH1 += Muv[i][k] * D[l][i][k] * (a - b);
          /* Do H2uv */
        }
        H1uv[l][i][j] = sumH1;
        
        sumH2 = 0;
        for (k = 0; k < n1; ++k) {
          for (p = 0; p < n2; ++p) {
            a = min(1.0, min(U[l][k] + h1, V[l][p] + h2));
            b = max(u, max(U[l][k] - h1, V[l][p] - h2));
            if (a > b)
              sumH2 += Muv[k][p] * D[l][k][p] * (a - b);
          }
          
        }
        H2uv[l][i][j] = sumH2;
        
        sumH1 = 0;
        
        /* Do H1vu */
        for (k = 0; k < n1; ++k) {
          a = min(1.0, U[l][k] + h1);
          b = max(u, U[l][k] - h1);
          if (a > b)
            sumH1 += Muv[k][j] * D[l][k][j] * (a - b);
        }
        H1vu[l][j][i] = sumH1;
      }
    }
    
    /* Handle H3uv */
    for (m = 0; m < d; ++m) {
      if(l==m) 
        continue;
      for (i = 0; i < n1; ++i) {
        for (j = 0; j < n2; ++j) {
          sumH3 = 0;
          for (k = 0; k < n1; ++k) {
            for (p = 0; p < n2; ++p) {
              a = min(1.0, U[l][k] + h1);
              b = max(max(U[l][i], V[l][p]), U[l][k] - h1);
              c = min(1.0, V[m][p] + h2);
              e = max(max(V[m][j], U[m][k]), V[m][p] - h2);
              if ((a > b) && (c > e))
                sumH3 +=
                  Muv[k][p] * D[l][k][p] * D[m][k][p] * (a - b) * (c - e);
            }
          }
          H3uv[l][m][i][j] = sumH3;
        }
      }
    }
  }
  for (l = 0; l < d; l++) {
    for (i = 0; i < n1; i++)
      free(D[l][i]);
    free(D[l]);
  }
  free(D);
  
}

void makeMmats(double **Muu, double **Mvv, double **Muv,
               double **Mvu, double **U, double **V, int n1, int n2, int d)
{
  int i, j, s;
  double prod;
  for (i = 0; i < n1; ++i) {
    for (j = 0; j <= i; ++j) {
      prod = 1.0;
      for (s = 0; s < d; ++s) {
        prod *= (1.0 - max(U[s][i], U[s][j]));
      }
      Muu[i][j] = prod;
      Muu[j][i] = prod;
    }
  }
  
  for (i = 0; i < n2; ++i) {
    for (j = 0; j <= i; ++j) {
      prod = 1.0;
      for (s = 0; s < d; ++s) {
        prod *= (1.0 - max(V[s][i], V[s][j]));
      }
      Mvv[i][j] = prod;
      Mvv[j][i] = prod;
    }
  }
  
  for (i = 0; i < n1; ++i) {
    for (j = 0; j < n2; ++j) {
      prod = 1.0;
      for (s = 0; s < d; ++s) {
        prod *= (1.0 - max(U[s][i], V[s][j]));
      }
      Muv[i][j] = prod;
      Mvu[j][i] = prod;
    }
  }
}
double calcS(double **Muu, double **Mvv, double **Muv, int n1, int n2)
{
  int i, j;
  double sumUu = 0;
  double sumVv = 0;
  double sumUv = 0;
  for (i = 0; i < n1; ++i) {
    for (j = 0; j < n1; ++j)
      sumUu += Muu[i][j];
    
    for (j = 0; j < n2; ++j)
      sumUv += Muv[i][j];
  }
  
  for (i = 0; i < n2; ++i) {
    for (j = 0; j < n2; ++j)
      sumVv += Mvv[i][j];
  }
  return (sumUu / (n1 * n1) + sumVv / (n2 * n2) -
          2 * sumUv / (n1 * n2)) / (1. / n1 + 1. / n2);
}
void makeDuu(double ***Duu, double **U, int n, int d)
{
  int i, j, l;
  for (l = 0; l < d; l++) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
        Duu[l][i][j] = 1. / (1 - max(U[l][i], U[l][j]));
    }
  }
}

void makeDuv(double ***Duv, double **U, double **V, int n1, int n2, int d)
{
  int i, j, l;
  for (l = 0; l < d; ++l) {
    for (i = 0; i < n1; ++i) {
      for (j = 0; j < n2; ++j)
        Duv[l][i][j] = 1. / (1 - max(U[l][i], V[l][j]));
    }
  }
}

void makeG1(double ***G1, double **U, int n, int d, double h)
{
  int i, j, l;
  double a, b;
  for (l = 0; l < d; l++) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        a = min(1., U[l][j] + h);
        b = max(U[l][i], U[l][j] - h);
        if (a > b) {
          G1[l][i][j] = a - b;
        } 
        else {
          G1[l][i][j] = 0;
        }
      }
    }
  }
}


double calcC(double ***H1uu, double ***H2uu, double ***H3uu, double **Muu,
             double *xi, int n, int d, double h)
{
  int i, j, l, m;
  double sumA, sumI, sumJd, sumJod;
  double X;
  
  sumA = 0;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j)
      sumA += xi[i] * xi[j] * Muu[i][j];
  }
  sumI = 0;
  sumJd = 0;
  sumJod = 0;
  for (l = 0; l < d; ++l) {
    for (i = 0; i < n; ++i) {
      for (j = 0; j < n; ++j) {
        X = xi[i] * xi[j];
        sumI += X * H1uu[l][i][j];
        sumJd += X * H2uu[l][i][j];
      }
    }
    for (m = 0; m < l; ++m) {
      for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
          sumJod += xi[i] * xi[j] * H3uu[XTOIJ(l, m)][i][j];
      }
    }
  }
  return sumA / ((double) n) - (sumI / (n * n * h))
    + .25 * (sumJd + 2 * sumJod) / (n * n * n * h * h);
}
double calcCD(double ***H1uv, double ***H1vu, double ***H2uv,
              double ****H3uv, double **Muv, double *xi, double *eta, int n1,
              int n2, int d, double h1, double h2)
{
  int i, j, l, m;
  double sumAG, sumIvu, sumIuv, sumJd, sumJod;
  double XY;
  
  sumAG = 0;
  for (i = 0; i < n1; ++i) {
    for (j = 0; j < n2; ++j)
      sumAG += xi[i] * eta[j] * Muv[i][j];
  }
  sumIuv = sumIvu = sumJd = sumJod = 0;
  for (l = 0; l < d; ++l) {
    for (i = 0; i < n1; ++i) {
      for (j = 0; j < n2; ++j) {
        XY = xi[i] * eta[j];
        sumIuv += XY * H1uv[l][i][j];
        sumIvu += XY * H1vu[l][j][i];
        sumJd += XY * H2uv[l][i][j];
      }
    }
    
    for (m = 0; m < d; ++m) {
      if(l ==m) 
        continue;
      for (i = 0; i < n1; ++i) {
        for (j = 0; j < n2; ++j)
          sumJod += xi[i] * eta[j] * H3uv[l][m][i][j];
      }
    }
  }
  return 2 * sumAG - sumIuv / (n2 * h2) - sumIvu / (n1 * h1)
    + .5 * (sumJd + sumJod) / (n1 * n2 * h1 * h2);
}




void pvalue2(double *X, double *Y, int *n, int *m, int *d, int *Nsim, int *paired, double *S, double *Shat, double *pval)
{
  
  double **XX, **YY;
  
  int i,j,k;
  
  XX = (double **) malloc(n[0] * sizeof(double *));
  for (i = 0; i < n[0]; i++)
    XX[i] = (double *) malloc(d[0]* sizeof(double));
  
  YY = (double **) malloc(m[0] * sizeof(double *));
  for (i = 0; i < m[0]; i++)
    YY[i] = (double *) malloc(d[0]* sizeof(double));
  
  k =0;
  for(j=0;j<d[0];j++)
  {
    for (i = 0; i < n[0]; i++)
    {  
      XX[i][j] = X[k];
      k++;
      
    }
    
    
  }
  
  k =0;
  for(j=0;j<d[0];j++)
  {
    for (i = 0; i < m[0]; i++)
    {
      YY[i][j] = Y[k];
      k++;
      
    }
    
  }
  
  pvalue(XX, YY, n[0], m[0], d[0], Nsim[0], paired[0], S, Shat, pval);
  
  
  for(i=0;i<n[0];i++)
    free(XX[i]);
  
  for(i=0;i<m[0];i++)
    free(YY[i]);
  
  free(XX); free(YY);
}
