/** @file libtwocop.h twoCopula test library headers.
 * This file contains the main libraries to compute the two copula equality
 * test.
 */


 /*! \mainpage Two Copula Equality tester
  * \section intro_sec Introduction
  * This document is describes the functions used to compute the P-values
  * for the two copula tester.
  * 
  * 
  * Reference: Remillard, B. and Scaillet, O. (2006). Two copula test. 
  * Working paper
  * \section com_use_sec Compiling and Usage
  * I've simplified the coding in this version to only have one library file,
  * it is still below 1000 lines which is reasonable for a source file.
  * to compile the library no extra library calls are needed since only
  * the tester uses GSL and the GSL libraries are included with the 
  * utils functions. I have to filter out functions used in the libtwocop
  * library that are in utils in order to make the library standalone
  * 
  * \section gsl_tester_sec Using the testing file with GSL.
  * The twoCop.c file is a tester program, for now I just apply the test
  * two to uniform matrices with independant colums, of course one can use
  * some copula simulators to implement some numerical applications
  */
#ifndef LIBTWOCOP_H_
#define LIBTWOCOP_H_

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/** Just a little max MACRO. */
#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

/** My min Macro. */
#ifndef min
#define min(a,b) ((a) > (b) ? (b) : (a))
#endif

/** */
#ifndef XTOIJ
#define XTOIJ(a,b) ((a*(a-1))/2 + b)
#endif


#define DEBUG 0									 
/** Function to get the bandwidth.
 * @param n the size of the sample
 */
   double geth(int n);
   void twoCopTest(double **U, double **V, int n1, int n2, int d,
                    double **xi_c, double **eta_c, int N, double *S, double *Shat, double *pval);  
   void makeMmats(double **Muu, double **Mvv, double **Muv, double **Mvu,
                 double **U, double **V, int n1, int n2, int d);
   double calcS(double **Muu, double **Mvv, double **Muv, int n1, int n2);
   void makeG1(double ***G1, double **U, int n, int d, double h);
/** Preparation of the Matrices for U.
 * @param H1uu
 * @param H2uu
 * @param H3uu
 * @param Muu
 * @param U
 * @param n
 * @param d
 * @param h
 */
   void makeHuu(double ***H1uu, double ***H2uu, double ***H3uu, double **Muu,
               double **U, int n, int d, double h);
/** Preparation of the Matrices for U and V.
 * @param H1uv
 * @param H1vu
 * @param H2uv
 * @param H3uv
 * @param Muv
 * @param U
 * @param V
 * @param n1
 * @param n2
 * @param d
 * @param h1
 * @param h2
 */
   void makeHuv(double ***H1uv, double ***H1vu, double ***H2uv, double ****H3uv,
               double **Muv, double **U, double **V, int n1, int n2, int d,
               double h1, double h2);

/** Preparation of the Denominator Matrix for U.
 * @param Duu
 * @param U
 * @param n
 * @param d
 */
   void makeDuu(double ***Duu, double **U, int n, int d);
/** Preparation of the Denominator Matrix for UV.
 * @param Duv
 * @param U
 * @param V
 * @param n1
 * @param n2
 * @param d
 */
   void makeDuv(double ***Duv, double **U, double **V, int n1, int n2, int d);

   double calcC(double ***H1uu, double ***H2uu, double ***H3uu, double **Muu,
               double *xi, int n, int d, double h);


   double calcCD(double ***H1uv, double ***H1vu, double ***H2uv, double ****H3uv,
                double **Muv, double *xi, double *eta, int n1, int n2, int d,
                double h1, double h2);
#endif /*LIBTWOCOP_H_ */
