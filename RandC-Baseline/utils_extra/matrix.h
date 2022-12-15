#ifndef MATRIX_H

#include "point.h"

#define	MATRIX_H

void linearSolve(int n, double A[n][n], double *b, double *c);

void gaussElimination(double *A, double *B, double *X, int n);

double* matToVector(int n, double mat[n][n]);

double* choleskyDecomp(double *A, int n, double *u);

double* transpose(double *A, int n);

double* forwardSubstitution(double *L, int n);

double* matProduct(double *A, double *B, int n);

void matxvec(double *C, double *A, double *b, int n);

void invertTriangular(double *A, double *x, int n);

void LUDecomposition(double *A, double *L, double *U, int n);

void LUInversion(double *L, double *U, int n, double *Li, double *Ui);

void covarianceMat(Lpoint **points, int numPoints, double convarianceMat[3][3]);

void covarianceMatInt(Lpoint **points, int numPoints, double convarianceMat[3][3]);

void zeroTriMat(int n, double mat[n][n]);

void sumMats(int n, double mat1[n][n], double mat2[n][n]);

#endif	/* MATRIX_H */