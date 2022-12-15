#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "matrix.h"
#include "vector.h"

#define mat_elem(a, y, x, n) (a + ((y) * (n) + (x)))

void printMat(int n, double A[n][n])
{
    int i = 0, j = 0;
    
    for(i = 0; i < n; i++)
    {
       for(j = 0; j < n; j++)
       {
            printf("%lf\t", A[i][j]);
       }
          
       printf("\n");
    }
    
    printf("\n");
}

void printMat2(double *A, int cols)
{
    #define A(i, j) A[i * cols + j]
    
    int i = 0, j = 0;
    
    for(i = 0; i < cols; i++)
    {
       for(j = 0; j < cols; j++)
       {
            printf("%.4lf\t", A(i,j));
       }
          
       printf("\n");
    }
    
    printf("\n");
    
    #undef A
}


void sumMats(int n, double mat1[n][n], double mat2[n][n])
{
    int i = 0, j = 0;

    for(i = 0; i < n; i++)
    {
	for(j = 0; j < n; j++)
	{
	    mat2[i][j] += mat1[i][j];
	}
    }
}

void divideMat(int n, float mat[n][n], double denom)
{
    int i = 0, j = 0;

    for(i = 0; i < n; i++)
    {
	for(j = 0; j < n; j++)
	{
	    mat[i][j] /= denom;
	}
    }
}

double* transpose(double *A, int n)
{
    int i = 0, j = 0;
    double *T = NULL;
    
    T = (double*)malloc(n * n * sizeof(double));

    for(i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            T[j * n + i] = A[i * n + j];
        }
    }
    
    return T;
}

// Performs a forward and backward substitution on a triangular matrix column
void fbSubstitution(double *L, double *U, double *b, double *y, double *x, int n)
{
    int i = 0, j = 0;
    
    for(i = 0; i < n; i++)      // Forward solve Ly = b
    {
        y[i] = b[i];
        
        for(j = 0; j < i; j++) 
        {
            y[i] -= L[i*n+j] * y[j];
        }
        
        y[i] /= L[i*n+i];
    }
    
    for(i = n - 1; i >= 0; i--) // Backward solve Ux = y
    {
        x[i] = b[i];
      
        for(j = i + 1; j < n; j++)
        {
            x[i] -= U[i*n+j] * x[j];
        }
        
        x[i] /= U[i*n+i];
    }
}

// Calculate the inverse of triangular matrix L and U
void LUInversion(double *L, double *U, int n, double *Li, double *Ui)
{
    int i = 0;
    double *x, *y, *b;
    
    b = (double*)calloc(n, sizeof(double));
    y = (double*)calloc(n, sizeof(double));
    x = (double*)calloc(n, sizeof(double));
    
    for(i = 0; i < n; i++)
    {
        b[i] = 1;
        fbSubstitution(L, U, b, y, x, n); 
        memcpy(&Li[i*n], y, n * sizeof(double));
        memcpy(&Ui[i*n], x, n * sizeof(double)); 
        b[i] = 0;   
    }
    
    memcpy(Li, transpose(Li, n), n*n * sizeof(double));
    memcpy(Ui, transpose(Ui, n), n*n * sizeof(double));
    
    free(b);
    free(y);
    free(x);
}

void inverseMatrix(double mat[3][3], double inv[3][3])
{  
   int i = 0, j = 0;

    inv[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    inv[0][1] = mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2];
    inv[0][2] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
    inv[1][0] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
    inv[1][1] = mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
    inv[1][2] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];
    inv[2][0] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
    inv[2][1] = mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1];
    inv[2][2] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

   double  determinant = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) 
                      -  mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]) 
                      +  mat[0][2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]);

   
   for(i = 0; i < 3; i++)
   {
	for(j = 0; j < 3; j++)
	{
	    inv[i][j] /= determinant;
	}
    } 
}

double* matToVector(int n, double mat[n][n])
{  
    int i = 0, j = 0, idx = 0;
    double *vector = NULL;

    vector = (double*)malloc((n * n) * sizeof(double));
    
    for(i = 0; i < n; i++)
    {
	for(j = 0; j < n; j++)
        {
            vector[idx++] = mat[i][j];
        }
    }

    return vector;
}

void swapRow(double *A, double *B, int n, int rowA, int rowB)
{
    int i = 0;
    double aux = 0, *p1 = NULL, *p2 = NULL;

    if(rowA == rowB)
    {
        return;
    }

    for(i = 0; i < n; i++) 
    {
        p1 = mat_elem(A, rowA, i, n);
        p2 = mat_elem(A, rowB, i, n);
        aux = *p1;
        *p1 = *p2;
        *p2 = aux;
    }

    aux = B[rowA]; 
    B[rowA] = B[rowB]; 
    B[rowB] = aux;
}
 
void gaussElimination(double *A, double *B, double *X, int n)
{
    #define A(y, x) (*mat_elem(A, y, x, n))
    
    int j = 0, col = 0, row = 0, maxRow = 0, dia = 0;
    double max = 0.0, aux = 0.0;
    
    for(dia = 0; dia < n; dia++) 
    {
        maxRow = dia;
        max = A(dia, dia);

        for(row = dia + 1; row < n; row++)
        {
            aux = fabs(A(row, dia));
            
            if(aux > max)
            {
                maxRow = row, max = aux;
            }
        }

        swapRow(A, B, n, dia, maxRow);

        for(row = dia + 1; row < n; row++) 
        {
            aux = A(row, dia) / A(dia, dia);
            
            for (col = dia+1; col < n; col++)
            {       
                A(row, col) -= aux * A(dia, col);
            }
            
            A(row, dia) = 0;
            B[row] -= aux * B[dia];
        }
    }
    
    for(row = n - 1; row >= 0; row--) 
    {
        aux = B[row];
        
        for (j = n - 1; j > row; j--)
        {
            aux -= X[j] * A(row, j);
        }
        
        X[row] = aux / A(row, row);
    }
    
    #undef A
}


// No row exchanges - Square matrix
void LUDecomposition(double *A, double *L, double *U, int n)
{	
    int i, j, k;

    for(k = 0; k< n ; k++) 
    {
        L[k*n+k] = 1;

        for(i = k + 1; i < n; i++) 
        {
            L[i*n+k] = A[i*n+k] / A[k*n+k];

            for(j = k + 1; j < n; j++) 
            {
                A[i*n+j] = A[i*n+j] - L[i*n+k] * A[k*n+j];
            }

        }

        for(j = k; j<n; j++) 
        {
            U[k*n+j] = A[k*n+j];
        }
    }
}


void matxvec(double *C, double *A, double *b, int n)
{
    int i = 0, j = 0; 
    
    memset(C, 0, n * sizeof(double));    
    
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            C[i] += A[i*n+j] * b[j];
        }
    }
}

double* matProduct(double *A, double *B, int n)
{
    int i = 0,j = 0,k = 0;
    double sum = 0, *C = NULL;
 
    C = (double*)calloc(n * n, sizeof(double));

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            sum = 0;
            
            for(k = 0; k < n; k++)
            {
                sum += A[i*n+k] * B[k*n+j];
            }
            
            C[i*n+j] = sum;
        }
    }
    
    return C;
}

// Calculate the covariance matrix of one point and their neighbors
void covarianceMat(Lpoint **points, int numPoints, double convarianceMat[3][3])
{
    int i = 0;
    Lpoint meanVector;
    double variance[3] = {0, 0, 0}, mat[3][3];

    meanVector = avgPointP(points, numPoints); 
   
    for(i = 0; i < numPoints; i++)
    {
	variance[0] = points[i]->x - meanVector.x;
	variance[1] = points[i]->y - meanVector.y;
	variance[2] = points[i]->z - meanVector.z;

	selfProduct(3, variance, mat); 
	sumMats(3, mat, convarianceMat);
    }
}

void covarianceMatInt(Lpoint **points, int numPoints, double convarianceMat[3][3])
{
    int i = 0;
    Lpoint meanVector;
    double variance[3] = {0, 0, 0}, mat[3][3];

    meanVector = avgPointInt(points, numPoints); 
   
    for(i = 0; i < numPoints; i++)
    {
	variance[0] = points[i]->x - meanVector.x;
	variance[1] = points[i]->y - meanVector.y;
	variance[2] = points[i]->intensity - meanVector.z;

	selfProduct(3, variance, mat); 
	sumMats(3, mat, convarianceMat);
    }
}

void zeroTriMat(int n, double mat[n][n])
{
    mat[0][0] = 0; mat[0][1] = 0; mat[0][2] = 0;
    mat[1][0] = 0; mat[1][1] = 0; mat[1][2] = 0;
    mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 0;
}