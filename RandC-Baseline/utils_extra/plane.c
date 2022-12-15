#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include "plane.h"
#include "matrix.h"
#include "util.h"


#define NEIGHBORHOOD_RADIUS 1.5 // Radius of the voxel for normal calculation
#define MAX_DEGREE_DIFF 7.5
#define REF_VOX_RAD  1.5
#define REFINE_MAX_PERCENT 100

int nanNormals = 0;

// Calculate the normals of the point by estimating the best fitted plane
double* leastSquares(Lpoint **points, int numPoints)
{   
    int i = 0;
    double xSum = 0.0, ySum = 0.0, xxSum = 0.0, yySum = 0.0, xySum = 0.0, zSum = 0.0, zxSum = 0.0, zySum = 0.0, meanX = 0.0, meanY = 0.0, meanZ = 0.0;
    double *matVector = NULL, *coeficients = NULL;
   
    coeficients = (double*)malloc(4 * sizeof(double));
  
    for(i = 0; i < numPoints; i++)
    { 
	meanX += points[i]->x;
        meanY += points[i]->y;
        meanZ += points[i]->z;
    }
    
    meanX /= numPoints;
    meanY /= numPoints;
    meanZ /= numPoints;
    
    for(i = 0; i < numPoints; i++)
    { 
	xSum += points[i]->x - meanX;
	ySum += points[i]->y - meanY;
	zSum += points[i]->z - meanZ;
	xxSum += (points[i]->x - meanX) * (points[i]->x - meanX);
	yySum += (points[i]->y - meanY) * (points[i]->y - meanY);
	xySum += (points[i]->x - meanX) * (points[i]->y - meanY);
	zxSum += (points[i]->z - meanZ) * (points[i]->x - meanX);
	zySum += (points[i]->z - meanZ) * (points[i]->y - meanY);
    }

    double matrix[3][3] = 
    {  
        {xxSum, xySum, xSum},  
        {xySum, yySum, ySum},  
        {xSum,  ySum, numPoints} 
    };

    double lvector[3] = {zxSum, zySum, zSum}; 
    
    matVector = matToVector(3, matrix); 
    gaussElimination(matVector, lvector, coeficients, 3);
    coeficients[3] = -(coeficients[0] * meanX + coeficients[1] * meanY + coeficients[2] * meanZ);
    
    free(matVector);
    
    return coeficients;
}

// Same as above but with intensity values instead of z
double* leastSquaresInt(Lpoint **points, int numPoints)
{   
    int i = 0;
    double xSum = 0.0, ySum = 0.0, xxSum = 0.0, yySum = 0.0, xySum = 0.0, iSum = 0.0, ixSum = 0.0, izSum = 0.0, meanX = 0.0, meanY = 0.0, meanI = 0.0;
    double *matVector = NULL, *coeficients = NULL;
   
    coeficients = (double*)malloc(4 * sizeof(double));
  
    for(i = 0; i < numPoints; i++)
    { 
	meanX += points[i]->x;
        meanY += points[i]->y;
        meanI += points[i]->intensity;
    }
    
    meanX /= numPoints;
    meanY /= numPoints;
    meanI /= numPoints;
    
    for(i = 0; i < numPoints; i++)
    { 
	xSum += points[i]->x - meanX;
	ySum += points[i]->y - meanY;
	iSum += points[i]->intensity - meanI;
	xxSum += (points[i]->x - meanX) * (points[i]->x - meanX);
	yySum += (points[i]->y - meanY) * (points[i]->y - meanY);
	xySum += (points[i]->x - meanX) * (points[i]->y - meanY);
	ixSum += (points[i]->intensity - meanI) * (points[i]->x - meanX);
	izSum += (points[i]->intensity - meanI) * (points[i]->y - meanY);
    }

    double matrix[3][3] = 
    {  
        {xxSum, xySum, xSum},  
        {xySum, yySum, ySum},  
        {xSum,  ySum, numPoints} 
    };

    double lvector[3] = {ixSum, izSum, iSum}; 
    
    matVector = matToVector(3, matrix); 
    gaussElimination(matVector, lvector, coeficients, 3);
    coeficients[3] = -(coeficients[0] * meanX + coeficients[1] * meanY + coeficients[2] * meanI);
    
    free(matVector);
    
    return coeficients;
}

// Same as above but with another variable "var" values
double* leastSquaresVar(Lpoint **points, double *var, int numPoints)
{   
    int i = 0;
    double xSum = 0.0, ySum = 0.0, xxSum = 0.0, yySum = 0.0, xySum = 0.0, vSum = 0.0, vxSum = 0.0, vzSum = 0.0, meanX = 0.0, meanY = 0.0, meanV = 0.0;
    double *matVector = NULL, *coeficients = NULL;
   
    coeficients = (double*)malloc(4 * sizeof(double));
  
    for(i = 0; i < numPoints; i++)
    { 
	meanX += points[i]->x;
        meanY += points[i]->y;
        meanV += var[i];
    }
    
    meanX /= numPoints;
    meanY /= numPoints;
    meanV /= numPoints;
    
    for(i = 0; i < numPoints; i++)
    { 
	xSum += points[i]->x - meanX;
	ySum += points[i]->y - meanY;
	vSum += var[i] - meanV;
	xxSum += (points[i]->x - meanX) * (points[i]->x - meanX);
	yySum += (points[i]->y - meanY) * (points[i]->y - meanY);
	xySum += (points[i]->x - meanX) * (points[i]->y - meanY);
	vxSum += (var[i] - meanV) * (points[i]->x - meanX);
	vzSum += (var[i] - meanV) * (points[i]->y - meanY);
    }

    double matrix[3][3] = 
    {  
        {xxSum, xySum, xSum},  
        {xySum, yySum, ySum},  
        {xSum,  ySum, numPoints} 
    };

    double lvector[3] = {vxSum, vzSum, vSum}; 
    
    matVector = matToVector(3, matrix); 
    gaussElimination(matVector, lvector, coeficients, 3);
    coeficients[3] = -(coeficients[0] * meanX + coeficients[1] * meanY + coeficients[2] * meanV);
    
    free(matVector);
    
    return coeficients;
}

double* fitCurve(Lpoint **points, int num)
{   
    int i = 0;
    double x = 0.0, y = 0.0, xSum = 0.0, ySum = 0.0, x2Sum = 0.0, x3Sum = 0.0, xySum = 0.0, x4Sum = 0.0, x2ySum = 0.0, meanX = 0.0, meanY = 0.0;
    double *matVector = NULL, *coeficients = NULL;
   
    coeficients = (double*)malloc(3 * sizeof(double));
  
    for(i = 0; i < num; i++)
    { 
	meanX += points[i]->x;
        meanY += points[i]->y;
    }
    
    meanX /= num;
    meanY /= num;
    
    for(i = 0; i < num; i++)
    { 
        x = (points[i]->x - meanX);
        y = (points[i]->y - meanY);
	xSum += x;
	ySum += y;
        xySum += x * y;
	x2Sum += x * x;      
        x3Sum += x * x * x;
        x4Sum += x * x * x * x;
        x2ySum +=  x * x * y;
    }

    double matrix[3][3] = 
    {  
        {num,  xSum,  x2Sum},  
        {xSum, x2Sum, x3Sum},  
        {x2Sum, x3Sum, x4Sum} 
    };

    double lvector[3] = {ySum, xySum, x2ySum}; 
    
    matVector = matToVector(3, matrix); 
    gaussElimination(matVector, lvector, coeficients, 3);
    
    free(matVector);
    
    return coeficients;
}

// Estime the best 2D fitted line of a set of points (no mean subtraction -> few points causes accuracy erros)
int fitLineP(Lpoint **points, int num, double **x)
{   
    int i = 0;
    double xSum = 0.0, ySum = 0.0, xxSum = 0.0, xySum = 0.0;
    double *Av = NULL;
   
    *x = calloc(2, sizeof(double)); 
 
    for(i = 0; i < num; i++)
    { 
	xSum += points[i]->x;
	ySum += points[i]->y;	
	xxSum += points[i]->x * points[i]->x;
	xySum += points[i]->x * points[i]->y;
    }

    double A[2][2] = 
    {  
        {xxSum, xSum},  
        { xSum,  num}  
    };
    
    double b[2] = {xySum, ySum}; 
    
    Av = matToVector(2, A);
    
    gaussElimination(Av, b, *x, 2); 
    
    free(Av);
    
    return number((*x)[0]) && number((*x)[1]); 
}

int fitLine(Lpoint *points, int num, double **x)
{   
    int i = 0;
    double xSum = 0.0, ySum = 0.0, xxSum = 0.0, xySum = 0.0;
    double *Av = NULL;
   
    *x = calloc(2, sizeof(double)); 
 
    for(i = 0; i < num; i++)
    { 
	xSum += points[i].x;
	ySum += points[i].y;	
	xxSum += points[i].x * points[i].x;
	xySum += points[i].x * points[i].y;
    }

    double A[2][2] = 
    {  
        {xxSum, xSum},  
        { xSum,  num}  
    };
    
    double b[2] = {xySum, ySum}; 
    
    Av = matToVector(2, A);
    
    gaussElimination(Av, b, *x, 2); 
    
    free(Av);
    
    return number((*x)[0]) && number((*x)[1]); 
}

Lvector planeNormal(Lpoint *point, double *coeficients) // TODO: Handle NaN
{
    Lvector v1, v2, normal; 
    double z = 0.0;
      
    z = coeficients[0] * point->x + coeficients[1] * point->y + coeficients[2]; // z = Ax + By + C
    v1.i = 1; v1.j = 0; v1.k = coeficients[0];
    v2.i = 0; v2.j = 1; v2.k = coeficients[1];
    normal = crossProduct(&v1, &v2);
    normalizeV3D(&normal);
        
    if(normal.i != normal.i || normal.j != normal.j || normal.k != normal.k)
    {	      
        normal.i = 0;
        normal.j = 0;
        normal.k = 0;
    }     
       
    return normal;
}

float* vectorAngles(float normal[3])
{
    float denom = 0, aux[3], *out = NULL;
    
    out = (float*)malloc(3 * sizeof(float));

    denom = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));
    
    if(denom != 0)
    {
        aux[0] = normal[0] / denom;
        aux[1] = normal[1] / denom;
        aux[2] = normal[2] / denom;
        out[0] = acos(aux[0]) * 180.0 / M_PI;
        out[1] = acos(aux[1]) * 180.0 / M_PI;
        out[2] = acos(aux[2]) * 180.0 / M_PI;
    }
    else
    {
        out[0] = 0;
        out[1] = 0;
        out[2] = 0;
    }
    
    return out;
}




double planeDistance(Lpoint *point, double *coefs) 
{
    double distance = 0.0;

    distance = coefs[0] * point->x + coefs[1] * point->y + coefs[2] * point->z + coefs[3] ;
    distance /= sqrt(coefs[0] * coefs[0] + coefs[1] * coefs[1] + coefs[2] * coefs[2]);

    return fabs(distance);
}

