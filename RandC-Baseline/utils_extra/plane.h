#ifndef PLANE_H

#include "octree.h"
#include "util.h"
#include "matrix.h"

#define	PLANE_H

double* leastSquares(Lpoint **points, int numNeighs);

double* leastSquaresInt(Lpoint **points, int numPoints);

double* leastSquaresVar(Lpoint **points, double *var, int numPoints);

double* leastSquaresP(Lpoint *points, int numPoints, FILE *file);

int fitLineP(Lpoint **points, int num, double **x);

int fitLine(Lpoint *points, int num, double **x);

double* fitCurve(Lpoint **points, int num);

double fitCurveRAW(Lpoint *points, int num, double **coefs); double* fitCurveRAW3(Lpoint *points, int num);

float* vectorAngles(float normal[3]);

int clusterNormals(Lvector *normals, int size, FILE *file);

Lvector planeNormal(Lpoint *point, double *coeficients);

double planeDistance(Lpoint *point, double *plane);


#endif  /* PLANE_H */

