#ifndef VECTOR_H

#include "point.h"

#define	VECTOR_H

typedef struct 
{
    double i;
    double j;
    double k;

} Lvector;  // Just a vector

void checkVecNAN(Lvector *v);

Lvector calcVector(Lpoint *p1, Lpoint *p2);

Lvector perpenVector(Lvector *v);

void multVector(Lvector *v, double number);

void subVectors(Lvector *v1, Lvector v2);

void addVectors(Lvector *v1, Lvector v2);

double dotProduct(Lvector *v1, Lvector *v2);

Lvector crossProduct(Lvector *v1, Lvector *v2);

void normalizeV2D(Lvector *v);

void normalizeV3D(Lvector *v);

double angleV2D(Lvector *v);

double angleV2DN(Lvector *v);

double angleVS2D(Lvector *v1, Lvector *v2);

double norm2D(Lvector *v);

double norm3D(Lvector *v);

Lpoint movedPoint(Lpoint *orig ,Lpoint *dest, Lvector v, double dist);

Lpoint getDest(Lpoint *orig, Lvector *v, double dist);

Lpoint* discretizeVector(Lpoint *orig, Lpoint *dest, double step, int *numPoints);

void segmentEnds(Lpoint *p, Lvector *v, Lpoint *e1, Lpoint *e2, int centered);

void selfProduct(int n, double vector[n], double product[n][n]);

void vecmult(double *dest, double number, double *v, int n);

void subVec(double number, double *v, int size);

void vecssub(double *sub, double *a, double *b, int n);

double* vecsadd(double *a, double *b, int n);

void addVecs(double *a, double *b, int size);

double* vecsProdMat(double *a, int n);

int nanVec(Lvector *v);

void printVec(double *v, int n);

#endif	/* VECTOR_H */
