#ifndef UTIL_H

#include "octree.h"
#include "vector.h"
#include <time.h>

#define	UTIL_H

void stclock(struct timespec *start);

void declock(struct timespec *start, struct timespec *end); void declock2(struct timespec *start, struct timespec end);

void declocksav(struct timespec *start, struct timespec *end, double *time);


int onRange(double value, double offset, double center);

int onDegree(float normal[3], float interval, float epicenter[3]);

int onIntensity(float int1, float int2, float interval);

float euclideanDistance3D(Lpoint point1, Lpoint point2);

float euclideanDistance2D(Lpoint *point1, Lpoint *point2);


int insideCircle(Lpoint *point, Lpoint *center, float radius);

double polygonPerimeter(Lpoint *points, int numPoints);

Lpoint polygonCentroid(Lpoint *points, int numPoints);


double polygonAreaP(Lpoint *points, int numPoints);


Lpoint* pointReposition(Lpoint *array, long int *orderIndices, int size);


void resetOrigin(Lpoint *points, int numPoints, Lpoint *center);


int hasValue(int *array, int size, int value);

int* removeValue(int *array, int size, int value);

double radians(double degrees);

double degrees(double radians);

double stdev(double *values, int num);

int number(double x);


void pointToVec(Vector3D *dest, Lpoint *source);

Lpoint nearestLinePoint(Lpoint *p, double A, double B, double C);

void convexHull(Lpoint* points, int num, Lpoint** out_hull, int* hullSize);

double* calcEigenvalues(double mat[3][3]);

int pointInPolygon(Lpoint *p, Lpoint *polygon, int polSize);


int inHull(Lpoint *point, Lpoint* hull, int num);

int rmvIntDuplicates(int **array, int n);

int* rmvIntDuplicates2(int *array, int n, int *newSize);

int removeDuplicates(int *a, int array_size);

Lvector eigenNormal(Lpoint *point, Octree octree, double radius);

Lvector eigenNormalInt(Lpoint *point, Octree octree, double radius);

void checkNAN(double x);

void addPoint(int pos, Lpoint *p, Lpoint **points, int *num);


int inHull2(Lpoint *point, Lpoint* hull, int num);

int pnpoly(int n, Lpoint *vert, Lpoint *test);

void sortBySize(int *idx, int *val, int n);

int diffOnRange(double a, double b, double factor);

#endif	/* UTIL_H */
