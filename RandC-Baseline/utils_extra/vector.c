#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vector.h"
#include "util.h"

void checkVecNAN(Lvector *v)
{
    if(nanVec(v))
    {
        v->i = 0.0;
        v->j = 0.0;
        v->k = 0.0;
    }
}

// Calculate the vector between two points
Lvector calcVector(Lpoint *p1, Lpoint *p2)
{
    Lvector v;

    v.i = p2->x - p1->x;
    v.j = p2->y - p1->y;
    v.k = p2->z - p1->z;
    
    return v;
}

// Left-wise perpendicular vector
Lvector perpenVector(Lvector *v)
{
    Lvector p;
    
    p.i = -v->j;
    p.j =  v->i;
    p.k =  v->k;
    
    return p;
}

void multVector(Lvector *v, double number)
{ 
    v->i *= number;
    v->j *= number;
    v->k *= number;
}

void subVectors(Lvector *v1, Lvector v2)
{ 
    v1->i -= v2.i;
    v1->j -= v2.j;
    v1->k -= v2.k;
}

void addVectors(Lvector *v1, Lvector v2)
{ 
    v1->i += v2.i;
    v1->j += v2.j;
    v1->k += v2.k;
}

Lvector crossProduct(Lvector *v1, Lvector *v2)
{
    Lvector r;

    r.i = v1->j * v2->k - v1->k * v2->j;
    r.j = v1->k * v2->i - v1->i * v2->k;
    r.k = v1->i * v2->j - v1->j * v2->i;

    return r;
}

double dotProduct(Lvector *v1, Lvector *v2)
{
    return v1->i * v2->i + v1->j * v2->j + v1->k * v2->k;
}

double norm2D(Lvector *v)
{
    return sqrt(v->i * v->i + v->j * v->j);
}

double norm3D(Lvector *v)
{
    return sqrt(v->i * v->i + v->j * v->j + v->k * v->k);
}

void normalizeV2D(Lvector *v)
{
    double norm = 0.0;
    
    norm = norm2D(v);
    v->i /= norm;
    v->j /= norm;
    v->k /= norm;
}

void normalizeV3D(Lvector *v)
{
    double norm = 0.0;
    
    norm = norm3D(v);
    v->i /= norm;
    v->j /= norm;
    v->k /= norm;
}

// Angle with the X axis
double angleV2D(Lvector *v)
{
    return atan2(v->j, v->i);
}

// Angle with the X axis 1-quadrant normalized
double angleV2DN(Lvector *v)
{
    double angle = 0.0;
    
    angle = atan2(v->j, v->i);
    
    if(angle < 0)
    {
        angle += M_PI;
    }
    
    return angle;
}

// Angle between two vectors quadrant mod
double angleVS2D(Lvector *v1, Lvector *v2)
{
    float angle = 0.0;
   
    angle = atan2(v2->j, v2->i) - atan2(v1->j, v1->i);
    
    if(angle < 0) angle += 2 * M_PI;
    
    return angle;
}

// Get a point moved from origin through the vector formed with another point
Lpoint movedPoint(Lpoint *orig, Lpoint *dest, Lvector v, double dist)
{
    Lpoint p;
    
    multVector(&v, dist);
    copyCoords(&p, orig);
    p.x += v.i;
    p.y += v.j;
    p.z = (orig->z + dest->z) / 2;
    
    return p;
}

// Get the destination point of moving a point form a vector a given distance
Lpoint getDest(Lpoint *orig, Lvector *v, double dist)
{
    Lpoint dest;
    
    normalizeV2D(v);
    dest.x = orig->x + v->i * dist;
    dest.y = orig->y + v->j * dist;
    dest.z = orig->z;
    
    return dest;
}

// Get points form a vector with a given step
Lpoint* discretizeVector(Lpoint *orig, Lpoint *dest, double step, int *numPoints)
{
    int i = 0, iters = 0;
    double dist = 0.0;
    Lpoint *points = NULL;
    Lvector dir;
      
    *numPoints = 0;
    dir = calcVector(orig, dest);
    dist = norm2D(&dir);
    iters = dist / step;
    
    for(i = 0; i < iters; i++)
    {
        points = (Lpoint*)realloc(points, ++*numPoints * sizeof(Lpoint));
        points[*numPoints-1] = getDest(orig, &dir, step * (i+1));
    }
    
    return points;
}

// Get the end points of a segment formed by a point and a vector
void segmentEnds(Lpoint *p, Lvector *v, Lpoint *e1, Lpoint *e2, int centered)
{
    copyCoords(e1, p);
    copyCoords(e2, p);
    
    if(centered)
    {    
        e1->x -= v->i;
        e1->y -= v->j;
        e2->x += v->i;
        e2->y += v->j;
    }  
    else
    {    
        e2->x += 2.0 * v->i;
        e2->y += 2.0 * v->j;
    }    
}

// Returns the matrix of multipliying a vector by itself
void selfProduct(int n, double vector[n], double product[n][n])
{
    int i = 0, j = 0;

    for(i = 0; i < n; i++)
    {
	for(j = 0; j < n; j++)
	{
	    product[i][j] = vector[i] * vector[j];
	}
    }
}

void vecmult(double *dest, double number, double *v, int n)
{ 
    int i = 0;
      
    for(i = 0; i < n; i++)
    {
        dest[i] = v[i] * number;
    }
}

void vecssub(double *sub, double *a, double *b, int n)
{ 
    int i = 0;
    
    for(i = 0; i < n; i++)
    {
        sub[i] = a[i] - b[i];
    }
}

double* vecsadd(double *a, double *b, int n)
{ 
    int i = 0;
    double *r = NULL;
    
    r = (double*)malloc(n * sizeof(double));
    
    for(i = 0; i < n; i++)
    {
        r[i] = a[i] + b[i];
    }
    
    return r;
}


void addVecs(double *a, double *b, int size)
{ 
    int i = 0;
    
    for(i = 0; i < size; i++)
    {
        a[i] += b[i];
    }
}

int nanVec(Lvector *v)
{
    return !number(v->i) || !number(v->j) || !number(v->k);   
}

void printVec(double *v, int n)
{
    int i = 0;
    
    for(i = 0; i < n; i++)
    {
        printf("%.4lf\t", v[i]);
    }
    
    printf("\n");
}



