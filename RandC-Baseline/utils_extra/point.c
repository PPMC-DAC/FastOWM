#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "point.h"
#include "util.h"

Lpoint createPoint(unsigned int id, double x, double y, double z, unsigned int intensity, double normal[3])
{
    Lpoint p;

    p.id = id;
    p.x = x;
    p.y = y;
    p.z = z;
    p.intensity = intensity;
    p.normal[0] = normal[0];
    p.normal[1] = normal[1];
    p.normal[2] = normal[2];

    return p;
}

// Copy points coordinates between points
void copyCoords(Lpoint *dest, Lpoint *source)
{
    dest->x = source->x;
    dest->y = source->y;
    dest->z = source->z;
}

// Calculate the mid point of two points
Lpoint midPoint(Lpoint *p1, Lpoint *p2)
{
    Lpoint mid;
    
    mid.x = (p1->x + p2->x) / 2.0;
    mid.y = (p1->y + p2->y) / 2.0;
    mid.z = (p1->z + p2->z) / 2.0;

    return mid;
}

// Calculate the average point of a set of points
Lpoint avgPointP(Lpoint **points, int num)
{
    int i = 0;
    Lpoint avg;
    
    avg.x = 0;
    avg.y = 0;
    avg.z = 0;
    
    for(i = 0; i < num; i++)
    {
        avg.x += points[i]->x;
        avg.y += points[i]->y;
        avg.z += points[i]->z;
    }
    
    avg.x /= num;
    avg.y /= num;
    avg.z /= num;
    
    return avg;
}

Lpoint avgPoint(Lpoint *points, int num)
{
    int i = 0;
    Lpoint avg;
    
    avg.x = 0;
    avg.y = 0;
    avg.z = 0;
    
    for(i = 0; i < num; i++)
    {
        avg.x += points[i].x;
        avg.y += points[i].y;
        avg.z += points[i].z;
    }
    
    avg.x /= num;
    avg.y /= num;
    avg.z /= num;
    
    return avg;
}

Lpoint avgPointInt(Lpoint **points, int num)
{
    int i = 0;
    Lpoint avg;
    
    avg.x = 0;
    avg.y = 0;
    avg.z = 0;
    
    for(i = 0; i < num; i++)
    {
        avg.x += points[i]->x;
        avg.y += points[i]->y;
        avg.z += points[i]->intensity;
    }
    
    avg.x /= num;
    avg.y /= num;
    avg.z /= num;
    
    return avg;
}

// Join two arrays with points in one array
Lpoint* mergePoints(Lpoint *pts1, int num1, Lpoint *pts2, int num2)
{    
    Lpoint *merged = NULL;
    
    merged = (Lpoint*)malloc((num1 + num2) * sizeof(Lpoint));   
    memcpy(merged, pts1, num1 * sizeof(Lpoint));
    memcpy(&merged[num1], pts2, num2 * sizeof(Lpoint));
    
    return merged;
}

// Return the different points of two sets of points
Lpoint** distinctPoints(Lpoint **pts1, int numPts1, Lpoint **pts2, int numPts2, int *distNum) // TODO: It should be a better way
{
    int i = 0, j = 0, duplicate = 0, numDist = 0;
    Lpoint **distinct = NULL;
    
    for(i = 0; i < numPts1; i++)
    {
        duplicate = 0;
        
        for(j = 0; j < numPts2; j++)
        {
            if(pts1[i]->id == pts2[j]->id)
            {
                duplicate = 1; break;
            }
        }
        
        if(!duplicate)
        {
            distinct = realloc(distinct, ++numDist * sizeof(Lpoint*));
            distinct[numDist-1] = pts1[i];
        }
    }
    
    *distNum = numDist;
    
    return distinct;
}

// Sort the points in inverse order
void invertPoints(Lpoint *points, int num)
{
    int i = 0;
    Lpoint *aux = NULL;
    
    aux = (Lpoint*)malloc(num * sizeof(Lpoint));
    memcpy(aux, points, num * sizeof(Lpoint));
    
    for(i = 0; i < num; i++)
    {
        points[i] = aux[num-1-i];
    }
    
    free(aux);
}

// Sort the points pointers in inverse order
void invertPointsP(Lpoint **points, int num)
{
    int i = 0;
    Lpoint **aux = NULL;
    
    aux = malloc(num * sizeof(Lpoint*));
    memcpy(aux, points, num * sizeof(Lpoint*));
    
    for(i = 0; i < num; i++)
    {
        points[i] = aux[num-1-i];
    }
    
    free(aux);
}

// Rotate a point from origin
void rotatePoint(Lpoint *p, double sin, double cos, Lpoint *orig)
{
    double x = 0.0, y = 0.0;
    
    x = cos * (p->x - orig->x) - sin * (p->y - orig->y) + orig->x;
    y = sin * (p->x - orig->x) + cos * (p->y - orig->y) + orig->y;
    p->x = x;
    p->y = y;
}

// Calculate the distance between a point and a line
double pointLineDistance(Lpoint *p, Lpoint* l1, Lpoint* l2)
{
    double dist = 0.0, y2_y1 = 0.0, x2_x1 = 0.0;
    
    y2_y1 = l2->y - l1->y;
    x2_x1 = l2->x - l1->x;
    
    dist = fabs(y2_y1 * p->x - x2_x1 * p->y + l2->x * l1->y - l2->y * l1->x);
    dist /= sqrt(y2_y1 * y2_y1 + x2_x1 * x2_x1);
    
    return dist;
}

// Get two points within a line
void linePoints(Lpoint *lp, double *coefs, Lpoint *p1, Lpoint *p2)
{
    copyCoords(p1, lp);
    copyCoords(p2, lp);
    p1->x += 3.0;
    p2->x -= 3.0;
    p1->y = coefs[0] * p1->x + coefs[1]; // y = Ax + B
    p2->y = coefs[0] * p2->x + coefs[1];
}

// Get the point where two lines intersect
int lineIntersection(Lpoint *p1, Lpoint *p2, Lpoint *p3, Lpoint *p4, Lpoint *intersect)
{
    double x1y2_y1x2 = 0.0, x3y4_y3x4 = 0.0, denom = 0.0;
    
    x1y2_y1x2 = p1->x * p2->y - p1->y * p2->x;
    x3y4_y3x4 = p3->x * p4->y - p3->y * p4->x;
    denom = (p1->x - p2->x) * (p3->y - p4->y) - (p1->y - p2->y) * (p3->x - p4->x);
    
    if(denom == 0.0) // Lines are parallel 
    {
        return 0;
    }
    
    intersect->x = (x1y2_y1x2 * (p3->x - p4->x) - (p1->x - p2->x) * x3y4_y3x4) / denom;
    intersect->y = (x1y2_y1x2 * (p3->y - p4->y) - (p1->y - p2->y) * x3y4_y3x4) / denom;
    intersect->z = (p1->z + p2->z) / 2;
    
    return 1;
}

// Get points of a curve with 1 meter step
Lpoint* curvePoints(Lpoint *points, int numPoints, double *coefs, int *steps)
{
    int i = 0, minMeters = 10; 
    double x = 0.0, y = 0.0, step = 0.5, minX = DBL_MAX, maxX = -DBL_MAX, meanX = 0.0, meanZ = 0.0;
    Lpoint *curvePts = NULL;
    
    for(i = 0; i < numPoints; i++)
    {
        meanX += points[i].x;
        meanZ += points[i].z;
        if(points[i].x < minX) minX = points[i].x;
        if(points[i].x > maxX) maxX = points[i].x;
    }
    
    meanX /= numPoints;
    meanZ /= numPoints;
    *steps = (maxX - minX) / step;
    
    if(*steps * step < minMeters)
    {
        minX = meanX - (minMeters / 2.0);
        *steps = minMeters / step;
    }
    
    curvePts = (Lpoint*)malloc(*steps * sizeof(Lpoint));  //printf("GET %lf %lf %lf\t %lf\n",  coefs[0],  coefs[1],  coefs[2], x); 
   
    for(i = 0; i < *steps; i++)
    {
        x = minX + step * i; 
        y = coefs[2] * x * x + coefs[1] * x + coefs[0]; // y = Axx + Bx + C
        curvePts[i].x = x;
        curvePts[i].y = y;
        curvePts[i].z = meanZ;
    }

    return curvePts;
}

// Get the slope-intercept equation
double* lineCoeficients(Lpoint *p1, Lpoint *p2)
{
    double *coefs = NULL;
    
    coefs = (double*)malloc(2 * sizeof(double));
    
    coefs[0] = (p2->y - p1->y) / (p2->x - p1->x);  
    coefs[1] = p1->y - coefs[0] * p1->x;
    
    return coefs;
}

double boundingCubeWidth(Lpoint *points, int numPoints, Lpoint center)
{
    int i = 0;
    double centerToMinDist = 0, centerToMaxDist = 0, distX = 0, distY = 0, minY = DBL_MAX, maxY = -DBL_MAX, minX = DBL_MAX, maxX = -DBL_MAX;
  
    for(i = 0; i < numPoints; i++)
    {               
        if(points[i].y < minY)   minY = points[i].y;      
        if(points[i].y > maxY)   maxY = points[i].y;
        if(points[i].x < minX)   minX = points[i].x;      
        if(points[i].x > maxX)   maxX = points[i].x;               				
    } 
    
    centerToMinDist = fabs(center.x - minX);
    centerToMaxDist = fabs(center.x - maxX);
    if(centerToMinDist > centerToMaxDist)   
        distX = centerToMinDist;
    else                                    
        distX = centerToMaxDist;  
    centerToMinDist = fabs(center.y - minY);
    centerToMaxDist = fabs(center.y - maxY);
    if(centerToMinDist > centerToMaxDist)   
        distY = centerToMinDist;
    else                                    
        distY = centerToMaxDist;
    
    if(distX > distY)   return distX;
    else                return distY;
}