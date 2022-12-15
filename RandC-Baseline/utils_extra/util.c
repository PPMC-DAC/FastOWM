#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "plane.h"
#include "util.h"
#include "time.h"
#include "matrix.h"
#include "file.h"


typedef struct coordinate 
{
    double value;
    long int index;
    
} Coord;

// Time
void stclock(struct timespec *start)
{
    clock_gettime(CLOCK_MONOTONIC, start);    
}

void declock(struct timespec *start, struct timespec *end)
{
    double time = 0.0;
    
    clock_gettime(CLOCK_MONOTONIC, end);
    time = ((double)end->tv_sec + 1.0e-9*end->tv_nsec) - ((double)start->tv_sec + 1.0e-9*start->tv_nsec);
    printf("%.2lf s (%.2lf min)\n", time, time / 60.0);
}

void declock2(struct timespec *start, struct timespec end)
{
    double time = 0.0;

    time = ((double)end.tv_sec + 1.0e-9*end.tv_nsec) - ((double)start->tv_sec + 1.0e-9*start->tv_nsec);
    printf("%.2lf s (%.2lf min)\n", time, time / 60.0);
}

void declocksav(struct timespec *start, struct timespec *end, double *time)
{ 
    clock_gettime(CLOCK_MONOTONIC, end);
    *time = ((double)end->tv_sec + 1.0e-9*end->tv_nsec) - ((double)start->tv_sec + 1.0e-9*start->tv_nsec);
}

int onRange(double value, double offset, double center)
{
    if(value >= (center - offset) && value <= (center + offset))
    {
        return 1;
    }

    return 0;
}

int onDegree(float normal[3], float interval, float epicenter[3])
{
    float *degrees = NULL, *epiDegrees = NULL, diffs[3];

    degrees = vectorAngles(normal);
    epiDegrees = vectorAngles(epicenter); 
    diffs[0] = fabs(epiDegrees[0] - degrees[0]);
    diffs[1] = fabs(epiDegrees[1] - degrees[1]);
    diffs[2] = fabs(epiDegrees[2] - degrees[2]);
 
    free(degrees);
    free(epiDegrees);
    
    if(diffs[0] < interval && diffs[1] < interval && diffs[2] < interval)
    {
        return 1;
    } 

    return 0;
}

int onIntensity(float int1, float int2, float interval)
{
    float offset = 0.0;
    
    if(int1 > int2)
    {
        offset = int1 * interval;
    }
    else
    {
        offset = int2 * interval;
    }

    return onRange(int1, offset, int2);
}



// Geometry
float euclideanDistance2D(Lpoint *point1, Lpoint *point2)
{
    float diffX, diffY;
    
    diffX = (point1->x - point2->x) * (point1->x - point2->x);
    diffY = (point1->y - point2->y) * (point1->y - point2->y);

    return sqrt(diffX + diffY);
}

int insideCircle(Lpoint *point, Lpoint *center, float radius)
{
    return (point->x - center->x) * (point->x - center->x) + (point->y - center->y) * (point->y - center->y) < radius * radius;
}




Lpoint nearestLinePoint(Lpoint *p, double A, double B, double C) 
{      
    Lpoint r;
    
    r.x = (B * ( B * p->x - A * p->y) - A * C) / A * A + B * B;
    r.y = (A * (-B * p->x + A * p->y) - B * C) / A * A + B * B;
    
    return r;
}


double polygonAreaP(Lpoint *points, int numPoints) // Shoelace method
{
    int i = 0;
    float area = 0.0;
    double xySum = 0.0, yxSum = 0.0;
    
    for(i = 0; i < numPoints - 1; i++)
    {
        xySum += points[i].x * points[i+1].y;
        yxSum += points[i].y * points[i+1].x;
    }
    
    xySum += points[i].x * points[0].y;
    yxSum += points[i].y * points[0].x;

    area = fabs((xySum - yxSum)) * 0.5;
    
    return area;
}

double polygonPerimeter(Lpoint *points, int numPoints) 
{
    int i = 0;
    double perimeter = 0.0;
    
    for(i = 0; i < numPoints - 1; i++)
    {   
        perimeter += euclideanDistance2D(&points[i], &points[i+1]);       
    }
    
    perimeter += euclideanDistance2D(&points[i], &points[0]);
    
    return perimeter;
}

// Calculate the centroid of a polygon formed by ordered pairs points
Lpoint polygonCentroid(Lpoint *points, int numPoints)
{
    int i = 0;
    double partialArea = 0.0, dobleArea = 0.0;
    Lpoint centroid;
    
    centroid.x = 0;
    centroid.y = 0;

    for(i = 0; i < numPoints - 1; i++)
    {
        partialArea = points[i].x * points[i+1].y - points[i+1].x * points[i].y;
        centroid.x += (points[i].x + points[i+1].x) * partialArea;
        centroid.y += (points[i].y + points[i+1].y) * partialArea;
        dobleArea += partialArea;
    }
    
    partialArea = points[i].x * points[0].y - points[0].x * points[i].y;
    centroid.x += (points[i].x + points[0].x) * partialArea;
    centroid.y += (points[i].y + points[0].y) * partialArea;
    dobleArea += partialArea;

    centroid.x /= 3 * dobleArea;
    centroid.y /= 3 * dobleArea;
    centroid.z = dobleArea / 2.0;
    
    return centroid;
}

Lpoint* pointReposition(Lpoint *array, long int *orderIndices, int size)
{
    int i = 0, position = 0;
    Lpoint *auxArray = NULL;
    
    auxArray = (Lpoint*)malloc(size * sizeof(Lpoint));
      
    for(i = 0; i < size; i++)
    {      
        position = orderIndices[i];
        auxArray[i] = array[position];
    }
    
    return auxArray;
}

int compareX (const void * a, const void * b)
{
    double v1 = ((struct coordinate *)a)->value; 
    double v2 = ((struct coordinate *)b)->value; 
   
    return (v1 > v2) - (v1 < v2);
}
  


void sortBySize(int *idx, int *val, int n)
{
    int i = 0;
    Coord coord, *coords = NULL;
    
    coords = (Coord*)malloc(n * sizeof(Coord));
    
    for(i = 0; i < n; i++)
    {      
        coord.index = idx[i];
        coord.value = val[i];
        coords[i] = coord;   
    }   
    
    qsort(coords, n, sizeof(struct coordinate), compareX); 
    
    for(i = 0; i < n; i++)
    {
       idx[i] = coords[n-1-i].index; // Bigger to smaller
    }   
    
    free(coords);
}


void resetOrigin(Lpoint *points, int numPoints, Lpoint *center)
{
    int i = 0;
    
    for(i = 0; i < numPoints; i++)
    {
        points[i].x -= center->x;
        points[i].y -= center->y;
    }
}

void planing(Lpoint p1, Lpoint p2, Lvector normal, FILE *file)
{
    Lvector tan, bitan, v1, v2, vtx1, vtx2, vtx3, vtx4;
    float D = 7.5;
    v1.i = p1.x; v1.j = p1.y; v1.k = p1.z;
    v2.i = p2.x; v2.j = p2.y; v2.k = p2.z;
    vtx1.i = v1.i; vtx1.j = v1.j; vtx1.k = v1.k;
    vtx2.i = v1.i; vtx2.j = v1.j; vtx2.k = v1.k;
    vtx3.i = v1.i; vtx3.j = v1.j; vtx3.k = v1.k;
    vtx4.i = v1.i; vtx4.j = v1.j; vtx4.k = v1.k;
    
    tan = calcVector(&p1, &p2);
    normalizeV2D(&tan);
    bitan = crossProduct(&tan, &normal);
    multVector(&tan, D);
    multVector(&bitan, D);
    
    subVectors(&vtx1, tan);  subVectors(&vtx1, bitan);
    addVectors(&vtx2, tan);  subVectors(&vtx2, bitan);
    addVectors(&vtx3, tan);  addVectors(&vtx3, bitan);
    subVectors(&vtx4, tan);  addVectors(&vtx4, bitan);

    fprintf(file, "%f\t%f\t%f\n", vtx1.i, vtx1.j, vtx1.k);
    fprintf(file, "%f\t%f\t%f\n", vtx2.i, vtx2.j, vtx2.k);
    fprintf(file, "%f\t%f\t%f\n", vtx3.i, vtx3.j, vtx3.k);
    fprintf(file, "%f\t%f\t%f\n", vtx4.i, vtx4.j, vtx4.k);
}

Lpoint* polarPoints(Lpoint center, float distance, int angleIdx, double *sins, double *coss)
{
    double x0 = 0, y0 = 0;
    Lpoint *points = NULL;
    
    points = (Lpoint*)malloc(2 * sizeof(Lpoint));
     
    x0 = coss[angleIdx] * distance;
    y0 = sins[angleIdx] * distance;

    points[0].x = center.x + (x0 + 1.0 * -sins[angleIdx]);
    points[0].y = center.y + (y0 + 1.0 *  coss[angleIdx]);
    points[1].x = center.x + (x0 - 1.0 * -sins[angleIdx]);
    points[1].y = center.y + (y0 - 1.0 *  coss[angleIdx]);
    
    return points;
}

int hasValue(int *array, int size, int value)
{
    int i = 0;

    for(i = 0; i < size; i++)
    {   
        if(array[i] == value)
	{ 
	    return i;
	}
    }

    return -1;
}

int* removeValue(int *array, int size, int value) // Segfaults prone (good for revealing neighboring erros): value to remove isn't in the array
{
    int i = 0, count = 0, *newArray = NULL;

    newArray = (int*)malloc((size - 1) * sizeof(int));
    
    for(i = 0; i < size; i++)
    {
        if(array[i] != value)
	{  
	    newArray[count] = array[i];
            
            count++;
	}
    }
    
    return newArray;
}

double radians(double degrees)
{
    return degrees * M_PI / 180.0;
}

double degrees(double radians)
{
    return radians * 180 / M_PI;
}

double stdev(double *values, int num)   
{
    int i = 0;
    double mean = 0.0, diff = 0.0;
 
    for(i = 0; i < num; i++)
    {
        mean += values[i];
    }
    
    mean /= (double)num;
    
    for(i = 0; i < num; i++)
    {
        diff += pow(values[i] - mean, 2);      
    }
    
    return sqrt(diff / (double)num);  
}

int number(double x)
{
    return !isnan(x) && !isinf(x);
}

void checkNAN(double x)
{
    if(!number(x))
    {
        x = 0.0; 
    }
}


void pointToVec(Vector3D *dest, Lpoint *source)
{
    dest->x = source->x;
    dest->y = source->y;
    dest->z = source->z;
}

//Lpoint *groupPts(Group g)
//{
//    int i = 0;
//    Node n = NULL;
//    Lpoint *points = NULL;
//
//    points = (Lpoint*)malloc(g->size * sizeof(Lpoint));
//    
//    for(n = g->start; n != g->end; n = n->next)
//    {
//        points[i++] = n->elem;
//    }
//    
//    return points;
//}

// Shuffle Vx Vy
void shuffle(double *array, double *array2, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          double t = array[j];
          double t2 = array2[j];
          array[j] = array[i];
          array[i] = t;
          array2[j] = array2[i];
          array2[i] = t2;
        }
    }
}

//// Normalize the intensities of all the segmented groups between 0 and 100
//double* normalizedIntensities(Group *groups, int numGroups)
//{
//    int i = 0;
//    float intensity = 0.0, min = FLT_MAX, max = FLT_MIN, range = 0.0;
//    double *normedIntensities = NULL;
//    
//    normedIntensities = (double*)malloc(numGroups * sizeof(double));
//    
//    for(i = 0; i < numGroups; i++)
//    {
//        if(groups[i]->size > 0)
//        {               
//            intensity = avgI(groups[i]); 
//
//            if(intensity > max)
//            {
//                max = intensity;
//            }
//            if(intensity < min)
//            {
//                min = intensity;
//            }
//        }
//    }
//    
//    range = max - min; 
//    
//    for(i = 0; i < numGroups; i++)
//    {
//        intensity = avgI(groups[i]);
//        normedIntensities[i] = ((intensity - min) / range) * 100; 
//    } 
//    
//    //printf("Min: %f Max: %f Range: %f\n", min, max, range);
//    
//    return normedIntensities;
//}

// Normalize the intensities of all the segmented groups between 0 and 100
//double* normalizedSizes(Group *groups, int numGroups)
//{
//    int i = 0;
//    float intensity = 0.0, min = FLT_MAX, max = FLT_MIN, range = 0.0;
//    double *normedIntensities = NULL;
//    
//    normedIntensities = (double*)malloc(numGroups * sizeof(double));
//    
//    for(i = 0; i < numGroups; i++)
//    {
//        if(groups[i]->size > 0)
//        {               
//            intensity = groups[i]->size; 
//
//            if(intensity > max)
//            {
//                max = intensity;
//            }
//            if(intensity < min)
//            {
//                min = intensity;
//            }
//        }
//    }
//    
//    range = max - min; 
//    
//    for(i = 0; i < numGroups; i++)
//    {
//        normedIntensities[i] = ((groups[i]->size - min) / range) * 100; 
//    } 
//    
//    //printf("Min: %f Max: %f Range: %f\n", min, max, range);
//    
//    return normedIntensities;
//}


// Counter-clockwise situation of 3 points (ccw > 0, cw < 0, colinear = 0)
double ccw(Lpoint* p1, Lpoint* p2, Lpoint* p3)
{
  return (p2->x - p1->x) * (p3->y - p1->y) - (p2->y - p1->y) * (p3->x - p1->x);
}

// Check if a point is inside the ccw convex hull
int inHull(Lpoint *point, Lpoint* hull, int num)
{
    int i = 0;
    
    for(i = 0; i < num - 1; i++)
    {
        if(ccw(&hull[i], &hull[i+1], point) < 0)
        {
            return 0;
        }
    }
    
    return 1;
}

int inHull2(Lpoint *point, Lpoint* hull, int num)
{
    int i = 0;
    
    for(i = 0; i < num - 1; i++)
    {
        if(ccw(&hull[i], &hull[i+1], point) > 0)
        {
            return 0;
        }
    }
    
    return 1;
}

// Monotone Chain Convex Hull (A. M. Andrew, 1979). O(n) when points sorted O(n log n) otherwise
void convexHull(Lpoint* points, int num, Lpoint** hull, int* hullSize)
{
    int i = 0, t = 0, k = 0;
    Lpoint* myHull  = NULL;

    myHull = *hull;

    for(i = 0; i < num; ++i) 
    {
        while(k >= 2 && ccw(&myHull[k-2], &myHull[k-1], &points[i]) <= 0)
        {
            --k;
        }

        myHull[k++] = points[i];
    }


    for(i = num-2, t = k+1; i >= 0; --i)  
    {
        while(k >= t && ccw(&myHull[k-2], &myHull[k-1], &points[i]) <= 0)
        {
            --k;
        }
        
        myHull[k++] = points[i];
    }

    *hull = myHull;
    *hullSize = k -1; // -1 to avoid last = first
   // *hull = realloc(*hull, *hullSize * sizeof(Lpoint)); // TODO: Iterative alloc instead of cut hull size at the end
}


// Smith and Oliver method (1961). WARNING: Modifies input matrix!
double* calcEigenvalues(double mat[3][3])
{
    int i = 0, j = 0;
    double *eigenvalues = NULL, trace = 0.0, p = 0.0, p1 = 0.0, p2 = 0.0, B[3][3], phi = 0.0, q = 0.0, r = 0.0, determinant = 0.0;

    eigenvalues = (double*)malloc(3 * sizeof(double));
    p1 = pow(mat[0][1], 2) + pow(mat[0][2], 2) + pow(mat[1][2], 2); 
    
    if (p1 == 0) // Diagonal mat
    {
	eigenvalues[0] = mat[0][0];
  	eigenvalues[1] = mat[1][1];
   	eigenvalues[2] = mat[2][2];
    }
    else
    {
	trace = mat[0][0] + mat[1][1] + mat[2][2]; 
   	q = trace / 3; 
   	p2 = pow(mat[0][0] - q, 2) + pow(mat[1][1] - q, 2) + pow(mat[2][2] - q, 2) + 2 * p1; 
   	p = sqrt(p2 / 6); 
	mat[0][0] -= q;  // A - Iq
        mat[1][1] -= q; 
        mat[2][2] -= q; 
   	
	for(i = 0; i < 3; i++) // B = (1 / p) * (A)  
	{
	    for(j = 0; j < 3; j++)
	    {
		B[i][j] = 1 / p * mat[i][j];
	    }
	}  
	
 	determinant = B[0][0] * (B[1][1] * B[2][2] - B[2][1] * B[1][2]) 
	      	   -  B[0][1] * (B[1][0] * B[2][2] - B[2][0] * B[1][2]) 
              	   +  B[0][2] * (B[1][0] * B[2][1] - B[2][0] * B[1][1]);

   	r = determinant / 2;  
   	
   	if(r <= -1) // Treat computation error
	{
      	    phi = M_PI / 3;
	}
   	else if(r >= 1)
	{
      	    phi = 0;
	}
   	else
	{
      	    phi = acos(r) / 3;
   	}

   	eigenvalues[0] = q + 2 * p * cos(phi);
   	eigenvalues[2] = q + 2 * p * cos(phi + ( 2 * M_PI / 3));
   	eigenvalues[1] = 3 * q - eigenvalues[0] - eigenvalues[2]; 
    }

    return eigenvalues;
}

// Eigen decomposition for symmetric 3x3 matrices (public domain Java Matrix library JAMA) //
#define MAX(a, b) ((a)>(b)?(a):(b))
#define n 3

static double hypot2(double x, double y) {
    return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form (Bowdler et al, Fortran subroutine EISPACK)
static void tred2(double V[n][n], double d[n], double e[n]) {
    
    int i,j,k;
    for (j = 0; j < n; j++) {
       d[j] = V[n-1][j];
     }

     // Householder reduction to tridiagonal form.
     for (i = n-1; i > 0; i--) {
       // Scale to avoid under/overflow.   
       double scale = 0.0;
       double h = 0.0;
       for (k = 0; k < i; k++) {
         scale = scale + fabs(d[k]);
       }
       if (scale == 0.0) {
         e[i] = d[i-1];
         for (j = 0; j < i; j++) {
           d[j] = V[i-1][j];
           V[i][j] = 0.0;
           V[j][i] = 0.0;
         }
       } else {

         // Generate Householder vector.
         for (k = 0; k < i; k++) {
           d[k] /= scale;
           h += d[k] * d[k];
         }
         double f = d[i-1];
         double g = sqrt(h);
         if (f > 0) {
           g = -g;
         }
         e[i] = scale * g;
         h = h - f * g;
         d[i-1] = f - g;
         for (j = 0; j < i; j++) {
           e[j] = 0.0;
         }

         // Apply similarity transformation to remaining columns.
         for (j = 0; j < i; j++) {
           f = d[j];
           V[j][i] = f;
           g = e[j] + V[j][j] * f;
           for (k = j+1; k <= i-1; k++) {
             g += V[k][j] * d[k];
             e[k] += V[k][j] * f;
           }
           e[j] = g;
         }
         f = 0.0;
         for (j = 0; j < i; j++) {
           e[j] /= h;
           f += e[j] * d[j];
         }
         double hh = f / (h + h);
         for (j = 0; j < i; j++) {
           e[j] -= hh * d[j];
         }
         for (j = 0; j < i; j++) {
           f = d[j];
           g = e[j];
           for (k = j; k <= i-1; k++) {
             V[k][j] -= (f * e[k] + g * d[k]);
           }
           d[j] = V[i-1][j];
           V[i][j] = 0.0;
         }
       }
       d[i] = h;
     }

     // Accumulate transformations.
     for (i = 0; i < n-1; i++) {
       V[n-1][i] = V[i][i];
       V[i][i] = 1.0;
       double h = d[i+1];
       if (h != 0.0) {
         for (k = 0; k <= i; k++) {
           d[k] = V[k][i+1] / h;
         }
         for (j = 0; j <= i; j++) {
           double g = 0.0;
           for (k = 0; k <= i; k++) {
             g += V[k][i+1] * V[k][j];
           }
           for (k = 0; k <= i; k++) {
             V[k][j] -= g * d[k];
           }
         }
       }
       for (k = 0; k <= i; k++) {
         V[k][i+1] = 0.0;
       }
     }
     for (j = 0; j < n; j++) {
       d[j] = V[n-1][j];
       V[n-1][j] = 0.0;
     }
     V[n-1][n-1] = 1.0;
     e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm (Bowdler et al, Fortran subroutine EISPACK)
static void tql2(double V[n][n], double d[n], double e[n]) {

    int i = 0, j = 0, k = 0, l =0;

    for (i = 1; i < n; i++) {
      e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    for (l = 0; l < n; l++) {

      // Find small subdiagonal element
      tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
      int m = l;
      while (m < n) {
        if (fabs(e[m]) <= eps*tst1) {
          break;
        }
        m++;
      }

      // If m == l, d[l] is an eigenvalue, otherwise, iterate.
      if (m > l) {
        int iter = 0;
        do {
          iter = iter + 1;  // (Could check iteration count here.)

          // Compute implicit shift
          double g = d[l];
          double p = (d[l+1] - g) / (2.0 * e[l]);
          double r = hypot2(p,1.0);
          if (p < 0) {
            r = -r;
          }
          d[l] = e[l] / (p + r);
          d[l+1] = e[l] * (p + r);
          double dl1 = d[l+1];
          double h = g - d[l];
          for (i = l+2; i < n; i++) {
            d[i] -= h;
          }
          f = f + h;

          // Implicit QL transformation.
          p = d[m];
          double c = 1.0;
          double c2 = c;
          double c3 = c;
          double el1 = e[l+1];
          double s = 0.0;
          double s2 = 0.0;
          for (i = m-1; i >= l; i--) {
            c3 = c2;
            c2 = c;
            s2 = s;
            g = c * e[i];
            h = c * p;
            r = hypot2(p,e[i]);
            e[i+1] = s * r;
            s = e[i] / r;
            c = p / r;
            p = c * d[i] - s * g;
            d[i+1] = h + s * (c * g + s * d[i]);

            // Accumulate transformation.
            for (k = 0; k < n; k++) {
              h = V[k][i+1];
              V[k][i+1] = s * V[k][i] + c * h;
              V[k][i] = c * V[k][i] - s * h;
            }
          }
          p = -s * s2 * c3 * el1 * e[l] / dl1;
          e[l] = s * p;
          d[l] = c * p;

          // Check for convergence.
        } while (fabs(e[l]) > eps*tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.
    for (i = 0; i < n-1; i++) {
      int k = i;
      double p = d[i];
      for (j = i+1; j < n; j++) {
        if (d[j] < p) {
          k = j;
          p = d[j];
        }
      }
      if (k != i) {
        d[k] = d[i];
        d[i] = p;
        for (j = 0; j < n; j++) {
          p = V[j][i];
          V[j][i] = V[j][k];
          V[j][k] = p;
        }
      }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////

void eigenDecomposition(double A[n][n], double V[n][n], double d[n]) 
{
    double e[n];
    int i = 0, j = 0;

    for (i = 0; i < n; i++) 
    {
        for (j = 0; j < n; j++)
        {
          V[i][j] = A[i][j];
        }
    } 

    tred2(V, d, e); // Reduces a real symmetric matrix to a symmetric tridiagonal matrix using and accumulating orthogonal similarity transformations.
    tql2(V, d, e);  // Finds the eigenvalues and eigenvectors of a symmetric tridiagonal matrix by the ql method.
}

#undef n

// Calculate the normal of a point from his eigenvectors
Lvector eigenNormal(Lpoint *point, Octree octree, double radius)
{
    int numNeighs = 0;
    double dot = 0, eigenvalues[3], covarMat[3][3], eigenvectors[3][3];
    Lpoint **neighbors = NULL, viewpoint;
    Lvector normal, dir;
    
    zeroTriMat(3, covarMat);                    
    neighbors = searchNeighbors3D(point, octree, radius, &numNeighs);  // TODO: KNN?
    
    if(numNeighs > 0)
    {
        covarianceMat(neighbors, numNeighs, covarMat); 
        eigenDecomposition(covarMat, eigenvectors, eigenvalues); 
        normal.i = eigenvectors[0][0];   // Normal = Eigenvector of smallest eigenvalue
        normal.j = eigenvectors[1][0];   
        normal.k = eigenvectors[2][0];            
        checkVecNAN(&normal);     
        free(neighbors);
    }   
    
    // Consistent normal orientations
    copyCoords(&viewpoint, point);
    viewpoint.z += 1000;                // Airborne LiDAR only (2.5D not real 3D)
    dir = calcVector(point, &viewpoint);
    dot = dotProduct(&normal, &dir);
    if(dot < 0) multVector(&normal, -1); // Flip normal vector
    
    return normal;
}

// Calculate the normal of a point from his eigensystem
Lvector eigenNormalInt(Lpoint *point, Octree octree, double radius)
{
    int numNeighs = 0;
    double aux = 0, dot = 0, eigenvalues[3], covarMat[3][3], eigenvectors[3][3];
    Lpoint **neighbors = NULL, viewpoint;
    Lvector normal, dir;
    
    zeroTriMat(3, covarMat);                     
    neighbors = searchNeighbors2D(point, octree, radius, &numNeighs);  // TODO: KNN?
    
    if(numNeighs > 0)
    {
        covarianceMatInt(neighbors, numNeighs, covarMat); 
        eigenDecomposition(covarMat, eigenvectors, eigenvalues); 
        normal.i = eigenvectors[0][0];   // Normal = Eigenvector of smallest eigenvalue
        normal.j = eigenvectors[1][0];   
        normal.k = eigenvectors[2][0];            
        checkVecNAN(&normal);     
        free(neighbors);
    }  
    
/*
    // Consistent normal orientations
    aux = point->z; 
    point->z = point->intensity;
    copyCoords(&viewpoint, point);
    viewpoint.z += 10;                // Airborne LiDAR only (2.5D not real 3D)
    dir = calcVector(point, &viewpoint);
    dot = dotProduct(&normal, &dir); //printf("Pi %f Vi %f dot %lf\n", point->z, viewpoint.z, dot);
    if(dot < 0) multVector(&normal, -1); // Flip normal vector
    point->z = aux;
*/
     
    return normal;
}


int cmpint (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

void printIntArray(int *array, int n)
{
    int i = 0;
    
    for(i = 0; i < n; i++)
    {
        printf("%d\t", array[i]);
    }
    
    printf("\n");
}

int removeDuplicates(int *v, int n)
{    
    int i = 0, j = 0;

    qsort(v, n, sizeof(int), cmpint); 

    for(i = 1; i < n; i++)
    {
        if(v[i] != v[j])
        {
            j++;
            v[j] = v[i]; 
        }
    }

    n = (j + 1);

    return(j + 1);
}

// Check if point is inside the polygon with ray casting method
int pnpoly(int n, Lpoint *pol, Lpoint *p)
{
  int i, j, c = 0;
  
  for(i = 0, j = n-1; i < n; j = i++)   
    if(((pol[i].y>p->y) != (pol[j].y>p->y)) && (p->x < (pol[j].x-pol[i].x) * (p->y-pol[i].y) / (pol[j].y-pol[i].y) + pol[i].x))
       c = !c;
  
  return c;
}

int diffOnRange(double a, double b, double factor)
{
    double avg = 0, diff = 0, maxDiff = 0;
    
    avg = (a + b) / 2;
    maxDiff = avg * factor;
    diff = fabs(a - b);

    return diff < maxDiff;
}
