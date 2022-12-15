#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "octree.h"
#include "file.h"
#include "util.h"

#define MIN_RADIUS 1.5  // Stop dividing the octree when childs have this radius
//#define MIN_RADIUS 0.10 //For the discretization
#define MAX_POINTS 100  // Stop dividing the octree when childs have this number of points
//#define MAX_POINTS 10000  // Stop dividing the octree when childs have this number of points

// Calculate the radius of the bounding box
Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius)
{
    Vector3D radius;
     
    radius.x = (max.x - min.x) / 2;
    radius.y = (max.y - min.y) / 2;
    radius.z = (max.z - min.z) / 2;
    
    if(radius.x > radius.y && radius.x > radius.z)  
    {
        *maxRadius = radius.x;
    }
    else if(radius.y > radius.x && radius.y > radius.z)
    {        
        *maxRadius = radius.y;
    }
    else       
    {
        *maxRadius = radius.z;
    }   
    
    return radius;
}

// Calculate the center of the bounding box
Vector3D getCenter(Vector3D min, Vector3D radius)
{
    Vector3D center;
    
    center.x = min.x + radius.x;
    center.y = min.y + radius.y;
    center.z = min.z + radius.z;
    
    return center;
}

// Examine the points to build the bounding box
Vector3D boundingBox(FILE *file, long numPoints, float *maxRadius)
{
    long i = 0;
    Vector3D center, min, max, radius;
    double x = 0, y = 0, z = 0;
    char ignoreLine[1024]; 

    min.x = DBL_MAX, min.y = DBL_MAX, min.z = DBL_MAX, max.x = -DBL_MAX, max.y = -DBL_MAX, max.z = -DBL_MAX;    
    fgets(ignoreLine, sizeof(ignoreLine), file);
    
    for(i = 0; i < numPoints; i++)
    {
        fscanf(file, "%lf\t%lf\t%lf\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\n", &x, &y, &z);
        if(x < min.x) min.x = x; if(x > max.x) max.x = x;
        if(y < min.y) min.y = y; if(y > max.y) max.y = y;
        if(z < min.z) min.z = z; if(z > max.z) max.z = z;
    }
 
    radius = getRadius(min, max, maxRadius);
    center = getCenter(min, radius);
    
    return center;
}
    
// Create a octree with the given center and radius
Octree createOctree(Vector3D center, float radius)
{
    int i = 0;
    Octree octree = NULL;
    
    //if(file != NULL) fprintf(file, "%.2f\t%.2f\t%.2f\t%.2f\n", center.x, center.y, center.z, radius);
    
    octree = malloc(sizeof(struct Octree));    
    octree->center = center;
    octree->radius = radius;
    octree->points = NULL;
    octree->numPoints = 0;
    octree->removed = 0;
    
    for(i = 0; i < 8; i++)
    {
        octree->childs[i] = NULL;
    }
    
    return octree;
}

// Find the child corresponding a given point
int getChild(Lpoint *point, Octree octree)
{
    int child = 0;
    
    if(point->x >= octree->center.x) child |= 4;
    if(point->y >= octree->center.y) child |= 2;
    if(point->z >= octree->center.z) child |= 1;

    return child;
}

int isLeaf(Octree octree)
{
    return octree->childs[0] == NULL;
}

int isEmpty(Octree octree)
{
    return octree->numPoints == 0;
}

// Insert a point in the octree creating the appropiate childs
void insertPoint(Lpoint *point, Octree octree)
{
    int i = 0;
    
    if(isLeaf(octree))
    {
        if(isEmpty(octree))             // Empty leaf -> insert point            
        {
            octree->points = malloc(sizeof(Lpoint*));
            octree->points[0] = point;
            octree->numPoints = 1;  
        }
        else                            // Not empty but still divisible -> divide     
        {
            if(octree->radius / 2.0 > MIN_RADIUS) 
            //if(octree->numPoints > MAX_POINTS)
            {          
                for(i = 0; i < 8; i++)
                {
                    Vector3D newCenter = octree->center;
                    newCenter.x += octree->radius * (i&4 ? 0.5f : -0.5f);
                    newCenter.y += octree->radius * (i&2 ? 0.5f : -0.5f);
                    newCenter.z += octree->radius * (i&1 ? 0.5f : -0.5f);
                    octree->childs[i] = createOctree(newCenter,  octree->radius * 0.5);
                }

                i = getChild(octree->points[0], octree);         
                octree->childs[i]->points = octree->points;  
                octree->childs[i]->numPoints++;
                octree->numPoints = 0;
                octree->points = NULL;                
                i = getChild(point, octree);
                insertPoint(point, octree->childs[i]);
            }
            else                         // Not empty and isn't divisible -> insert point
            {
                octree->points = realloc(octree->points, ++octree->numPoints * sizeof(Lpoint*));
                octree->points[octree->numPoints-1] = point; 
            }
        }
    }
    else                                // No leaf -> search the correct one          
    {
        i = getChild(point, octree);
        insertPoint(point, octree->childs[i]);
    }
}

int insideBox2D(Lpoint *point, Vector3D min, Vector3D max)
{   
    if(point->x > min.x && point->y > min.y)
    {
        if(point->x < max.x && point->y < max.y)
        {  
            return 1;
        }
    }

    return 0;
}

int insideBox3D(Lpoint *point, Vector3D min, Vector3D max)
{   
    if(point->x > min.x && point->y > min.y  && point->z > min.z)
    {
        if(point->x < max.x && point->y < max.y && point->z < max.z)
        {  
            return 1;
        }
    }

    return 0;
}

Vector3D subVector(Vector3D v1, float number)
{
    Vector3D result;
    
    result.x = v1.x - number;
    result.y = v1.y - number;
    result.z = v1.z - number;
    
    return result;  
}

Vector3D addVector(Vector3D v1, float number)
{
    Vector3D result;
    
    result.x = v1.x + number;
    result.y = v1.y + number;
    result.z = v1.z + number;
    
    return result;  
}

int childTouches3D(Vector3D boxMin, Vector3D boxMax, Octree octree)
{
    Vector3D min = subVector(octree->center, octree->radius);
    Vector3D max = addVector(octree->center, octree->radius);
    
    if(max.x < boxMin.x || max.y < boxMin.y || max.z < boxMin.z)  return 0;
    if(min.x > boxMax.x || min.y > boxMax.y || min.z > boxMax.z)  return 0;
    
    return 1;
}

int childTouches2D(Vector3D boxMin, Vector3D boxMax, Octree octree)
{
    Vector3D min = subVector(octree->center, octree->radius);
    Vector3D max = addVector(octree->center, octree->radius);
    
    if(max.x < boxMin.x || max.y < boxMin.y)  return 0;
    if(min.x > boxMax.x || min.y > boxMax.y)  return 0;
  
    return 1;
}

Lpoint** neighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {                    
            for(i = 0; i < octree->numPoints; i++)
            {   
                if(insideBox2D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id) 
                {   
                    ptsInside = realloc(ptsInside, ++*numInside * sizeof(Lpoint*));
                    ptsInside[*numInside-1] = octree->points[i];

                }
            }
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            if(!childTouches2D(boxMin, boxMax, octree->childs[i]))
            {
                continue;
            }
            else
            {           
                ptsInside = neighbors2D(point, boxMin, boxMax, octree->childs[i], ptsInside, numInside);
            }
        }
    }
    
    return ptsInside;
}

Lpoint** neighbors2DReset(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {                    
            for(i = 0; i < octree->numPoints; i++)
            {   
                if(insideBox2D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id) 
                {   
                    ptsInside = realloc(ptsInside, ++*numInside * sizeof(Lpoint*));
                    octree->points[i]->id  = *numInside - 1;
                    ptsInside[*numInside-1] = octree->points[i];
                }
            }
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            if(!childTouches2D(boxMin, boxMax, octree->childs[i]))
            {
                continue;
            }
            else
            {           
                ptsInside = neighbors2DReset(point, boxMin, boxMax, octree->childs[i], ptsInside, numInside);
            }
        }
    }
    
    return ptsInside;
}

Lpoint** neighbors3D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {                    
            for(i = 0; i < octree->numPoints; i++)
            {
                if(insideBox3D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id) 
                {
                    ptsInside = realloc(ptsInside, ++*numInside * sizeof(Lpoint*));
                    ptsInside[*numInside-1] = octree->points[i];

                }
            }
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            if(!childTouches3D(boxMin, boxMax, octree->childs[i]))
            {
                continue;
            }
            else
            {           
                ptsInside = neighbors3D(point, boxMin, boxMax, octree->childs[i], ptsInside, numInside);
            }
        }
    }
    
    return ptsInside;
}

Lpoint** circleNeighbors(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *numInside, float circleRadius)
{
    int i = 0;

    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {                    
            for(i = 0; i < octree->numPoints; i++)
            {
                if(insideCircle(octree->points[i], point, circleRadius) & point->id != octree->points[i]->id) 
                {
                    ptsInside = realloc(ptsInside, ++*numInside * sizeof(Lpoint*));
                    ptsInside[*numInside-1] = octree->points[i];

                }
            }
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            if(!childTouches3D(boxMin, boxMax, octree->childs[i]))
            {
                continue;
            }
            else
            {           
                ptsInside = circleNeighbors(point, boxMin, boxMax, octree->childs[i], ptsInside, numInside, circleRadius);
            }
        }
    }
    
    return ptsInside;
}

int numNeighbors3D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {                    
            for(i = 0; i < octree->numPoints; i++)
            {
                if(insideBox3D(octree->points[i], boxMin, boxMax) && octree->points[i]->intensity != -1.0) 
                {
                    ++*numInside;         
                }
            }
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            if(!childTouches3D(boxMin, boxMax, octree->childs[i]))
            {
                continue;
            }
            else
            {           
                *numInside = numNeighbors3D(point, boxMin, boxMax, octree->childs[i], numInside);
            }
        }
    }
    
    return *numInside;
}

int numNeighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {                    
            for(i = 0; i < octree->numPoints; i++)
            {
                if(insideBox2D(octree->points[i], boxMin, boxMax)) 
                {
                    ++*numInside;         
                }
            }
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            if(!childTouches2D(boxMin, boxMax, octree->childs[i]))
            {
                continue;
            }
            else
            {           
                *numInside = numNeighbors2D(point, boxMin, boxMax, octree->childs[i], numInside);
            }
        }
    }
    
    return *numInside;
}


// Return also the octrees containing inside points
Lpoint** neighborsNBH(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint** ptsInside, Octree **octrees, int *numInside)
{
    int i = 0;
    
    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {                    
            for(i = 0; i < octree->numPoints; i++)
            {
                if(insideBox3D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id && octree->points[i]->intensity != -1) 
                {
                    *numInside += 1;
                    ptsInside = realloc(ptsInside, *numInside * sizeof(Lpoint*));
                    ptsInside[*numInside-1] = octree->points[i];
                    *octrees = realloc(*octrees, *numInside * sizeof(Octree));
                    (*octrees)[*numInside-1] = octree;
                }
            }
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            if(!childTouches3D(boxMin, boxMax, octree->childs[i]))
            {
                continue;
            }
            else
            {           
                ptsInside = neighborsNBH(point, boxMin, boxMax, octree->childs[i], ptsInside, octrees, numInside);
            }
        }
    }
       
    return ptsInside;
}

// Make a box with center the point and the specified radius
void makeBox(Lpoint *point, float radius, Vector3D *min, Vector3D *max)
{
    min->x = point->x - radius;
    min->y = point->y - radius;
    min->z = point->z - radius;
    max->x = point->x + radius;
    max->y = point->y + radius;
    max->z = point->z + radius;  
}

// Make the box that represents an octree
void ownBox(Octree octree, Vector3D *min, Vector3D *max)
{
    min->x = octree->center.x - octree->radius;
    min->y = octree->center.y - octree->radius;
    min->z = octree->center.z - octree->radius;
    max->x = octree->center.x + octree->radius;
    max->y = octree->center.y + octree->radius;
    max->z = octree->center.z + octree->radius;  
}

int numOctreeNeighbors(Lpoint *point, Octree octree, float radius)
{
    Vector3D boxMin, boxMax;
    int numInside = 0;

    makeBox(point, radius, &boxMin, &boxMax); 
    numInside = numNeighbors3D(point, boxMin, boxMax, octree, &numInside);
    
    return numInside; 
}

int numOctreeNeighbors2D(Lpoint *point, Octree octree, float radius)
{
    Vector3D boxMin, boxMax;
    int numInside = 0;

    makeBox(point, radius, &boxMin, &boxMax); 
    numInside = numNeighbors2D(point, boxMin, boxMax, octree, &numInside);
    
    return numInside; 
}

Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    
    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax); 
    ptsInside = neighbors2D(point, boxMin, boxMax, octree, ptsInside, numNeighs);
   
    return ptsInside;
}

Lpoint** searchNeighbors2DReset(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    
    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax); 
    ptsInside = neighbors2DReset(point, boxMin, boxMax, octree, ptsInside, numNeighs);
   
    return ptsInside;
}

Lpoint** searchNeighbors3D(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    
    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax); 
    ptsInside = neighbors3D(point, boxMin, boxMax, octree, ptsInside, numNeighs);
   
    return ptsInside;
}

Lpoint** searchCircleNeighbors(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int numInside = 0;

    makeBox(point, radius, &boxMin, &boxMax); 
    ptsInside = circleNeighbors(point, boxMin, boxMax, octree, ptsInside, &numInside, radius);
    *numNeighs = numInside;
    
    return ptsInside;
}

Neighborhood searchNeighborsNBH(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    Octree *octrees = NULL;
    Neighborhood neigh;
    int numInside = 0;

    makeBox(point, radius, &boxMin, &boxMax); 
    ptsInside = neighborsNBH(point, boxMin, boxMax, octree, ptsInside, &octrees, &numInside);
    *numNeighs = numInside;
    neigh.points = ptsInside;
    neigh.octrees = octrees;
    
    return neigh;
}

void removeOcpoint(Lpoint *point, Octree octree)
{
    Lpoint aux, *last = NULL;
   
    last = octree->points[octree->numPoints-1];
    aux = *point;
    *point = *last;
    *last = aux;
    octree->numPoints--;
   
    octree->removed++;
}

void recoverOctreePoint(Lpoint *point, Lpoint *original, Octree octree)
{
    Lpoint aux;
    
    aux = *point;
    *point = *original;
    *original = aux;
    octree->numPoints++;
}

Octree searchOctree(Lpoint *point, Octree octree, Octree *leaf)
{
    int i = 0;
    Vector3D boxMin, boxMax;

    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {   
            for(i = 0; i < octree->numPoints; i++)
            {
                if(octree->points[i]->id == point->id) 
                {
                    *leaf = octree;     
                }
            }  
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            ownBox(octree->childs[i], &boxMin, &boxMax); 
            
            if(!insideBox3D(point, boxMin, boxMax))
            {
                continue;
            }
            else
            {           
                *leaf = searchOctree(point, octree->childs[i], leaf);
            }
        }
    }
    
    return *leaf;
}

void printOctreeCenters(Octree octree, FILE *file)
{
    int i = 0;
    
    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {   
            fprintf(file, "%.2f\t%.2f\t%.2f\t%.2f\n", octree->center.x, octree->center.y, octree->center.z, octree->radius);
        }
    }
    else
    {      
        fprintf(file, "%.2f\t%.2f\t%.2f\t%.2f\n", octree->center.x, octree->center.y, octree->center.z, octree->radius);
        
        for(i = 0; i < 8; i++)
        {
            printOctreeCenters(octree->childs[i], file);
        }
    }
}

void printOctreePoints(Octree octree, FILE *file)
{
    int i = 0;
    
    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {   
            for(i = 0; i < octree->numPoints; i++)
            {
                if(file == NULL) printf("%d\t%.2f\t%.2f\t%.2f\n", octree->points[i]->id, octree->points[i]->x, octree->points[i]->y, octree->points[i]->z);
                else fprintf(file, "%d\t%.2f\t%.2f\t%.2f\n", octree->points[i]->id, octree->points[i]->x, octree->points[i]->y, octree->points[i]->z);
            }
        }
    }
    else
    {      
        for(i = 0; i < 8; i++)
        {
            printOctreePoints(octree->childs[i], file);
        }
    }
}

int numOctreePoints(Octree octree, int *num)
{
    int i = 0;
    
    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {   
            for(i = 0; i < octree->numPoints; i++)
            {
                ++*num;
            }
        }
    }
    else
    {      
        for(i = 0; i < 8; i++)
        {
            *num = numOctreePoints(octree->childs[i], num);
        }
    }
    
    return *num;
}

void searchMinMaxOct(Octree octree, double *min, double *max)
{
    int i = 0;
    
    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {   
            for(i = 0; i < octree->numPoints; i++)
            {  
                if(octree->points[i]->intensity > *max)
                {
                    *max = octree->points[i]->intensity;
                }
                if(octree->points[i]->intensity < *min)
                {
                    *min = octree->points[i]->intensity;
                }   
            }
        }
    }
    else
    {              
        for(i = 0; i < 8; i++)
        {
            searchMinMaxOct(octree->childs[i], min ,max);
        }
    }
}

void normalizeIntOct(Octree octree, double min, double max)
{
    int i = 0;
    double range = 0;
    
    range = max - min; //printf("Range %.2lf Min %.2lf Max %.2lf\n", range, min, max);
    
    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {   
           for(i = 0; i < octree->numPoints; i++)
           {
               octree->points[i]->intensity = ((octree->points[i]->intensity - min) / range); 
           }
        }
    }
    else
    {              
        for(i = 0; i < 8; i++)
        {
            normalizeIntOct(octree->childs[i], min ,max);
        }
    }
}

void normalizeIntensity(Octree octree)
{
    double min = DBL_MAX, max = -DBL_MAX;
    
    searchMinMaxOct(octree, &min, &max); printf("Min %f max %f\n", min, max);
    normalizeIntOct(octree, min, max);
}

//
//void octreeDiscretization(Octree octree, FILE *file)
//{
//    int i = 0;
//
//    if(isLeaf(octree))
//    {
//        if(!isEmpty(octree))
//        {
//            writePointRaw(file, octree->points[0]);
//        }
//    }
//    else
//    {
//        for(i = 0; i < 8; i++)
//        {
//            octreeDiscretization(octree->octants[i], file);
//        }
//    }
//}

/*
 * OSCAR
*/

Lpoint** donutNeighbors(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *numInside, float extRadius, float intRadius)
{
    int i = 0;

    if(isLeaf(octree))
    { 
        if(!isEmpty(octree))
        {                    
            for(i = 0; i < octree->numPoints; i++)
            {
                if(insideCircle(octree->points[i], point, extRadius) & point->id != octree->points[i]->id) 
                {
			if(insideCircle(octree->points[i], point, intRadius) & point->id != octree->points[i]->id) 
			{
				//do nothing
			}else{
	                    ptsInside = realloc(ptsInside, ++*numInside * sizeof(Lpoint*));
        	            ptsInside[*numInside-1] = octree->points[i];
			}
                }
            }
        }
    }
    else
    {            
        for(i = 0; i < 8; i++)
        {
            if(!childTouches3D(boxMin, boxMax, octree->childs[i]))
            {
                continue;
            }
            else
            {           
                ptsInside = donutNeighbors(point, boxMin, boxMax, octree->childs[i], ptsInside, numInside, extRadius, intRadius);
            }
        }
    }
    
    return ptsInside;
}

Lpoint** searchDonutNeighbors(Lpoint *point, Octree octree, float extRadius, float intRadius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int numInside = 0;

    makeBox(point, extRadius, &boxMin, &boxMax); 
    ptsInside = donutNeighbors(point, boxMin, boxMax, octree, ptsInside, &numInside, extRadius, intRadius);
    *numNeighs = numInside;
    
    return ptsInside;
}

