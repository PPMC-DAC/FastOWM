#ifndef OCTREE_H

#include <stdio.h>
#include "point.h"

#define	OCTREE_H

typedef struct 
{
    double x;
    double y;
    double z;
    
} Vector3D;


typedef struct Octree *Octree;

struct Octree
{   
    Octree childs[8];
    Vector3D center;
    Lpoint **points;
    float radius;
    int numPoints; 
    int removed;
    
}; // A tree with eight childs

typedef struct 
{
    Lpoint **points;
    Octree *octrees;   
    
} Neighborhood;


Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

Vector3D boundingBox(FILE *file, long numPoints, float *radius); 

Octree createOctree(Vector3D center, float radius);

void insertPoint(Lpoint *point, Octree octree);


int numOctreeNeighbors(Lpoint *point, Octree octree, float radius);

Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchNeighbors2DReset(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchNeighbors3D(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchCircleNeighbors(Lpoint *point, Octree octree, float radius, int *numNeighs);

Neighborhood searchNeighborsNBH(Lpoint *point, Octree octree, float radius, int *numNeighs);

void removeOcpoint(Lpoint *point, Octree octree);

void recoverOctreePoint(Lpoint *point, Lpoint *original, Octree octree);

Octree searchOctree(Lpoint *point, Octree octree, Octree *leaf);


void printOctreeCenters(Octree octree, FILE *file);

void printOctreePoints(Octree octree, FILE *file);


void makeBox(Lpoint *point, float radius, Vector3D *min, Vector3D *max);

int insideBox2D(Lpoint *point, Vector3D min, Vector3D max);

int numOctreePoints(Octree octree, int *num);

void normalizeIntensity(Octree octree);

int numOctreeNeighbors2D(Lpoint *point, Octree octree, float radius);

Lpoint** searchNeighbors2D_withcutoff(Lpoint *point, Octree octree, float radius, int *numNeighs, float cutoff, double *plane);

int isLeaf(Octree octree);

int isEmpty(Octree octree);

//void octreeDiscretization(Octree octree, FILE *file);

/*
 * OSCAR
*/
Lpoint** searchDonutNeighbors(Lpoint *point, Octree octree, float extRadius, float intRadius, int *numNeighs);

#endif  /* OCTREE_H */

