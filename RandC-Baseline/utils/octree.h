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

/**
 * @struct Octree
 * A tree with eight childs
 */
struct Octree
{
    Octree octants[8];
    Vector3D center;
    Lpoint **points;
    int numPts;
    float radius;
    // double plane[5];
};

typedef struct
{
    Lpoint **points;
    Octree *octrees;

} Neighborhood; // Save for each point pointer also the octree where it is

int isLeaf(Octree oct);

int isEmpty(Octree oct);

Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

Octree createOctree(Vector3D center, float radius);

void insertPoints(Lpoint *points, int numPts, Octree oct);

Octree buildOctree(Lpoint *points, int numPts, Vector3D *octCenter, float *octRadius, int label);

Octree* searchOctNeighbors(Octree oct, Octree gOctree, int *numNeighs);

Octree* searchOctNeighborsFelipe(Lpoint *center, Octree octree, float radius, int *numNeighs);

void insertPoint(Lpoint *point, Octree octree);

void insertPointMinRadius(Lpoint *point, Octree octree, float minRadius);

void insertPointMaxPoints(Lpoint *point, Octree octree, int maxPoints);

//void insertPoint(Lpoint *point, Octree octree);

int getPointGroup(Lpoint point, Octree octree);

void mergeOctants(Octree dst, Octree src);

int numOctreeNeighbors3D(Lpoint *point, Octree octree, float radius);

Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchNeighbors2DResetId(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchCatNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** catNeighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size, int *numInside);


Lpoint** searchNeighbors2DWrap(Lpoint *point, Octree octree, float radiusX, float radiusY, int *numNeighs);


Lpoint** searchNeighbors2DReset(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchNeighbors3D(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchNeighbors2DNoGroupLowerZ(Lpoint *point, Octree octree, float radius, int gId, int *numNeighs);

Lpoint** searchNeighbors3DNoGroup(Lpoint *point, Octree octree, float radius, int gId, int *numNeighs);

//Lpoint** searchSphereNeighbors(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchSphereNeighborhood(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchNeighborsBorder(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint** searchCircleNeighbors(Lpoint *point, Octree octree, float radius, int *numNeighs);


void hideOctPoint(Lpoint *point, Octree octree);

void hideOctPointByIndex(int index, Octree oct);

void recoverOctPoint(Lpoint *point, Lpoint *original, Octree octree);

Octree searchOctree(Lpoint *point, Octree octree, Octree *leaf);


void printOctreeCenters(Octree octree, FILE *file, int alsoEmpty);


void printOctreePoints(Octree octree, FILE *file);


void makeBox(Lpoint *point, float radius, Vector3D *min, Vector3D *max);

void makeBoxWrap(Lpoint *point, float radiusX, float radiusY, Vector3D *min, Vector3D *max);

int insideBox2D(Lpoint *point, Vector3D min, Vector3D max);

int numOctreePoints(Octree octree, int *num);

void searchMinFelipe(Octree octree, double* min, unsigned int* id);

// void normalizeOctIntensity(Octree octree);

int numOctreeNeighbors2D(Lpoint *point, Octree octree, float radius);

int numOctreeNeighbors2DDiffGroup(Lpoint *point, Octree octree, float radius);

int numNeighbors2DDistincGroup(Lpoint *point, Octree octree, float radius);

Lpoint** searchNeighbors2DGroupDistinct(Lpoint *point, Octree octree, float radius, int *numNeighs);

void createOctants(Octree octree);

void insertPointNoSplit(Lpoint *point, Octree octree);

void writeOctreeNormals(Octree octree, FILE *file);

void populateChilds(Octree octree);

void cpyOctNoDiv(Octree dst, Octree src);

void sortOctsByFeature(Octree *octrees, int num, int* sortIdxs);

Lpoint octCentroid (Octree oct);

void octreeSpacing(Octree oct, int *numPoints, double *area, FILE *file);

void removeOctree(Octree oct);

//Octree buildOctreeUnsplit(Vector3D center, float radius, Lpoint *points, int numPoints, char label);

void gatherLeafs(Octree oct, Octree **fitLeafs, int *numFit, Octree **unFitLeafs, int *numUnFit);

void octPadding();

int insideOct(Lpoint *p, Octree oct);

int insideBox3D(Lpoint *point, Vector3D min, Vector3D max);


/*
 * Landing
*/
Lpoint** searchDonutNeighbors(Lpoint *point, Octree octree, float extRadius, float intRadius, int *numNeighs);

void octreeDiscretization(Octree octree, Octree octreeDisc, int *totalPoints, int *newPoints, float minRadius);

void getLeafMinRadius(Octree octree, float *radius);
void printLeafRadius(Octree octree);

#endif  /* OCTREE_H */
