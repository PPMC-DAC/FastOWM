#ifndef ENVIRONMENT_H

#define ENVIRONMENT_H

typedef struct
{
    double x;
    double y;
    double z;

} Vector3D;

typedef struct
{
    unsigned int id;
    double x;
    double y;
    double z;

} Lpoint;

typedef struct Octree *Octree;

struct Octree
{
    Octree octants[8];
    Vector3D center;
    Lpoint **points;
    int numPts;
    float radius;
};

typedef struct{
    unsigned int id;
    double z;
} zSearch;

unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min, const unsigned chunk);

unsigned int stage1t(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min, const unsigned chunk);

unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

unsigned int stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

// unsigned int findMin(Lpoint** neighbors, double zzmin, unsigned int cellPoints);

double round2d(double z);

unsigned int findMin(Lpoint** neighbors, unsigned int cellPoints);

unsigned int findMin2(Lpoint** neighbors, unsigned int cellPoints);


// Lpoint createPoint(unsigned int id, double x, double y, double z);

Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

Octree createOctree(Vector3D center, float radius);

int isLeaf(Octree oct);

int isEmpty(Octree oct);

void insertPoint(Lpoint *point, Octree octree);

void insertPoint2(Lpoint *point, Octree octree);

void insertPointMinRadius(Lpoint *point, Octree octree, float minRadius);

Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs);


#endif // ENVIRONMENT_H
