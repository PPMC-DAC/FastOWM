#include "../include/envi.h"
#include "tbb/parallel_for.h"
#include <numeric>

double round2d(double z){
  return round(z*100)/100;
}

// Make a box with center the point and the specified radius
void makeBox(Lpoint *point, float radius, Vector3D *min, Vector3D *max)
{
    // printf("Radio: %.2f\n", radius);
    // printf("Centro: [ %.2lf, %.2lf]\n",point->x,point->y );
    min->x = point->x - radius;
    min->y = point->y - radius;
    min->z = point->z - radius;
    max->x = point->x + radius;
    max->y = point->y + radius;
    max->z = point->z + radius;
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

int boxOverlap2D(Vector3D boxMin, Vector3D boxMax, Octree_t* oct)
{
    if(oct->center.x + oct->radius < boxMin.x ||
       oct->center.y + oct->radius < boxMin.y)
        return 0;

    if(oct->center.x - oct->radius > boxMax.x ||
       oct->center.y - oct->radius > boxMax.y)
        return 0;

    return 1;
}

int isLeaf(Octree oct)
{
    return oct->octants[0] == NULL;
}

int isEmpty(Octree oct)
{
    return oct->numPts == 0;
}

/** Wrapper to handle malloc() nicely */
void* mallocWrap(size_t size)
{
    void *ptr = malloc(size);
    if(!ptr)
    {
        fprintf(stderr, "Not enough memory\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void* reallocWrap(void *ptr, size_t size)
{
    void *tmp = realloc(ptr, size);
    if(tmp == NULL)
    {
        fprintf(stderr, "Error in realloc() of size %zu\n", size);
        exit(EXIT_FAILURE);
    }
    else
    {
        ptr = tmp;
    }
    return ptr;
}


unsigned int stage1cpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min)
{


    double initX = min.x + Wsize/2 - Wsize*Overlap;
    double initY = min.y + Wsize/2 - Wsize*Overlap;

    double Displace = round2d(Wsize*(1-Overlap));

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    // Lpoint** neighbors = NULL;
    Lpoint newmin;

    unsigned int countMin;

    int cellPoints;

    std::vector<int> vcount(Ccol,0);

    #pragma omp parallel for private(cellCenter,cellPoints,newmin,countMin) schedule(dynamic)
    for(int jj = 0 ; jj < Ccol ; jj++ ){

        countMin=0; // counts mins per row

        cellCenter.y = initY + jj*Displace;
        // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);

        for(int ii=0 ; ii < Crow ; ii++ ){

            // std::vector<std::function<bool(Vector3D,float)>> conditions = {
            //     [&diag,&cellCenter] (Vector3D c, float r) {
            //         double d = sqrt(pow(cellCenter.x - c.x,2) + pow(cellCenter.y - c.y,2));
            //         return d < (diag+r*sqrt(2)); // suma de las diagonales
            //     },
            // };

            cellCenter.x = initX + ii*Displace;
            // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
            // printf("Busco los vecinos\n");
            newmin = searchNeighborsMin(&cellCenter, octreeIn, Wsize/2, &cellPoints);
            // printf("Numero de elementos de la celda: %d\n", cellPoints );
            if(cellPoints >= minNumPoints ){

                minIDs[jj*Crow+ii] = newmin.id; //los voy guardando por posición de celda; esto hace que tenga que modificar stage2
                countMin++;
            }
        }
        vcount[jj] = countMin;
        // countMin=0;
    }

    countMin = std::accumulate(vcount.begin(), vcount.end(), 0);

    return countMin;
}


unsigned int stage2(unsigned int countMin, int* minIDs){

    unsigned int index = 0;

    int ii,jj,id;

    for( ii=0 ; ii<countMin ; ii=jj ){
        id = minIDs[ii];

        for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

        /* esta es la forma de descartar las posiciones no utilizadas, inicializadas a -1 */
        if((id >= 0) && (jj-ii > 1)){
            // if(jj-ii > 1){
            minIDs[index]=id;
            index++;
            // }
        }
    }

    return index;
}


unsigned int stage3cpp(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min){

    unsigned int addMin = 0;

    int cellPoints;

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    // Lpoint** neighbors = NULL;

    // Lpoint** minimos = NULL;

    Lpoint newmin;

    int ii,jj;

    #pragma omp parallel for private(cellCenter,cellPoints,newmin) schedule(static)
    for(int jj = 0 ; jj < Ccol ; jj++ ){

        cellCenter.y = min.y + Bsize/2 + jj*Bsize;

        for(int ii = 0 ; ii < Crow ; ii++ ){
            cellCenter.x = min.x + Bsize/2 + ii*Bsize;

            // Candidata a VOID porque no necesito el mínimo, solo el número de puntos encontrados.
            // newmin = searchNeighborsMin(&cellCenter, grid, Bsize/2, &cellPoints);
            countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
            // printf("Numero de elementos de la celda: %d\n", cellPoints );
            // Tengo que hacerlo porque el algoritmo no lo hace
            if(cellPoints == 0){
                // printf("    Tengo una celda vacia en la malla\n" );
                newmin = searchNeighborsMin(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                if(cellPoints>0){

                    // idmin = findMin(neighbors, cellPoints);
                    #pragma omp critical
                    {
                        minGridIDs[addMin] = newmin.id;
                        addMin++;
                    }
                    // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
                }
            }
        }
    }

    return addMin;

}


/** Calculate the radius in each axis and save the max radius of the bounding box */
Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius)
{
    Vector3D radii;

    radii.x = (max.x - min.x) / 2.0;
    radii.y = (max.y - min.y) / 2.0;
    radii.z = (max.z - min.z) / 2.0;

    if(radii.x >= radii.y && radii.x >= radii.z)
    {
        *maxRadius = radii.x;
    }
    else if(radii.y >= radii.x && radii.y >= radii.z)
    {
        *maxRadius = radii.y;
    }
    else
    {
        *maxRadius = radii.z;
    }

    return radii;
}

/** Calculate the center of the bounding box */
Vector3D getCenter(Vector3D min, Vector3D radius)
{
    Vector3D center;

    center.x = min.x + radius.x;
    center.y = min.y + radius.y;
    center.z = min.z + radius.z;

    return center;
}

// Find the child corresponding a given point
int octantIdx(Lpoint *point, Octree octree)
{
    int child = 0;

    if(point->x >= octree->center.x) child |= 4;
    if(point->y >= octree->center.y) child |= 2;
    if(point->z >= octree->center.z) child |= 1;

    // printf("child: %d\n",child);

    return child;
}


// Lpoint findValidMin(Octree octree, const std::vector<std::function<bool(Vector3D, float)>>& conditions)
Lpoint findValidMin(Octree octree, Vector3D* boxMin, Vector3D* boxMax, int* numInside)
{
    // Lpoint tmp, min = nomin;
    Lpoint tmp, min = {0,0.0,0.0,99999.0};
    int i;

    if(isLeaf(octree))
    {
        // printf("octree with %d points\n", octree->numPts);
        // printf("nodo hoja con radio %g\n", octree->radius/2.0);
        if(!isEmpty(octree))
        {
            // // Check conditions are fulfilled by the value or store MAX in (0,0) otherwise
            // for (auto& condition : conditions) {
            //     if (!condition(octree->center, octree->radius)) {
            //         return nomin
            //     }
            // }

            // // Don't need to check again because I've checked before
            // if(!boxOverlap2D(boxMin, boxMax, octree->octants[i]))
            //     return nomin

            // min = *(octree->points[0]);
            // *numNeighs += octree->numPts;

            for (i = 0; i < octree->numPts; i++) {
                if(insideBox2D(octree->points[i], *boxMin, *boxMax))
                {
                    tmp = *(octree->points[i]);
                    if (tmp.z < min.z) {
                        min = tmp;
                    }
                    (*numInside)++;
                }
            }

        }

    } else {

        for(i = 0; i < 8; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
                continue;
            else {
                tmp = findValidMin(octree->octants[i], boxMin, boxMax, numInside);
                if (tmp.z < min.z) {
                    min = tmp;
                }
            }

        }

    }

    return min;
}

Lpoint searchNeighborsMin(Lpoint* point, Octree octree, float radius, int* numNeighs)
{
    Vector3D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMin(octree, &boxMin, &boxMax, numNeighs);
}

void countNeighbors(Octree octree, Vector3D* boxMin, Vector3D* boxMax, int* numInside)
{
    int i;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {

            for (i = 0; i < octree->numPts; i++) {
                if(insideBox2D(octree->points[i], *boxMin, *boxMax))
                {
                  (*numInside)++;
                }
            }

        }

    } else {

        for(i = 0; i < 8; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
              continue;
            else {
              countNeighbors(octree->octants[i], boxMin, boxMax, numInside);
            }

        }

    }

    return;
}

void countNeighbors2D(Lpoint* point, Octree octree, float radius, int* numNeighs)
{
    Vector3D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    countNeighbors(octree, &boxMin, &boxMax, numNeighs);

    return;
}

// Lpoint findValidMinPath(Octree octree, Vector3D* boxMin, Vector3D* boxMax,
//                     int* numInside, std::vector<int>& path)
// {
//     // Lpoint tmp, min = nomin;
//     Lpoint tmp, min = {0,0.0,0.0,99999.0};
//     // std::vector<int> local_path, best_path;
//     int i;
//
//     if(isLeaf(octree))
//     {
//         if(!isEmpty(octree))
//         {
//             for (i = 0; i < octree->numPts; i++) {
//                 if(insideBox2D(octree->points[i], *boxMin, *boxMax))
//                 {
//                     tmp = *(octree->points[i]);
//                     if (tmp.z < min.z) {
//                         min = tmp;
//                     }
//                     (*numInside)++;
//                 }
//             }
//
//         }
//
//     } else {
//
//         std::vector<int> local_path, best_path;
//
//         for(i = 0; i < 8; i++) {
//             // Check
//             if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
//                 continue;
//             else {
//                 local_path.push_back(i);
//                 tmp = findValidMinPath(octree->octants[i], boxMin, boxMax, numInside, local_path);
//                 if (tmp.z < min.z) {
//                     min = tmp;
//                     best_path = local_path;
//                 }
//                 local_path.clear();
//             }
//         }
//         path.insert(path.end(), best_path.begin(), best_path.end());
//
//     }
//
//     return min;
// }
//
// Lpoint searchNeighborsMinPath(Lpoint* point, Octree octree, float radius, int* numNeighs, std::vector<int>& path)
// {
//     Vector3D boxMin, boxMax;
//
//     *numNeighs = 0;
//     makeBox(point, radius, &boxMin, &boxMax);
//
//     return findValidMinPath(octree, &boxMin, &boxMax, numNeighs, path);
// }

/** Create a octree with the given center and radius */
Octree createOctreeI(Vector3D center, float radius)
{
    int i = 0;
    Octree oct = NULL;

    oct = (Octree_t*)mallocWrap(sizeof(Octree_t));
    // oct = std::unique_ptr<Octree_t>{ (Octree_t*)mallocWrap(sizeof(Octree_t)) };
    oct->center = center;
    oct->radius = radius;
    // oct->points = NULL;
    oct->numPts = 0;
    for(i = 0; i < 8; i++)
      oct->octants[i] = NULL;

    return oct;
}

Octree createOctreeF(Vector3D center, float radius)
{
    int i = 0;
    // std::unique_ptr<Octree_t> oct = std::unique_ptr<Octree_t>{(Octree_t*)malloc(sizeof(Octree_t))};
    // std::unique_ptr<Octree_t> oct( new Octree_t );

    Octree oct = (Octree_t*)mallocWrap(sizeof(Octree_t));
    // oct = std::unique_ptr<Octree_t>{ (Octree_t*)mallocWrap(sizeof(Octree_t)) };
    // oct = std::unique_ptr<Octree_t>{ new Octree_t };
    oct->center = center;
    oct->radius = radius;
    oct->points = NULL;
    oct->numPts = 0;
    for(i = 0; i < 8; i++)
      oct->octants[i] = NULL;

    return std::move(oct);
}


void createOctantsF(Octree_t* oct)
{
    int i = 0;
    Vector3D newCenter;
    float newRadius = oct->radius * 0.5;
    // if(newRadius / 2.0 > MIN_RADIUS){
    //   for(i = 0; i < 8; i++)
    //   {
    //       newCenter = oct->center;
    //       newCenter.x += oct->radius * (i&4 ? 0.5f : -0.5f);
    //       newCenter.y += oct->radius * (i&2 ? 0.5f : -0.5f);
    //       newCenter.z += oct->radius * (i&1 ? 0.5f : -0.5f);
    //       oct->octants[i] = createOctreeI(newCenter, oct->radius * 0.5);
    //   }
    //   // printf("octante FI\n");
    // }
    // else
    // {
      for(i = 0; i < 8; i++)
      {
          newCenter = oct->center;
          newCenter.x += oct->radius * (i&4 ? 0.5f : -0.5f);
          newCenter.y += oct->radius * (i&2 ? 0.5f : -0.5f);
          newCenter.z += oct->radius * (i&1 ? 0.5f : -0.5f);
          oct->octants[i] = createOctreeF(newCenter, oct->radius * 0.5);
      }
    // }
}


// Insert a point in the octree creating the appropiate childs
void insertPointF(Lpoint *point, Octree octree, float minRadius, int nivel)
{
    int idx = 0;

    if(isLeaf(octree))
    {
        // printf("octante hoja nivel %d\n",nivel);
        if(octree->radius / 2.0 > minRadius)    // still divisible -> divide
        {
            createOctantsF(octree);
            // fillOctants(octree);
            idx = octantIdx(point, octree);
            insertPointF(point, octree->octants[idx], minRadius, nivel+1);

        } else {
            // printf("\toctree con %d puntos\n", octree->numPts);
            if(isEmpty(octree))             // Empty leaf -> insert point
            {
                // printf("\ta nodo vacio con radio %g\n", octree->radius/2.0);
                octree->points = (Lpoint**)mallocWrap(sizeof(Lpoint*));
                octree->points[0] = point;
                octree->numPts = 1;
            }
            else                         // Not empty and isn't divisible -> insert point
            {
                // printf("\t\ten uno lleno\n");
                octree->points = (Lpoint**)reallocWrap(octree->points, ++octree->numPts * sizeof(Lpoint*));
                octree->points[octree->numPts-1] = point;
            }

        }
        // else                            // Not empty but still divisible -> divide
        // {
        //     //if(octree->numPts > MAX_POINTS)
        //     if(octree->radius / 2.0 > MIN_RADIUS)
        //     {
        //         createOctants(octree);
        //         fillOctants(octree);
        //         idx = octantIdx(point, octree);
        //         insertPoint(point, octree->octants[idx]);

        //     }
        //     else                         // Not empty and isn't divisible -> insert point
        //     {
        //         octree->points = (Lpoint**)reallocWrap(octree->points, ++octree->numPts * sizeof(Lpoint*));
        //         octree->points[octree->numPts-1] = point;
        //     }
        // }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
      idx = octantIdx(point, octree);
      insertPointF(point, octree->octants[idx], minRadius, nivel+1);
      // insertPoint(point, octree->octants[idx]);
    }
}
