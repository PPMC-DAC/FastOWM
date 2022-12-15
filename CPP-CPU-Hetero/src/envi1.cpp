#include "../include/envi1.h"
#include "tbb/parallel_for.h"
#include <numeric>

double round2d(double z){
  return round(z*100)/100;
}

unsigned int findMin(Lpoint** neighbors, unsigned int cellPoints) {
  unsigned int idmin=0;
  double zzmin =round2d(neighbors[0]->z);
  for(int i=1; i<cellPoints; i++)
    if(neighbors[i]->z < zzmin){
      zzmin = round2d(neighbors[i]->z);
      idmin=i;
    }
  return idmin;
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

int boxOverlap2D(Vector3D boxMin, Vector3D boxMax, Octree oct)
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


Lpoint** neighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox2D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id)
                {
                    if (*numInside >= *ptsInside_size) {
                        (*ptsInside_size) += REALLOC_INCREMENT;
                        ptsInside = (Lpoint**)reallocWrap(ptsInside, *ptsInside_size * sizeof(Lpoint*));
                    }
                    ptsInside[(*numInside)++] = octree->points[i];
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            if(!boxOverlap2D(boxMin, boxMax, octree->octants[i]))
            {
                continue;
            }
            else
            {
                ptsInside = neighbors2D(point, boxMin, boxMax, octree->octants[i], ptsInside, ptsInside_size, numInside);
            }
        }
    }

    return ptsInside;
}

Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int ptsInside_size = 0;


    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors2D(point, boxMin, boxMax, octree, ptsInside, &ptsInside_size, numNeighs);

    return ptsInside;
}

unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min)
{

    // Tamaño del bloque del scheduler de OMP
    // unsigned short block_size = 1;

    double initX = min.x + Wsize/2 - Wsize*Overlap;
    double initY = min.y + Wsize/2 - Wsize*Overlap;

    double Displace = round2d(Wsize*(1-Overlap));

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    // Lpoint** neighbors = NULL;
    Lpoint newmin;

    unsigned int countMin = 0, idmin = 0;

    int cellPoints = 0;

    int ii,jj;

    #pragma omp parallel shared(countMin) firstprivate(Ccol,Crow,minIDs,octreeIn,Wsize,Displace,initX,initY) \
                            private(ii,jj,cellCenter,cellPoints,idmin,newmin)
    {
        // printf("Thread: %d\n", omp_get_thread_num());
        // printf("Thread num: %d\n", omp_get_num_threads());
        for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=omp_get_num_threads() ){
        // for( jj = 1 ; jj < Ccol ; jj++ ){

            cellCenter.y = initY + jj*Displace;
            // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);

            #pragma omp parallel for shared(countMin) schedule(dynamic)
            for( ii=0 ; ii < Crow ; ii++ ){

                cellCenter.x = initX + ii*Displace;
                // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
                // printf("Busco los vecinos\n");
                // neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
                newmin = searchNeighborsMin(&cellCenter, octreeIn, Wsize/2, &cellPoints);
                // printf("Numero de elementos de la celda: %d\n", cellPoints );
                if(cellPoints >= minNumPoints ){
                    // printf("Numero de elementos de la celda: %d\n", cellPoints );
                    // idmin = findMin(neighbors, cellPoints);
                    // idmin = min.id;
                    #pragma omp critical
                    {
                        minIDs[countMin] = newmin.id;
                        countMin++;
                    }
                    // printf("El mínimo %d de la celda es: %.2f\n",countMin ,neighbors[idmin]->z);

                }

                // free(neighbors);
                // neighbors = NULL;
            }
        }
    }

  return countMin;
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

int stage1tbb(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min)
{


    double initX = min.x + Wsize/2 - Wsize*Overlap;
    double initY = min.y + Wsize/2 - Wsize*Overlap;

    double Displace = round2d(Wsize*(1-Overlap));

    // Lpoint cellCenter = {0,0.0,0.0,0.0};

    // Lpoint** neighbors = NULL;
    // Lpoint newmin;

    // unsigned int finalcount;

    // int cellPoints = 0;

    // std::vector<int> vcount(Ccol,0);

    tbb::parallel_for(0, static_cast<int>(Ccol),
        [&](int jj){

        Lpoint newmin;
        Lpoint cellCenter;
        int cellPoints;
        // unsigned int countMin=0;

        cellCenter.y = initY + jj*Displace;

        for(int ii=0 ; ii < Crow ; ii++ ){

            cellCenter.x = initX + ii*Displace;

            newmin = searchNeighborsMin(&cellCenter, octreeIn, Wsize/2, &cellPoints);

            if(cellPoints >= minNumPoints ){
                //los voy guardando por posición de celda; esto hace que tenga que modificar stage2
                minIDs[jj*Crow+ii] = newmin.id;
                // countMin++;
            }
        }
        // mejor aquí para evitar false-sharing
        // vcount[jj] = countMin;
    });

    // finalcount = std::accumulate(vcount.begin(), vcount.end(), 0);

    return 0;
}



unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min){

  // Tamaño del bloque del scheduler de OMP
  // unsigned short block_size = 1;

  double initX = min.x + Wsize/2 - Wsize*Overlap;
  double initY = min.y + Wsize/2 - Wsize*Overlap;

  double Displace = round2d(Wsize*(1-Overlap));

  Lpoint cellCenter = {0,0.0,0.0,0.0};

  Lpoint** neighbors = NULL;

  unsigned int countMin = 0, idmin = 999;

  int cellPoints = 0;

  int ii,jj;

      for( jj = 1 ; jj < Ccol ; jj++ ){

          cellCenter.y = initY + jj*Displace;
          // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);

              for( ii=0 ; ii < Crow ; ii++ ){

                  cellCenter.x = initX + ii*Displace;
                  // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
                  // printf("Busco los vecinos\n");
                  neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
                  // printf("Numero de elementos de la celda: %d\n", cellPoints );
                  if(cellPoints >= minNumPoints ){
                      // printf("Numero de elementos de la celda: %d\n", cellPoints );
                      idmin = findMin(neighbors, cellPoints);
                      // #pragma omp critical
                      // {
                          minIDs[countMin] = neighbors[idmin]->id;
                          countMin++;
                      // }
                      // printf("El mínimo %d de la celda es: %.2f\n",countMin ,neighbors[idmin]->z);

                  }

                  free(neighbors);
                  neighbors = NULL;
              }
          }
      // }
  // }

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

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min){

    unsigned int addMin = 0, idmin = 999;

    int cellPoints = 0;

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    Lpoint** neighbors = NULL;

    Lpoint** minimos = NULL;

    int ii,jj;

    #pragma omp parallel shared(addMin) firstprivate(Ccol,Crow,minGridIDs,octreeIn,grid,Bsize) \
                              private(ii,jj,cellCenter,neighbors,minimos,cellPoints,idmin)
    {

        for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=omp_get_num_threads() ){
        // for( jj = 0 ; jj < Ccol ; jj++ ){
           cellCenter.y = min.y + Bsize/2 + jj*Bsize;

           for( ii = 0 ; ii < Crow ; ii++ ){
               cellCenter.x = min.x + Bsize/2 + ii*Bsize;

               minimos = searchNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
               // printf("Numero de elementos de la celda: %d\n", cellPoints );
               //Tengo que hacerlo porque el algoritmo no lo hace
               if(cellPoints == 0){
                   // printf("    Tengo una celda vacia en la malla\n" );
                   neighbors = searchNeighbors2D(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                   if(cellPoints>0){

                       idmin = findMin(neighbors, cellPoints);
                       #pragma omp critical
                       {
                          minGridIDs[addMin] = neighbors[idmin]->id;
                          addMin++;
                       }
                       // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
                   }
                   free(neighbors);
                   neighbors = NULL;

               }
               free(minimos);
               minimos=NULL;
           }
        }
    }

    return addMin;

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

unsigned int stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min){

    unsigned int addMin = 0, idmin = 999;

    int cellPoints = 0;

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    Lpoint** neighbors = NULL;

    Lpoint** minimos = NULL;

    int ii,jj;

        for( jj = 0 ; jj < Ccol ; jj++ ){
           cellCenter.y = min.y + Bsize/2 + jj*Bsize;

           for( ii = 0 ; ii < Crow ; ii++ ){
               cellCenter.x = min.x + Bsize/2 + ii*Bsize;

               minimos = searchNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
               // printf("Numero de elementos de la celda: %d\n", cellPoints );
               //Tengo que hacerlo porque el algoritmo no lo hace
               if(cellPoints == 0){
                   // printf("    Tengo una celda vacia en la malla\n" );
                   neighbors = searchNeighbors2D(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                   if(cellPoints>0){

                       idmin = findMin(neighbors, cellPoints);
                       // #pragma omp critical
                       // {
                          minGridIDs[addMin] = neighbors[idmin]->id;
                          addMin++;
                       // }
                       // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
                   }
                   free(neighbors);
                   neighbors = NULL;

               }
               free(minimos);
               minimos=NULL;
           }
        }
    // }

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

/** Create a octree with the given center and radius */
Octree createOctree(Vector3D center, float radius)
{
    int i = 0;
    Octree oct = NULL;

    oct = (Octree_t*)mallocWrap(sizeof(Octree_t));
    oct->center = center;
    oct->radius = radius;
    oct->points = NULL;
    oct->numPts = 0;
    // oct->plane[0] = 0.0;
    // oct->plane[1] = 0.0;
    // oct->plane[2] = 0.0;
    // oct->plane[3] = 0.0;
    for(i = 0; i < 8; i++)
        oct->octants[i] = NULL;

    return oct;
}

void createOctants(Octree oct)
{
    int i = 0;
    Vector3D newCenter;

    for(i = 0; i < 8; i++)
    {
        newCenter = oct->center;
        newCenter.x += oct->radius * (i&4 ? 0.5f : -0.5f);
        newCenter.y += oct->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.z += oct->radius * (i&1 ? 0.5f : -0.5f);
        oct->octants[i] = createOctree(newCenter, oct->radius * 0.5);
    }
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

/** Move the points inside the octree to the corresponding octants */
void fillOctants(Octree octree)
{
    int i = 0, idx = 0;

    for(i = 0; i < octree->numPts; i++)
    {
        idx = octantIdx(octree->points[i], octree);
        insertPoint(octree->points[i], octree->octants[idx]);
    }
    octree->numPts = 0;
    octree->points = NULL;
}

// Insert a point in the octree creating the appropiate childs
void insertPoint(Lpoint *point, Octree octree)
{
    int idx = 0;

    if(isLeaf(octree))
    {
        if(isEmpty(octree))             // Empty leaf -> insert point
        {
            octree->points = (Lpoint**)mallocWrap(sizeof(Lpoint*));
            octree->points[0] = point;
            octree->numPts = 1;
        }
        else                            // Not empty but still divisible -> divide
        {
            //if(octree->numPts > MAX_POINTS)
            if(octree->radius / 2.0 > MIN_RADIUS)
            {
                createOctants(octree);
                fillOctants(octree);
                idx = octantIdx(point, octree);
                insertPoint(point, octree->octants[idx]);

            }
            else                         // Not empty and isn't divisible -> insert point
            {
                octree->points = (Lpoint**)reallocWrap(octree->points, ++octree->numPts * sizeof(Lpoint*));
                octree->points[octree->numPts-1] = point;
            }
        }
    }
    else                                // No leaf -> search the correct one
    {
        idx = octantIdx(point, octree);
        insertPoint(point, octree->octants[idx]);
    }
}

// Insert a point in the octree creating the appropiate childs. Keep dividing until reaching radius
void insertPointMinRadius(Lpoint *point, Octree octree, float minRadius)
{
    int idx = 0;

    if(isLeaf(octree))
    {
        if(isEmpty(octree))             // Empty leaf -> insert point
        {
            octree->points = (Lpoint**)mallocWrap(sizeof(Lpoint*));
            octree->points[0] = point;
            octree->numPts = 1;
        }
        else                            // Not empty but still divisible -> divide
        {
            if(octree->radius / 2.0 > minRadius)
            {
                createOctants(octree);
                fillOctants(octree);
                idx = octantIdx(point, octree);
                insertPointMinRadius(point, octree->octants[idx],minRadius);

            }
            else                         // Not empty and isn't divisible -> insert point
            {
                octree->points = (Lpoint**)reallocWrap(octree->points, ++octree->numPts * sizeof(Lpoint*));
                octree->points[octree->numPts-1] = point;
            }
        }
    }
    else                                // No leaf -> search the correct one
    {
        idx = octantIdx(point, octree);
        insertPointMinRadius(point, octree->octants[idx],minRadius);
    }
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

Lpoint findValidMinPath(Octree octree, Vector3D* boxMin, Vector3D* boxMax,
                    int* numInside, std::vector<int>& path)
{
    // Lpoint tmp, min = nomin;
    Lpoint tmp, min = {0,0.0,0.0,99999.0};
    // std::vector<int> local_path, best_path;
    int i;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
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

        std::vector<int> local_path, best_path;

        for(i = 0; i < 8; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
                continue;
            else {
                local_path.push_back(i);
                tmp = findValidMinPath(octree->octants[i], boxMin, boxMax, numInside, local_path);
                if (tmp.z < min.z) {
                    min = tmp;
                    best_path = local_path;
                }
                local_path.clear();
            }
        }
        path.insert(path.end(), best_path.begin(), best_path.end());

    }

    return min;
}

Lpoint searchNeighborsMinPath(Lpoint* point, Octree octree, float radius, int* numNeighs, std::vector<int>& path)
{
    Vector3D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMinPath(octree, &boxMin, &boxMax, numNeighs, path);
}

/** Create a octree with the given center and radius */
Octree createOctreeFI(Vector3D center, float radius)
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
    Octree oct = NULL;

    oct = (Octree_t*)mallocWrap(sizeof(Octree_t));
    // oct = std::unique_ptr<Octree_t>{ (Octree_t*)mallocWrap(sizeof(Octree_t)) };
    oct->center = center;
    oct->radius = radius;
    oct->points = NULL;
    oct->numPts = 0;
    for(i = 0; i < 8; i++)
      oct->octants[i] = NULL;

    return oct;
}


void createOctantsF(Octree oct)
{
    int i = 0;
    Vector3D newCenter;
    float newRadius = oct->radius * 0.5;
    if(newRadius / 2.0 > MIN_RADIUS){
      for(i = 0; i < 8; i++)
      {
          newCenter = oct->center;
          newCenter.x += oct->radius * (i&4 ? 0.5f : -0.5f);
          newCenter.y += oct->radius * (i&2 ? 0.5f : -0.5f);
          newCenter.z += oct->radius * (i&1 ? 0.5f : -0.5f);
          oct->octants[i] = createOctreeFI(newCenter, oct->radius * 0.5);
      }
      // printf("octante FI\n");
    }
    else
    {
      for(i = 0; i < 8; i++)
      {
          newCenter = oct->center;
          newCenter.x += oct->radius * (i&4 ? 0.5f : -0.5f);
          newCenter.y += oct->radius * (i&2 ? 0.5f : -0.5f);
          newCenter.z += oct->radius * (i&1 ? 0.5f : -0.5f);
          oct->octants[i] = createOctreeF(newCenter, oct->radius * 0.5);
      }
    }
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
