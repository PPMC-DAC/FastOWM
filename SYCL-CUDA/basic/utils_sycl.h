#ifndef UTILS_H

#define UTILS_H

#include <unistd.h>
#include <math.h>
#include <string>
#include <fstream>
#include <bitset>
#include <cstring>
// mmap
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include "tbb/global_control.h"
#include "tbb/parallel_for.h"

void handle_error(const char* msg) {
  perror(msg); 
  exit(255);
}

// inline cudaError_t checkCuda(cudaError_t result)
// {
//   if (result != cudaSuccess) {
//     fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
//     assert(result == cudaSuccess);
//   }
//   return result;
// }

double round2d(double z){
  return round(z*100.0)/100.0;
}

void readHeader( std::string inputTXT, whole_t& BBox, uint32_t& N)
{
  if( inputTXT.find("INAER_2011_Alcoy.xyz") != std::string::npos ){
    N = 2772832;
    BBox.lower.x   = 715244.96;
    BBox.lower.y   = 4286623.63;
    BBox.lower.z   = 836.424;
    BBox.upper.x   = 716057.75;
    BBox.upper.y   = 4287447.70;
    BBox.upper.z   = 976.790;

  } else if( inputTXT.find("INAER_2011_Alcoy_Core.xyz") != std::string::npos ){ // Alcoy
    N = 20380212;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    BBox.upper.z = 991.516;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX2.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*2;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    BBox.upper.z = 991.516;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX4.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*4;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    BBox.upper.z = 991.516;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX6.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*6;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    BBox.upper.z = 991.516;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX8.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*8;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    BBox.upper.z = 991.516;

  } else if( inputTXT.find("BABCOCK_2017_Arzua_3B.xyz") != std::string::npos ){ //Arzua
    N = 40706503;
    BBox.lower.x = 568000.00;
    BBox.lower.y = 4752320.00;
    BBox.lower.z = 331.620;
    BBox.upper.x = 568999.99;
    BBox.upper.y = 4753319.99;
    BBox.upper.z = 495.630;

  } else if( inputTXT.find("V21_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion forestal
    N = 42384876;
    BBox.lower.x = 526964.093;
    BBox.lower.y = 4742610.292;
    BBox.lower.z = 38.656;
    BBox.upper.x = 527664.647;
    BBox.upper.y = 4743115.738;
    BBox.upper.z = 112.269;

  } else if( inputTXT.find("V19_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion urban
    N = 48024480;
    BBox.lower.x = 526955.908;
    BBox.lower.y = 4742586.025;
    BBox.lower.z = 38.150;
    BBox.upper.x = 527686.445;
    BBox.upper.y = 4743124.373;
    BBox.upper.z = 119.833;

  } else {
    handle_error("No header data");
  }
}

void map_file_atof_tbb_th(std::string fname, point_t* cloud, uint32_t num_objects)
{
    size_t length = 0u;
    int fd = open(fname.c_str(), O_RDONLY);
    if (fd < 0)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) < 0)
        handle_error("fstat");

    length = sb.st_size;

    // creates a new mapping in the virtual address space of the calling process
    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");

    // fd can be closed immediately without invalidating the mapping
    if(close(fd) < 0)
        handle_error("close");

    int stride = int(num_objects/NUM_PROCS) + 1;

    tbb::parallel_for( 0, NUM_PROCS, 1,
        [&](int id){

            int myline = stride * id;
            int fin = (stride*(id+1) <= num_objects)? stride*(id+1) : num_objects;
            auto mylist = addr;
            uint32_t m_numLines = 0;
            while(m_numLines < fin) {
                if(m_numLines == myline){
                    for(int i = myline; i<fin; i++){
                        cloud[i].x = atof(mylist);
                        cloud[i].y = atof(mylist+10);
                        cloud[i].z = atof(mylist+22);

                        mylist = static_cast<const char*>(memchr(mylist, '\n', 64)) + 1;
                    }
                    m_numLines += stride; //fin
                } 
                else {
                    mylist = static_cast<const char*>(memchr(mylist, '\n', 64)) + 1;
                    m_numLines++;
                }
                
            }
    });


    // unmap
    if(munmap((void*)addr, length) < 0)
        handle_error("munmap");

    return;
}

inline aabb_t merge(const aabb_t& lhs, const aabb_t& rhs) noexcept
{
    aabb_t merged;
    merged.upper.x = sycl::fmax(lhs.upper.x, rhs.upper.x);
    merged.upper.y = sycl::fmax(lhs.upper.y, rhs.upper.y);
    merged.lower.x = sycl::fmin(lhs.lower.x, rhs.lower.x);
    merged.lower.y = sycl::fmin(lhs.lower.y, rhs.lower.y);
    return merged;
}

inline real_v getCenter(const aabb_t& box) noexcept
{
    real_v c;
    c.x = (box.upper.x + box.lower.x) * 0.5;
    c.y = (box.upper.y + box.lower.y) * 0.5;
    return c;
}

inline void swap(uint32_t& left, uint32_t& right)
{
  uint32_t aux = left;
  left = right;
  right = aux;
}

inline int boxOverlap2D(const aabb_t& global, const aabb_t& local)
{

  if((local.upper.x - global.lower.x < TOL) ||
    (local.upper.y - global.lower.y < TOL))
      return 0;

  if((local.lower.x - global.upper.x > TOL) ||
    (local.lower.y - global.upper.y > TOL))
      return 0;

  return 1;

}

/* Checks if the search area is completely overlapped by the node */
inline int isAllOverlaped(const aabb_t& global, const aabb_t& local)
{
  if((local.upper.x - global.upper.x > TOL) ||
    (local.upper.y - global.upper.y > TOL))
      return 0;

  if((local.lower.x - global.lower.x < TOL) ||
    (local.lower.y - global.lower.y < TOL))
      return 0;

  return 1;
}

/* Checks if a point is contained in the search area */
inline int insideBox(const aabb_t& global, const point_t& point)
{
    if((point.x - global.lower.x > TOL) && (point.y - global.lower.y > TOL))
    {
        if((point.x - global.upper.x < TOL) && (point.y - global.upper.y < TOL))
        {
            return 1;
        }
    }

    return 0;
}

struct Split_task
{
  Split_task() {}
  Split_task(const uint32_t id, const uint32_t in, const uint32_t n, const uint32_t min) : 
    m_node( id ), m_input( in ), m_numPts(n), m_minIdx(min) {}

  uint32_t m_node; // nuevo nodo del octree
  uint32_t m_input; // antiguo nodo del bintree
  uint32_t m_numPts; // numPts
  uint32_t m_minIdx; // indice del minimo
};

// void prefetch()
// {
//   return;
// }

// template <typename pointer, typename size, typename ...Rest>
// void prefetch( pointer p, size s, int device, Rest... args )
// {
//   cudaMemPrefetchAsync(p, s, device);
//   prefetch( args... );
//   return;
// }

uint32_t stage2CPU(uint32_t countMin, uint32_t* minIDs){

  uint32_t index = 0, id;

  int ii,jj;

  /* esta es la forma de descartar las posiciones no utilizadas, inicializadas a 0 */
  for( jj=0 ;  minIDs[jj] == 0 ; ++jj );

  for( ii=jj ; ii<countMin ; ii=jj ){

    id = minIDs[ii];

    for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

    if(jj-ii > 1){
        // if(jj-ii > 1){
        minIDs[index]=id;
        index++;
        // }
    }
  }

    return index;
}


#endif
