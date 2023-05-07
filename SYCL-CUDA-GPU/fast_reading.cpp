#include <algorithm>
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <chrono>
#include <sstream>
#include <fstream>
#include <array>
#include <omp.h>
#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include "tbb/global_control.h"
#include <tbb/parallel_for.h>
// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

using tempo_t = std::chrono::steady_clock;

using cast_t = std::chrono::duration<double, std::milli>;

struct double4{
    double x;
    double y;
    double z;
    double w;
};

const char* map_file(const char* fname, size_t& length);
void map_file_stringsstream(const char* fname, uint32_t num_objects);
void map_file_atof(const char* fname, uint32_t num_objects, double4* cloud);
void map_file_atof_line(const char* fname, uint32_t num_objects, double4* cloud);
void map_file_atof_tbb(const char* fname, uint32_t num_objects, double4* cloud);
void map_file_atof_tbb_th(const char* fname, uint32_t num_objects, double4* cloud);
static uintmax_t wc(char const *fname);
static uintmax_t lectura(const char *fname);
void readHeader( std::string inputTXT, uint32_t& N);

int NUM_PROCS = 8;

int main( int argc, char* argv[])
{

    size_t length;
    uintmax_t m_numLines = 0;
    std::string inputTXT = (argc>1)? argv[1] : "data/INAER_2011_Alcoy_Core.xyz";
    NUM_PROCS = (argc>2)? atoi(argv[2]) : 8;
    uint32_t num_objects;
    readHeader(inputTXT, num_objects);
    
    tbb::global_control c(tbb::global_control::max_allowed_parallelism, NUM_PROCS);

    double4* cloud = new double4[num_objects];

    std::chrono::time_point<tempo_t> i_start = tempo_t::now();
    
    map_file_atof_tbb_th(inputTXT.c_str(), num_objects, cloud);
    // map_file_atof_tbb(inputTXT.c_str(), num_objects, cloud);
    // map_file_atof_line(inputTXT.c_str(), num_objects, cloud);
    // map_file_atof(inputTXT.c_str(), num_objects, cloud);
    // map_file_stringsstream(inputTXT.c_str(), length);
    // auto f = map_file(inputTXT.c_str(), length);
    // auto start = f;
    // auto l = f + length;
    // // printf("length: %zu\n", length);
    // const char* list[8];
    // list[0] = f;
    // int stride = int(num_objects/8) + 1;

    // std::istringstream strValue(f);
    // f = static_cast<const char*>(memchr(f, '\n', l-f));
    // int line_length = f-start+1;

    // int x, y, z;
    // double xd, yd, zd;
    // while(!strValue.eof()){
    //     strValue >> xd >> yd >> zd;
    // }
    // for (char a[100]; strValue.getline(&a[0], 100, '\n'); ) {
    //     xd = atof(a);
    //     yd = atof(a);
    //     zd = atof(a);
    // }
    // if(std::ifstream is{inputTXT.c_str(), std::ios::binary | std::ios::ate}) {
    //     auto size = is.tellg();
    //     std::string str(size, '\0'); // construct string to stream size
    //     is.seekg(0);
    //     if(is.read(&str[0], size)){}
    //         // std::cout << str << '\n';
    // }

    // FILE* fileTXT;
    // if((fileTXT = fopen(inputTXT.c_str(),"rb")) == NULL){
    //     printf("Unable to open file!\n");
    //     return -1;
    // }

    // if(fscanf(fileTXT, "%lf %lf %lf", &xd, &yd, &zd) < 3){
    //   printf("Imposible to obtain values\n");
    //   return -1;
    // }
    // printf("%lf %lf %lf\n", xd, yd, zd);

    // if(fclose(fileTXT)){
    //     printf("Cannot close the file\n");
    //     return -1;
    // } 

    // sscanf(f, "%d %d %d", &x,&y,&z);

    // strValue >> xd >> yd >> zd;
    // printf("%d %d %d\n",int(xd),int(yd),int(zd));
    // printf("%lf %lf %lf\n", xd, yd, zd);
    // x = atoi(start);
    // y = atoi(start+10);
    // z = atoi(start+22);
    // printf("%d %d %d\n",x,y,z);
    // strValue >> xd >> yd >> zd;
    // f = static_cast<const char*>(memchr(f, '\n', l-f));
    // int len = f-start+1;
    // // // f+=len;
    // start+=len;
    // printf("%d %d %d\n",int(xd),int(yd),int(zd));
    // printf("%lf %lf %lf\n", xd, yd, zd);
    // x = atoi(start);
    // y = atoi(start+10);
    // z = atoi(start+22);
    // printf("%d %d %d\n",x,y,z);



    // int p_count = 0;
    // int chunkStart = stride;
    // // f = static_cast<const char*>(memchr(f, '\n', l-f));
    // // int len = f-start+1;
    // // printf("stride: %d\n", chunkStart);
    // while (f && f!=l)
    //     if ((f = static_cast<const char*>(memchr(f, '\n', l-f)))){
    //         m_numLines++, f++;
    //         if(m_numLines == chunkStart){
    //             p_count++;
    //             chunkStart += stride-1;
    //             list[p_count] = f;
    //             // printf("stride: %d\n", chunkStart);
    //         }
    //         // strValue >> x >> y >> z;
    //         // int x1 = atoi(f);
    //         // int y1 = atoi(f+10);
    //         // int z1 = atoi(f+22);
    //         // if(int(x)!=x1 || int(y)!=y1 || int(z)!=z1){
    //         //     printf("ERROR!!\n");
    //         //     printf("%d %d %d\n",int(x),int(y),int(z));
    //         //     printf("%d %d %d\n",int(x1),int(y1),int(z1));
    //         // }
    //         // if(m_numLines==1) len = f-start;
    //     }
    // #pragma omp parallel shared(list,stride)
    // {
    // int id = omp_get_thread_num();
    // int fin = stride * (id+1);
    // int myline = stride * id;
    // std::istringstream strValue(list[id]);
    // double xd, yd, zd;
    // while(!strValue.eof() && myline < fin){
    //     strValue >> xd >> yd >> zd;
    //     myline++;
    // }
    // // printf("[%d] line: %d\n", id, myline);

    // }

    // auto f = map_file(inputTXT.c_str(), length);
    // auto l = f + length;
    // while(f!=l)
    // {
    //     // strValue >> x >> y >> z;
    //     double x1 = atof(f);
    //     double y1 = atof(f+10);
    //     double z1 = atof(f+22);
    //     // int x1 = atoi(f);
    //     // int y1 = atoi(f+10);
    //     // int z1 = atoi(f+22);
    //     if(cloud[m_numLines].x!=x1 || cloud[m_numLines].y!=y1 || cloud[m_numLines].z!=z1){
    //         printf("ERROR!! line: %zu\n", m_numLines);
    //         // printf("%d %d %d\n",int(x),int(y),int(z));
    //         printf("%.3f %.3f %.3f\n",cloud[m_numLines].x,cloud[m_numLines].y,cloud[m_numLines].z);
    //         printf("%.3f %.3f %.3f\n",x1,y1,z1);
    //         // printf("%d %d %d\n",int(x1),int(y1),int(z1));
    //         break;
    //     }
    //     f = static_cast<const char*>(memchr(f, '\n', l-f)) + 1;
    //     // start += f-start;
    //     // start += line_length;
    //     m_numLines++;
    // }

    std::cout << "m_numLines = " << m_numLines << "\n";
    double time = cast_t(tempo_t::now() - i_start).count();
    std::cout << " map_file time: " << time << " ms\n";

    i_start = tempo_t::now();

    m_numLines = wc(inputTXT.c_str());
    
    std::cout << "m_numLines = " << m_numLines << "\n";
    time = cast_t(tempo_t::now() - i_start).count();
    std::cout << " wc time: " << time << " ms\n";

    delete(cloud);

    return 0;
}

void handle_error(const char* msg) {
    perror(msg); 
    exit(255);
}

const char* map_file(const char* fname, size_t& length)
{
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");

    length = sb.st_size;

    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");

    // TODO close fd at some point in time, call munmap(...)
    return addr;
}

void map_file_stringsstream(const char* fname, uint32_t num_objects)
{
    size_t length = 0u;
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");

    length = sb.st_size;

    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");


    auto f = addr;
    auto l = f + length;
    // printf("length: %zu\n", length);
    const char* list[8];
    list[0] = f;
    int stride = int(num_objects/8) + 1;
    int p_count = 0;
    int chunkStart = stride;
    uintmax_t m_numLines = 0;

    while (f && f!=l)
        if ((f = static_cast<const char*>(memchr(f, '\n', l-f)))){
            m_numLines++, f++;
            if(m_numLines == chunkStart && m_numLines < num_objects){
                p_count++;
                chunkStart = (chunkStart+stride <= num_objects)? chunkStart+stride : num_objects;
                list[p_count] = f;
                printf("stride: %d, p_count: %d\n", chunkStart, p_count);
            }
        }

    #pragma omp parallel shared(list,stride)
    {
    int id = omp_get_thread_num();
    int fin = stride * (id+1);
    int myline = stride * id;
    std::istringstream strValue(list[id]);
    double xd, yd, zd;
    while(!strValue.eof() && myline < fin){
        strValue >> xd >> yd >> zd;
        myline++;
    }

    }

    // TODO close fd at some point in time, call munmap(...)
    if(munmap((void*)addr, length) < 0)
        handle_error("munmap");

    return;
}

void map_file_atof(const char* fname, uint32_t num_objects, double4* cloud)
{
    size_t length = 0u;
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");

    length = sb.st_size;

    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");


    auto f = addr;
    auto l = f + length;
    // printf("length: %zu\n", length);
    const char* list[8];
    list[0] = f;
    int stride = int(num_objects/8) + 1;
    int p_count = 0;
    int chunkStart = stride;
    uintmax_t m_numLines = 0;
    // printf("stride: %d\n", chunkStart);

    while (f && f!=l)
        if ((f = static_cast<const char*>(memchr(f, '\n', l-f)))){
            m_numLines++, f++;
            if(m_numLines == chunkStart && m_numLines < num_objects){
                p_count++;
                chunkStart = (chunkStart+stride <= num_objects)? chunkStart+stride : num_objects;
                list[p_count] = f;
                // printf("stride: %d, p_count: %d\n", chunkStart, p_count);
            }
        }
    printf("m_numlines = %zu\n", m_numLines);

    #pragma omp parallel shared(list,stride,cloud)
    {
    int id = omp_get_thread_num();
    int fin = (stride*(id+1) <= num_objects)? stride*(id+1) : num_objects;
    int myline = stride * id;
    auto myf = list[id];
    // double x, y, z;
    for(int i=myline; i<fin; i++){
        cloud[i].x = atof(myf);
        cloud[i].y = atof(myf+10);
        cloud[i].z = atof(myf+22);
        // x = atof(f);
        // y = atof(f+10);
        // z = atof(f+22);
        myf = static_cast<const char*>(memchr(myf, '\n', l-myf)) + 1;
    }
    // printf("[%d] fin-myline = %d\n", id, fin-myline);
    // printf("[%d] myline = %d\n", id, myline);
    // int i=stride*id;
    // int end = i+100;
    // while(i<end){
    //     printf("%.3f %.3f %.3f\n",cloud[i].x,cloud[i].y,cloud[i].z);    
    //     i++;
    // }

    }
    // int i=0;
    // while(i<100){
    //     printf("%.3f %.3f %.3f\n",cloud[i].x,cloud[i].y,cloud[i].z);    
    //     i++;
    // }

    // TODO close fd at some point in time, call munmap(...)
    if(munmap((void*)addr, length) < 0)
        handle_error("munmap");

    return;
}

void map_file_atof_tbb(const char* fname, uint32_t num_objects, double4* cloud)
{
    size_t length = 0u;
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");

    length = sb.st_size;

    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");


    auto f = addr;
    auto l = f + length;
    // printf("length: %zu\n", length);
    // const char* list[8];
    const char** list = new const char*[NUM_PROCS];
    list[0] = f;
    int stride = int(num_objects/NUM_PROCS) + 1;
    int p_count = 0;
    int chunkStart = stride;
    uintmax_t m_numLines = 0;
    // printf("stride: %d\n", chunkStart);

    while (f && f!=l)
        if ((f = static_cast<const char*>(memchr(f, '\n', l-f)))){
            m_numLines++, f++;
            if(m_numLines == chunkStart && m_numLines < num_objects){
                p_count++;
                chunkStart = (chunkStart+stride <= num_objects)? chunkStart+stride : num_objects;
                list[p_count] = f;
                // printf("stride: %d, p_count: %d\n", chunkStart, p_count);
            }
        }
    printf("m_numlines = %zu\n", m_numLines);

    tbb::parallel_for( 0, NUM_PROCS, 1,
        [&](int id){
            int fin = (stride*(id+1) <= num_objects)? stride*(id+1) : num_objects;
            // int myline = stride * id;
            auto mylist = list[id];
            // double x, y, z;
            for(int i = stride*id; i<fin; i++){
                cloud[i].x = atof(mylist);
                cloud[i].y = atof(mylist+10);
                cloud[i].z = atof(mylist+22);

                mylist = static_cast<const char*>(memchr(mylist, '\n', l-mylist)) + 1;
            }
    });


    // TODO close fd at some point in time, call munmap(...)
    if(munmap((void*)addr, length) < 0)
        handle_error("munmap");

    delete(list);

    return;
}

void map_file_atof_tbb_th(const char* fname, uint32_t num_objects, double4* cloud)
{
    size_t length = 0u;
    int fd = open(fname, O_RDONLY);
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

void map_file_atof_line(const char* fname, uint32_t num_objects, double4* cloud)
{
    size_t length = 0u;
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");

    length = sb.st_size;

    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");


    auto f = addr;
    auto l = f + length;
    // printf("length: %zu\n", length);
    const char* list[8];
    list[0] = f;
    int stride = int(num_objects/8) + 1;

    f = static_cast<const char*>(memchr(f, '\n', l-f));
    int line_length = f-addr+1;
    int bytes_stride = line_length*stride;

    for(int i=1; i<8; i++){
        list[i] = i*bytes_stride + addr;
    }

#pragma omp parallel shared(list,stride,cloud)
{
    int id = omp_get_thread_num();
    int fin = (stride*(id+1) <= num_objects)? stride*(id+1) : num_objects;
    int myline = stride * id;
    auto myf = list[id];
    double x, y, z;
    while(myline < fin){
        cloud[myline].x = atof(myf);
        cloud[myline].y = atof(myf+10);
        cloud[myline].z = atof(myf+22);
        // x = atof(f);
        // y = atof(f+10);
        // z = atof(f+22);
        myf = static_cast<const char*>(memchr(myf, '\n', l-myf)) + 1;
        myline++;
    }

}

    // TODO close fd at some point in time, call munmap(...)
    if(munmap((void*)addr, length) < 0)
        handle_error("munmap");

    return;
}


static uintmax_t wc(char const *fname)
{
    static const auto BUFFER_SIZE = 16*1024;
    int fd = open(fname, O_RDONLY);
    if(fd == -1)
        handle_error("open");

    /* Advise the kernel of our access pattern.  */
    posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL

    char buf[BUFFER_SIZE + 1];
    uintmax_t lines = 0;

    while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
    {
        if(bytes_read == (size_t)-1)
            handle_error("read failed");
        if (!bytes_read)
            break;

        for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
            ++lines;
    }

    return lines;
}


static uintmax_t lectura(const char *fname)
{
    std::string Text1;
    double Number1 = 0;
    double Number2 = 0;	
    double Number3 = 0;	
    char Space;
    uintmax_t nlines=0;
    	std::ifstream LargeFile(fname);
    	while( getline(LargeFile, Text1, ' ') )		
    	{						
    	     LargeFile >> Number1;		        
    	    //  LargeFile >> Space;                   
    	     LargeFile >> Number2;			
    	    //  LargeFile >> Space;                   
    	     LargeFile >> Number3;
             nlines++;			
    	     LargeFile.get();	
    	}
    return nlines;
}

void readHeader( std::string inputTXT, uint32_t& N)
{
  if( inputTXT.find("INAER_2011_Alcoy.xyz") != std::string::npos ){
    N = 2772832;

  } else if( inputTXT.find("INAER_2011_Alcoy_Core.xyz") != std::string::npos ){ // Alcoy
    N = 20380212;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX2.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*2;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX4.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*4;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX6.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*6;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX8.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*8;

  } else if( inputTXT.find("BABCOCK_2017_Arzua_3B.xyz") != std::string::npos ){ //Arzua
    N = 40706503;

  } else if( inputTXT.find("V21_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion forestal
    N = 42384876;

  } else if( inputTXT.find("V19_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion urban
    N = 48024480;

  } else {
    printf("No header data!\n");
    exit(-1);
  }
}
