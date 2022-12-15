#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/** Wrapper to handle callocWrap() nicely */
void* callocWrap(int nmemb, size_t  size)
{
    void *ptr = calloc(nmemb, size);
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
