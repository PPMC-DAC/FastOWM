#ifndef WRAPPER_H
#include <stdlib.h>
#define	WRAPPER_H

void* mallocWrap(size_t size);
void* callocWrap(int nmemb, size_t  size);
void* reallocWrap(void *ptr, size_t size);

#endif	/* WRAPPER_H */
