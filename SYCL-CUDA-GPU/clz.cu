#include <iostream>
#include <string>


__global__
void clz_op(const uint64_t l, const uint64_t r)
{
	printf("resultado: %d\n", ::__clzll(l ^ r));
	return;
}

int main(int argc, char* argv[])
{
	int resultado = 0;
	uint64_t left = atol(argv[1]);
	uint64_t right = atol(argv[2]);
	printf("entrada %lu %lu\n", left, right);
	clz_op<<<1,1>>>( left, right);
	cudaDeviceSynchronize();
	return 0;
}
