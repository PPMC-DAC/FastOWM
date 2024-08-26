#ifndef OCTREE_TEST_H

#define OCTREE_TEST_H

#pragma once

void octree_test(std::string inputTXT, const uint32_t chunkDim);

#include <sycl_octree/octree_test.inl>

void octree_traverse(std::string inputTXT, const uint32_t chunkDim, const uint32_t Wsize, const real_t Overlap);
void octree_traverse_heter(std::string inputTXT, const uint32_t chunkDim, const float factor, const uint32_t Wsize, const real_t Overlap);
void octree_traverse_scheduler(std::string inputTXT, const uint32_t chunkDim, const int chunkGPU, const uint32_t _Wsize, const real_t _Overlap);

#include <sycl_octree/octree_traverse.inl>

#endif // OCTREE_TEST_H