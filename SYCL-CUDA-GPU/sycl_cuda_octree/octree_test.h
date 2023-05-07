#ifndef OCTREE_TEST_H

#define OCTREE_TEST_H

#pragma once

void octree_test(std::string inputTXT, const uint32_t chunkDim);

#include <sycl_cuda_octree/octree_test.inl>

void octree_traverse(std::string inputTXT, const uint32_t chunkDim);

void octree_traverse_heter(std::string inputTXT, const uint32_t chunkDim, const float factor);

#include <sycl_cuda_octree/octree_traverse.inl>

#endif // OCTREE_TEST_H