#ifndef OCTREE_TEST_H

#define OCTREE_TEST_H

#pragma once

#include<string>

void octree_test(std::string inputTXT, const uint32_t chunkDim);

#include <octree/octree_test.inl>

void octree_traverse(std::string inputTXT, const uint32_t chunkDim);

#include <octree/octree_traverse.inl>

#endif // OCTREE_TEST_H