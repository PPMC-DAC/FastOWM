#ifndef SYCL_TEST_H

#define SYCL_TEST_H

#pragma once

void sycl_test(std::string inputTXT, const uint32_t chunkDim);

#include <sycl_ordered/sycl_test.inl>

void sycl_traverse(std::string inputTXT, const uint32_t chunkDim);

#include <sycl_ordered/sycl_traverse.inl>

#endif // SYCL_TEST_H