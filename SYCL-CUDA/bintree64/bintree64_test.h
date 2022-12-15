#ifndef BINTREE64_TEST_H

#define BINTREE64_TEST_H

#pragma once

void bintree64_test(std::string inputTXT, const uint32_t chunkDim);

#include <bintree64/bintree64_test.inl>

void bintree64_traverse(std::string inputTXT, const uint32_t chunkDim);

#include <bintree64/bintree64_traverse.inl>

#endif // BINTREE64_TEST_H