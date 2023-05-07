#ifndef ORDERED_TEST_H

#define ORDERED_TEST_H

#pragma once

void ordered_test(std::string inputTXT, const uint32_t chunkDim);

#include <ordered/ordered_test.inl>

void ordered_traverse(std::string inputTXT, const uint32_t chunkDim);

#include <ordered/ordered_traverse.inl>

#endif // ORDERED_TEST_H