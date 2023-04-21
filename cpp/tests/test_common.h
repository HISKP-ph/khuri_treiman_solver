#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include "type_aliases.h"

#include "gtest/gtest.h"

using type_aliases::Complex;

void expect_near(const Complex& z1, const Complex& z2, double tolerance=1e-16)
{
    EXPECT_NEAR(z1.real(), z2.real(), tolerance);
    EXPECT_NEAR(z1.imag(), z2.imag(), tolerance);
}


#endif
