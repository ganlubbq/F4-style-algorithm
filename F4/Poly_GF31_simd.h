#pragma once

#include "Poly.h"
//windows
#ifdef _MSC_VER
#include <intrin.h>
//linux
#else
#include <x86intrin.h>
#endif

class Poly_GF31_simd : public Poly<Degree>
{
	using Poly::Poly;
public:
	virtual void operator+(Poly &poly);
	virtual void operator*(Poly &monomial);
	virtual void operator*(unsigned char &coeff);

};