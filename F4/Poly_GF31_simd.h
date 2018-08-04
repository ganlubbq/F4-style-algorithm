#pragma once
#include "Poly.h"

class Poly_GF31_simd : public Poly<Degree>
{
	using Poly::Poly;
public:
	virtual void operator+(Poly &poly);
	virtual void operator*(Poly &monomial);
	virtual void operator*(unsigned char &coeff);

};