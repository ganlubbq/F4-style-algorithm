#include "Poly_GF31_simd.h"

//PolyŠÖ˜A‰Šú‰»
#ifdef D1
Degree_table d(variables);
Degree_table Poly::_Degree = d;
#endif //D1
int Poly::_Max_degree = 7;

int main()
{
	vector<unsigned char> a(2);
	Poly_GF31_simd p(a);

	return 0;
}