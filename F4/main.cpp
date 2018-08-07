#include <iostream>
#include "F4.h"

static const int variables = 10;
//PolyŠÖ˜A‰Šú‰»
#ifdef D1
Degree_table d(variables);
Degree_table Poly::_Degree = d;
#endif //D1

int Poly::_Max_degree = 7;

#ifdef _GF31
int F4<GF31>::_Variables = variables;
#endif //_GF31

void printvec(vector<unsigned char> vec)
{
	cout << "[ ";
	for (int i = 0; i < vec.size(); i++)
	{
		cout << (int)vec[i] << " ";
	}
	cout << "]" << endl;
}

int main()
{
	vector<unsigned char> a = {2,7,4,5,6,2,8,9,5,14,17,25,28,30,5,3,7,9,12,4};
	vector<unsigned char> b = {2,7,4,5,6,2,8,9,5,14,17,25,28,30,5,3,7,9,12,4};
	GF31 f(a),g(b);

	printvec(f._Coeff);
	printvec(g._Coeff);

	f + g;
	printvec(f._Coeff);
	unsigned char e = 5;
	g * e;
	printvec(g._Coeff);

	system("pause");
	return 0;
}