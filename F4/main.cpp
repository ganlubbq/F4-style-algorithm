#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "F4.h"

//‘½€®include
#define _GF31_
#ifdef _GF31_
#include "GF31.h"
#endif //_GF31_

//”»’èinclude
#define _decision_
#ifdef _decision_
#include "Decision.h"
#endif //_decision_

//Spoly include
#include "Spoly.h"

static const int variables = 2;
string filename = "GF31_2-0.txt";
//PolyŠÖ˜A‰Šú‰»
#ifdef D1
Degree_table d(variables);
Degree_table Poly::_Degree = d;
#endif //D1

int Poly::_Max_degree = 7;

#ifdef _GF31_ 
Decision<GF31> dd;
int F4<GF31,Decision<GF31>>::_Variables = variables;
Decision<GF31> F4<GF31,Decision<GF31>>::_Decision = dd;
#endif //_GF31_

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
	/*vector<unsigned char> a = {2,7,4,5,6,2,8,9,5,14,17,25,28,30,5,3,7,9,12,4};
	vector<unsigned char> b = {2,7,4,5,6,2,8,9,5,14,17,25,28,30,5,3,7,9,12,4};
	GF31 f(a),g(b);

	printvec(f._Coeff);
	printvec(g._Coeff);

	f + g;
	printvec(f._Coeff);
	unsigned char e = 5;
	g * e;
	printvec(g._Coeff);*/

	F4<GF31,Decision<GF31>> f4(filename);
	f4.F4_style();

	system("pause");
	return 0;
}