#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "F4.h"

//ëΩçÄéÆinclude
#define _GF31_
#ifdef _GF31_
#include "GF31.h"
vector<unsigned char> GF31::_Inverse = { 0,1,16,21,8,25,26,9,4,7,28,17,13,12,20,29,2,11,19,18,14,3,24,27,22,5,6,23,10,15,30 };
vector<unsigned char> GF31::_Add_inverse = { 0,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1 };
string GF31::_DX = "d";
#endif //_GF31_

//îªíËinclude
#define _decision_
#ifdef _decision_
#include "Decision.h"
#endif //_decision_

//Spoly include
#include "Spoly.h"

//Red include
#include "Red.h"

//LB include
#include "LB.h"

static const int variables = 12;
string filename = "GF31_12-0.txt";
//Polyä÷òAèâä˙âª
#ifdef D1
Degree_table d(variables);
Degree_table Poly::_Degree = d;
#endif //D1

int Poly::_Max_degree = 5;

#ifdef _GF31_ 
Decision<GF31> dd;
Spoly<GF31> ss;
Red<GF31> rr;
LB<GF31> ll;
int F4<GF31,Decision<GF31>,Spoly<GF31>, Red<GF31>, LB<GF31>>::_Variables = variables;
Decision<GF31> F4<GF31,Decision<GF31>,Spoly<GF31>, Red<GF31>, LB<GF31>>::_Decision = dd;
Spoly<GF31> F4<GF31,Decision<GF31>,Spoly<GF31>, Red<GF31>, LB<GF31>>::_Spoly = ss;
Red<GF31> F4<GF31,Decision<GF31>,Spoly<GF31>, Red<GF31>, LB<GF31>>::_Red = rr;
LB<GF31> F4<GF31,Decision<GF31>,Spoly<GF31>, Red<GF31>, LB<GF31>>::_LB = ll;
int F4<GF31, Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>>::_Parallel_div = 256;
#endif //_GF31_

int main()
{
	F4<GF31,Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>> f4(filename);
	f4.F4_style();

	//system("pause");
	return 0;
}