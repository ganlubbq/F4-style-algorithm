#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <time.h>
#include <fstream>
#include <string>

#include "F4.h"

//多項式include
#define _GF31_
#ifdef _GF31_
#include "GF31.h"
vector<unsigned char> GF31::_Inverse = { 0,1,16,21,8,25,26,9,4,7,28,17,13,12,20,29,2,11,19,18,14,3,24,27,22,5,6,23,10,15,30 };
vector<unsigned char> GF31::_Add_inverse = { 0,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1 };
string GF31::_DX = "d";
#endif //_GF31_

//判定include
#define _decision_
#ifdef _decision_
#include "Decision.h"
#endif //_decision_

//Poly関連初期化
#ifdef D1
Degree_table d(3);
Degree_table Poly::_Degree = d;
Degree_table Decision<GF31>::_Degree_deci = d;
#endif //D1

//Spoly include
#include "Spoly.h"

//Red include
#include "Red.h"

//LB include
#include "LB.h"

int ctoi(char *argv[]);

int Poly::_Max_degree = 5;

#ifdef _GF31_ 
Decision<GF31> dd;
Spoly<GF31> ss;
Red<GF31> rr;
LB<GF31> ll;
Decision<GF31> F4<GF31, Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>>::_Decision = dd;
Spoly<GF31> F4<GF31, Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>>::_Spoly = ss;
Red<GF31> F4<GF31, Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>>::_Red = rr;
LB<GF31> F4<GF31, Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>>::_LB = ll;
int F4<GF31, Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>>::_Parallel_div = 256;
///int F4<GF31, Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>>::_Seiki = 0;//0->normal 1->gauss 2->reduct 3 modify reduct
#endif //_GF31_

//a.out inputfile variables resultfile seiki alltimefile gauss_time_file(方程式の個数含む) LB_time_size_file red_time_size_file  deicision_file
int main(int argc, char *argv[])
{
	//for (int i = 0; i < argc; i++) cout << argv[i] << endl;
	//system("pause");
	int variables = ctoi(argv);

	//file
	string filename = argv[1];
	string writing_file = argv[3];
	string all_file = argv[5];
	string gauss_file = argv[6];
	string LB_file = argv[7];
	string red_file = argv[8];
	string decision_file = argv[9];

	int seiki;
	//正規化なし
	if (argv[4][0] == '0') seiki = 0;
	//ガウス掃き出し
	else if (argv[4][0] == '1') seiki = 1;
	//正規化あり
	else if (argv[4][0] == '2') seiki = 2;
	//メラー改良 正規化なし
	else if (argv[4][0] == '3') seiki = 3;
	//normal改良 多項式順序整理
	else if (argv[4][0] == '4') seiki = 4;

	std::cout << seiki << endl;

	F4<GF31, Decision<GF31>, Spoly<GF31>, Red<GF31>, LB<GF31>> f4(filename, variables, writing_file, seiki, all_file,gauss_file,LB_file,red_file,decision_file);

	f4.F4_style();
	//system("pause");
	return 0;
}

//くそ仕様注意　7~39以外バグる
int ctoi(char *argv[]) {
	return 3;
	switch (argv[2][0]) {
	case '1':
		switch (argv[2][1])
		{
		case '0': return 10;
		case '1': return 11;
		case '2': return 12;
		case '3': return 13;
		case '4': return 14;
		case '5': return 15;
		case '6': return 16;
		case '7': return 17;
		case '8': return 18;
		case '9': return 19;
		default: return 0;
		}
	case '2':
		switch (argv[2][1])
		{
		case '0': return 20;
		case '1': return 21;
		case '2': return 22;
		case '3': return 23;
		case '4': return 24;
		case '5': return 25;
		case '6': return 26;
		case '7': return 27;
		case '8': return 28;
		case '9': return 29;
		default: return 0;
		}
	case '3':
		switch (argv[2][1])
		{
		case '0': return 30;
		case '1': return 31;
		case '2': return 32;
		case '3': return 33;
		case '4': return 34;
		case '5': return 35;
		case '6': return 36;
		case '7': return 37;
		case '8': return 38;
		case '9': return 39;
		default: return 0;
		}
	case '7':return 7;
	case '8':return 8;
	case '9':return 9;
	default: return 0;
	}
}
