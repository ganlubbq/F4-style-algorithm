#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <sstream>
#include <string>

#define _GF31
#ifdef _GF31
#include "GF31.h"
#endif //GF31

using namespace std;

//F4�A���S���Y���@�t�@�C���ǂݍ��݂͕W�������@��{�I�ɂ�Spoly class�����C���N���[�h�������@�傫�ȏ��������͌p�����ׂ�
template <class GF>
class F4
{
public:
	F4(string filename);

	//variables
	string file_name;
	static int _Variables;

	vector<GF> _Equations;

	//function
	void file_read();
	int var_deg_comb(int n, int r);
};

template <class GF>
F4<GF>::F4(string filename)
{
	file_name = filename;
}

template <class GF>
void F4<GF>::file_read()
{
	FILE *fp;
	fp = fopen(file_name.c_str(), "r");

	unsigned int coeff;

	fscanf(fp, "%d", &coeff);

	if (coeff != _Variables) {
		puts("VAL_NUM is not correct");
		exit(0);
	}

	cout << "init" << endl;
	unsigned char coeff_char;
	int coeff_int;

	//MQcharennge�^�����ǂݍ��߂Ȃ��ύX
	for (int i = 0; i < 2 * _Variables + 1; i++)
	{
		vector<unsigned char> temp;
		for (int j = 0; j < var_deg_comb(_Variables, 2); j++) {
			if(GF._DX == "d")	fscanf(fp, "%d", &coeff_int);
			else if (GF._DX == "x") fscanf(fp, "%x", &coeff_int);
			temp.push_back((unsigned char)coeff_int);
		}
		reverse(temp.begin(), temp.end());
		_Equations.push_back(GF31(temp));
	}
}

//n�ϐ�r���������܂ł̑S����
template <class GF>
inline int F4<GF>::var_deg_comb(int n, int r) {//n�ϐ�r���������܂ł̑S����
	int ans = 1;
	n += r;
	for (int div = 1; div <= r; ++div, --n) {
		ans = ans * n / div;
	}
	return ans;
}