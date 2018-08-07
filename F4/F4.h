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
class F4
{
public:
	F4(string filename);

	//variables
	string file_name;
	static int _Variables;

#ifdef _GF31
	vector<GF31> _Equations;
#endif

	//function
	void file_read();
	int var_deg_comb(int n, int r);
};