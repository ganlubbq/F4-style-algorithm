#pragma once
#define D1
#include <vector>
#ifdef D1
#include "Degree_table.h"
#endif //D1

int variables = 10;

using namespace std;

/*��ƂȂ�N���X�@������p�����邱�Ƃŋ@�\����������\��
Degree�͈�ʓI�Ȃ��̂̂ݑΉ��@�ߓx�ȍ������͔ėp��������
*/
class Poly
{
public:
	Poly(vector<unsigned char> &coeff);

	//variable
	int _LM;
	vector<unsigned char> _LMdeg;
	int _LMdeg_index;
	vector<unsigned char> _Coeff;
#ifdef D1
	static Degree_table _Degree;
#endif //D1

	//function
	void* operator new(size_t size);
	void operator delete(void* pv);

	virtual void set_LM();
	virtual void set_LMdeg();
	//inline virtual void set_LMdeg_index();
};

#ifdef D1
Degree_table d(variables);
Degree_table Poly::_Degree = d;
#endif //D1