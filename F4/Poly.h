#pragma once
#include <vector>
#define D1
#ifdef D1
#include "Degree_table.h"
#endif //D1

using namespace std;

static const int variables = 10;

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

//main�Œ�`
	static int _Max_degree;
#ifdef D1
	static Degree_table _Degree;
#endif //D1
//

	//function
	void* operator new(size_t size);
	void operator delete(void* pv);

	virtual void set_LM();
	virtual void set_LMdeg();
	//inline virtual void set_LMdeg_index();
};

