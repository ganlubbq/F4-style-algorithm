#pragma once
#include <vector>
//degree�؂�ւ�
#define D1
#ifdef D1
#include "Degree_table.h"
#endif //D1

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

//main�Œ�`
	static int _Max_degree;
#ifdef D1
	static Degree_table _Degree;
#endif //D1
//

	//operator
	void* operator new(size_t size);
	void operator delete(void* pv);

	//function
	virtual void set_LM();
	virtual void set_LMdeg();
};

