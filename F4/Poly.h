#pragma once
#include <vector>

using namespace std;

/*��ƂȂ�N���X�@������p�����邱�Ƃŋ@�\����������\��
Degree�͈�ʓI�Ȃ��̂̂ݑΉ��@�ߓx�ȍ������͔ėp��������
*/
template<class Degree>
class Poly
{
public:
	Poly(vector<unsigned char> &coeff);

	//variable
	int LM;
	vector<unsigned char> LMdeg;
	int LMdeg_index;
	vector<unsigned char> Coeff;

	//function
	virtual void set_LM();
	virtual void set_LMdeg();
	virtual void set_LMdeg_index();
};

Poly<class Degree>::Poly(vector<unsigned char> &coeff)
{
	Coeff = coeff;
	set_LM();
	set_LMdeg();
	set_LMdeg_index();
}

