#pragma once
#include <vector>

using namespace std;

//��ƂȂ�N���X�@������p�����邱�Ƃŋ@�\����������\��
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
	virtual void operator+(Poly &poly);
	virtual void operator*(Poly &monomial);
	virtual void operator*(unsigned char &coeff);

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

void Poly<class Degree>::operator+(Poly &poly)
{

}