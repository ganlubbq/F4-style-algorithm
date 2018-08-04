#pragma once
#include <vector>

using namespace std;

/*基準となるクラス　これを継承することで機能を実現する予定
Degreeは一般的なもののみ対応　過度な高速化は汎用性を失う
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

