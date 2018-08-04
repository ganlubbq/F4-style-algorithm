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
	void* operator new(size_t size);
	void operator delete(void* pv);

	inline virtual void set_LM();
	inline virtual void set_LMdeg();
	//inline virtual void set_LMdeg_index();
};

Poly<class Degree>::Poly(vector<unsigned char> &coeff)
{
	Coeff = coeff;
	set_LM();
	set_LMdeg();
}

//アライメント　simdに必要なほか高速化にも寄与
void* Poly<Degree>::operator new(size_t size) {
	return _mm_malloc(size, 32);
}

//アライメント
void Poly<Degree>::operator delete(void* pv) {
	_mm_free(pv);
}

//LMとLMのindexを代入
inline void Poly<Degree>::set_LM()
{
	for (int i = Coeff.size() - 1; i >= 0; i--)
	{
		if (Coeff[i] != 0)
		{
			LM = Coeff[i];
			LMdeg_index = i;
			break;
		}
	}
}

//set_LMの後じゃないと動かない
inline void Poly<Degree>::set_LMdeg()
{
	LMdeg = Degree.index_to_deg(LMdeg_index);
}