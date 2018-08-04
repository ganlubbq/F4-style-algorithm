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
	Poly(vector<unsigned char> &coeff,Degree degree);

	//variable
	int _LM;
	vector<unsigned char> _LMdeg;
	int _LMdeg_index;
	vector<unsigned char> _Coeff;
	Degree _degree;

	//function
	void* operator new(size_t size);
	void operator delete(void* pv);

	inline virtual void set_LM();
	inline virtual void set_LMdeg();
	//inline virtual void set_LMdeg_index();
};

Poly<class Degree>::Poly(vector<unsigned char> &coeff)
{
	_Coeff = coeff;
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
	for (int i = _Coeff.size() - 1; i >= 0; i--)
	{
		if (_Coeff[i] != 0)
		{
			_LM = _Coeff[i];
			_LMdeg_index = i;
			break;
		}
	}
}

//set_LMの後じゃないと動かない
inline void Poly<Degree>::set_LMdeg()
{
	_LMdeg = _degree.index_to_deg(_LMdeg_index);
}