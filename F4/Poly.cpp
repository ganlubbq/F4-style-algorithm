#include "Poly.h"

Poly::Poly(vector<unsigned char> &coeff)
{
	_Coeff = coeff;
	set_LM();
	set_LMdeg();
}

//アライメント　simdに必要なほか高速化にも寄与
void* Poly::operator new(size_t size) {
	return _mm_malloc(size, 32);
}

//アライメント
void Poly::operator delete(void* pv) {
	_mm_free(pv);
}

//LMとLMのindexを代入 init
inline void Poly::set_LM()
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

//set_LMの後じゃないと動かない init
inline void Poly::set_LMdeg()
{
	_LMdeg = _Degree.index_to_degree(_LMdeg_index);
}