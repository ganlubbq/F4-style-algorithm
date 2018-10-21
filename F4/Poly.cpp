#include "Poly.h"

Poly::Poly(vector<unsigned char> &coeff)
{
	_Coeff = coeff;
	set_LM();
	set_LMdeg();
}

//0多項式生成
Poly::Poly()
{
	_LM = 0;
	_LMdeg_index = -1;
	_LMdeg.resize(0);
	//1次元　やばいかも
	vector<unsigned char> coeff = { 0 };
	_Coeff = coeff;
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
	bool flag = true;
	for (int i = _Coeff.size() - 1; i >= 0; i--)
	{
		if (_Coeff[i] != 0)
		{
			flag = false;
			_LM = _Coeff[i];
			_LMdeg_index = i;
			break;
		}
	}

	if (flag)
	{
		_LM = 0;
		_LMdeg_index = -1;
	}
}

//set_LMの後じゃないと動かない init
inline void Poly::set_LMdeg()
{
	if (_LMdeg_index == -1)
	{
		_LMdeg.resize(0);
	}
	else
	{
		_LMdeg = _Degree.index_to_degree(_LMdeg_index);
	}
}