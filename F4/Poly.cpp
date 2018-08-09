#include "Poly.h"

Poly::Poly(vector<unsigned char> &coeff)
{
	_Coeff = coeff;
	set_LM();
	set_LMdeg();
}

//�A���C�����g�@simd�ɕK�v�Ȃق��������ɂ���^
void* Poly::operator new(size_t size) {
	return _mm_malloc(size, 32);
}

//�A���C�����g
void Poly::operator delete(void* pv) {
	_mm_free(pv);
}

//LM��LM��index���� init
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

//set_LM�̌ザ��Ȃ��Ɠ����Ȃ� init
inline void Poly::set_LMdeg()
{
	_LMdeg = _Degree.index_to_degree(_LMdeg_index);
}