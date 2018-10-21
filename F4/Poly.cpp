#include "Poly.h"

Poly::Poly(vector<unsigned char> &coeff)
{
	_Coeff = coeff;
	set_LM();
	set_LMdeg();
}

//0����������
Poly::Poly()
{
	_LM = 0;
	_LMdeg_index = -1;
	_LMdeg.resize(0);
	//1�����@��΂�����
	vector<unsigned char> coeff = { 0 };
	_Coeff = coeff;
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

//set_LM�̌ザ��Ȃ��Ɠ����Ȃ� init
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