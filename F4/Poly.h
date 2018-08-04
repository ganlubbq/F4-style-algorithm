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

//�A���C�����g�@simd�ɕK�v�Ȃق��������ɂ���^
void* Poly<Degree>::operator new(size_t size) {
	return _mm_malloc(size, 32);
}

//�A���C�����g
void Poly<Degree>::operator delete(void* pv) {
	_mm_free(pv);
}

//LM��LM��index����
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

//set_LM�̌ザ��Ȃ��Ɠ����Ȃ�
inline void Poly<Degree>::set_LMdeg()
{
	_LMdeg = _degree.index_to_deg(_LMdeg_index);
}