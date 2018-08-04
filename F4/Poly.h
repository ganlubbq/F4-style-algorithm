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

//set_LM�̌ザ��Ȃ��Ɠ����Ȃ�
inline void Poly<Degree>::set_LMdeg()
{
	LMdeg = Degree.index_to_deg(LMdeg_index);
}