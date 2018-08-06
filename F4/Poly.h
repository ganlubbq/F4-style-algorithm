#pragma once
#include <vector>
#define D1
#ifdef D1
#include "Degree_table.h"
#endif //D1

using namespace std;

static const int variables = 10;

/*基準となるクラス　これを継承することで機能を実現する予定
Degreeは一般的なもののみ対応　過度な高速化は汎用性を失う
*/
class Poly
{
public:
	Poly(vector<unsigned char> &coeff);

	//variable
	int _LM;
	vector<unsigned char> _LMdeg;
	int _LMdeg_index;
	vector<unsigned char> _Coeff;

//mainで定義
	static int _Max_degree;
#ifdef D1
	static Degree_table _Degree;
#endif //D1
//

	//function
	void* operator new(size_t size);
	void operator delete(void* pv);

	virtual void set_LM();
	virtual void set_LMdeg();
	//inline virtual void set_LMdeg_index();
};

