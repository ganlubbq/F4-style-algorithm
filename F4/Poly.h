#pragma once
#include <vector>
//degree切り替え
#define D1
#ifdef D1
#include "Degree_table.h"
#endif //D1

using namespace std;

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

	//operator
	void* operator new(size_t size);
	void operator delete(void* pv);

	//function
	virtual void set_LM();
	virtual void set_LMdeg();
};

