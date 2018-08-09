#pragma once

#include <sstream>
#include <string>
#include<vector>


using namespace std;

//F4アルゴリズム　ファイル読み込みは標準実装　大きな書き換えは継承すべし
template <class GF, class Deci>
class F4
{
public:
	F4(string filename);

	//variables
	string file_name;
	static int _Variables;
	static Deci _Decision;
	GF _GFf;

	vector<GF> _Spolies;
	vector<GF> _Equations;

	//function
	template<typename T>
	void printvec(T vec)
	{
		{
			cout << "[ ";
			for (int i = 0; i < vec.size(); i++)
			{
				cout << (int)vec[i] << " ";
			}
			cout << "]" << endl;
		}
	}

	void file_read();
	int var_deg_comb(int n, int r);
	void F4_style();
};

template <class GF, class Deci>
F4<GF,Deci>::F4(string filename)
{
	file_name = filename;
	file_read();
}

template <class GF, class Deci>
void F4<GF,Deci>::file_read()
{
	FILE *fp;
	fp = fopen(file_name.c_str(), "r");
	unsigned int coeff;

	fscanf(fp, "%d", &coeff);

	if (coeff != _Variables) {
		puts("VAL_NUM is not correct");
		system("pause");
	}

	unsigned char coeff_char;
	int coeff_int;

	//MQcharennge型しか読み込めない変更
	for (int i = 0; i < 2 * _Variables; i++)
	{
		vector<unsigned char> temp;
		for (int j = 0; j < var_deg_comb(_Variables, 2); j++) {
			if(_GFf._DX == "d")	fscanf(fp, "%d", &coeff_int);
			else if (_GFf._DX == "x") fscanf(fp, "%x", &coeff_int);
			temp.push_back((unsigned char)coeff_int);
		}
		reverse(temp.begin(), temp.end());
		_Equations.push_back(GF(temp));
	}
}

//n変数r次多項式までの全項数
template <class GF, class Deci>
inline int F4<GF,Deci>::var_deg_comb(int n, int r) {//n変数r次多項式までの全項数
	int ans = 1;
	n += r;
	for (int div = 1; div <= r; ++div, --n) {
		ans = ans * n / div;
	}
	return ans;
}

template <class GF, class Deci>
inline void F4<GF, Deci>::F4_style()
{
	for (auto itr = _Equations.begin(); itr != _Equations.end(); itr++)
	{
		printvec((*itr)._Coeff);
	}
	//init
	_Decision.decision(_Equations);

	for (int i = 0; i < _Decision._D_sort.size(); i++)
	{
		for (auto itr = _Decision._D_sort[i].begin(); itr != _Decision._D_sort[i].end(); itr++)
		{
			printvec(*itr);
		}
	}
}