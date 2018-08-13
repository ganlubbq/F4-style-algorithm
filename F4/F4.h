#pragma once

#include <sstream>
#include <string>
#include<vector>

using namespace std;

//F4アルゴリズム　ファイル読み込みは標準実装　大きな書き換えは継承すべし
template <class GF, class Deci,class Spol,class Red,class LB>
class F4
{
public:
	F4(string filename);

	//variables
	string file_name;
	static int _Variables;
	static Deci _Decision;
	static Spol _Spoly;
	static Red _Red;
	static LB _LB;
	static int _Parallel_div;
	GF _GFf;
	vector<GF> _Answer;
	vector<GF> _Equations;

	//function
	void printvec(vector<unsigned char> a);
	void printvec(vector<int> a);
	void file_read();
	int var_deg_comb(int n, int r);
	void F4_style();
};

template <class GF, class Deci, class Spol, class Red, class LB>
F4<GF,Deci, Spol,Red,LB>::F4(string filename)
{
	file_name = filename;
	file_read();
}

template <class GF, class Deci, class Spol, class Red, class LB>
void F4<GF, Deci, Spol, Red, LB>::printvec(vector<unsigned char> vec)
{
	cout << "[ ";
	for (int i = 0; i < vec.size(); i++)
	{
		cout << (int)vec[i] << " ";
	}
	cout << "]" << endl;
}

template <class GF, class Deci, class Spol, class Red, class LB>
void F4<GF, Deci, Spol, Red, LB>::printvec(vector<int> vec)
{
	cout << "[ ";
	for (int i = 0; i < vec.size(); i++)
	{
		cout << (int)vec[i] << " ";
	}
	cout << "]" << endl;
}

template <class GF, class Deci, class Spol, class Red, class LB>
void F4<GF,Deci, Spol, Red, LB>::file_read()
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
template <class GF, class Deci, class Spol, class Red, class LB>
inline int F4<GF,Deci, Spol, Red, LB>::var_deg_comb(int n, int r) {//n変数r次多項式までの全項数
	int ans = 1;
	n += r;
	for (int div = 1; div <= r; ++div, --n) {
		ans = ans * n / div;
	}
	return ans;
}

//F4アルゴリズム
template <class GF, class Deci, class Spol, class Red, class LB>
inline void F4<GF, Deci, Spol, Red, LB>::F4_style()
{
	_Answer.resize(_Variables + 1);
	int count = 0;
	/*for (int i = 0; i < _Equations.size(); i++)
	{
		printvec(_Equations[i]._Coeff);
	}*/
	
	//Spoly絞り込み init
	_Decision.decision(_Equations);
	/*for (int i = 0; i < _Decision._D.size(); i++)
	{
		printvec(_Decision._D[i]);
	}*/

	while (_Decision._D.size() > 0)
	{
		_Spoly.spoly_erase();
		_Spoly.calc_Spoly(_Equations, _Decision._D);
		/*for (int i = 0;i < _Spoly._Spolies.size();i++)
		{
			printvec(_Spoly._Spolies[i]._Coeff);
		}*/
		_Decision.d_erase();
		_Red.calc_red(_Spoly._Spolies,_Equations);
		//cout << _Red._Reds.size() << "REd" << endl;
		//printvec(_Red._Reds[0]._Coeff);

		//ここSpolyうまく使えば消せる?
		_Spoly._Spolies.insert(_Spoly._Spolies.end(), _Red._Reds.begin(), _Red._Reds.end()); // 連結 S = S or Red
		/*for (int i = 0; i < _Spoly._Spolies.size(); i++)
		{
			printvec(_Spoly._Spolies[i]._Coeff);
		}*/
		_LB.calc_LB(_Spoly._Spolies);

		/*cout << "LB" << endl;
		for (int i = 0; i < _Spoly._Spolies.size(); i++)
		{
			printvec(_Spoly._Spolies[i]._Coeff);
		}*/
		for (int i = 0; i < _Spoly._Spolies.size(); i++)
		{
			//0多項式判定
			if (_Spoly._Spolies[i]._LMdeg_index != -1)
			{
				bool flag = true;
				for (int j = 0;j <_Red._Reds.size();j++)
				{
					if (_Decision.veceq(_Spoly._Spolies[i]._LMdeg, _Red._Reds[j]._LMdeg))
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					_Equations.push_back(_Spoly._Spolies[i]);
					if (_Spoly._Spolies[i]._LMdeg_index <= _Variables)
					{
						_Answer[_Spoly._Spolies[i]._LMdeg_index] = _Spoly._Spolies[i];
						count++;
						if (count == _Variables) break;
					}
					//_LB.Gauss_rev(_Equations);
					_Decision.Gebauer_Moller_mono(_Equations);
				}
			}
			if (count == _Variables) break;
		}
		if (count == _Variables) break;
		/*cout << "add" << endl;
		for (int i = 0; i < _Equations.size(); i++)
		{
			printvec(_Equations[i]._Coeff);
		}*/
		_Decision.Buchberger(_Equations);
		/*for (int i = 0; i < _Decision._D.size(); i++)
		{
			printvec(_Decision._D[i]);
		}*/
		cout << _Decision._D.size() << endl;
	}
	for (int i = 1; i < _Answer.size(); i++)
	{
		printvec(_Answer[i]._Coeff);
	}
}