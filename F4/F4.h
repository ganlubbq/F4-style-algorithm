#pragma once

#include <sstream>
#include <string>
#include<vector>
#include<time.h>
#include <algorithm>
#include "memory.h"

using namespace std;

//F4アルゴリズム　ファイル読み込みは標準実装　大きな書き換えは継承すべし
template <class GF, class Deci, class Spol, class Red, class LB>
class F4
{
public:
	F4(string filename, int variables, string writing_file, int seiki, string all,string gauss,string LB_file,string red_file,string decision_file);

	//file
#pragma region
	string file_name;
	string w_file_name;
	string all_file_name;
	string gauss_file_time;
	string LB_file_time_size;
	string red_file_time_size;
	string decision_file_time_size;
#pragma endregion
	//variables
	int _Variables;
	static Deci _Decision;
	static Spol _Spoly;
	static Red _Red;
	static LB _LB;
	static int _Parallel_div;
	int _Seiki;
	int count;
	GF _GFf;
	vector<GF> _Answer;
	vector<GF> _Equations;
	vector<int> _LMplace;
	vector<int> _LMplace_new;
	vector<int> de_i;//gauss掃き出しに必要なindex
	//clock_t all_time;

	//function
	void printvec(vector<unsigned char> a);
	void printvec(vector<int> a);
	void file_read();
	int var_deg_comb(int n, int r);
	void F4_style();
	vector<int> seikika(vector<GF> &G);
	void Equation_reduction(GF &Spoly);
	bool one_erace_check3();
};

template <class GF, class Deci, class Spol, class Red, class LB>
F4<GF, Deci, Spol, Red, LB>::F4(string filename, int variables, string writing_file, int seikia, string all,string gauss, string LB_file, string red_file,string decision_file)
{
	_Variables = variables;
	_Seiki = seikia;
	file_name = filename;
	w_file_name = writing_file;
	all_file_name = all;
	gauss_file_time = gauss;
	decision_file_time_size = decision_file;
	LB_file_time_size = LB_file;
	red_file_time_size = red_file;

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
void F4<GF, Deci, Spol, Red, LB>::file_read()
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
			if (_GFf._DX == "d")	fscanf(fp, "%d", &coeff_int);
			else if (_GFf._DX == "x") fscanf(fp, "%x", &coeff_int);
			temp.push_back((unsigned char)coeff_int);
		}
		reverse(temp.begin(), temp.end());
		_Equations.push_back(GF(temp));
	}
}

//n変数r次多項式までの全項数
template <class GF, class Deci, class Spol, class Red, class LB>
inline int F4<GF, Deci, Spol, Red, LB>::var_deg_comb(int n, int r) {//n変数r次多項式までの全項数
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
//#define DEBUG
//#define ONE_ERASE
	//writing file
#pragma region
	std::ofstream writing_file;
	writing_file.open(w_file_name, std::ios::out);

	std::ofstream LB_file_time_size_;
	LB_file_time_size_.open(LB_file_time_size, std::ios::app);

	std::ofstream gauss_file_time_;
	gauss_file_time_.open(gauss_file_time, std::ios::app);

	ofstream red_file_time_size_;
	red_file_time_size_.open(red_file_time_size, ios::app);

	ofstream decision_file_time_size_;
	decision_file_time_size_.open(decision_file_time_size, ios::app);
#pragma endregion

	auto start = clock();
	//all_time = 0;

	_Answer.resize(_Variables + 1);
	count = 0;
	int a_count = 0;
	bool reset = false;
	bool flag_3 = false;


#ifdef DEBUG
	cout << "init" << endl;
	/*for (int x = 0; x < _Equations.size(); x++)
	{
		printvec(_Equations[x]._Coeff);
	}
	cout << endl;*/
#endif // DEBUG

	//old
	if(_Seiki != 3 && _Seiki != 4) _LB.Gauss_rev(_Equations);
	//only
	else  if(_Seiki == 3)
	{
	//Equations resize 4じ　あんど　gaussrev and lm配置
		_LB.Gauss_rev_only_init(_Equations,_Variables,_LMplace);
	}
	else  if (_Seiki == 4)
	{
		//Equations resize 5じ　あんど　gaussrev and lm配置
		_LB.Gauss_rev_only_init_4(_Equations, _Variables, _LMplace);
	}

	auto decision_start = clock();
	if(_Seiki != 3 && _Seiki != 4) _Decision.decision(_Equations);
	else if(_Seiki == 3)
	{
		_Decision.decision_kai(_LMplace);
	}
	else
	{
		_Decision.decision_0_kai(_LMplace,_Equations);
	}

	auto decision_end =  clock();
	decision_file_time_size_ << _Equations.size() << "\t" << decision_end - decision_start << endl;

#ifdef DEBUG
	cout << "init2" << endl;
/*	for (int x = 0; x < _Equations.size(); x++)
	{
		printvec(_Equations[x]._Coeff);
	}
	cout << endl;*/
#endif // DEBUG

	//次数ごと　pが次数　2次から
	for (int p = 1; p < _Decision._D_sort.size(); p++)
	{
		reset = false;
		
		while (_Decision._D_sort[p].size() > 0)
		{
			cout << "p" << p << endl;
			cout << "size" << _Decision._D_sort[p].size() << endl;
			writing_file << "p" << p << endl;
			writing_file << "size" << _Decision._D_sort[p].size() << endl;
			LB_file_time_size_ << p << "\t" << _Decision._D_sort[p].size() << "\t";

			_Spoly.spoly_erase();
			vector<vector<int>> DD;
			if (_Decision._D_sort[p].size() > _Parallel_div)
			{
				//DD.resize(_Parallel_div);
				DD.insert(DD.end(), _Decision._D_sort[p].begin(), _Decision._D_sort[p].begin() + _Parallel_div);
				_Decision._D_sort[p].erase(_Decision._D_sort[p].begin(), _Decision._D_sort[p].begin() + _Parallel_div);
#ifdef DEBUG
				cout << "DD" << endl;
			/*	for (int x = 0; x < DD.size(); x++)
				{
					printvec(DD[x]);
				}
				cout << endl;
				system("pause");*/
#endif // DEBUG
				_Spoly.calc_Spoly(_Equations, DD);
			}
			else
			{
#ifdef DEBUG
				cout << "D" << endl;
				/*for (int x = 0; x < _Decision._D_sort[p].size(); x++)
				{
					printvec(_Decision._D_sort[p][x]);
				}
				cout << endl;*/
#endif // DEBUG
				_Spoly.calc_Spoly(_Equations, _Decision._D_sort[p]);
				_Decision.d_sort_erase(p);
			}

			red_file_time_size_ << _Spoly._Spolies.size() << "\t" << _Equations.size() << "\t";
			
#ifdef DEBUG
			cout << "Sp" << endl;
		/*	for (int x = 0; x < _Spoly._Spolies.size(); x++)
			{
				printvec(_Spoly._Spolies[x]._Coeff);
			}
			cout << endl;*/
#endif // DEBUG

			auto red_start = clock();
			_Red.calc_red(_Spoly._Spolies, _Equations);
			auto red_end = clock();
			red_file_time_size_ << red_end - red_start << endl;
#ifdef DEBUG
			cout << "Red" << endl;
			/*for (int x = 0; x < _Red._Reds.size(); x++)
			{
				printvec(_Red._Reds[x]._Coeff);
			}
			cout << endl;*/
#endif // DEBUG

			//ここSpolyうまく使えば消せる?
			_Spoly._Spolies.insert(_Spoly._Spolies.end(), _Red._Reds.begin(), _Red._Reds.end()); // 連結 S = S or Red

#ifdef DEBUG
			cout << "Sp or Red" << endl;
			/*for (int x = 0; x < _Spoly._Spolies.size(); x++)
			{
				printvec(_Spoly._Spolies[x]._Coeff);
			}
			cout << endl;*/
#endif // DEBUG

			a_count += 1;
			LB_file_time_size_ << _Spoly._Spolies.size() << "\t";
			auto LB_start = clock();
			_LB.calc_LB(_Spoly._Spolies);
			auto LB_end = clock();
			LB_file_time_size_ << LB_end - LB_start << endl;

#ifdef DEBUG
			cout << "LB" << endl;
			/*for (int x = 0; x < _Spoly._Spolies.size(); x++)
			{
				printvec(_Spoly._Spolies[x]._Coeff);
			}
			cout << endl;*/

#endif // DEBUG

			for (int i = 0; i < _Spoly._Spolies.size(); i++)
			{
				//0多項式判定
				if (_Spoly._Spolies[i]._LMdeg_index != -1)
				{
					bool flag = true;
#ifdef ONE_ERASE
					if (_Spoly._Spolies[i]._LMdeg_index <= _Variables)
					{
						cout << "ONE_Erace" << endl;
						if (_Answer[_Spoly._Spolies[i]._LMdeg_index]._Coeff.size() == 0)
						{
							_Answer[_Spoly._Spolies[i]._LMdeg_index] = _Spoly._Spolies[i];
							count++;
						}
						if (count == _Variables) break;
					}
					if (count == _Variables) break;
#endif //ONE_ERASE
					for (int j = 0; j < _Red._Reds.size(); j++)
					{
						if (_Decision.veceq(_Spoly._Spolies[i]._LMdeg, _Red._Reds[j]._LMdeg))
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						if(_Seiki != 3 && _Seiki != 4)_Equations.push_back(_Spoly._Spolies[i]);
						else
						{
							//LMplace_newにindex追加まで _Spolie[i]は割られて変形されることあり
							Equation_reduction(_Spoly._Spolies[i]);
							flag_3 = one_erace_check3();
						}
						//_Seiki == 3も大丈夫

						if (_Seiki == 0)
						{
							_Decision.Gebauer_Moller_mono(_Equations);
						}
						else if (_Seiki == 1 || _Seiki == 2)
						{
							de_i.push_back(_Equations.size() - 1);
						}
					}
				}
				if (count == _Variables) break;
				if (flag_3) break;
			}
			if (count == _Variables) break;
			if (flag_3) break;
			if (_Seiki == 1 || _Seiki == 2)
			{
				if (_Seiki == 1)
				{
					vector<int> de_temp;
					auto gauss_start = clock();
					
					de_temp = _LB.Gauss_rev_Eq(_Equations);
					auto gauss_end = clock();
					gauss_file_time_  << _Equations.size() << "\t" << gauss_end - gauss_start << endl;
					for (int sss = 0; sss < de_temp.size(); sss++)
					{
						de_i.push_back(de_temp[sss]);
					}
				}
				else if (_Seiki == 2)
				{
					vector<int> de_temp;
					auto gauss_start = clock();
					de_temp = seikika(_Equations);
					auto gauss_end = clock();
					gauss_file_time_ << _Equations.size() << "\t" << gauss_end - gauss_start << endl;
					for (int sss = 0; sss < de_temp.size(); sss++)
					{
						de_i.push_back(de_temp[sss]);
					}
				}
#ifdef ONE_ERASE
				for (int n = 0; n < _Equations.size(); n++)
				{
					if (_Equations[n]._LMdeg_index != -1)
					{
						if (_Equations[n]._LMdeg_index <= _Variables)
						{
							if (_Answer[_Equations[n]._LMdeg_index]._Coeff.size() == 0)
							{
								cout << "ONE_ERACE" << endl;
								_Answer[_Equations[n]._LMdeg_index] = _Equations[n];
								count++;
							}
							if (count == _Variables)
							{
								break;
							}
						}
					}
				}
				if (count == _Variables) break;
#endif //ONE_ERASE
				reset = true;
			}
			//////////////////////////Spoly関連ここまで/////////////////////////////////////////////////
			if (_Seiki == 1 || _Seiki == 2)
			{
				for (auto itr = de_i.begin(); itr != de_i.end(); itr++)
				{
					cout << (*itr) << endl;
					printvec(_Equations[*itr]._Coeff);
					printvec(_Equations[0]._Coeff);
					_Decision.Gebauer_Moller_num(_Equations, *itr);
				}
				
				de_i.resize(0);
			}
			//げばうわめらー & ぶっふばーがー
			else if (_Seiki == 3)
			{
				std::sort(_LMplace_new.begin(), _LMplace_new.end());
				//改良可能　newを分ければ　ブッフバーガーの判定削れる
				_Decision.decision_kai_2(_LMplace_new,_LMplace);
				//LMplace oldにnewをつけ加える
				_LMplace.insert(_LMplace.end(), _LMplace_new.begin(), _LMplace_new.end());
				std::sort(_LMplace.begin(), _LMplace.end());
				_LMplace_new.resize(0);
			}
			else if (_Seiki == 4)
			{
				std::sort(_LMplace_new.begin(), _LMplace_new.end());
				//改良可能　newを分ければ　ブッフバーガーの判定削れる
				_Decision.decision_0_kai_2(_LMplace_new, _LMplace,_Equations);
				//LMplace oldにnewをつけ加える
				_LMplace.insert(_LMplace.end(), _LMplace_new.begin(), _LMplace_new.end());
				std::sort(_LMplace.begin(), _LMplace.end());
				_LMplace_new.resize(0);
			}

			if(_Seiki != 3 && _Seiki != 4)_Decision.Buchberger(_Equations);

#ifdef DEBUG
			cout << "init3" << endl;
			for (int x = 0; x < _Equations.size(); x++)
			{
				printvec(_Equations[x]._Coeff);
			}
			cout << endl;
#endif // DEBUG
		}
		if (reset == true) p = 1;
		if (count == _Variables) break;
		if (flag_3) break;
	}

	auto end = clock();

	std::ofstream writing_all;
	writing_all.open(all_file_name, std::ios::app);

	writing_file << endl;
	writing_file << end - start << endl;
	writing_all << end - start << endl;
	for (int i = 0; i <= a_count;i++)
	{
		writing_all << endl;
	}
	gauss_file_time_ << endl << endl << endl;
	decision_file_time_size_ << endl << endl;
	LB_file_time_size_ << endl << endl;
	red_file_time_size_ << endl << endl;

#ifdef ONE_ERASE
	_LB.Gauss_rev_fin(_Answer);
	for (int i = 1; i < _Answer.size(); i++)
	{
		_Answer[i]._Coeff.resize(_Variables + 1);

		writing_file << "[";
		for (int j = 0; j < _Answer[i]._Coeff.size(); j++)
		{
			writing_file << (int)_Answer[i]._Coeff[j] << " ";
		}
		writing_file << "]" << endl;
	}
#endif // !ONE_ERASE

#ifndef ONE_ERASE
	if (_Seiki != 3 && _Seiki != 4) {
		for (int i = 0; i < _Equations.size(); i++)
		{
			writing_file << "[";
			for (int j = 0; j < _Equations[i]._Coeff.size(); j++)
			{
				writing_file << (int)_Equations[i]._Coeff[j] << " ";
			}
			writing_file << "]" << endl;
		}
	}
	else if(_Seiki == 3){
		for (int i = 1; i <= _Variables; i++)
		{
			writing_file << "[";
			_Equations[i]._Coeff.resize(_Variables + 1);
			for (int j = 0; j < _Equations[i]._Coeff.size(); j++)
			{
				writing_file << (int)_Equations[i]._Coeff[j] << " ";
			}
			writing_file << "]" << endl;
		}
	}
	else if (_Seiki == 4)
	{
		//vector size調整
		int max_length = 0;
		for (int i = 1; i <= _Variables; i++)
		{
			if (max_length < _Equations[i]._LMdeg_index + 1) max_length = _Equations[i]._LMdeg_index + 1;
		}
		for (int i = 1; i <= _Variables; i++)
		{
			_Equations[i]._Coeff.resize(max_length);
			_Equations[i]._Coeff_size = max_length;
			_Equations[i]._Div_single_size = _Equations[i]._Coeff_size / single_size;
		}

		//gauss
		for (int i = 1; i <= _Variables; i++)
		{
			int index = _Equations[i]._LMdeg_index;
			for (int j = 1; j <= _Variables; j++)
			{
				if (i == j) continue;
				if (_Equations[j]._Coeff[index] != 0)
				{
					GF temp = _Equations[i];
					temp * (_GFf._Add_inverse[_Equations[j]._Coeff[index]]);
					_Equations[j] + temp;
				}
			}
		}

		for (int i = 1; i <= _Variables; i++)
		{
			writing_file << "[";
			_Equations[i]._Coeff.resize(_Variables + 1);
			for (int j = 0; j < _Equations[i]._Coeff.size(); j++)
			{
				writing_file << (int)_Equations[i]._Coeff[j] << " ";
			}
			writing_file << "]" << endl;
		}
	}
#endif // ONE_ERASE

	writing_file.close();
	writing_all.close();
	gauss_file_time_.close();
	LB_file_time_size_.close();
	red_file_time_size_.close();
	decision_file_time_size_.close();
}

//LC = 1前提
template <class GF, class Deci, class Spol, class Red, class LB>
inline vector<int> F4<GF, Deci, Spol, Red, LB>::seikika(vector<GF> &G)
{
	vector<int> result;
	bool flag = true;
	while (flag)
	{
		flag = false;
		for (int i = 0; i < G.size(); i++)
		{
			if (G[i]._LMdeg_index == -1) continue;
			vector<unsigned char> i_LMdeg = G[i]._LMdeg;
#pragma omp parallel for reduction(||:flag)
			for (int j = 0; j < G.size(); j++)
			{
				if (i == j) continue;
				if (G[j]._LMdeg_index == -1) continue;

				for (int k = G[j]._LMdeg_index; k > 0; k--)
				{

					if (G[j]._Coeff[k] != 0)
					{
						//	cout << (int)G[j]._Coeff[k] << endl;
						GF temp_i = G[i];
						if (_GFf._Degree.reducible(i_LMdeg, _GFf._Degree.index_to_degree(k)))
						{
							if (_GFf._Degree.calc_total_deg(i_LMdeg) < _GFf._Degree.calc_total_deg(_GFf._Degree.index_to_degree(k)))
							{
								flag = true;
							}
							unsigned char real_temp = (G[j]._Coeff[k] * 30) % 31;
							//	cout << (int)real_temp << endl;

							temp_i * real_temp;
							vector<unsigned char> temp_deg = _GFf._Degree.vec_sub(_GFf._Degree.index_to_degree(k), i_LMdeg);
							temp_i * temp_deg;

							if (G[j]._Coeff_size < temp_i._Coeff_size)
							{
								G[j]._Coeff.resize(temp_i._Coeff_size);
								G[j]._Coeff_size = temp_i._Coeff_size;
								G[j]._Div_single_size = G[j]._Coeff_size / single_size;
							}
							G[j] + temp_i;
							result.push_back(j);
						}
					}
				}
				if (G[j]._LMdeg_index == -1) continue;

				G[j] * (_GFf._Inverse[G[j]._LM]);
				if (G[j]._LMdeg_index <= _Variables)
				{
					if (_Answer[G[j]._LMdeg_index]._Coeff.size() == 0)
					{
						_Answer[G[j]._LMdeg_index] = G[j];
						count++;
					}
					if (count == _Variables) break;
				}
				if (count == _Variables) break;
			}
			if (count == _Variables) break;
		}
		if (count == _Variables) break;
	}
	return result;
}

template <class GF, class Deci, class Spol, class Red, class LB>
inline void F4<GF, Deci, Spol, Red, LB>::Equation_reduction(GF &Spoly)
{
	while (true)
	{
		if (Spoly._LMdeg_index == -1) break;
		//該当箇所存在しない
		else if (_Equations[Spoly._LMdeg_index]._LMdeg_index == -1)
		{
			_Equations[Spoly._LMdeg_index] = Spoly;
			_LMplace_new.push_back(Spoly._LMdeg_index);
			break;
		}
		//reduction
		else
		{
			unsigned char inv = 30;
			GF temp = _Equations[Spoly._LMdeg_index];
			temp * inv;
			Spoly + temp;
		}
	}
}

template <class GF, class Deci, class Spol, class Red, class LB>
inline bool F4<GF, Deci, Spol, Red, LB>::one_erace_check3()
{
	for (int i = 1; i <= _Variables; i++)
	{
		if (_Equations[i]._LMdeg_index == -1) return false;
	}
	return true;
}