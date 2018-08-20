#pragma once

#include <sstream>
#include <string>
#include<vector>
#include<time.h>

using namespace std;

//F4アルゴリズム　ファイル読み込みは標準実装　大きな書き換えは継承すべし
template <class GF, class Deci, class Spol, class Red, class LB>
class F4
{
public:
	F4(string filename, int variables, string writing_file, int seiki, string all,string gauss,string LB_file,string red_file,string decision_file);

	//variables
	string file_name;
	string w_file_name;
	string all_file_name;
	string gauss_file_time;
	string LB_file_time_size;
	string red_file_time_size;
	string decision_file_time_size;
	int _Variables;
	static Deci _Decision;
	static Spol _Spoly;
	static Red _Red;
	static LB _LB;
	static int _Parallel_div;
	int _Seiki;
	GF _GFf;
	vector<GF> _Answer;
	vector<GF> _Equations;
	vector<int> de_i;//gauss掃き出しに必要なindex
	//clock_t all_time;

	//function
	void printvec(vector<unsigned char> a);
	void printvec(vector<int> a);
	void file_read();
	int var_deg_comb(int n, int r);
	void F4_style();
	vector<int> seikika(vector<GF> &G);
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

	auto start = clock();
	//all_time = 0;

	_Answer.resize(_Variables + 1);
	int count = 0;
	int a_count = 0;
	bool reset = false;

	if (_Seiki != 0) _LB.Gauss_rev(_Equations);

	auto decision_start = clock();
	//Spoly絞り込み init
	_Decision.decision(_Equations);
	auto decision_end =  clock();
	decision_file_time_size_ << _Equations.size() << "\t" << decision_end - decision_start << endl;

	for (int p = 1; p < _Decision._D_sort.size(); p++)
	{
		reset = false;
		
		while (_Decision._D_sort[p].size() > 0)
		{
			//writing_file << _Decision._D_sort[p].size() << endl;
			//cout << _Decision._D_sort[p].size() << endl;

			cout << "p" << p << endl;
			cout << "size" << _Decision._D_sort[p].size() << endl;

			writing_file << "p" << p << endl;
			writing_file << "size" << _Decision._D_sort[p].size() << endl;

			LB_file_time_size_ << p << "\t" << _Decision._D_sort[p].size() << "\t";

			_Spoly.spoly_erase();
			/*vector<vector<int>> DD;
			if (_Decision._D_sort[p].size() > _Parallel_div)
			{
				//DD.resize(_Parallel_div);
				DD.insert(DD.end(), _Decision._D_sort[p].begin(), _Decision._D_sort[p].begin() + _Parallel_div);
				_Decision._D_sort[p].erase(_Decision._D_sort[p].begin(), _Decision._D_sort[p].begin() + _Parallel_div);
				_Spoly.calc_Spoly(_Equations, DD);
			}
			else
			{*/
				_Spoly.calc_Spoly(_Equations, _Decision._D_sort[p]);
				_Decision.d_sort_erase(p);
			//}

			red_file_time_size_ << _Spoly._Spolies.size() << "\t" << _Equations.size() << "\t";
			auto red_start = clock();
			_Red.calc_red(_Spoly._Spolies, _Equations);
			auto red_end = clock();
			red_file_time_size_ << red_end - red_start << endl;

			//ここSpolyうまく使えば消せる?
			_Spoly._Spolies.insert(_Spoly._Spolies.end(), _Red._Reds.begin(), _Red._Reds.end()); // 連結 S = S or Red

			a_count += 1;
			LB_file_time_size_ << _Spoly._Spolies.size() << "\t";
			auto LB_start = clock();
			_LB.calc_LB(_Spoly._Spolies);
			auto LB_end = clock();
			LB_file_time_size_ << LB_end - LB_start << endl;

			for (int i = 0; i < _Spoly._Spolies.size(); i++)
			{
				//0多項式判定
				if (_Spoly._Spolies[i]._LMdeg_index != -1)
				{
					bool flag = true;
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
						_Equations.push_back(_Spoly._Spolies[i]);
						//if (_Seiki == false)
						//{
						/*if (_Spoly._Spolies[i]._LMdeg_index <= _Variables)
						{
							_Answer[_Spoly._Spolies[i]._LMdeg_index] = _Spoly._Spolies[i];
							count++;
							if (count == _Variables) break;
						}*/
						if (_Seiki == 0 )
						{
							_Decision.Gebauer_Moller_mono(_Equations);
						}
						else if (_Seiki != 0)
						{
							de_i.push_back(_Equations.size() - 1);
						}
						//}
					}
				}
				//if (count == _Variables) break;
			}
		//	if (count == _Variables) break;
			if (_Seiki != 0)
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
				else if (_Seiki == 3);

				/*for (int n = 0; n < _Equations.size(); n++)
				{
					if (_Equations[n]._LMdeg_index != -1)
					{
						if (_Equations[n]._LMdeg_index <= _Variables)
						{
							_Answer[_Equations[n]._LMdeg_index] = _Equations[n];
							count++;
							if (count == _Variables) break;
						}
					}
				}*/
			//	if (count == _Variables) break;
				reset = true;
			}
			if (_Seiki != 0)
			{
				for (auto itr = de_i.begin(); itr != de_i.end(); itr++)
				{
					_Decision.Gebauer_Moller_num(_Equations, *itr);
				}
				
				de_i.resize(0);
			}
			_Decision.Buchberger(_Equations);
		}
		if (reset == true) p = 1;
		//if (count == _Variables) break;
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
	for (int i = 1; i < _Answer.size(); i++)
	{
		_Answer[i].set_LM();
		_Answer[i].set_LMdeg();
		_Answer[i] * _GFf._Inverse[_Answer[i]._LM];
		_Answer[i]._Coeff.resize(_Variables + 1);
		writing_file << "[";
		for (int j = 0; j < _Answer[i]._Coeff.size(); j++)
		{
			writing_file << (int)_Answer[i]._Coeff[j] << " ";
		}
		writing_file << "]" << endl;
	}
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
			}
		}
	}
	return result;
}