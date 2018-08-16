#pragma once
#include<vector>

template<class GF>
class Decision
{
public:
	Decision() {};
	//variables
	//vector<vector<int>> _D;
	vector<vector<vector<int>>> _D_sort;//total degree順
	GF _GFd;

	//function
	void decision(vector<GF> &G);
	void Gebauer_Moller(vector<GF> &G);
	void Gebauer_Moller_mono(vector<GF> &G);
	void Buchberger(vector<GF> &G);
	//void sort_D(vector<GF> &G);
	//void d_erase();
	void d_sort_erase(int n);

	//&使うとバグる　正確には計算結果をそのまま引数にとれない
	virtual bool veceq(vector<unsigned char> f, vector<unsigned char> g);

};

template <class GF>
inline void Decision<GF>::decision(vector<GF> &G)
{
	_D_sort.resize(_GFd._Max_degree + 1);
	Gebauer_Moller(G);
	Buchberger(G);
	//sort_D(G);
}


template<class GF>
inline void Decision<GF>::Gebauer_Moller(vector<GF> &G)
{
	for (int i = 0; i < G.size(); i++)
	{
		for (int j = i + 1; j < G.size(); j++)
		{
			//F
			int flag = 0;
			int k = 0;
			vector<unsigned char> T_ij = _GFd._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg);
			while (k < i)
			{
				if (veceq(_GFd._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), T_ij))
				{
					flag++;
					int temp_T_ij_deg = _GFd._Degree.calc_total_deg(T_ij);
					if (_GFd._Max_degree < temp_T_ij_deg)
					{
						_GFd._Max_degree = temp_T_ij_deg;
						_GFd._Degree.update_degree(temp_T_ij_deg);
						_D_sort.resize(_GFd._Max_degree + 1);
					}
					_D_sort[_GFd._Degree.calc_total_deg(T_ij)].push_back({ i,j });
					break;
				}
				k++;
			}

			//M
			if (flag == 0)
			{
				int k = 0;
				while (k < j)
				{
					if (_GFd._Degree.reducible(G[k]._LMdeg, T_ij))
					{
						if (!veceq(_GFd._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), T_ij))
						{
							flag++;
							int temp_T_ij_deg = _GFd._Degree.calc_total_deg(T_ij);
							if (_GFd._Max_degree < temp_T_ij_deg)
							{
								_GFd._Max_degree = temp_T_ij_deg;
								_GFd._Degree.update_degree(temp_T_ij_deg);
								_D_sort.resize(_GFd._Max_degree + 1);
							}
							_D_sort[_GFd._Degree.calc_total_deg(T_ij)].push_back({ i,j });
							break;
						}
					}
					k++;
				}
			}

			//B
			if (flag == 0)
			{
				int k = G.size() - 1;
				while (k > j)
				{
					if (_GFd._Degree.reducible(G[k]._LMdeg, T_ij))
					{
						if (!veceq(_GFd._Degree.LCM(G[i]._LMdeg, G[k]._LMdeg), T_ij))
						{
							if (!veceq(_GFd._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), T_ij))
							{
								int temp_T_ij_deg = _GFd._Degree.calc_total_deg(T_ij);
								if (_GFd._Max_degree < temp_T_ij_deg)
								{
									_GFd._Max_degree = temp_T_ij_deg;
									_GFd._Degree.update_degree(temp_T_ij_deg);
									_D_sort.resize(_GFd._Max_degree + 1);
								}
								_D_sort[_GFd._Degree.calc_total_deg(T_ij)].push_back({ i,j });
								break;
							}
						}
					}
					k--;
				}
			}
		}
	}
}

//多項式一個追加したときのみ動く
template<class GF>
inline void Decision<GF>::Gebauer_Moller_mono(vector<GF> &G)
{
	for (int i = 0; i < G.size() - 1; i++)
	{
		for (int j = G.size() - 1; j < G.size(); j++)
		{
			//F
			int flag = 0;
			int k = 0;
			vector<unsigned char> T_ij = _GFd._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg);
			while (k < i)
			{
				if (veceq(_GFd._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), T_ij))
				{
					flag++;
					int temp_T_ij_deg = _GFd._Degree.calc_total_deg(T_ij);
					if (_GFd._Max_degree < temp_T_ij_deg)
					{
						_GFd._Max_degree = temp_T_ij_deg;
						_GFd._Degree.update_degree(temp_T_ij_deg);
						_D_sort.resize(_GFd._Max_degree + 1);
					}
					_D_sort[_GFd._Degree.calc_total_deg(T_ij)].push_back({ i,j });
					break;
				}
				k++;
			}

			//M
			if (flag == 0)
			{
				int k = 0;
				while (k < j)
				{
					if (_GFd._Degree.reducible(G[k]._LMdeg, T_ij))
					{
						if (!veceq(_GFd._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), T_ij))
						{
							flag++;
							int temp_T_ij_deg = _GFd._Degree.calc_total_deg(T_ij);
							if (_GFd._Max_degree < temp_T_ij_deg)
							{
								_GFd._Max_degree = temp_T_ij_deg;
								_GFd._Degree.update_degree(temp_T_ij_deg);
								_D_sort.resize(_GFd._Max_degree + 1);
							}
							_D_sort[_GFd._Degree.calc_total_deg(T_ij)].push_back({ i,j });
							break;
						}
					}
					k++;
				}
			}
		}
	}
}

//ゲバウアメラーの後じゃないと動かない　書き換えれば先こっちで絞り込みも可能
template<class GF>
inline void Decision<GF>::Buchberger(vector<GF> &G)
{
	vector<vector<int>> temp;
	for (int j = 0; j < _D_sort.size(); j++)
	{
		for (int i = 0; i < _D_sort[j].size(); i++)
		{
			if (!(_GFd._Degree.gcd_1(G[_D_sort[j][i][0]]._LMdeg, G[_D_sort[j][i][1]]._LMdeg)))
			{
				temp.push_back(_D_sort[j][i]);
			}
		}
		_D_sort[j].resize(temp.size());
		_D_sort[j] = temp;
	}
}

template<class GF>
void Decision<GF>::d_sort_erase(int n)
{
	_D_sort[n].resize(0);
}

//veq イコール判定機
template<class GF>
inline bool Decision<GF>::veceq(vector<unsigned char> f, vector<unsigned char> g)
{
	for (int i = 0; i < f.size(); i++)
	{
		if (f[i] != g[i]) return false;
	}
	return true;
}