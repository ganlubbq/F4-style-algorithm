#pragma once
#include<vector>
#include "Degree_table.h"

template<class GF>
class Decision
{
public:
	Decision() {};
	//variables
	//vector<vector<int>> _D;
	vector<vector<vector<int>>> _D_sort;//total degree順
	GF _GFd;
	static Degree_table _Degree_deci;

	//function
	void decision(vector<GF> &G);
	void decision_0_kai(vector<int> &LMplace, vector<GF> &G);
	void decision_0_kai_2(vector<int> &LMplace, vector<int> &LMplace_old, vector<GF> &G);
	void decision_kai(vector<int> &LMplace);
	void decision_kai_2(vector<int> &LMplace, vector<int> &LMplace_old);
	void Gebauer_Moller(vector<GF> &G);
	void Gebauer_Moller_0_kai(vector<int> &LMplace, vector<GF> &G);
	void Gebauer_Moller_0_kai_2(vector<int> &LMplace, vector<int> &LMplace_old,vector<GF> &G);
	void Gebauer_Moller_kai(vector<int> &LMplace);
	void Gebauer_Moller_kai_2(vector<int> &LMplace, vector<int> &LMplace_old);
	void Gebauer_Moller_mono(vector<GF> &G);
	void Gebauer_Moller_num(vector<GF> &G,int num);
	void Buchberger(vector<GF> &G);
	void Buchberger_kai();
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

template <class GF>
inline void Decision<GF>::decision_0_kai(vector<int> &LMplace, vector<GF> &G)
{
	_D_sort.resize(_GFd._Max_degree + 1);
	Gebauer_Moller_0_kai(LMplace,G);
	Buchberger(G);
}

template <class GF>
inline void Decision<GF>::decision_0_kai_2(vector<int> &LMplace, vector<int> &LMplace_old,vector<GF> &G)
{
	_D_sort.resize(_GFd._Max_degree + 1);
	Gebauer_Moller_0_kai_2(LMplace,LMplace_old,G);
	Buchberger(G);
}

template <class GF>
inline void Decision<GF>::decision_kai(vector<int> &LMplace)
{
	_D_sort.resize(_GFd._Max_degree + 1);
	Gebauer_Moller_kai(LMplace);
	Buchberger_kai();
}

template <class GF>
inline void Decision<GF>::decision_kai_2(vector<int> &LMplace, vector<int> &LMplace_old)
{
	_D_sort.resize(_GFd._Max_degree + 1);
	Gebauer_Moller_kai_2(LMplace,LMplace_old);
	Buchberger_kai();
}

//LMplaceはsort済み
template<class GF>
inline void Decision<GF>::Gebauer_Moller_kai(vector<int> &LMplace)
{
	//{i,j}消してもよいか判断 消せないときはd_sortへ
	for (int i = 0; i < LMplace.size(); i++)
	{
		vector<unsigned char> i_deg = _Degree_deci.index_to_degree(LMplace[i]);
		for (int j = i + 1; j < LMplace.size(); j++)
		{
			int temp_T_ij_deg = 0;
			bool flag = false;
			vector<unsigned char> temp_j = _Degree_deci.index_to_degree(LMplace[j]);
			vector<unsigned char> T_ij = _Degree_deci.LCM(i_deg,temp_j);

			for (int k = 0; k < i; k++)
			{
				if (_Degree_deci.reducible(_Degree_deci.index_to_degree(LMplace[k]), T_ij))
				{
					flag = true;
					break;
				}
			}
			if (flag == false)
			{
				temp_T_ij_deg = _Degree_deci.calc_total_deg(T_ij);
				if (_GFd._Max_degree < temp_T_ij_deg)
				{
					_GFd._Max_degree = temp_T_ij_deg;
					_GFd._Degree.update_degree(temp_T_ij_deg);
					_D_sort.resize(_GFd._Max_degree + 1);
				}
				_D_sort[temp_T_ij_deg].push_back({ LMplace[i],LMplace[j] });
			}
		}
	}
}

//LMplaceはsort済み
template<class GF>
inline void Decision<GF>::Gebauer_Moller_kai_2(vector<int> &LMplace, vector<int> &LMplace_old)
{
	//{i,j}消してもよいか判断 消せないときはd_sortへ まずnew同士のi,j
	for (int i = 0; i < LMplace.size(); i++)
	{
		vector<unsigned char> i_deg = _Degree_deci.index_to_degree(LMplace[i]);
		for (int j = i + 1; j < LMplace.size(); j++)
		{
			int temp_T_ij_deg = 0;
			bool flag = false;
			vector<unsigned char> temp_j = _Degree_deci.index_to_degree(LMplace[j]);
			vector<unsigned char> T_ij = _Degree_deci.LCM(i_deg, temp_j);

			//new　めらー判定
			for (int k = 0; k < i; k++)
			{
				if (_Degree_deci.reducible(_Degree_deci.index_to_degree(LMplace[k]), T_ij))
				{
					flag = true;
					break;
				}
			}
			//old　めらー判定
			for (int k = 0; k < LMplace_old.size(); k++)
			{
				if (LMplace_old[k] > LMplace[i]) break;
				if (_Degree_deci.reducible(_Degree_deci.index_to_degree(LMplace_old[k]), T_ij))
				{
					flag = true;
					break;
				}
			}
			if (flag == false)
			{
				temp_T_ij_deg = _Degree_deci.calc_total_deg(T_ij);
				if (_GFd._Max_degree < temp_T_ij_deg)
				{
					_GFd._Max_degree = temp_T_ij_deg;
					_GFd._Degree.update_degree(temp_T_ij_deg);
					_D_sort.resize(_GFd._Max_degree + 1);
				}
				_D_sort[temp_T_ij_deg].push_back({ LMplace[i],LMplace[j] });
			}
		}
	}

	//oldとnew
	for (int i = 0; i < LMplace.size(); i++)
	{
		vector<unsigned char> i_deg = _Degree_deci.index_to_degree(LMplace[i]);
		for (int j = 0; j < LMplace_old.size(); j++)
		{
			if (LMplace[i] > LMplace_old[j]) break;
			int temp_T_ij_deg = 0;
			bool flag = false;
			vector<unsigned char> temp_j = _Degree_deci.index_to_degree(LMplace_old[j]);
			vector<unsigned char> T_ij = _Degree_deci.LCM(i_deg, temp_j);

			//new　めらー判定
			for (int k = 0; k < i; k++)
			{
				if (_Degree_deci.reducible(_Degree_deci.index_to_degree(LMplace[k]), T_ij))
				{
					flag = true;
					break;
				}
			}
			//old　めらー判定
			for (int k = 0; k < LMplace_old.size(); k++)
			{
				if (LMplace_old[k] > LMplace[i]) break;
				if (_Degree_deci.reducible(_Degree_deci.index_to_degree(LMplace_old[k]), T_ij))
				{
					flag = true;
					break;
				}
			}
			if (flag == false)
			{
				temp_T_ij_deg = _Degree_deci.calc_total_deg(T_ij);
				if (_GFd._Max_degree < temp_T_ij_deg)
				{
					_GFd._Max_degree = temp_T_ij_deg;
					_GFd._Degree.update_degree(temp_T_ij_deg);
					_D_sort.resize(_GFd._Max_degree + 1);
				}
				_D_sort[temp_T_ij_deg].push_back({ LMplace[i],LMplace[j] });
			}
		}
	}
}

template<class GF>
inline void Decision<GF>::Gebauer_Moller_0_kai(vector<int> &LMplace, vector<GF> &G)
{
	for (int i = 0; i < LMplace.size(); i++)
	{
		for (int j = i + 1; j < LMplace.size(); j++)
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

template<class GF>
inline void Decision<GF>::Gebauer_Moller_0_kai_2(vector<int> &LMplace,vector<int> &LMplace_old ,vector<GF> &G)
{
	//判定されるやつnewの中の一人
	for (int i = LMplace[0]; i < LMplace.size() - 1; i++)
	{
		//判定するやつold
		for (int j = LMplace_old[0]; j < LMplace_old.size(); j++)
		{
			//F
			int flag = 0;
			int k = LMplace_old[0];
			vector<unsigned char> T_ij = _GFd._Degree.LCM(G[LMplace[i]]._LMdeg, G[LMplace_old[j]]._LMdeg);
			while (LMplace_old[k] < LMplace[i])
			{
				if (veceq(_GFd._Degree.LCM(G[LMplace_old[j]]._LMdeg, G[LMplace_old[k]]._LMdeg), T_ij))
				{
					flag++;
					int temp_T_ij_deg = _GFd._Degree.calc_total_deg(T_ij);
					if (_GFd._Max_degree < temp_T_ij_deg)
					{
						_GFd._Max_degree = temp_T_ij_deg;
						_GFd._Degree.update_degree(temp_T_ij_deg);
						_D_sort.resize(_GFd._Max_degree + 1);
					}
					_D_sort[_GFd._Degree.calc_total_deg(T_ij)].push_back({ LMplace[i],LMplace_old[j] });
					break;
				}
				k++;
			}

			//M
			if (flag == 0)
			{
				int k = LMplace_old[0];
				while (LMplace_old[k] < LMplace[j])
				{
					if (_GFd._Degree.reducible(G[LMplace_old[k]]._LMdeg, T_ij))
					{
						if (!veceq(_GFd._Degree.LCM(G[LMplace_old[j]]._LMdeg, G[LMplace_old[k]]._LMdeg), T_ij))
						{
							flag++;
							int temp_T_ij_deg = _GFd._Degree.calc_total_deg(T_ij);
							if (_GFd._Max_degree < temp_T_ij_deg)
							{
								_GFd._Max_degree = temp_T_ij_deg;
								_GFd._Degree.update_degree(temp_T_ij_deg);
								_D_sort.resize(_GFd._Max_degree + 1);
							}
							_D_sort[_GFd._Degree.calc_total_deg(T_ij)].push_back({ LMplace[i],LMplace_old[j] });
							break;
						}
					}
					k++;
				}
			}
		}
	}
	for (int i = 0; i < LMplace.size(); i++)
	{
		for (int j = i + 1; j < LMplace.size(); j++)
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

//numの場所に対して他全員でチェック　 Gebauer_Moller_monoの上位互換
template<class GF>
inline void Decision<GF>::Gebauer_Moller_num(vector<GF> &G, int num)
{
	cout << "Gebauer_Moller" << endl;
	cout << G.size() << endl;
	cout << num << endl;

	//i < num = j
	for (int i = 0; i < num; i++)
	{
		if (G[i]._LMdeg_index != -1)
		{
			int j = num;
			//F
			int flag = 0;
			int k = 0;
			cout << "eeee" << endl;
			if (G[j]._LMdeg_index != -1)
			{
				vector<unsigned char> T_ij = _GFd._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg);
				cout << "mazika" << endl;
				while (k < i)
				{
					cout << "don!!" << endl;
					if (G[k]._LMdeg_index != -1)
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
					}
					k++;
				}


				//M
				if (flag == 0)
				{
					int k = 0;
					cout << "seyana" << endl;
					while (k < j)
					{
						if (G[k]._LMdeg_index != -1)
						{
							cout << "one two three" << endl;
							if (_GFd._Degree.reducible(G[k]._LMdeg, T_ij))
							{
								cout << "one" << endl;
								if (!veceq(_GFd._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), T_ij))
								{
									cout << "two" << endl;
									flag++;
									int temp_T_ij_deg = _GFd._Degree.calc_total_deg(T_ij);
									cout << "three" << endl;
									if (_GFd._Max_degree < temp_T_ij_deg)
									{
										cout << "four" << endl;
										_GFd._Max_degree = temp_T_ij_deg;
										_GFd._Degree.update_degree(temp_T_ij_deg);
										_D_sort.resize(_GFd._Max_degree + 1);
									}
									cout << "$$" << endl;
									_D_sort[_GFd._Degree.calc_total_deg(T_ij)].push_back({ i,j });
									break;
								}
								cout << "sa-senn" << endl;
							}
						}
						k++;
					}
					cout << "oissu" << endl;
				}

				//B
				if (flag == 0)
				{
					int k = G.size() - 1;
					while (k > j)
					{
						if (G[k]._LMdeg_index != -1)
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
						}
						k--;
					}
				}
			}
		}
	}

	cout << "complete" << endl;

	//i = num < j
	int i = num;
	if (G[i]._LMdeg_index != -1)
	{
		cout << G.size() << endl;
		cout << num << endl;
		for (int j = num + 1; j < G.size(); j++)
		{
			if (G[j]._LMdeg_index != -1)
			{
				cout << "not in" << endl;
				//F
				int flag = 0;
				int k = 0;
				vector<unsigned char> T_ij = _GFd._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg);
				while (k < i)
				{
					if (G[k]._LMdeg_index != -1)
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
					}
					k++;
				}

				//M
				if (flag == 0)
				{
					int k = 0;
					while (k < j)
					{
						if (G[k]._LMdeg_index != -1)
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
						if (G[k]._LMdeg_index != -1)
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
						}
						k--;
					}
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
		if (temp.size() > 0)
		{
			_D_sort[j].resize(temp.size());
			//ここバグっぽい　3次から3次出たとき　かつ　パラレルディブ以上の時　前のやつ上書きされ消える
			_D_sort[j] = temp;
		}
		temp.resize(0);
	}
}

template<class GF>
inline void Decision<GF>::Buchberger_kai()
{
	vector<vector<int>> temp;
	//次数
	for (int j = 0; j < _D_sort.size(); j++)
	{
		for (int i = 0; i < _D_sort[j].size(); i++)
		{
			vector <unsigned char> temp_0 = _Degree_deci.index_to_degree(_D_sort[j][i][0]);
			vector <unsigned char> temp_1 = _Degree_deci.index_to_degree(_D_sort[j][i][1]);
			if (!(_Degree_deci.gcd_1(temp_0, temp_1)))
			{
				temp.push_back(_D_sort[j][i]);
			}
		}
		if (temp.size() > 0)
		{
			_D_sort[j].resize(temp.size());
			_D_sort[j] = temp;
		}
		temp.resize(0);
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
	if (f.size() != g.size()) return false;
	for (int i = 0; i < f.size(); i++)
	{
		if (f[i] != g[i]) return false;
	}
	return true;
}