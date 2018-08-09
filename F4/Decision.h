#pragma once
#include<vector>

template<class GF>
class Decision
{
public:
	//variables
	vector<vector<int>> _D;
	vector<vector<vector<int>>> _D_sort;//total degree順

	//function
	vector<int> decision(vector<GF> &G);
	void Gebauer_Moller(vector<GF> &G);
	void Buchberger(vector<GF> &G);
	void sort_D(vector<GF> &G);

	virtual bool veceq(vector<unsigned char> &f, vector<unsigned char> &g);

};

template <class GF>
vector<int> Decision<GF>::decision(vector<GF> &G)
{
	Gebauer_Moller(G);
	Buchberger(G);
	sort_D(G);
}

template<class GF>
void Decision<GF>::Gebauer_Moller(vector<GF> &G)
{
	for (int i = 0; i < G.size(); i++)
	{
		for (int j = i + 1; j < G.size(); j++)
		{
			//F
			int flag = 0;
			int k = 0;
			while (k < i)
			{
				if (veceq(GF._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), GF._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg)))
				{
					flag++;
					_D.push_back({ i,j });
					break;
				}
				k ++ ;
			}

			//M
			if (flag == 0)
			{
				int k = 0;
				while (k < j)
				{
					if (GF._Degree.reducible(GF._LMdeg, GF._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg)))
					{
						if (!veceq(GF._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), GF._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg)))
						{
							flag++;
							_D.push_back({ i,j });
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
					if (GF._Degree.reducible(GF._LMdeg, GF._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg)))
					{
						if (!veceq(GF._Degree.LCM(G[i]._LMdeg, G[k]._LMdeg), GF._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg)))
						{
							if (!veceq(GF._Degree.LCM(G[j]._LMdeg, G[k]._LMdeg), GF._Degree.LCM(G[i]._LMdeg, G[j]._LMdeg)))
							{
								_D.push_back({ i,j });
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

//ゲバウアメラーの後じゃないと動かない　書き換えれば先こっちで絞り込みも可能
template<class GF>
void Decision<GF>::Buchberger(vector<GF> &G)
{
	vector<int> temp;
	for (iny i = 0;i < _D.size();i++)
	{
		if (GF._Degree.gcd_1(G[_D[i][0]]._LMdeg,G[_D[i][1]]._LMdeg))
		{
			temp.push_back(_D[i]);
		}
	}
	_D.resize(temp.size());
	_D = temp;
}

//Spolyは左と右分けることが前提　sortする Max_degreeの更新も
template<class GF>
void Decision<GF>::sort_D(vector<GF> &G)
{
	_D_sort.resize(GF._Max_degree);

	for (int i = 0; i < _D.size(); i++)
	{
		int lcm_deg = GF._Degree.calc_total_deg(GF._Degree.LCM(G[_D[i][0]].LMdeg, G[_D[i][1]].LMdeg));
		if (GF._Max_degree < lcm_deg)
		{
			GF._Max_degree = lcm_deg;
			GF._Degree.update_degree(_Max_degree);
			_D_sort.resize(GF._Max_degree);
		}
		_D_sort[GF._Max_degree].push_back(_D[i]);
	}
}


//veq イコール判定機
template<class GF>
bool Decision<GF>::veceq(vector<unsigned char> &f,vector<unsigned char> &g)
{
	for (inti = 0; i < f.size(); i++)
	{
		if (f[i] != g[i]) return false;
	}
	return true;
}