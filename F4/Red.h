#pragma once

#include<vector>
#include <set>
#include <algorithm>

template<class GF>
class Red
{
public:
	//variables
	GF _GFr;
	vector<GF> _Red;

	//function
	void calc_red(vector<GF> &Spolies,vector<GF> &G);

};

template<class GF>
void Red<GF>::calc_red(vector<GF> &Spolies, vector<GF> &G)
{
	//spoly正規化してればここはしょれる
	vector<int> M_index;
	for (int i = 0; i < Spolies.size(); i++)
	{
		for (int j = 0; j < Spolies[i].size(); j++)
		{
			if (Spolies[i][j] != 0)
			{
				bool flag_wao = true;
				for (int k = 0; k < M_index.size(); k++)
				{
					if (j == M_index[k])
					{
						flag_wao = false;
						break;
					}
				}
				if(flag_wao) M_index.push_back(j);
			}
		}
	}
	/*//変換　後でやると時間食いそう
	vector<vector<unsigned char>> M_deg;
	M_deg.resize(M_index);
	for (int i = 0; i < M_deg.size(); i++)
	{
		M_deg[i] = _Gfr._Degree.index_to_degree(M_index[i]);
	}*/

	//calc_red
	bool flag = true;
	while (flag)
	{
		flag = false;
		int size_before = M_index.size();
#pragma omp parallel for reduction(||:flag)
		for(int i = 0;i < M_index.size();i++)
		{
			for (int j = 0; j < G.size(); j++)
			{
				if (_Gfr._Degree.reducible(G[j]._LMdeg, _Gfr._Degree.index_to_degree(M_index[i])))
				{
					flag = true;
					vector<unsigned char> temp_deg = _Gfr._degree.vec_sub(_Gfr._Degree.index_to_degree(M_index[i]),G[j]._LMdeg);
					GF temp_g = G[j];
					temp_g * temp_deg;
					_Red.push_back(temp_g);
					temp_g = G[j];
					temp_g.LM_del();
					temp_g * temp_deg;
					for (int k = 0; k < temp_g.size(); k++)
					{
						if (temp_g._Coeff[k] != 0)
						{
							M_index.push_back(k);
						}
					}
				}
			}
		}
		//使ったやつ削除
		M_index.erase(M_index.begin(), M_index.begin() + size_before);
		std::sort(M_index.begin(), M_index.end());
		//重複削除
		M_index.erase(std::unique(M_index.begin(), M_index.end()), M_index.end());
	}
}