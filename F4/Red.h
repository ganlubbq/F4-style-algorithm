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
	vector<GF> _Reds;

	//function
	void calc_red(vector<GF> &Spolies,vector<GF> &G);

};

template<class GF>
void Red<GF>::calc_red(vector<GF> &Spolies, vector<GF> &G)
{
	//spoly正規化してればここはしょれる?
	vector<int> M_index;
	for (int i = 0; i < Spolies.size(); i++)
	{
		for (int j = 0; j < Spolies[i]._Coeff.size(); j++)
		{
			if (Spolies[i]._Coeff[j] != 0)
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
	for (int i = 0; i < M_index.size(); i++)
	{
		cout << M_index[i] << "M" << endl;
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
				if (_GFr._Degree.reducible(G[j]._LMdeg, _GFr._Degree.index_to_degree(M_index[i])))
				{
					flag = true;
					vector<unsigned char> temp_deg = _GFr._Degree.vec_sub(_GFr._Degree.index_to_degree(M_index[i]),G[j]._LMdeg);
					GF temp_g = G[j];
					temp_g * temp_deg;
					_Reds.push_back(temp_g);
					temp_g = G[j];
					temp_g.LM_del();
					temp_g * temp_deg;
					//Mに追加
					for (int k = 0; k < temp_g._Coeff.size(); k++)
					{
						if (temp_g._Coeff[k] != 0)
						{
							M_index.push_back(k);
						}
					}
					break;
				}
			}
		}
		for (int i = 0; i < M_index.size(); i++)
		{
			cout << M_index[i] << "MM" << endl;
		}
		//使ったやつ削除
		M_index.erase(M_index.begin(), M_index.begin() + size_before);
		for (int i = 0; i < M_index.size(); i++)
		{
			cout << M_index[i] << "MMM" << endl;
		}
		std::sort(M_index.begin(), M_index.end());
		//重複削除
		M_index.erase(std::unique(M_index.begin(), M_index.end()), M_index.end());
		for (int i = 0; i < M_index.size(); i++)
		{
			cout << M_index[i] << "MMMM" << endl;
		}
	}
}