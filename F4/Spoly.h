#pragma once
#include<vector>
#include <omp.h>

template <class GF>
class Spoly
{
public:

	//variables
	vector<GF> _Spolies;
	GF _GFs;

	//function
	void calc_Spoly(vector<GF> &G, vector<vector<int>> &D);
	void spoly_erase();

};

//SpolyåvéZÅ@
template<class GF>
inline void Spoly<GF>::calc_Spoly(vector<GF> &G, vector<vector<int>> &D)
{
	_Spolies.resize(D.size() * 2);

#pragma omp parallel for
	for (int i = 0; i < D.size(); i++)
	{
		vector<unsigned char> lcm_deg = _GFs._Degree.LCM(G[D[i][0]]._LMdeg, G[D[i][1]]._LMdeg);
		vector<unsigned char> temp = _GFs._Degree.vec_sub(lcm_deg, G[D[i][0]]._LMdeg);

		G[D[i][0]] * _GFs._Inverse[G[D[i][0]]._LM];
		G[D[i][0]] * temp;

		temp = _GFs._Degree.vec_sub(lcm_deg, G[D[i][1]]._LMdeg);
		G[D[i][1]] * _GFs._Inverse[G[D[i][1]]._LM];
		G[D[i][1]] * temp;

		_Spolies[2 * i] = G[D[i][0]];
		_Spolies[2 * i + 1] = G[D[i][1]];
	}
}

template<class GF>
inline void Spoly<GF>::spoly_erase()
{
	_Spolies.resize(0);
}