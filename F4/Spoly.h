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
	void calc_Spoly(vector<GF> &G, vector<int> &D);
	void spoly_erase();

};

//SpolyåvéZÅ@
template<class GF>
inline void Spoly<GF>::calc_Spoly(vector<GF> &G, vector<int> &D)
{
	_Spolies.resize(D.size() * 2);

#pragma omp parallel for
	for (int i = 0; i < D.size(); i++)
	{
		vector<unsigned char> lcm_deg = _GFs._Degree.LCM(G[D[0]]._LMdeg, G[D[1]]._LMdeg);

		GF temp_left(G[D[0]] * _GFs._Inverse(G[D[0]]._LM) * _GFs._Degree.vec_sub(lcm_deg, G[D[0]]._LMdeg));
		GF temp_right(G[D[1]] * _GFs._Inverse(G[D[1]]._LM) * _GFs._Degree.vec_sub(lcm_deg, G[D[1]]._LMdeg));
		_Spolies[2 * i] = temp_left;
		_Spolies[2 * i + 1] = temp_right;
	}
}

template<class GF>
inline void Spoly<GF>::spoly_erase()
{
	_Spolies.resize(0);
}