#pragma once
#include<vector>

template <class GF>
class Spoly
{
public:
	//function
	vector<GF> calc_Spoly(GF f, GF g);

};

template<class GF>
inline vector<GF> Spoly<GF>::calc_Spoly(GF f, GF g)
{
	vector<GF> temp;
	vector<unsigned char> lcm_deg = GF._Degree.LCM(f._LMdeg,g._LMdeg);

	GF temp_left(f * GF._Inverse(f._LM) * GF._Degree.vec_sub(lcm_deg,f._LMdeg));
	GF temp_right(g * GF._Inverse(g._LM) * GF._Degree.vec_sub(lcm_deg,g._LMdeg));
	temp.push_back(temp_left);
	temp.push_back(temp_right);

	return temp;
}