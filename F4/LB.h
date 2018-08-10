#pragma once

#include<vector>
#include<omp.h>

template<class GF>
class LB
{
public:
	//variables
	GF _GFl;
	//vector<GF> _LB;

	//function
	void calc_LB(vector<GF> &Sp_red);
	void Gauss_rev(vector<GF> &Sp_red);
};

//LBŒvZ@Sp_red‚ğ‚»‚Ì‚Ü‚Ü•ÏŠ·‚·‚é
template <class GF>
void LB<GF>::calc_LB(vector<GF> &Sp_red)
{
	Gauss_rev(vector<GF> &Sp_red);
	del_0poly(vector<GF> &Sp_red);
	
}

template <class GF>
void LB<GF>::Gauss_rev(vector<GF> &Sp_red)
{
	for (int i = 0; i < Sp_red.size(); i++)
	{
		//0‘½€®
		if (Sp_red[i]._LMdeg == -1) continue;
		//LC = 1‚É
		if (Sp_red[i]._LM != 1)
		{
			Sp_red[i] * _GFl._Inverse(Sp_red[i]._LM);
		}
		int index = Sp_red[i]._LMdeg_index;
#pragma omp parallel for
		for (int j = 0; j < Sp_red.size(); j++)
		{
			if (i == j) continue;
			if (Sp_red[j]._Coeff[index] != 0)
			{
				Sp_red[j] + (Sp_red[i] * (_GFl._Add_inverse(Sp_red[j]._Coeff[index])));
			}
		}
	}
}


