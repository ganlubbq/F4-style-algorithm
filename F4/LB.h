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

//LB計算　Sp_redをそのまま変換する
template <class GF>
void LB<GF>::calc_LB(vector<GF> &Sp_red)
{
	Gauss_rev(Sp_red);
}

template <class GF>
void LB<GF>::Gauss_rev(vector<GF> &Sp_red)
{
	for (int i = 0; i < Sp_red.size(); i++)
	{
		//0多項式
		if (Sp_red[i]._LMdeg_index == -1) continue;
		//LC = 1に
		if (Sp_red[i]._LM != 1)
		{
			Sp_red[i] * _GFl._Inverse[Sp_red[i]._LM];
		}
		int index = Sp_red[i]._LMdeg_index;
#pragma omp parallel for
		for (int j = 0; j < Sp_red.size(); j++)
		{
			if (i == j) continue;
			if (Sp_red[j]._Coeff[index] != 0)
			{
				GF temp = Sp_red[i];
				temp * (_GFl._Add_inverse[Sp_red[j]._Coeff[index]]);
				Sp_red[j] + temp;
			}
		}
	}
}


