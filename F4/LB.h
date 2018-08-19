#pragma once

#include<vector>
#include<omp.h>
#include <time.h>
#include <algorithm>

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
	vector<int> Gauss_rev_Eq(vector<GF> &Eq);
};

//LBåvéZÅ@Sp_redÇÇªÇÃÇ‹Ç‹ïœä∑Ç∑ÇÈ
template <class GF>
void LB<GF>::calc_LB(vector<GF> &Sp_red)
{
	Gauss_rev(Sp_red);
}

template <class GF>
void LB<GF>::Gauss_rev(vector<GF> &Sp_red)
{
	//vector sizeí≤êÆ
	int max_length = 0;
	for (int i = 0; i < Sp_red.size(); i++)
	{
		if (max_length < Sp_red[i]._LMdeg_index + 1) max_length = Sp_red[i]._LMdeg_index + 1;
	}
	for (int i = 0; i < Sp_red.size(); i++)
	{
		Sp_red[i]._Coeff.resize(max_length);
		Sp_red[i]._Coeff_size = max_length;
		Sp_red[i]._Div_single_size = Sp_red[i]._Coeff_size / single_size;
	}

	//cout << max_length << endl;

	for (int i = 0; i < Sp_red.size(); i++)
	{
		//0ëΩçÄéÆ
		if (Sp_red[i]._LMdeg_index == -1) continue;
		//LC = 1Ç…
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

template <class GF>
vector<int> LB<GF>::Gauss_rev_Eq(vector<GF> &Eq)
{
	vector<int> result;

	//vector sizeí≤êÆ
	int max_length = 0;
	for (int i = 0; i < Eq.size(); i++)
	{
		if (max_length < Eq[i]._LMdeg_index + 1) max_length = Eq[i]._LMdeg_index + 1;
	}
	for (int i = 0; i < Eq.size(); i++)
	{
		Eq[i]._Coeff.resize(max_length);
		Eq[i]._Coeff_size = max_length;
		Eq[i]._Div_single_size = Eq[i]._Coeff_size / single_size;
	}

	//cout << max_length << endl;

	for (int i = 0; i < Eq.size(); i++)
	{
		//0ëΩçÄéÆ
		if (Eq[i]._LMdeg_index == -1) continue;
		//LC = 1Ç…
		if (Eq[i]._LM != 1)
		{
			Eq[i] * _GFl._Inverse[Eq[i]._LM];
		}
		int index = Eq[i]._LMdeg_index;
#pragma omp parallel for
		for (int j = 0; j < Eq.size(); j++)
		{
			if (i == j) continue;
			if (Eq[j]._Coeff[index] != 0)
			{
				GF temp = Eq[i];
				temp * (_GFl._Add_inverse[Eq[j]._Coeff[index]]);
				Eq[j] + temp;
				result.push_back(j);
			}
		}
	}
	std::sort(result.begin(), result.end());
	//èdï°çÌèú
	result.erase(std::unique(result.begin(), result.end()), result.end());
	return result;
}
