#pragma once
#include <vector>
#include <iostream>

using namespace std;


/*indexとdegreeの対応を保存 
Variables:変数の個数
Degree:この次数までの対応を持つ
Degree_first_index:次数の最初のindexをもつ 0->0,1->1,2->Variables + 1,...
*/
class Degree_table
{
public:
	//variable
	vector<vector<unsigned char>> _Degree_table;
	int _Degree;
	int _Variables;
	vector<int> _Degree_first_index;

	//function
	Degree_table(int variables,int degree = 7);
	void degree_gene(int min,int max);
	void update_degree(int degree);
	void update_index(int index);
	vector<unsigned char> index_to_degree(int index);
	int degree_to_index(vector<unsigned char> degree);
	int calc_total_deg(vector<unsigned char> degree);

	virtual vector<unsigned char> LCM(vector<unsigned char> &f, vector<unsigned char> &g);
	//&使うとバグる　正確には計算結果をそのまま引数にとれない
	virtual bool reducible(vector<unsigned char> f, vector<unsigned char> g);
	virtual bool gcd_1(vector<unsigned char> &f, vector<unsigned char> &g);

	vector<unsigned char> vec_add(vector<unsigned char> a, vector<unsigned char> &b);
	vector<unsigned char> vec_sub(vector<unsigned char> a, vector<unsigned char> &b);

};