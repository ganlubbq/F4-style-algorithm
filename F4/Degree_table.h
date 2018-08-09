#pragma once
#include <vector>
#include <iostream>

using namespace std;


/*index��degree�̑Ή���ۑ� 
Variables:�ϐ��̌�
Degree:���̎����܂ł̑Ή�������
Degree_first_index:�����̍ŏ���index������ 0->0,1->1,2->Variables + 1,...
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
	//&�g���ƃo�O��@���m�ɂ͌v�Z���ʂ����̂܂܈����ɂƂ�Ȃ�
	virtual bool reducible(vector<unsigned char> f, vector<unsigned char> g);
	virtual bool gcd_1(vector<unsigned char> &f, vector<unsigned char> &g);

	vector<unsigned char> vec_add(vector<unsigned char> a, vector<unsigned char> &b);
	vector<unsigned char> vec_sub(vector<unsigned char> a, vector<unsigned char> &b);

};