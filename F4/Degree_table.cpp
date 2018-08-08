#include "Degree_table.h"
//
//degree��defalt�͂V
Degree_table::Degree_table(int variables,int degree)
{
	_Variables = variables;
	degree_gene(0,degree);
}

//min����max�܂ł�degree�Ή���Degree_table�ɒǉ����� Degree�̍X�V�����˂� _Degree_first_index�� �ŏ���2���ȏ���O��
void Degree_table::degree_gene(int min,int max)
{
	_Degree = max;
	vector<unsigned char> deg(_Variables, 0);
	//a���̂͂��߂ƏI����index counter
	int  degree_end, degree_start, degree_end_count;

	if (min == 0)
	{//const�����@�S��0
		_Degree_first_index.push_back(0);
		_Degree_table.push_back(deg);
	
		//1��
		for (int i = _Variables - 1; i >= 0; i--)
		{
			deg[i] += 1;
			_Degree_table.push_back(deg);
			deg[i] -= 1;
		}
		_Degree_first_index.push_back(1);
		min = 2;
	}
	if(max >= 2){
		degree_end_count = _Degree_table.size() - 1;
		degree_end = _Degree_table.size() - 1;
		degree_start = 1;

		//2���ȏ�
		for (int k = min; k <= max; k++)
		{
			_Degree_first_index.push_back(_Degree_table.size());
			for (int i = _Variables - 1; i >= 0; i--)
			{
				deg[i] += 1;
				//�ЂƂO�̕ϐ�search
				for (int j = degree_start; ; j++)
				{
					if (_Degree_table[j][i] == k - 1)
					{
						degree_start = j;
						break;
					}
				}
				//����
				for (int j = degree_start; j <= degree_end; j++)
				{
					deg = vec_add(deg, _Degree_table[j]);
					_Degree_table.push_back(deg);
					deg = vec_sub(deg, _Degree_table[j]);
					degree_end_count += 1;
				}
				deg[i] -= 1;
			}
			degree_start = degree_end + 1;
			degree_end = degree_end_count;
		}
	}

}

//degree�܂�degree_table��update����
void Degree_table::update_degree(int degree)
{
	if (_Degree < degree)
	{
		degree_gene(_Degree + 1, degree);
	}
}

//index���܂�degree�܂�Degree_table��update����
void Degree_table::update_index(int index)
{
	while (_Degree_table.size() - 1 < index)
	{
		degree_gene(_Degree + 1, _Degree + 1);
	}
}

vector<unsigned char>  Degree_table::index_to_degree(int index)
{
	return _Degree_table[index];
}

int Degree_table::degree_to_index(vector<unsigned char> degree)
{
	for (int i = _Degree_first_index[calc_total_deg(degree)]; i < _Degree_table.size(); i++)
	{
		if (_Degree_table[i] == degree) return i;
	}
	return -1;
}

int Degree_table::calc_total_deg(vector<unsigned char> degree)
{
	int temp = 0;
	for (int i = 0; i < degree.size(); i++)
	{
		temp += degree[i];
	}
	return temp;
}

//�ʏ��add,_Degree_table�Ɏg��
vector<unsigned char> Degree_table::vec_add(vector<unsigned char> a, vector<unsigned char> &b)
{
	for (int i = 0; i < a.size(); i++)
	{
		a[i] = a[i] + b[i];
	}
	return a;
}

//����
vector<unsigned char> Degree_table::vec_sub(vector<unsigned char> a, vector<unsigned char> &b)
{
	for (int i = 0; i < a.size(); i++)
	{
		a[i] = a[i] - b[i];
	}
	return a;
}