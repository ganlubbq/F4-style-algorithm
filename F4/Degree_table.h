#pragma once
#include <vector>

using namespace std;


/*index‚Ædegree‚Ì‘Î‰ž‚ð•Û‘¶ 
Variables:•Ï”‚ÌŒÂ”
Degree:‚±‚ÌŽŸ”‚Ü‚Å‚Ì‘Î‰ž‚ðŽ‚Â
Degree_first_index:ŽŸ”‚ÌÅ‰‚Ìindex‚ð‚à‚Â 0->0,1->1,2->Variables + 1,...
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

	vector<unsigned char> vec_add(vector<unsigned char> a, vector<unsigned char> &b);
	vector<unsigned char> vec_sub(vector<unsigned char> a, vector<unsigned char> &b);

};