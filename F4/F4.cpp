#include "F4.h"

F4::F4(string filename)
{
	file_name = filename;
}

void F4::file_read()
{
	FILE *fp;
	fp = fopen(file_name.c_str(), "r");

	unsigned int coeff;

	fscanf(fp, "%d", &coeff);

	if (coeff != _Variables) {
		puts("VAL_NUM is not correct");
		exit(0);
	}

	cout << "init" << endl;
	unsigned char coeff_char;
	int coeff_int;

	//MQcharenngeŒ^‚µ‚©“Ç‚Ýž‚ß‚È‚¢•ÏX
	for (int i = 0; i < 2 * _Variables + 1; i++)
	{
		vector<unsigned char> temp;
		for (int j = 0; j < var_deg_comb(_Variables, 2); j++) {
			fscanf(fp, "%x", &coeff_int);
			temp.push_back((unsigned char)coeff_int);
		}
		reverse(temp.begin(), temp.end());
		_Equations.push_back(GF31(temp));
	}
}

//n•Ï”rŽŸ‘½€Ž®‚Ü‚Å‚Ì‘S€”
inline int F4:: var_deg_comb(int n, int r) {//n•Ï”rŽŸ‘½€Ž®‚Ü‚Å‚Ì‘S€”
	int ans = 1;
	n += r;
	for (int div = 1; div <= r; ++div, --n) {
		ans = ans * n / div;
	}
	return ans;
}