#include <iostream>
#include <fstream>
#include <cstring>
#define B 3
#define L 200
#define NUM_INDIVIDUAL 10015

using namespace std;

int main() {
	ifstream in ("ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.simplified.txt");
	char c = ' ';
	int count = 0;
	while (c!='\n') {
		c = in.get();
		if (c != ' ')
			count++;
	}
	cout<<count;
	in.close();
	return 0;
}
