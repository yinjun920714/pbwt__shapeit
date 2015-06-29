#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;

int main() {
	ifstream in ("phase1.chr20.10-11Mb.haps");
	char c = ' ';
	int count = 0;
	while (c!=EOF) {
		c = in.get();
		if(c != '\n')
			cout<<c<<' ';
		else
			cout<<endl;
	}
	in.close();
	return 0;
}
