#include <iostream>
#include <fstream>
#include <cstring>
#define B 3
#define L 200
#define NUM_INDIVIDUAL 5008

using namespace std;

int main() {
	ifstream in ("ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.small.txt");
//	ofstream out ("test_output.txt");
	int data[L][NUM_INDIVIDUAL];
	for ( int i = 0; i < L; i++) {
			for ( int j = 0; j < NUM_INDIVIDUAL; j++) {
				in>>data[i][j];
			}
		}


	int begin = 0, end = 0;
	int ref[L];
	for ( int i = 0; i < L; ++i )
		ref[i] = data[i][0] + data[i][1];
	
	//first segment
	int tmp = 0;
	while (end < L) {
		if (ref[end++] == 1) {
			tmp++;
			if (tmp == B)
				break;
		}
	}

	int count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	int record[3];
	int record_idx;
	bool consistent_flag;
	//	compare
	for ( int j = 2; j < NUM_INDIVIDUAL; j++) {
		consistent_flag = true;
		record_idx = 0;
		for ( int i = begin; i < end; i++) {
			if (ref[i] != 1) {
				if (data[i][j] != ref[i] / 2) {
					consistent_flag = false;
					break;
				}
			} else {
				record[record_idx++] = data[i][j];
		    }
		}
		if (consistent_flag) {
			count[record[0]*4 + record[1]*2 + record[2]]++;
		}
	}
	
	cout<<count[0]<<endl;
	cout<<count[1]<<endl;
	cout<<count[2]<<endl;
	cout<<count[3]<<endl;
	cout<<count[4]<<endl;
	cout<<count[5]<<endl;
	cout<<count[6]<<endl;
	cout<<count[7]<<endl;
	in.close();
//	out.close();


	return 0;
}
