#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#define B 3
#define L 10000
#define NUM_INDIVIDUAL 5008

using namespace std;

int main() {
	ifstream in ("ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.simplified.txt");
//	ofstream out ("test_output.txt");
	int **data = new int*[L];
	for (int i = 0; i < L; i++)
		data[i] = new int[NUM_INDIVIDUAL];
	int *prev_record = new int[NUM_INDIVIDUAL];      //record the consistent haplotypes labels of previous segment
	int end = 0;
	int ref[L];
	int tmp = 0;
	//first segment
	for ( int i = 0; i < L; i++) {
		for ( int j = 0 ; j < NUM_INDIVIDUAL; j++) {
			in>>data[i][j];
		}
		ref[i] = data[i][0] + data[i][1];
		if (ref[i] == 1) {
			tmp++;	
			if (tmp == B) {
			end = i + 1;
			break;
			}				
		}
	}

	int count[64];
	for ( int i = 0; i < 64; i++)
		count[i] = 0;
	int record[3];
	int record_idx;
	bool consistent_flag;
	//	compare
	for ( int j = 2; j < NUM_INDIVIDUAL; j++) {
		consistent_flag = true;
		record_idx = 0;
		for ( int i = 0; i < end; i++) {
			if (ref[i] != 1) {
				if (data[i][j] != ref[i] / 2) {
					consistent_flag = false;
					prev_record[j] = 8;      // for non-consistent haplotypes
					break;
				}
			} else {
				record[record_idx++] = data[i][j];
		    }
		}
		if (consistent_flag) {
			count[record[0] * 4 + record[1] * 2 + record[2]]++;
			prev_record[j] = record[0] * 4 + record[1] * 2 + record[2];
		}
	}
	//print the first segment
	cout<<"the first segment"<<endl;
	cout<<"000 "<<count[0]<<endl;
	cout<<"001 "<<count[1]<<endl;
	cout<<"010 "<<count[2]<<endl;
	cout<<"011 "<<count[3]<<endl;
	cout<<"100 "<<count[4]<<endl;
	cout<<"101 "<<count[5]<<endl;
	cout<<"110 "<<count[6]<<endl;
	cout<<"111 "<<count[7]<<endl;

	bool last_segment = false;
	int seg_num = 2;
	while (true) {
		//from the second segment   
		tmp = 0;
		for ( int i = 0; i < 64; i++)
			count[i] = 0;
		for ( int i = 0; i < L; i++) {
			for ( int j = 0 ; j < NUM_INDIVIDUAL; j++) {
				if( !(in>>data[i][j])) {
					last_segment = true;
					break;				
				}
			}
			if (last_segment) {
				end = i;				
				break;
			}
			ref[i] = data[i][0] + data[i][1];
			if (ref[i] == 1) {
				tmp++;	
				if (tmp == B) {
				end = i + 1;
				break;
				}				
			}
			if (i == L - 1) {
				cout<<"exceed max length of segment"<<endl;
				return -1;							
			}
		}
		
		if (last_segment) {
			break;		
		}
		//	compare
		for ( int j = 2; j < NUM_INDIVIDUAL; j++) {
			consistent_flag = true;
			record_idx = 0;
			for ( int i = 0; i < end; i++) {
				if (ref[i] != 1) {
					if (data[i][j] != ref[i] / 2) {
						consistent_flag = false;
						prev_record[j] = 8;      // for non-consistent haplotypes
						break;
					}
				} else {
					record[record_idx++] = data[i][j];
			    }
			}
			if (consistent_flag) {
				if (prev_record[j] != 8) {
					count[prev_record[j] * 8 + record[0] * 4 + record[1] * 2 + record[2]]++;
				}			
				prev_record[j] = record[0] * 4 + record[1] * 2 + record[2];

			}
		}
		
		cout<<"\n\nsegment num\t"<<setw(4)<<seg_num++<<endl;
		cout<<setw(5)<<"index";
		for (int j = 0; j < 8; j++)
			cout<<setw(10)<<j;
		cout<<endl;
		for (int i = 0; i < 8; i++) {
			cout<<setw(5)<<i;
			for (int j = 0; j < 8; j++)
				cout<<setw(10)<<count[i * 8 + j];
			cout<<endl;
		}
	}
	
	//compare last segment
	int count_idx;
	if (!tmp) {
		for ( int j = 2; j < NUM_INDIVIDUAL; j++) {
			consistent_flag = true;
			record_idx = 0;
			for ( int i = 0; i < end; i++) {
				if (ref[i] != 1) {
					if (data[i][j] != ref[i] / 2) {
						consistent_flag = false;
						prev_record[j] = 8;      // for non-consistent haplotypes
						break;
					}
				} else {
					record[record_idx++] = data[i][j];
			    }
			}
			if (consistent_flag) {
				if (prev_record[j] != 8) {
					count_idx = record[0];					
					for (int k = 1; k < tmp; k++) {
						count_idx = count_idx * 2 + record[k];
					}
					count[prev_record[j] * 8 + count_idx]++;
				}			
			}
		}

		for (int i = 0; i < tmp * 8; i++)
			cout<<"index "<<i<<"     "<<count[i]<<endl;
	}
	in.close();
//	out.close();

	delete [] prev_record;
	for (int i = 0; i < L; i++) {
		delete [] data[i];
	}
	delete [] data;
	return 0;
}
