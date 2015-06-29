// command   ./a.out input_file M N ref_num 


#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#define B 3
//#define N 10000
//#define M 2184

using namespace std;

int main(int argc, char *argv[]) {
	ifstream in (argv[1]);
	int M = atoi(argv[2]);
	int N = atoi(argv[3]);
	int ref_num = atoi(argv[4]);
//	ofstream out ("test_output.txt");
	int **data = new int*[N];
	for (int i = 0; i < N; i++)
		data[i] = new int[M];
	int *prev_record = new int[M];      //record the consistent haplotypes labels of previous segment
	int end;
	int tmp = 0;
    //calculate the ref seq. and 1 position
    int *x;
    int *pos; /* for the 1 position for ref */
    int seg_num = 1; /* for the segment number  */
    int num_1 = 0; /* for the num of 1 */
    int k;  /* for the k th  site */
	//input data
	for ( int i = 0; i < N; i++ )
		for ( int j = 0; j < M; j++ )
			in>>data[i][j];

	struct timeval tstart, tend;
    gettimeofday( &tstart, NULL );


    x = new int[N];
    pos = new int[N];

    for ( int i = 0; i < N; ++i)
    	x[i] = data[i][ref_num * 2] + data[i][ref_num * 2 + 1];

    for (int i = 0, j = 0; i < N; ++i) {
    	if (x[i] == 1) {
    		++j;
    		pos[num_1++] = i;
    		if (j == 3) {
    			++seg_num;
    			j = 0;
    		}
    	}
    }


	//first segment
	int count[64];
	for ( int i = 0; i < 64; i++)
		count[i] = 0;
	int record[3];
	int record_idx;
	bool consistent_flag;
	//	compare
	for ( int j = 0; j < M; j++) {
		consistent_flag = true;
		record_idx = 0;
		for ( k = 0; k <= pos[2]; k++) {
			if (x[k] != 1) {
				if (data[k][j] != x[k] / 2) {
					consistent_flag = false;
					prev_record[j] = 8;      // for non-consistent haplotypes
					break;
				}
			} else {
				record[record_idx++] = data[k][j];
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

	int s = 2;
	while (true) {
		//from the second segment   
		for ( int i = 0; i < 64; i++)
			count[i] = 0;
		if ( s > seg_num - 1 )   //last segment 
			break;    
		end = pos[s * 3 - 1];
		for ( int j = 0; j < M; j++) {
			k = pos[s * 3 -4] + 1;
			consistent_flag = true;
			record_idx = 0;
			for ( ; k <= end; k++) {
				if (x[k] != 1) {
					if (data[k][j] != x[k] / 2) {
						consistent_flag = false;
						prev_record[j] = 8;      // for non-consistent haplotypes
						break;
					}
				} else {
					record[record_idx++] = data[k][j];
			    }
			}
			if (consistent_flag) {
				if (prev_record[j] != 8) {
					count[prev_record[j] * 8 + record[0] * 4 + record[1] * 2 + record[2]]++;
				}			
				prev_record[j] = record[0] * 4 + record[1] * 2 + record[2];

			}
		}
		
		cout<<"\n\nsegment num\t"<<setw(4)<<s++<<endl;
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
	tmp = num_1 - (s - 1) * 3;
	if (tmp > 0) {
		for ( int j = 0; j < M; j++) {
			consistent_flag = true;
			record_idx = 0;
			k = pos[(seg_num - 1)* 3 - 1] + 1;
			for ( ; k < N; k++) {
				if (x[k] != 1) {
					if (data[k][j] != x[k] / 2) {
						consistent_flag = false;
						prev_record[j] = 8;      // for non-consistent haplotypes
						break;
					}
				} else {
					record[record_idx++] = data[k][j];
			    }
			}
			if (consistent_flag) {
				if (prev_record[j] != 8) {
					count_idx = record[0];					
					for (int i = 1; i < tmp; i++) {
						count_idx = count_idx * 2 + record[i];
					}
					count[prev_record[j] * 8 + count_idx]++;
				}			
			}
		}

		int row = tmp == 1 ? 2 : 4;
		cout<<"\n\nsegment num\t"<<setw(4)<<s++<<endl;
		cout<<setw(5)<<"index";
		for (int j = 0; j < row; j++)
			cout<<setw(10)<<j;
		cout<<endl;
		for (int i = 0; i < 8; i++) {
			cout<<setw(5)<<i;
			for (int j = 0; j < row; j++)
				cout<<setw(10)<<count[i * row + j];
			cout<<endl;
		}
	}
	in.close();
//	out.close();

	delete [] prev_record;
	for (int i = 0; i < N; i++) {
		delete [] data[i];
	}
	delete [] data;
	delete [] pos;
	delete [] x;

	gettimeofday( &tend, NULL );
    int timeuse = 1000000 * ( tend.tv_sec - tstart.tv_sec ) + tend.tv_usec -tstart.tv_usec;
	printf("time: %d us\n", timeuse);
	return 0;
}
