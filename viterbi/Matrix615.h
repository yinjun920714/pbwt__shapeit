#ifndef __MATRIX_615_H // to avoid multiple inclusion of same headers #define __MATRIX_615_H
#define __MATRIX_615_H
#include <vector>
template <class T> 
class Matrix615 { 
public:
	std::vector< std::vector<T> > data; 
	Matrix615(int nrow, int ncol, T val = 0) {
		data.resize(nrow); // make n rows 
		for(int i=0; i < nrow; ++i) {
			data[i].resize(ncol,val); // make n cols with default value val 
		}
}
	int rowNums() { return (int)data.size(); }
	int colNums() { return ( data.size() == 0 ) ? 0 : (int)data[0].size(); }
};

#endif // __MATRIX_615_H