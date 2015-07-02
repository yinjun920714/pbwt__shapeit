#include <iostream>
#include <iomanip>
#include "HMM615.h"
int main(int argc, char** argv) {
	std::vector<int> toss; 
	std::string tok;
	while( std::cin >> tok ) {
		if ( tok == "H" ) toss.push_back(0);
		else if ( tok == "T" ) toss.push_back(1);
		else {
			std::cerr << "Cannot recognize input " << tok << std::endl;
			return -1;
		}
	}
	int T = toss.size();
	HMM615 hmm(2, 2, T);
	
	hmm.trans.data[0][0] = 0.95; hmm.trans.data[0][1] = 0.05;
	hmm.trans.data[1][0] = 0.2; hmm.trans.data[1][1] = 0.8;
	hmm.emis.data[0][0] = 0.5; hmm.emis.data[0][1] = 0.5;
	hmm.emis.data[1][0] = 0.9; hmm.emis.data[1][1] = 0.1;

	hmm.pis[0] = 0.5; hmm.pis[1] = 0.5;
	hmm.outs = toss;
	
	hmm.viterbi();
	
	std::cout << "TIME\tTOSS\tP(FAIR)\tP(BIAS)\tMLSTATE" << std::endl;
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(4);
	for(int t=0; t < T; ++t) {
		std::cout << t+1 << "\t" << (toss[t] == 0 ? "H" : "T") << "\t"
				<< hmm.deltas.data[t][0] << "\t" << hmm.deltas.data[t][1] << "\t" 
				<< (hmm.path[t] == 0 ? "FAIR" : "BIASED" ) << std::endl;
	}
	return 0;
}