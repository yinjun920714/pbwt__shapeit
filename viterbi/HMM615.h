#ifndef __HMM_615_H 
#define __HMM_615_H 
#include "Matrix615.h" 
class HMM615 {
public:
	// parameters
	int nStates; // n : number of possible states
	int nObs; // m : number of possible output values
	int nTimes; // t : number of time slots with observations 
	std::vector<double> pis; // initial states
	std::vector<int> outs; // observed outcomes
	Matrix615<double> trans; // trans[i][j] corresponds to A_{ij} 
	Matrix615<double> emis;

	// storages for dynamic programming
	Matrix615<double> deltas; 
	Matrix615<int> phis;
	std::vector<int> path;
	HMM615(int states, int obs, int times) : nStates(states), nObs(obs), 
		nTimes(times), trans(states, states, 0), emis(states, obs, 0), 
		deltas(times, states, 0),
		phis(times, states, 0) {
		pis.resize(nStates);
		path.resize(nTimes); }

	void viterbi() {
		for(int i=0; i < nStates; ++i) {
			deltas.data[0][i] = pis[i] * emis.data[i][ outs[0] ]; }
		for(int t=1; t < nTimes; ++t) { 
			for(int i=0; i < nStates; ++i) {
				int maxIdx = 0;
				double maxVal = deltas.data[t-1][0] * trans.data[0][i]
								* emis.data[i][ outs[t] ]; 
				for(int j=1; j < nStates; ++j) {
					double val = deltas.data[t-1][j] * trans.data[j][i] * emis.data[i][ outs[t] ];
					if ( val > maxVal ) { maxIdx = j; maxVal = val; } }
				deltas.data[t][i] = maxVal;
				phis.data[t][i] = maxIdx; 
			}
		}

		// backtrack viterbi path
		double maxDelta = deltas.data[nTimes-1][0];
		path[nTimes-1] = 0;
		for(int i=1; i < nStates; ++i) {
			if ( maxDelta < deltas.data[nTimes-1][i] ) {
				maxDelta = deltas.data[nTimes-1][i];
				path[nTimes-i] = i;
			}
		}
		for(int t=nTimes-2; t >= 0; --t) { 
			path[t] = phis.data[t+1][ path[t+1] ];
		} 
	}
};
#endif  // __HMM_615_H

