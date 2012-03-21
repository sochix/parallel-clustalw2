#ifndef __SWALGO_H__
#define __SWALGO_H__

#include "PairwiseAlignBase.h"
#include "I_ExtendData.h"

namespace clustalw {

class SWAlgo
{
	public:
		SWAlgo();
		~SWAlgo();
		void Pass(const vector<int>* seq1, const vector<int>* seq2, int n, int m, const int, const int);
		
		/*FIXME: it might be consts!*/
		int maxScore;
		int se1;
		int se2;
		int sb1;
		int sb2;
		
	private:
		void forwardPass(const vector<int>* seq1, const vector<int>* seq2, int n, int m, const int, const int);
		void reversePass(const vector<int>* ia, const vector<int>* ib, const int, const int);
		
		vector<int> HH;
    vector<int> DD;		
};
}
#endif /*__SWALGO_H__*/
