#ifndef __MMALGO_H__
#define __MMALGO_H__

#include "PairwiseAlignBase.h"
#include "I_ExtendData.h"

namespace clustalw {

class MMAlgo
{
public:
	MMAlgo();
	~MMAlgo();
	int Pass(int, int, int, int, int, int, const vector<int>*, const vector<int>*, const int, const int);

	//FIXME: must be in private
	int printPtr;
	vector<int> displ;
	
private:
	int diff(int, int, int, int, int, int);
      void addToDisplay(int);
      void delFromDisplay(int);
      int calcScore(int iat, int jat, int v1, int v2); 
      int tbgap(int k, int tb);
      int tegap(int k, int te);
      
      int lastPrint;
      vector<int> HH;
      vector<int> DD;
      vector<int> RR;
      vector<int> SS;
      
      const vector<int>* _ptrToSeq1;
      const vector<int>* _ptrToSeq2;
      
      int _gapOpen;
      int _gapExtend;

	//Moved from diff, where they were static
      int f;
      int hh;
      int e;
      int s;
      int t;
      
};

}
#endif /*__MMALGO_H__*/
