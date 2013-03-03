/*
*     Implementation of diff algorithm described in paper of Eugene W. Mylers & Webb Miller 
*     "Optimal alignments in linear space" //CABIOS Vol.4. no.1.1988 Pages 11-17
*      
*     Author: Mark Larkin
*     Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
*
*     Refactored by Ilya Pirozhenko. 2012
*/

#ifndef __MMALGO_H__
#define __MMALGO_H__

#include "PairwiseAlignBase.h"
#include "../parallel/ParallelAlgo.h"

namespace clustalw {

class MMAlgo
{
public:
	MMAlgo(ExternalData*);
	~MMAlgo();
	int Pass(int, int, int, int, int, int, const vector<int>*, const vector<int>*, const int, const int);

	//FIXME: must be in private
	int printPtr;
	vector<int> displ;
      const ExternalData* data;
	
private:
      enum type_t {TYPE_1, TYPE_2};

	int diff(int, int, int, int, int, int);
      void addToDisplay(int);
      void delFromDisplay(int);
      int calcScore(int iat, int jat, int v1, int v2); 
      int gapAffineFunction(int k, int t);
      void forwardPass(int, int, int, int, int, int);
      void backwardPass(int, int, int, int, int, int);

      //for legacy
      int tbgap(int k, int tb);
      int tegap(int k, int te);
      
      int lastPrint;
      vector<int> HH; // HH(j) - minimum cost of a conversion seq1(i, len1) to seq2(j,len2)
      vector<int> DD; // DD(j) - minimum cost of a conversion seq1(i, len1) to seq2(j, len2) that ends with a delete
      vector<int> RR; // RR(N-j) - minumum cost of a conversion of seq1(i, len1)^T to seq2(j, len2)^T
      vector<int> SS; // SS(N-j) - minumum cost of a conversion of seq1(i, len1)^T to seq2(j, len2)^T begins with a delete
      
      const vector<int>* _ptrToSeq1;
      const vector<int>* _ptrToSeq2;
      
      int _gapOpen;
      int _gapExtend;

	//Moved from diff, where they were static
   /*   int f;
      int hh;
      int e;
      int s;
      int t;*/
      
};

}
#endif /*__MMALGO_H__*/
