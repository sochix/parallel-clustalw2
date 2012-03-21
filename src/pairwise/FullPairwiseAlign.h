/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef FULLPAIRWISEALIGN_H
#define FULLPAIRWISEALIGN_H

#include "PairwiseAlignBase.h"
#include "I_SWAlgo.h"
#include "I_MMAlgo.h"
#include "I_ExtendData.h"

namespace clustalw
{

class FullPairwiseAlign : public PairwiseAlignBase
{
    public:
        /* Functions */
        FullPairwiseAlign();
				virtual ~FullPairwiseAlign(){};

        virtual void pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, 
                                   int iEnd, int jStart, int jEnd); 
        /* Attributes */

    private:
        /* Functions */
        //Algo MM
        /*
        int diff(int A, int B, int M, int N, int tb, int te);
        void add(int v);
        void del(int k);
        int calcScore(int iat, int jat, int v1, int v2); 
        int tbgap(int k, int tb);
        int tegap(int k, int te);\
        */
        //Other stuff
        float tracePath(int tsb1, int tsb2, vector<int>&, int );       
        int gap(int k);
        
        /* Attributes */
        // I have constant pointers to the data. This allows for the fastest access.
        const vector<int>* _ptrToSeq1;
        const vector<int>* _ptrToSeq2;
        float mmScore;
   /*
        int printPtr;
        int lastPrint;
        vector<int> displ;
      
        vector<int> HH;
        vector<int> DD;
        vector<int> RR;
        vector<int> SS;
		*/
        int _gapOpen; // scaled to be an integer, this is not a mistake
        int _gapExtend; // scaled to be an integer, not a mistake
     
        int seq1;
        int seq2;
     
        int maxScore;

};

}
#endif
