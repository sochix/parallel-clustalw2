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
      FullPairwiseAlign() {};
		  virtual ~FullPairwiseAlign(){};
      virtual void pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, 
                                   int iEnd, int jStart, int jEnd); 
      
    private:
      /* Functions */
      //MPI
      void BroadcastExtendData();
      void BroadcastSequencesAndBounds(Alignment*, int,int,int,int);
      void GatherDistMatrix(DistMatrix*);
      void CalculateNumbersOfSeqToAlignForeachIteration(const int numOfSeq, int jStart, int jEnd, std::vector<int>& vec, int& countOfSequences);
      void SchedulePortions(const vector<int>& vec, int);
    
      static int NUMBER_OF_SEQ;    
};

}
#endif
