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
      void BroadcastSequences(Alignment*, int,int,int,int);
      void GatherDistMatrix(DistMatrix*);
      void SchedulePortionOfSequencesForEachProc(const int, int*); //has side effect!
};

}
#endif
