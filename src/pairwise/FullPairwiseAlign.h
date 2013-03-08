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
        /* Attributes */

    private:
        /* Attributes */
        int isInteger;
        int portionPerProc;
        int lastProcPortion;        
        /* Functions */
       //MPI
       void broadcastExtendData();
       void sendSequences(Alignment*, int,int,int,int);
       void recieveDistMatrix(DistMatrix*);
       void scheduleSequences(int, int*); //has side effect!
};

}
#endif
