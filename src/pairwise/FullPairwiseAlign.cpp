/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "FullPairwiseAlign.h"
#include <math.h>
#include <omp.h>

namespace clustalw
{

//=========START FULLPAIRWISE ALIGN======================

FullPairwiseAlign::FullPairwiseAlign()
: mmScore(0),
  _gapOpen(0),
  _gapExtend(0),
  seq1(0),
  seq2(0),
  maxScore(0)
{

}

void FullPairwiseAlign::pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, int iEnd, int jStart, int jEnd)
{
    int si, sj, i;
    int n, m, len1, len2;
    int res;
    double _score;
    
    try
    {
        
        if(distMat->getSize() != alignPtr->getNumSeqs() + 1)
        {
            cerr << "The distance matrix is not the right size!\n"
                 << "Need to terminate program.\n";
            exit(1);
        }
        if((iStart < 0) || (iEnd < iStart) || (jStart < 0) || (jEnd < jStart))
        {
            cerr << "The range for pairwise Alignment is incorrect.\n"
                 << "Need to terminate program.\n";
            exit(1);
        }
        
        if(ExtendData::numSeqs == 0)
        {
            return;
        }
        
        if (ExtendData::maxRes == 0)
        {
            cerr << "Could not get the substitution matrix\n";
            return;
        }
           
       	const SeqArray* _ptrToSeqArray = alignPtr->getSeqArray(); //This is faster! 
    
		    SWAlgo swalgo;
		    MMAlgo mmalgo;
    
    		double startTime = omp_get_wtime();
    
        for (si = utilityObject->MAX(0, iStart); si < ExtendData::numSeqs && si < iEnd; si++)
        {
            n = alignPtr->getSeqLength(si + 1);
            len1 = 0;
            for (i = 1; i <= n; i++)
            {
                res = (*_ptrToSeqArray)[si + 1][i];
                if ((res != ExtendData::gapPos1) && (res != ExtendData::gapPos2))
                {
                    len1++;
                }
            }
						
						for (sj = utilityObject->MAX(si+1, jStart+1); sj < ExtendData::numSeqs && sj < jEnd; sj++)
            {
                m = alignPtr->getSeqLength(sj + 1);
                if (n == 0 || m == 0)
                {
                    distMat->SetAt(si + 1, sj + 1, 1.0);
                    distMat->SetAt(sj + 1, si + 1, 1.0);
                    continue;
                }
                len2 = 0;
                for (i = 1; i <= m; i++)
                {
                    res = (*_ptrToSeqArray)[sj + 1][i];
                    if ((res != ExtendData::gapPos1) && (res != ExtendData::gapPos2))
                    {
                        len2++;
                    }
                }
								ExtendData::UpdateGapOpenAndExtend(_gapOpen, _gapExtend, n, m);
                // align the sequences
            
                seq1 = si + 1;
                seq2 = sj + 1;

                _ptrToSeq1 = alignPtr->getSequence(seq1);
                _ptrToSeq2 = alignPtr->getSequence(seq2);
            
            		swalgo.Pass(_ptrToSeq1, _ptrToSeq2, n, m, _gapOpen, _gapExtend);
            
   	            // use Myers and Miller to align two sequences 

								
								maxScore = mmalgo.Pass(swalgo.sb1 - 1, swalgo.sb2 - 1, swalgo.se1 - swalgo.sb1 + 1, swalgo.se2 - swalgo.sb2 + 1,
                    (int)0, (int)0, _ptrToSeq1, _ptrToSeq2, _gapOpen, _gapExtend);
      
                 // calculate percentage residue identity

                mmScore = tracePath(swalgo.sb1, swalgo.sb2, mmalgo.displ, mmalgo.printPtr);

                if (len1 == 0 || len2 == 0)
                {
                    mmScore = 0;
                }
                else
                {
                    mmScore /= (float)utilityObject->MIN(len1, len2);
                }

                _score = ((float)100.0 - mmScore) / (float)100.0;
                distMat->SetAt(si + 1, sj + 1, _score);
                distMat->SetAt(sj + 1, si + 1, _score);
                
                if(userParameters->getDisplayInfo())
                {
                    utilityObject->info("Sequences (%d:%d) Aligned. Score:  %d",
                                        si+1, sj+1, (int)mmScore);     
                }
            }
        }
        double endTime = omp_get_wtime() - startTime;
        cout << endl << "[OMP] Elapsed time: " << endTime << " .sec" << endl;
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FullPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }
}

float FullPairwiseAlign::tracePath(int tsb1, int tsb2, vector<int>& displ, int printPtr)
{
    int res1, res2;
    int i1, i2;
    int i, k, pos, toDo;
    int count;
    float score;

    toDo = printPtr - 1;
    i1 = tsb1;
    i2 = tsb2;

    pos = 0;
    count = 0;
    for (i = 1; i <= toDo; ++i)
    {
        if (displ[i] == 0)
        {
            res1 = (*_ptrToSeq1)[i1];
            res2 = (*_ptrToSeq2)[i2];

            if ((res1 != userParameters->getGapPos1()) && 
                (res2 != userParameters->getGapPos2()) && (res1 == res2))
            {
                count++;
            }
            ++i1;
            ++i2;
            ++pos;
        }
        else
        {
            if ((k = displ[i]) > 0)
            {
                i2 += k;
                pos += k;
            }
            else
            {
                i1 -= k;
                pos -= k;
            }
        }
    }
    
    score = 100.0 *(float)count;
    return (score);
}




int FullPairwiseAlign::gap(int k)
{
    if(k <= 0)
    {
        return 0;
    }
    else
    {
        return _gapOpen + _gapExtend * k;
    }
}


}
