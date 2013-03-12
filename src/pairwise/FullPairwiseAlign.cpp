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
#include <mpi.h>
#include <memory>
#include <numeric>

namespace clustalw
{

void FullPairwiseAlign::pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, int iEnd, int jStart, int jEnd)
{
  try
    {
      if(distMat->getSize() != alignPtr->getNumSeqs() + 1) {
        cerr << "The distance matrix is not the right size!\n"
             << "Need to terminate program.\n";
        exit(1);
      }

      if((iStart < 0) || (iEnd < iStart) || (jStart < 0) || (jEnd < jStart)) {
        cerr << "The range for pairwise Alignment is incorrect.\n"
             << "Need to terminate program.\n";
        exit(1);
      }
      
      if(ExtendData::numSeqs == 0) {
        return;
      }
      
      if (ExtendData::maxRes == 0) {
        cerr << "Could not get the substitution matrix\n";
        return;
      }

      double startTime = MPI_Wtime();

      broadcastExtendData();
      sendSequences(alignPtr, iStart, iEnd, jStart, jEnd);
      recieveDistMatrix(distMat);   
      
      double endTime = MPI_Wtime();
      cout << "Elapsed time for FullPairwiseAlign: " << setprecision(10) << endTime - startTime << " sec" << endl;       
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FullPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }
}

void FullPairwiseAlign::recieveDistMatrix(DistMatrix* distMat){
  int procNum;
  MPI_Comm_size(MPI_COMM_WORLD, &procNum);

  vector<dmRecord> emptyBuf(0);
  int emptySize = 0;
  vector<int> recieveCounts(procNum);
  vector<int> displacements(procNum);

  MPI_Gather(&emptySize, 1, MPI_INT, recieveCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  #ifdef DEBUG
    for (int i=0; i<recieveCounts.size(); i++) {
      cout << "Proc #" << i << " will send " << recieveCounts[i] << "sequences" << endl;
    }
  #endif
  displacements[0] = 0;
  displacements[1] = 0;
  for (int i=2; i<recieveCounts.size(); i++) {
    displacements[i] = displacements[i-1] + recieveCounts[i-1] ;
  }

  int sum = recieveCounts[recieveCounts.size()-1] + displacements[displacements.size()-1]; //std::accumulate(recieveCounts.begin(),recieveCounts.end(),0);
  vector<dmRecord> recvBuf(sum);

  cout << "Sum: " << sum << endl;
  // for (int i=0; i<displacements.size(); i++) {
  //   cout << "Displacement["<< i << "]=" << displacements[i] << endl;
  // }

  MPI_Gatherv(MPI_IN_PLACE, 0, ExtendData::mpi_dmRecord_type,
              recvBuf.data(), recieveCounts.data(), displacements.data(), ExtendData::mpi_dmRecord_type, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  
  cout << "After: \tRow: " << recvBuf[0].col << "\tCol: " << recvBuf[0].row << "\tVal: " << recvBuf[0].val << endl;
  
  for (int i=0; i<recvBuf.size(); i++) {
      int si = recvBuf[i].col,
          sj = recvBuf[i].row;
      float _score = recvBuf[i].val;

     //cout << "\tRow: " << si << "\tCol: " << sj << "\tVal: " << _score << endl;

      if ((si != -1) && (sj != -1)) { 
    //    cout << "DistMat(" << si << ", " << sj << ") = " << _score << endl;
        
        distMat->SetAt(si, sj, _score);    
        distMat->SetAt(sj, si, _score);
      
      // if(userParameters->getDisplayInfo()) {
      //     utilityObject->info("[%d]Sequences (%d:%d) Aligned. Score:  %d", proc, si, sj, (int)(100.f - (_score*100.f)));     
      // }        
      }      
    }  
    return;
  }
  

void FullPairwiseAlign::scheduleSequences(int numOfSeq, int* bounds) {
  //TODO: should be improved!
  int procNum;
  MPI_Comm_size( MPI_COMM_WORLD, &procNum ) ;
  --procNum; //#0 proc is master
  
  isInteger = true;

  int jStart = bounds[2],
      jEnd = bounds[3];
  int* seqHist = new int[numOfSeq];
  
  int sum = 0;
  //fill number of subseq to be aligned for aligning of one seq in array
  //      &&
  //calculate the sum of seq to be aligned 
  cout << "Histogram of seq" <<  endl;
  for (int i=0; i<numOfSeq; i++) {
    int initSj = utilityObject->MAX(i+1, jStart+1),
        boundSj = utilityObject->MIN(numOfSeq,jEnd);
    seqHist[i] = boundSj - initSj;
    sum += seqHist[i];
  }
  cout << endl;
  
  //TODO: should be moved soon
  cout << "Count of sequences: " << numOfSeq << endl;
  cout << "Count of sequences to be aligned: " << sum << endl;
  //calculate the average number of seq per proc

  int averageNumOfSeqPerProc = sum / procNum;
  cout << "Average per proc: " << averageNumOfSeqPerProc << endl;

  //simply increase averageNumOfSeqPerProc for 25%
  maxSeqCount = averageNumOfSeqPerProc + averageNumOfSeqPerProc*0.25;
  //TODO: should be refactored
  MPI_Bcast(&maxSeqCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

  cout << "Max sequences count: " << maxSeqCount << endl;
  int* portionPerProc = new int[procNum+1]; //+1 for saving initial index
  portionPerProc[0] = 0; //start index
  //calculate portion per proc
  int temporalSum = 0;
  int procIdx = 1;

  cout << "Portion per proc" << endl;
  for (int i=0; i<numOfSeq; i++) {
    temporalSum += seqHist[i];
    if (temporalSum >= averageNumOfSeqPerProc) {
      portionPerProc[procIdx++] = i;
      cout << "Proc#" << procIdx-1 << " should align " 
           << portionPerProc[procIdx-1] - portionPerProc[procIdx-2] 
           << " temporalSum = " << temporalSum << endl;
      temporalSum = 0;
    }
  }

  //keep in mind that procNum also have one minus
  if (procIdx - procNum == 0) {
    portionPerProc[procNum] = numOfSeq;
    cout  << "Proc#" << procNum << " should align " 
          << portionPerProc[procNum] - portionPerProc[procNum-1] << endl;
  }

  MPI_Bcast(portionPerProc, procNum+1, MPI_INT, 0, MPI_COMM_WORLD);
}

void FullPairwiseAlign::sendSequences(Alignment* alignPtr, int iStart, int iEnd, int jStart, int jEnd) {
  const SeqArray* _ptrToSeqArray = alignPtr->getSeqArray(); //This is faster! 

  int bounds[4] = {
     iStart,
     iEnd,
     jStart,
     jEnd
  };

  int initSi = utilityObject->MAX(0, iStart),
      boundSi = utilityObject->MIN(ExtendData::numSeqs,iEnd);
  
  int procNum;
  MPI_Comm_size( MPI_COMM_WORLD, &procNum ) ;
  //send sequences
  const int NUMBER_OF_SEQ = ExtendData::numSeqs - initSi;
  scheduleSequences(NUMBER_OF_SEQ, bounds);
  
  MPI_Bcast(&bounds, 4, MPI_INT, 0, MPI_COMM_WORLD);
  
  for (int si=initSi; si<initSi+NUMBER_OF_SEQ; si++) {
    std::vector<int> seq = (*_ptrToSeqArray)[si + 1];
    int* data = seq.data();
    int size = seq.size();

  #ifdef DEBUG
    cout << "Size: " << size << endl;        
  #endif
    
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data, seq.size(), MPI_INT, 0, MPI_COMM_WORLD);              
  }   

  return;
}

void FullPairwiseAlign::broadcastExtendData(){
  int intBuf[8] = {
    ExtendData::intScale,
    ExtendData::matAvgScore,
    ExtendData::maxRes,
    (int)ExtendData::DNAFlag,
    ExtendData::gapPos1,
    ExtendData::gapPos2,
    ExtendData::maxAlnLength,
    ExtendData::numSeqs
  };

  MPI_Bcast(&intBuf, 8, MPI_INT, 0, MPI_COMM_WORLD);

  float floatBuf[4] = {
    ExtendData::gapOpenScale,
    ExtendData::gapExtendScale, 
    ExtendData::pwGapOpen,
    ExtendData::pwGapExtend
  };


  MPI_Bcast(&floatBuf, 4, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  int* matrix = new int[clustalw::NUMRES*clustalw::NUMRES]; // linearization of 2d matrix
  for (int i=0; i<clustalw::NUMRES; i++) 
    for (int j=0; j<clustalw::NUMRES; j++) {
      matrix[i*clustalw::NUMRES+j] = ExtendData::matrix[i][j];
    }
  
  MPI_Bcast(matrix, clustalw::NUMRES*clustalw::NUMRES, MPI_INT, 0, MPI_COMM_WORLD);
  
  delete[] matrix;

  return;
}
}
