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

      broadcastExtendData();
      sendSequences(alignPtr, iStart, iEnd, jStart, jEnd);
      recieveDistMatrix(distMat);     

      
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FullPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }
}

void FullPairwiseAlign::recieveDistMatrix(DistMatrix* distMat){
  int size;
  MPI_Recv(&size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  #ifdef DEBUG
    cout << "Size of distMat: " << size << endl;
  #endif

  float* unwindedDistMat = new float[size];
  MPI_Recv(unwindedDistMat, size, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  for (int i=0; i<size; i+=3) {
    int si = (int)unwindedDistMat[i],
        sj = (int)unwindedDistMat[i+1];

    float _score = unwindedDistMat[i+2];
    
    distMat->SetAt(si, sj, _score);    
    distMat->SetAt(sj, si, _score);
    
    // if(userParameters->getDisplayInfo()) {
    //    utilityObject->info("Sequences (%d:%d) Aligned. Score:  %d", si, sj, (int)(100.f - (_score*100.f)));     
    // }
  }

  delete[] unwindedDistMat;             
  return;
}

void FullPairwiseAlign::scheduleSequences(int numOfSeq) {
  //TODO: should be improved!
  int procNum;
  MPI_Comm_size( MPI_COMM_WORLD, &procNum ) ;
  --procNum; //#0 proc is master
  
  isInteger = false;

  if (numOfSeq % procNum == 0)  {
    portionPerProc = numOfSeq/procNum;
  } else {
    isInteger = true;
    portionPerProc = numOfSeq/procNum;
    lastProcPortion = numOfSeq - portionPerProc*(procNum-1);
  }

  cout << "Scheduler information" << endl;
  cout << "numOfSeq: " << numOfSeq << endl;
  cout << "Num of procs: " << procNum << endl;
  cout << "isInteger: " << isInteger << endl;
  cout << "portion of seq per proc: " << portionPerProc << endl;
  cout << "last proc portion: " << lastProcPortion << endl;
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

  //send sequences
  const int NUMBER_OF_SEQ = ExtendData::numSeqs - initSi;
  scheduleSequences(NUMBER_OF_SEQ);

  //send init data
  
  MPI_Send(&bounds, 4, MPI_INT, 1, 0, MPI_COMM_WORLD); 

  for (int si=initSi; si<initSi+NUMBER_OF_SEQ; si++) {
    std::vector<int> seq = (*_ptrToSeqArray)[si + 1];
    int* data = seq.data();
    int size = seq.size();

  #ifdef DEBUG
    cout << "Size: " << size << endl;        
  #endif
    
    MPI_Send(&size, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Send(data, seq.size(), MPI_INT, 1, 0, MPI_COMM_WORLD);              
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
