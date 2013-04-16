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

int FullPairwiseAlign::NUMBER_OF_SEQ;

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

      BroadcastExtendData();
      BroadcastSequencesAndBounds(alignPtr, iStart, iEnd, jStart, jEnd);

      vector<int> numberOfSeqToAlignPerIteration(NUMBER_OF_SEQ);
      int countOfSeq = 0;

      CalculateNumbersOfSeqToAlignForeachIteration(NUMBER_OF_SEQ, jStart, jEnd, numberOfSeqToAlignPerIteration, countOfSeq);
      SchedulePortions(numberOfSeqToAlignPerIteration, countOfSeq);
      GatherDistMatrix(distMat);   
      
      double endTime = MPI_Wtime();

      cerr << "\tElapsed time for FullPairwiseAlign: " << setprecision(10) << endTime - startTime << " sec" << endl;       
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FullPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }
}

void FullPairwiseAlign::GatherDistMatrix(DistMatrix* distMat){
  
  int procNum;
  MPI_Comm_size(MPI_COMM_WORLD, &procNum);
  
  //Recieve the size of distMatrix for each proc
  vector<int> recieveCounts(procNum);
  MPI_Gather(MPI_IN_PLACE, 0, MPI_INT, recieveCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
 
  vector<int> displacements(procNum);
  displacements[0] = 0;
  displacements[1] = 0;
  for (int i=2; i<recieveCounts.size(); i++) {
    displacements[i] = displacements[i-1] + recieveCounts[i-1] ;
  }

  //Receive the distMatrix from each proc
  int sum = recieveCounts[recieveCounts.size()-1] + displacements[displacements.size()-1];
  vector<dmRecord> recvBuf(sum);
  MPI_Gatherv(MPI_IN_PLACE, 0, ExtendData::mpi_dmRecord_type,
              recvBuf.data(), recieveCounts.data(), displacements.data(), ExtendData::mpi_dmRecord_type, 0, MPI_COMM_WORLD);
 
  for (int i=0; i<recvBuf.size(); i++) {
      int si = recvBuf[i].col,
          sj = recvBuf[i].row;
      float _score = recvBuf[i].val;

      #ifdef DEBUG
        cout << "\tRow: " << si << "\tCol: " << sj << "\tVal: " << _score << endl;
      #endif

      if ((si != -1) && (sj != -1)) { 
        distMat->SetAt(si, sj, _score);    
        distMat->SetAt(sj, si, _score);
        
        //TODO: don't forget to uncomment
        // if(userParameters->getDisplayInfo()) {
        //     utilityObject->info("[%d]Sequences (%d:%d) Aligned. Score:  %d", proc, si, sj, (int)(100.f - (_score*100.f)));     
        // }        
      }      
    }  
    return;
  }
  

void FullPairwiseAlign::CalculateNumbersOfSeqToAlignForeachIteration(const int numOfSeq, int jStart, int jEnd, vector<int>& vec, int& countOfSequences) {

  const int boundSj = utilityObject->MIN(numOfSeq,jEnd);
  for (int i=0; i<numOfSeq; i++) {
    const int initSj = utilityObject->MAX(i+1, jStart+1);
    vec[i] = boundSj - initSj;
    countOfSequences += vec[i];
  }
    
  cout << endl << "Scheduler data" << endl;
  cout << "\tCount of sequences: " << numOfSeq << endl;
  cout << "\tCount of sequences to be aligned:" << countOfSequences << endl;
}

void FullPairwiseAlign::BroadcastSequencesAndBounds(Alignment* alignPtr, int iStart, int iEnd, int jStart, int jEnd) {

  //broadcast initial bounds
  int bounds[4] = {
     iStart,
     iEnd,
     jStart,
     jEnd
  };

  MPI_Bcast(&bounds, 4, MPI_INT, 0, MPI_COMM_WORLD);

  int initSi = utilityObject->MAX(0, iStart);
  NUMBER_OF_SEQ = ExtendData::numSeqs - initSi; //actual number of sequences to align for all procs

  //Broadcast sequences
  const SeqArray* _ptrToSeqArray = alignPtr->getSeqArray(); //This is faster! 
  
  //TODO: think how to improve
  for (int si=initSi; si<initSi+NUMBER_OF_SEQ; si++) {
    vector<int> seq = (*_ptrToSeqArray)[si + 1];
    int size = seq.size();
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(seq.data(), seq.size(), MPI_INT, 0, MPI_COMM_WORLD);              
  }   

  return;
}

void FullPairwiseAlign::BroadcastExtendData(){
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
  
  vector<int> linearizedMatrix(clustalw::NUMRES*clustalw::NUMRES);

  for (int i=0; i<clustalw::NUMRES; i++) 
    for (int j=0; j<clustalw::NUMRES; j++) {
      linearizedMatrix[i*clustalw::NUMRES+j] = ExtendData::matrix[i][j];
    }
  
  MPI_Bcast(linearizedMatrix.data(), clustalw::NUMRES*clustalw::NUMRES, MPI_INT, 0, MPI_COMM_WORLD);
  
  return;
}

void FullPairwiseAlign::SchedulePortions(const vector<int>& vec, int count) {
  
  int procNum;
  MPI_Comm_size(MPI_COMM_WORLD, &procNum);
  procNum--;

  const int k = 16; //some coefficient, power of 2
  
  bool hasNextPortion = true;
  int i = 2;
  int currentSeqPortion = count / (procNum * k),
      start = 0,
      length = 0;  
  int curProc = 0;
  int rank = -1;
  int iterNum = 0;

  while (true) {
    if (hasNextPortion) {
      
      int sum = 0,
          idx = start;

      //decrease portion size
      if (iterNum == procNum) {
        if (i < k) {
          currentSeqPortion = (count * i) / (procNum * k);
          i+=2;
          iterNum = 0;
          cout << "currentSeqPortion is "<<currentSeqPortion << endl;  
        }
        continue;
      }
      

      //finding correct iterations to perform for align neede number of seq
      while ((sum < currentSeqPortion) && (idx<vec.size())) {
        sum += vec[idx++];
      }

      //means that all seq are aligned
      if (idx == vec.size()) 
        hasNextPortion = false;


      length = idx - start;

      int sendBuf[2] = {start,idx};

      start = idx;

      MPI_Recv(&rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      #ifdef DEBUG
        cout << "Master process recieved request form proc #" << rank << endl;
      #endif
      MPI_Send(&sendBuf, 2, MPI_INT, rank, 0, MPI_COMM_WORLD);  
      iterNum++;

    } else {
      int errorBuf[2] = {-1, -1};
      MPI_Recv(&rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&errorBuf, 2, MPI_INT, rank, 0, MPI_COMM_WORLD);
      curProc++;

      if (curProc == procNum)
        return;
    }    
  }
}

}
