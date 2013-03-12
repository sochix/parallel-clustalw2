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

      BroadcastExtendData();
      BroadcastSequences(alignPtr, iStart, iEnd, jStart, jEnd);
      GatherDistMatrix(distMat);   
      
      double endTime = MPI_Wtime();

      cout << "\tElapsed time for FullPairwiseAlign: " << setprecision(10) << endTime - startTime << " sec" << endl;       
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
  MPI_Barrier(MPI_COMM_WORLD);

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
  MPI_Barrier(MPI_COMM_WORLD);
  
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
  

void FullPairwiseAlign::SchedulePortionOfSequencesForEachProc(const int numOfSeq, int* bounds) {
  //TODO: should be improved!
  int countOfAllProcs;
  MPI_Comm_size( MPI_COMM_WORLD, &countOfAllProcs);
  int countOfWorkerProcs = countOfAllProcs - 1; //#0 proc is master
  
  int jStart = bounds[2],
      jEnd = bounds[3];
  
  vector<int> numberOfSequencesToAlign(numOfSeq); //number of sequences to align for getting an align of index sequence
  int countOfSequences = 0;
  
  const int boundSj = utilityObject->MIN(numOfSeq,jEnd);
  for (int i=0; i<numOfSeq; i++) {
    const int initSj = utilityObject->MAX(i+1, jStart+1);
    numberOfSequencesToAlign[i] = boundSj - initSj;
    countOfSequences += numberOfSequencesToAlign[i];
  }
    
  cout << endl << "Scheduler data" << endl;
  cout << "\tCount of sequences: " << numOfSeq << endl;
  cout << "\tCount of sequences to be aligned:" << countOfSequences << endl;

  //calculate the average number of seq per proc
  int averageNumOfSeqPerProc = countOfSequences / countOfWorkerProcs;
  cout << "\tAverage per proc: " << averageNumOfSeqPerProc << endl;

  //calculate portion per proc
  vector<int> portionPerProc(countOfAllProcs, 0);
  
  int procIdx,
      temporalSum;

  //TODO: need new method for approximation, maybe recursive binary search?
  //iterates to find best portion for each proc
  
  while (procIdx != countOfAllProcs) { //procIdx != countOfAllProcs means that not all procs scheduled
    procIdx = 1;
    temporalSum = 0;  

    for (int i=0; i<numOfSeq; i++) {
      temporalSum += numberOfSequencesToAlign[i];
      if (temporalSum >= averageNumOfSeqPerProc) {
        portionPerProc[procIdx++] = i;
        #ifdef DEBUG
          cout << "\tProc#" << procIdx-1 << " should align " << portionPerProc[procIdx-1] - portionPerProc[procIdx-2] << " sequences" 
                << "Temporal sum: " << temporalSum << endl;  
        #endif
        temporalSum = 0;
      }
    }
    averageNumOfSeqPerProc--; //TODO: think!!!    
  }
  
  //HACK: expand last proc
  portionPerProc[countOfAllProcs-1] = numOfSeq - 1;
  averageNumOfSeqPerProc++; //TODO: think!!!    

  cout << "\tAverage per proc after iterating: " << averageNumOfSeqPerProc << endl;
  for (int i=1; i<countOfAllProcs; i++) {
    cout << "\tProc#" << i << " should align " << portionPerProc[i] - portionPerProc[i-1] << " sequences" << endl;  
  }
  
  MPI_Bcast(portionPerProc.data(), countOfAllProcs, MPI_INT, 0, MPI_COMM_WORLD);
}

void FullPairwiseAlign::BroadcastSequences(Alignment* alignPtr, int iStart, int iEnd, int jStart, int jEnd) {

  //broadcast initial bounds
  int bounds[4] = {
     iStart,
     iEnd,
     jStart,
     jEnd
  };

  MPI_Bcast(&bounds, 4, MPI_INT, 0, MPI_COMM_WORLD);

  int initSi = utilityObject->MAX(0, iStart);
  const int NUMBER_OF_SEQ = ExtendData::numSeqs - initSi; //actual number of sequences to align for all procs

  SchedulePortionOfSequencesForEachProc(NUMBER_OF_SEQ, bounds);
  
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
}
