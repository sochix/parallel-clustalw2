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
#include <mpi.h>

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
    
    const vector<int>* _ptrToSeq1 = NULL;
    const vector<int>* _ptrToSeq2 = NULL;
    float mmScore = mmScore;

    int _gapOpen = _gapOpen; // scaled to be an integer, this is not a mistake
    int _gapExtend = _gapExtend; // scaled to be an integer, not a mistake
 
    int seq1 = seq1;
    int seq2 = seq2;
 	  int maxScore = maxScore;
    
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
         
     	const SeqArray* _ptrToSeqArray = alignPtr->getSeqArray(); //This is faster! 
      
      //MPI         
      int myrank, ntasks;

      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

      sendExtendData();
      

      //start sending sequences
      int bounds[4] = {
        iStart,
        iEnd,
        jStart,
        jEnd
      };

      MPI_Send(&bounds, 4, MPI_INT, 1, 0, MPI_COMM_WORLD);

      int initSi = utilityObject->MAX(0, iStart),
          boundSi = utilityObject->MIN(ExtendData::numSeqs,iEnd);

      //send sequences
      const int NUMBER_OF_SEQ = ExtendData::numSeqs - initSi;

      for (si=initSi; si<initSi+NUMBER_OF_SEQ; si++) {
        std::vector<int> seq = (*_ptrToSeqArray)[si + 1];
        int* data = seq.data();

        int size = seq.size();

      #ifdef DEBUG
        cout << "Size: " << size << endl;        
      #endif
        
        MPI_Send(&size, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Send(data, seq.size(), MPI_INT, 1, 0, MPI_COMM_WORLD);              
      }   

      int size;
      MPI_Recv(&size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      cout << "Size of distMat: " << size << endl;

      float* unwindedDistMat = new float[size];
      MPI_Recv(unwindedDistMat, size, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for (int i=0; i<size; i+=3) {
        distMat->SetAt((int)unwindedDistMat[i], (int)unwindedDistMat[i+1], unwindedDistMat[i+2]);        
      }             
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FullPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }
}

void FullPairwiseAlign::sendExtendData(){
  //Broadcast ExtendData to workers
  //TODO: change from send to broadcast
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

  MPI_Send(&intBuf, 8, MPI_INT, 1, 0, MPI_COMM_WORLD);

  float floatBuf[4] = {
    ExtendData::gapOpenScale,
    ExtendData::gapExtendScale, 
    ExtendData::pwGapOpen,
    ExtendData::pwGapExtend
  };

  MPI_Send(&floatBuf, 4, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);

  int* matrix = new int[clustalw::NUMRES*clustalw::NUMRES]; // linearization of 2d matrix
  for (int i=0; i<clustalw::NUMRES; i++) 
    for (int j=0; j<clustalw::NUMRES; j++) {
      matrix[i*clustalw::NUMRES+j] = ExtendData::matrix[i][j];
    }
  
  MPI_Send(matrix, clustalw::NUMRES*clustalw::NUMRES, MPI_INT, 1, 0, MPI_COMM_WORLD);
  delete[] matrix;

  return;
}

float FullPairwiseAlign::tracePath(int tsb1, int tsb2, const vector<int>& displ, int printPtr,   const vector<int>* _ptrToSeq1,  const vector<int>* _ptrToSeq2)
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

}
