#include <mpi.h>
#include <iostream>
#include "../alignment/Alignment.h"
#include "ParallelAlgo.h"
#include "../pairwise/I_MMAlgo.h"
#include "../pairwise/I_SWAlgo.h"

//#define DEBUG

using namespace std;
using namespace clustalw;

ExternalData ParallelAlgo::data;
SeqArray ParallelAlgo::seqArray;
int ParallelAlgo::iStart;
int ParallelAlgo::iEnd;
int ParallelAlgo::jStart;
int ParallelAlgo::jEnd;


void ParallelAlgo::DoFullPairwiseAlignment() {

	int r;
	MPI_Comm_rank(MPI_COMM_WORLD, &r);
	
	//init steps
	recieveExtendData();
	recieveSequences();

	int initSi = utilityObject->MAX(0, iStart),
        boundSi = utilityObject->MIN(data.numSeqs,iEnd),
        len1,
        i,
        si,
        sj,
        n,
        m,
        res, 
        seq1,
        seq2;
    double mmScore,
        	_score;
        

    int _gapOpen = 0; // scaled to be an integer, this is not a mistake
    int _gapExtend = 0; // scaled to be an integer, not a mistake
    int maxScore = 0;

    const vector<int>* _ptrToSeq1 = NULL;
    const vector<int>* _ptrToSeq2 = NULL;
    vector<distMatrixRecord> distMat;

	for (si = initSi; si < boundSi; si++) {
	    //n = alignPtr->getSeqLength(si + 1);
	    n = seqArray[translateIndex(si+1, initSi+1)].size()-1;
	    len1 = 0;
	    for (i = 1; i <= n; i++) {
	    //res = (*_ptrToSeqArray)[si + 1][i];	
	      res = seqArray[translateIndex(si + 1, initSi+1)][i];
	      if ((res != data.gapPos1) && (res != data.gapPos2)) {
	        len1++;
	      }
	    }
					
		// Simplify loop vars for OMP
		int initSj = utilityObject->MAX(si+1, jStart+1),
		    boundSj = utilityObject->MIN(data.numSeqs,jEnd),
		    len2;
																
		for (sj = initSj; sj < boundSj; sj++) {
			//m = alignPtr->getSeqLength(sj + 1);
			m = seqArray[translateIndex(sj + 1, initSj+1)].size()-1;
		          
		  	if (n == 0 || m == 0) {
		  		//THINK!
		 		distMat.push_back(distMatrixRecord(si+1, sj+1,1.0));
		 		distMat.push_back(distMatrixRecord(sj+1, si+1,1.0));
		 		//->SetAt(si + 1, sj + 1, 1.0);
		   		//distMat->SetAt(sj + 1, si + 1, 1.0);
		    	continue;
  			}

	      	len2 = 0;
	      	for (i = 1; i <= m; i++) {
	        	//res = (*_ptrToSeqArray)[sj + 1][i];
	        	res = seqArray[translateIndex(sj + 1, initSj+1)][i];
	        	if ((res != data.gapPos1) && (res != data.gapPos2)) {
	          		len2++;
	        	}
	      	}
              
     		data.UpdateGapOpenAndExtend(_gapOpen, _gapExtend, n, m);
            
      		// align the sequences
    		seq1 = translateIndex(si + 1, initSi+1);
      		seq2 = translateIndex(sj + 1, initSj+1);

		    _ptrToSeq1 = &seqArray[seq1];
		    _ptrToSeq2 = &seqArray[seq2];

			SWAlgo swalgo(&data);
			MMAlgo mmalgo(&data);
			
			#ifdef DEBUG
				cout << "Starting align: " << seq1 << ":" <<seq2 << endl;
			#endif

    		swalgo.Pass(_ptrToSeq1, _ptrToSeq2, n, m, _gapOpen, _gapExtend);
                               
    		//use Myers and Miller to align two sequences 
			maxScore = mmalgo.Pass(swalgo.sb1 - 1, swalgo.sb2 - 1, swalgo.se1 - swalgo.sb1 + 1, swalgo.se2 - swalgo.sb2 + 1,
			    (int)0, (int)0, _ptrToSeq1, _ptrToSeq2, _gapOpen, _gapExtend);	
			#ifdef DEBUG
				cout << "maxScore of mmalgo: " << maxScore << endl;		               
			#endif
     		// calculate percentage residue identity
      		mmScore = tracePath(swalgo.sb1, swalgo.sb2, mmalgo.displ, mmalgo.printPtr, _ptrToSeq1, _ptrToSeq2);
			    
			#ifdef DEBUG  		
      			cout << "len1: " << len1 << "; len2: " << len2 << endl;
      		#endif
 			if (len1 == 0 || len2 == 0) {
 				mmScore = 0;
      		} 
			else {
        		mmScore /= (float)utilityObject->MIN(len1, len2);
      		}

      		_score = ((float)100.0 - mmScore) / (float)100.0; 
      		#ifdef DEBUG
	      		cout << "mmScore: " << mmScore << endl;		               
	      		cout << "_score: " << _score << endl;		               
      		#endif
		   	
		   	distMat.push_back(distMatrixRecord(si+1, sj+1, _score));
		   	distMat.push_back(distMatrixRecord(sj+1, si+1, _score));
		   	// distMat->SetAt(si + 1, sj + 1, _score);
    		// distMat->SetAt(sj + 1, si + 1, _score);
	      //TODO: please uncomment me later!	
	      //   #pragma omp critical
	      //   {
	      // if(userParameters->getDisplayInfo())
	      // {
	      //     utilityObject->info("Sequences (%d:%d) Aligned. Score:  %d",
	        //                             si+1, sj+1, (int)mmScore);     
	      // }
	      //   }             
		}
	}    
	
	sendDistMat(&distMat);
	return;
}

void ParallelAlgo::sendDistMat(std::vector<distMatrixRecord>* distMat) {
	int size = distMat->size() * 3;
    MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    //temporal solution, should bde MPI_Type
    float* unwindedMat  = new float[size]; 
    for (int i=0; i<distMat->size(); i++) {
    	unwindedMat[i*3+0] = (*distMat)[i].row;
    	unwindedMat[i*3+1] = (*distMat)[i].col;
    	unwindedMat[i*3+2] = (*distMat)[i].val;
    }

    MPI_Send(unwindedMat, size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

    delete[] unwindedMat;
	return;
}

int ParallelAlgo::translateIndex(int i, int delta)
{
	if (i-delta < 0) {
		cout << "Error during index translation!!!!" << endl;
		return 0;	
	} 

	return i - delta;
}

void ParallelAlgo::recieveSequences()
{
	int bounds[4];
    MPI_Recv(&bounds, 4, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    iStart = bounds[0];
    iEnd = bounds[1];
    jStart = bounds[2];
    jEnd = bounds[3];

    initSi = utilityObject->MAX(0, iStart); //TODO
    int  boundSi = utilityObject->MIN(data.numSeqs,iEnd);
  
  	// [0] = [initSi+1]
    // [1] = [initSi+1+1]
    // [2] = [initSi+1+1+1]
    //recieve sequences
    const int NUMBER_OF_SEQ = data.numSeqs - initSi;

    for (int si=initSi; si<initSi+NUMBER_OF_SEQ; si++) {
        
        int size;
        MPI_Recv(&size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	 	#ifdef DEBUG
	        cout << "[ALGO] Size: "<< size << endl;
	    #endif

        std::vector<int> seq(size);
        int* seqPtr = seq.data();

        MPI_Recv(seqPtr, size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //TODO: maybe memory leak
        seqArray.push_back(seq);

    }
}

void ParallelAlgo::recieveExtendData() {
	//get ExtendData
	int intBuf[8];
	MPI_Recv(&intBuf, 8, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    data.intScale = intBuf[0];
	data.matAvgScore = intBuf[1];
	data.maxRes = intBuf[2];
	data.DNAFlag = (bool)intBuf[3];
	data.gapPos1 = intBuf[4];
	data.gapPos2 = intBuf[5];
	data.maxAlnLength = intBuf[6];
	data.numSeqs = intBuf[7];

	float floatBuf[4];
	MPI_Recv(&floatBuf, 4, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	data.gapOpenScale = floatBuf[0];
	data.gapExtendScale = floatBuf[1];
	data.pwGapOpen = floatBuf[2];
	data.pwGapExtend = floatBuf[3];
    
    #ifdef DEBUG
    	cout << "clustalw::NUMRES: " << clustalw::NUMRES << endl;
    #endif

	int* matrix = new int[clustalw::NUMRES*clustalw::NUMRES]; 
	MPI_Recv(matrix, clustalw::NUMRES*clustalw::NUMRES, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    for (int i=0; i<clustalw::NUMRES; i++) 
      for (int j=0; j<clustalw::NUMRES; j++) {
        data.matrix[i][j] = matrix[i*clustalw::NUMRES+j];
      }    

   	delete[] matrix;    

    #ifdef DEBUG
    	cout << "data.intScale:       " << data.intScale << endl;
		cout << "data.matAvgScore:    " << data.matAvgScore << endl;
		cout << "data.maxRes:         " << data.maxRes << endl;
		cout << "data.DNAFlag:        " << data.DNAFlag << endl;
		cout << "data.gapPos1:        " << data.gapPos1 << endl;
		cout << "data.gapPos2:        " << data.gapPos2 << endl;
		cout << "data.maxAlnLength:   " << data.maxAlnLength << endl;
		cout << "data.numSeqs:        " << data.numSeqs << endl;
		cout << "data.gapOpenScale:   " << data.gapOpenScale << endl;
		cout << "data.gapExtendScale: " << data.gapExtendScale << endl;
		cout << "data.pwGapOpen:      " << data.pwGapOpen << endl;
		cout << "data.pwGapExtend:    " << data.pwGapExtend << endl;

		cout << "data.matrix: " << endl;
		for (int i=0; i<clustalw::NUMRES; i++) {
			cout << data.matrix[clustalw::NUMRES-1][i] << ", ";
		}
		cout << endl;
    #endif

	return;
}

void ExternalData::UpdateGapOpenAndExtend(int& _gapOpen, int& _gapExtend, int n, int m)
{
	 if (DNAFlag)
   {
        _gapOpen = static_cast<int>(2 * pwGapOpen * intScale *
                        gapOpenScale);
        _gapExtend = static_cast<int>(pwGapExtend * intScale * gapExtendScale);
   }
   else
   {
        if (matAvgScore <= 0)
        {
            _gapOpen = 2 * static_cast<int>((pwGapOpen +
                   log(static_cast<double>(utilityObject->MIN(n, m)))) * intScale);
        }
        else
        {
            _gapOpen = static_cast<int>(2 * matAvgScore * (pwGapOpen +
            log(static_cast<double>(utilityObject->MIN(n, m)))) * gapOpenScale);
        }
        _gapExtend = static_cast<int>(pwGapExtend * intScale);
    }
}

float ParallelAlgo::tracePath(int tsb1, int tsb2, const vector<int>& displ, int printPtr,   const vector<int>* _ptrToSeq1,  const vector<int>* _ptrToSeq2)
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

          //maybe uncorrect change from userParam to exData
          if ((res1 != data.gapPos1) && 
              (res2 != data.gapPos2) && (res1 == res2))
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
