#ifndef __PARALLELALGO_H__
#define __PARALLELALGO_H__

#include "../alignment/Alignment.h"
#include "../pairwise/I_ExtendData.h"

struct ExternalData {
		//SubMatrixParameters
	int matrix[clustalw::NUMRES][clustalw::NUMRES];
	int intScale;
	float gapOpenScale;
	float gapExtendScale;
	int matAvgScore;
	int maxRes;
	//UserParameters
	bool DNAFlag;
	float pwGapOpen;
	float pwGapExtend;
	int gapPos1;
	int gapPos2;
	//AlignmentParameters
	int maxAlnLength;
	int numSeqs;
	//Funcs
	void UpdateGapOpenAndExtend(int&, int&, int, int);	
};

class ParallelAlgo
{
	public:
		void DoFullPairwiseAlignment();	

	private:

		//MPI
		void recieveExtendData();
		void recieveSequences();
		void sendDistMat(std::vector<clustalw::dmRecord>*);

		//Helpers
		int translateIndex(int, int);
		float tracePath(int, int, const std::vector<int>&, int printPtr, const std::vector<int>*,  const std::vector<int>*);
		

		//data
		static clustalw::SeqArray seqArray;
		static int iStart;
    	static int iEnd;
    	static int jStart;
    	static int jEnd;
    	static ExternalData data;
    	static int* portionPerProc;
		//TODO: should be const
    	static int maxSeqCount;
    	

};


#endif
