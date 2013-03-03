#ifndef __PARALLELALGO_H__
#define __PARALLELALGO_H__

#include "../alignment/Alignment.h"

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

struct distMatrixRecord {
	int row;
	int col;
	double val;

	distMatrixRecord(int r, int c, double v):
	row(r),
	col(c),
	val(v) {		
	}
};

class ParallelAlgo
{
	public:
		void DoFullPairwiseAlignment();	

	private:

		//MPI
		void recieveExtendData();
		void recieveSequences();
		void sendDistMat(std::vector<distMatrixRecord>*);

		//Helpers
		int translateIndex(int, int);
		float tracePath(int, int, const std::vector<int>&, int printPtr, const std::vector<int>*,  const std::vector<int>*);
		

		//data
		int initSi;
		static clustalw::SeqArray seqArray;
		static int iStart;
    	static int iEnd;
    	static int jStart;
    	static int jEnd;
    	static ExternalData data;

};


#endif
