#ifndef __EXTENDDATA_H__
#define __EXTENDDATA_H__

#include "PairwiseAlignBase.h"
#include <mpi.h>

namespace clustalw {

class ExtendData
{
	/*FIXME: it might be consts!*/
	public:
		//SubMatrixParameters
		static int matrix[clustalw::NUMRES][clustalw::NUMRES];
		static int intScale;
		static float gapOpenScale;
		static float gapExtendScale;
		static int matAvgScore;
		static int maxRes;
		//UserParameters
		static bool DNAFlag;
		static float pwGapOpen;
		static float pwGapExtend;
		static int gapPos1;
		static int gapPos2;
		//AlignmentParameters
		static int maxAlnLength;
		static int numSeqs;
		static MPI_Datatype mpi_dmRecord_type;
		//Funcs
		static void InitSubMatrixParameters(clustalw::SubMatrix* subMat);
		static void InitUserParameters(clustalw::UserParameters* userParameters);
		static void InitAlignmentParameters(clustalw::Alignment* alignPtr);
		static void UpdateGapOpenAndExtend(int&, int&, int, int);	
	
	private:
		ExtendData(const ExtendData&) {};
		ExtendData() {};		
};

typedef struct distMatrixRecord {
	int row;
	int col;
	float val;

	distMatrixRecord(int r, int c, float v):
	row(r),
	col(c),
	val(v) {		
	}

	distMatrixRecord():
	row(-1),
	col(-1),
	val(0) {

	}

	distMatrixRecord operator=(const distMatrixRecord& other) {
		return distMatrixRecord(other.row, other.col, other.val);
	}

} dmRecord; //TODO: should be better naming

}


#endif /* __EXTENDDATA_H__ */