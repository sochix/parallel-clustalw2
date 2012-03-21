#ifndef __EXTENDDATA_H__
#define __EXTENDDATA_H__

#include "PairwiseAlignBase.h"

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
		//Funcs
		static void InitSubMatrixParameters(clustalw::SubMatrix* subMat);
		static void InitUserParameters(clustalw::UserParameters* userParameters);
		static void InitAlignmentParameters(clustalw::Alignment* alignPtr);
		static void UpdateGapOpenAndExtend(int&, int&, int, int);
	
	private:
		ExtendData(const ExtendData&) {};
		ExtendData() {};		
};

}

#endif /* __EXTENDDATA_H__ */