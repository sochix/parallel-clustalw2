#include "I_ExtendData.h"

using namespace clustalw;

//SubMatrixParameters
int 	ExtendData::matrix[clustalw::NUMRES][clustalw::NUMRES];
int 	ExtendData::intScale;
int 	ExtendData::matAvgScore;
int 	ExtendData::maxRes;
float ExtendData::gapOpenScale;
float ExtendData::gapExtendScale;
//UserParameters
bool	ExtendData::DNAFlag;
int 	ExtendData::gapPos1;
int 	ExtendData::gapPos2;
float ExtendData::pwGapOpen;
float ExtendData::pwGapExtend;
//AlignmentParameters
int 	ExtendData::maxAlnLength;
int		ExtendData::numSeqs;

void ExtendData::InitSubMatrixParameters(SubMatrix* subMat)
{
		PairScaleValues scaleValues;
				
	  maxRes = subMatrix->getPairwiseMatrix(matrix, scaleValues, matAvgScore);
	
	  intScale = scaleValues.intScale;
    gapOpenScale = scaleValues.gapOpenScale;
    gapExtendScale = scaleValues.gapExtendScale;
}

void ExtendData::InitUserParameters(clustalw::UserParameters* userParameters)
{
		DNAFlag = userParameters->getDNAFlag();
		pwGapOpen = userParameters->getPWGapOpen();
    pwGapExtend = userParameters->getPWGapExtend();
		gapPos1 = userParameters->getGapPos1();
		gapPos2 = userParameters->getGapPos2();    
}

void ExtendData::UpdateGapOpenAndExtend(int& _gapOpen, int& _gapExtend, int n, int m)
{
	 if (ExtendData::DNAFlag)
   {
        _gapOpen = static_cast<int>(2 * ExtendData::pwGapOpen * ExtendData::intScale *
                        ExtendData::gapOpenScale);
        _gapExtend = static_cast<int>(ExtendData::pwGapExtend * ExtendData::intScale * ExtendData::gapExtendScale);
   }
   else
   {
        if (ExtendData::matAvgScore <= 0)
        {
            _gapOpen = 2 * static_cast<int>((ExtendData::pwGapOpen +
                   log(static_cast<double>(utilityObject->MIN(n, m)))) * ExtendData::intScale);
        }
        else
        {
            _gapOpen = static_cast<int>(2 * ExtendData::matAvgScore * (ExtendData::pwGapOpen +
            log(static_cast<double>(utilityObject->MIN(n, m)))) * ExtendData::gapOpenScale);
        }
        _gapExtend = static_cast<int>(ExtendData::pwGapExtend * ExtendData::intScale);
    }
}

void ExtendData::InitAlignmentParameters(clustalw::Alignment* alignPtr)
{
	maxAlnLength = alignPtr->getMaxAlnLength();
	numSeqs = alignPtr->getNumSeqs();
	
}
