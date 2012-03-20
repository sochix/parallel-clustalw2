/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef FULLPAIRWISEALIGN_H
#define FULLPAIRWISEALIGN_H

#include "PairwiseAlignBase.h"
//include "I_ExtendData.h"

namespace clustalw
{


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

class SWAlgo
{
	public:
		SWAlgo();
		~SWAlgo();
		void Pass(const vector<int>* seq1, const vector<int>* seq2, int n, int m, const int, const int);
		
		/*FIXME: it might be consts!*/
		int maxScore;
		int se1;
		int se2;
		int sb1;
		int sb2;
		
	private:
		void forwardPass(const vector<int>* seq1, const vector<int>* seq2, int n, int m, const int, const int);
		void reversePass(const vector<int>* ia, const vector<int>* ib, const int, const int);
		
		/*
		int maxScore;
		int se1;
		int se2;
		int sb1;
		int sb2;
		*/
		vector<int> HH;
    vector<int> DD;
 
		
};


class FullPairwiseAlign : public PairwiseAlignBase
{
    public:
        /* Functions */
        FullPairwiseAlign();
	virtual ~FullPairwiseAlign(){};

        virtual void pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, 
                                   int iEnd, int jStart, int jEnd); 
        /* Attributes */

    private:
        /* Functions */
				//Algo SW
//				void forwardPass(const vector<int>* seq1, const vector<int>* seq2, int n, int m);
//        void reversePass(const vector<int>* ia, const vector<int>* ib);
        //Algo MM
        int diff(int A, int B, int M, int N, int tb, int te);
        void add(int v);
        void del(int k);
        int calcScore(int iat, int jat, int v1, int v2); 
        int tbgap(int k, int tb);
        int tegap(int k, int te);
        //Other stuff
        float tracePath(int tsb1, int tsb2);       
        int gap(int k);
        
        /* Attributes */
        // I have constant pointers to the data. This allows for the fastest access.
        const vector<int>* _ptrToSeq1;
        const vector<int>* _ptrToSeq2;
        float mmScore;
        int printPtr;
        int lastPrint;
        vector<int> displ;
        vector<int> HH;
        vector<int> DD;
        vector<int> RR;
        vector<int> SS;

        int _gapOpen; // scaled to be an integer, this is not a mistake
        int _gapExtend; // scaled to be an integer, not a mistake
     
        int seq1;
        int seq2;
     
        int maxScore;
     /* 
        int sb1;
        int sb2;
        int se1;
        int se2;
      */

};

}
#endif
