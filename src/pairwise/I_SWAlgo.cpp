#include "I_SWAlgo.h"
#include <iostream>

using namespace clustalw;
using namespace std;

SWAlgo::SWAlgo(ExternalData* d, const int* m):
	sb1(0),
	sb2(0),
	se1(0),
	se2(0),
	maxScore(0)
{
    data = d;
    matrix = m;
    HH.resize(data->maxAlnLength);
    DD.resize(data->maxAlnLength);
}

SWAlgo::~SWAlgo()
{
	HH.clear();
	DD.clear();
}

void SWAlgo::Pass(const vector<int>* seq1, const vector<int>* seq2, int n, int m, const int gapOpen, const int gapExtend)
{
    // cout << "Seq1: " <<endl;
    // for (int i=0; i<seq1->size(); i++)
    //     cout << (*seq1)[i];
    // cout << endl <<  "Seq2: " <<endl;
    // for (int i=0; i<seq2->size(); i++)
    //     cout << (*seq2)[i];
    // cout << endl;


	forwardPass(seq1, seq2, n,m, gapOpen, gapExtend);
	reversePass(seq1, seq2, gapOpen, gapExtend);
}

void SWAlgo::forwardPass(const vector<int>* seq1, const vector<int>* seq2, int n, int m, const int _gapOpen, const int _gapExtend)
{
    int i, j;
    int f, hh, p, t;

    maxScore = 0;
    se1 = se2 = 0;
    for (i = 0; i <= m; i++)
    {
        HH[i] = 0;
        DD[i] =  -_gapOpen;
    }

    for (i = 1; i <= n; i++)
    {
        hh = p = 0;
        f =  -_gapOpen;

        for (j = 1; j <= m; j++)
        {

            f -= _gapExtend;
            t = hh - _gapOpen - _gapExtend;
            if (f < t)
            {
                f = t;
            }

            DD[j] -= _gapExtend;
            t = HH[j] - _gapOpen - _gapExtend;
            if (DD[j] < t)
            {
                DD[j] = t;
            }

            int index = (*seq1)[i]*clustalw::NUMRES+(*seq2)[j];
            //hh  = p +index;
            hh = p + matrix[index];//data->matrix[(*seq1)[i]][(*seq2)[j]];
            
            if (hh < f)
            {
                hh = f;
            }
            if (hh < DD[j])
            {
                hh = DD[j];
            }
            if (hh < 0)
            {
                hh = 0;
            }

            p = HH[j];
            HH[j] = hh;

            if (hh > maxScore)
            {
                maxScore = hh;
                se1 = i;
                se2 = j;
            }
        }
    }

}

void SWAlgo::reversePass(const vector<int>* seq1, const vector<int>* seq2, const int _gapOpen, const int _gapExtend)
{
    int i, j;
    int f, hh, p, t;
    int cost;

    cost = 0;
    sb1 = sb2 = 1;
    for (i = se2; i > 0; i--)
    {
        HH[i] =  - 1;
        DD[i] =  - 1;
    }

    for (i = se1; i > 0; i--)
    {
        hh = f =  - 1;
        if (i == se1)
        {
            p = 0;
        }
        else
        {
            p =  - 1;
        }

        for (j = se2; j > 0; j--)
        {

            f -= _gapExtend;
            t = hh - _gapOpen - _gapExtend;
            if (f < t)
            {
                f = t;
            }

            DD[j] -= _gapExtend;
            t = HH[j] - _gapOpen - _gapExtend;
            if (DD[j] < t)
            {
                DD[j] = t;
            }

            hh = p + data->matrix[(*seq1)[i]][(*seq2)[j]];
            if (hh < f)
            {
                hh = f;
            }
            if (hh < DD[j])
            {
                hh = DD[j];
            }

            p = HH[j];
            HH[j] = hh;

            if (hh > cost)
            {
                cost = hh;
                sb1 = i;
                sb2 = j;
                if (cost >= maxScore)
                {
                    break;
                }
            }
        }
        if (cost >= maxScore)
        {
            break;
        }
    }

}
