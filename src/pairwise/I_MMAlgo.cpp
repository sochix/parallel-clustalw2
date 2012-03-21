#include "I_MMAlgo.h"

using namespace clustalw;

MMAlgo::MMAlgo(): 
	printPtr(0),
  lastPrint(0),
  _ptrToSeq1(NULL),
  _ptrToSeq2(NULL),
  _gapOpen(0),
  _gapExtend(0)
{
	displ.resize((2 * ExtendData::maxAlnLength) + 1);
  HH.resize(ExtendData::maxAlnLength);
  DD.resize(ExtendData::maxAlnLength);
  RR.resize(ExtendData::maxAlnLength);
  SS.resize(ExtendData::maxAlnLength);  
}

MMAlgo::~MMAlgo()
{
	  displ.clear();
    HH.clear();
    DD.clear();
    RR.clear();
    SS.clear();
}

int MMAlgo::Pass(int A, int B, int M, int N, int tb, int te, const vector<int>* seq1, const vector<int>* seq2, const int gapOpen, const int gapExtend)
{
		_ptrToSeq1 = seq1;
		_ptrToSeq2 = seq2;
		_gapExtend = gapExtend;
		_gapOpen = gapOpen;
		
		lastPrint = 0;
		printPtr = 1;
		
		return diff(A,B,M,N,tb,te);
		
}

inline int MMAlgo::calcScore(int iat, int jat, int v1, int v2)
{
    return ExtendData::matrix[(*_ptrToSeq1)[v1 + iat]][(*_ptrToSeq2)[v2 + jat]];
}

void MMAlgo::add(int v)
{
    if (lastPrint < 0)
    {
        displ[printPtr - 1] = v;
        displ[printPtr++] = lastPrint;
    }
    else
    {
        lastPrint = displ[printPtr++] = v;
    }
}

int MMAlgo::diff(int A, int B, int M, int N, int tb, int te)
{
    int type;
    int midi, midj, i, j;
    int midh;
    static int f, hh, e, s, t;

    if (N <= 0)
    {
        if (M > 0)
        {
            del(M);
        }

        return ( -(int)tbgap(M, tb));
    }

    if (M <= 1)
    {
        if (M <= 0)
        {
            add(N);
            return ( -(int)tbgap(N, tb));
        }

        midh =  - (tb + _gapExtend) - tegap(N, te);
        hh =  - (te + _gapExtend) - tbgap(N, tb);
        if (hh > midh)
        {
            midh = hh;
        }
        midj = 0;
        for (j = 1; j <= N; j++)
        {
            hh = calcScore(1, j, A, B) - tegap(N - j, te) - tbgap(j - 1, tb);
            if (hh > midh)
            {
                midh = hh;
                midj = j;
            }
        }

        if (midj == 0)
        {
            del(1);
            add(N);
        }
        else
        {
            if (midj > 1)
            {
                add(midj - 1);
            }
            displ[printPtr++] = lastPrint = 0;
            if (midj < N)
            {
                add(N - midj);
            }
        }
        return midh;
    }

    // Divide: Find optimum midpoint (midi,midj) of cost midh

    midi = M / 2;
    HH[0] = 0;
    t =  - tb;
    for (j = 1; j <= N; j++)
    {
        HH[j] = t = t - _gapExtend;
        DD[j] = t - _gapOpen;
    }

    t =  - tb;
    for (i = 1; i <= midi; i++)
    {
        s = HH[0];
        HH[0] = hh = t = t - _gapExtend;
        f = t - _gapOpen;
        for (j = 1; j <= N; j++)
        {
            if ((hh = hh - _gapOpen - _gapExtend) > (f = f - _gapExtend))
            {
                f = hh;
            }
            if ((hh = HH[j] - _gapOpen - _gapExtend) > (e = DD[j] - _gapExtend))
            {
                e = hh;
            }
            hh = s + calcScore(i, j, A, B);
            if (f > hh)
            {
                hh = f;
            }
            if (e > hh)
            {
                hh = e;
            }

            s = HH[j];
            HH[j] = hh;
            DD[j] = e;
        }
    }

    DD[0] = HH[0];

    RR[N] = 0;
    t =  - te;
    for (j = N - 1; j >= 0; j--)
    {
        RR[j] = t = t - _gapExtend;
        SS[j] = t - _gapOpen;
    }

    t =  - te;
    for (i = M - 1; i >= midi; i--)
    {
        s = RR[N];
        RR[N] = hh = t = t - _gapExtend;
        f = t - _gapOpen;

        for (j = N - 1; j >= 0; j--)
        {

            if ((hh = hh - _gapOpen - _gapExtend) > (f = f - _gapExtend))
            {
                f = hh;
            }
            if ((hh = RR[j] - _gapOpen - _gapExtend) > (e = SS[j] - _gapExtend))
            {
                e = hh;
            }
            hh = s + calcScore(i + 1, j + 1, A, B);
            if (f > hh)
            {
                hh = f;
            }
            if (e > hh)
            {
                hh = e;
            }

            s = RR[j];
            RR[j] = hh;
            SS[j] = e;

        }
    }

    SS[N] = RR[N];

    midh = HH[0] + RR[0];
    midj = 0;
    type = 1;
    for (j = 0; j <= N; j++)
    {
        hh = HH[j] + RR[j];
        if (hh >= midh)
        if (hh > midh || (HH[j] != DD[j] && RR[j] == SS[j]))
        {
            midh = hh;
            midj = j;
        }
    }

    for (j = N; j >= 0; j--)
    {
        hh = DD[j] + SS[j] + _gapOpen;
        if (hh > midh)
        {
            midh = hh;
            midj = j;
            type = 2;
        }
    }

    // Conquer recursively around midpoint 


    if (type == 1)
    {
        // Type 1 gaps
        diff(A, B, midi, midj, tb, _gapOpen);
        diff(A + midi, B + midj, M - midi, N - midj, _gapOpen, te);
    }
    else
    {
        diff(A, B, midi - 1, midj, tb, 0);
        del(2);
        diff(A + midi + 1, B + midj, M - midi - 1, N - midj, 0, te);
    }

    return midh; // Return the score of the best alignment
}

void MMAlgo::del(int k)
{
    if (lastPrint < 0)
    {
        lastPrint = displ[printPtr - 1] -= k;
    }
    else
    {
        lastPrint = displ[printPtr++] =  - (k);
    }
}

int MMAlgo::tbgap(int k, int tb)
{
    if(k <= 0)
    {
        return 0;
    }
    else
    {
        return tb + _gapExtend * k;
    }
}

int MMAlgo::tegap(int k, int te)
{
    if(k <= 0)
    {
        return 0;
    }
    else
    {
        return te + _gapExtend * k;
    }
}