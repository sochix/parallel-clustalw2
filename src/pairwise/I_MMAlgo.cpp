#include "I_MMAlgo.h"

using namespace clustalw;

MMAlgo::MMAlgo(ExternalData* d): 
    printPtr(0),
  lastPrint(0),
  _ptrToSeq1(NULL),
  _ptrToSeq2(NULL),
  _gapOpen(0),
  _gapExtend(0)
{
  data = d;
  displ.resize((2 * data->maxAlnLength) + 1);
  HH.resize(data->maxAlnLength);
  DD.resize(data->maxAlnLength);
  RR.resize(data->maxAlnLength);
  SS.resize(data->maxAlnLength);  
}

MMAlgo::~MMAlgo()
{
    displ.clear();
    HH.clear();
    DD.clear();
    RR.clear();
    SS.clear();
}

int MMAlgo::Pass(int sb1, int sb2, int len1, int len2, int tb, int te, const vector<int>* seq1, const vector<int>* seq2, const int gapOpen, const int gapExtend)
{
		_ptrToSeq1 = seq1;
		_ptrToSeq2 = seq2;
		_gapExtend = gapExtend;
		_gapOpen = gapOpen;
		
		lastPrint = 0;
		printPtr = 1;
		
		return diff(sb1,sb2,len1,len2,tb,te);
		
}

inline int MMAlgo::calcScore(int i, int j, int sb1, int sb2)
{
    return data->matrix[(*_ptrToSeq1)[sb1 + i]][(*_ptrToSeq2)[sb2 + j]];
}

void MMAlgo::addToDisplay(int v)
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

int MMAlgo::diff(int sb1, int sb2, int len1, int len2, int tb, int te)
{
    int midi, midj, i, j;
    int midh;
    int f, hh, e, s, t; 

    if (len2 <= 0)
    {
        if (len1 > 0)
        {
            delFromDisplay(len1);
        }

        return ( -(int)tbgap(len1, tb));
    }

    if (len1 <= 1)
    {
        if (len1 <= 0)
        {
            addToDisplay(len2);
            return ( -(int)tbgap(len2, tb));
        }

        midh =  - (tb + _gapExtend) - tegap(len2, te);
        hh =  - (te + _gapExtend) - tbgap(len2, tb);
        if (hh > midh)
        {
            midh = hh;
        }
        midj = 0;
        for (j = 1; j <= len2; j++)
        {
            hh = calcScore(1, j, sb1, sb2) - tegap(len2 - j, te) - tbgap(j - 1, tb);
            if (hh > midh)
            {
                midh = hh;
                midj = j;
            }
        }

        if (midj == 0)
        {
            delFromDisplay(1);
            addToDisplay(len2);
        }
        else
        {
            if (midj > 1)
            {
                addToDisplay(midj - 1);
            }
            displ[printPtr++] = lastPrint = 0;
            if (midj < len2)
            {
                addToDisplay(len2 - midj);
            }
        }
        return midh;
    }

    // Divide: Find optimum midpoint (midi,midj) of cost midh
    midi = len1 / 2;
    forwardPass(sb1, sb2, len1, len2, tb, midi);
    backwardPass(sb1, sb2, len1, len2, te, midi);

    midh = HH[0] + RR[0];
    midj = 0;
    type_t type = TYPE_1;
    for (j = 0; j <= len2; j++)
    {
        hh = HH[j] + RR[j];
        if (hh >= midh)
        if (hh > midh || (HH[j] != DD[j] && RR[j] == SS[j]))
        {
            midh = hh;
            midj = j;
        }
    }

    for (j = len2; j >= 0; j--)
    {
        hh = DD[j] + SS[j] + _gapOpen;
        if (hh > midh)
        {
            midh = hh;
            midj = j;
            type = TYPE_2;
        }
    }

    // Conquer recursively around midpoint
    if (type == TYPE_1)
    {
        // Type 1 gaps
        diff(sb1, sb2, midi, midj, tb, _gapOpen);
        diff(sb1 + midi, sb2 + midj, len1 - midi, len2 - midj, _gapOpen, te);
    }
    else
    {
        diff(sb1, sb2, midi - 1, midj, tb, 0);
        delFromDisplay(2);
        diff(sb1 + midi + 1, sb2 + midj, len1 - midi - 1, len2 - midj, 0, te);
    }

    return midh; // Return the score of the best alignment
}

void MMAlgo::delFromDisplay(int k)
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

int MMAlgo::gapAffineFunction(int k, int t) 
{
    if(k <= 0)
    {
        return 0;
    }
    else
    {
        return t + _gapExtend * k;
    }
}

int MMAlgo::tbgap(int k, int tb)
{
    return gapAffineFunction(k, tb);
}

int MMAlgo::tegap(int k, int te)
{
    return gapAffineFunction(k, te);
}

void MMAlgo::forwardPass(int sb1, int sb2, int len1, int len2, int tb, int midi) 
{
    int t = -tb, 
        s, hh, f, e;

    HH[0] = 0;
    
    for (int j = 1; j <= len2; j++)
    {
        HH[j] = t = t - _gapExtend;
        DD[j] = t - _gapOpen;
    }

    t =  - tb;
    for (int i = 1; i <= midi; i++)
    {
        s = HH[0];
        HH[0] = hh = t = t - _gapExtend;
        f = t - _gapOpen;
        for (int j = 1; j <= len2; j++)
        {
            if ((hh = hh - _gapOpen - _gapExtend) > (f = f - _gapExtend))
            {
                f = hh;
            }
            if ((hh = HH[j] - _gapOpen - _gapExtend) > (e = DD[j] - _gapExtend))
            {
                e = hh;
            }
            hh = s + calcScore(i, j, sb1, sb2);
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
    return;
}

void MMAlgo::backwardPass(int sb1, int sb2, int len1, int len2, int te, int midi) 
{
    int t = -te, 
        hh, f, e, s;

    RR[len2] = 0;

    for (int j = len2 - 1; j >= 0; j--)
    {
        RR[j] = t = t - _gapExtend;
        SS[j] = t - _gapOpen;
    }

    t =  - te;
    for (int i = len1 - 1; i >= midi; i--)
    {
        s = RR[len2];
        RR[len2] = hh = t = t - _gapExtend;
        f = t - _gapOpen;

        for (int j = len2 - 1; j >= 0; j--)
        {

            if ((hh = hh - _gapOpen - _gapExtend) > (f = f - _gapExtend))
            {
                f = hh;
            }
            if ((hh = RR[j] - _gapOpen - _gapExtend) > (e = SS[j] - _gapExtend))
            {
                e = hh;
            }
            hh = s + calcScore(i + 1, j + 1, sb1, sb2);
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

    SS[len2] = RR[len2];
    return;
}