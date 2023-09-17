const char* docstring=""
"fasta2freq input.fasta freq.txt\n"
"\n"
"Input:\n"
"    input.fasta  - input fasta\n"
"\n"
"Output:\n"
"    freq.txt     - frequency of ACGT\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void adjust_ACGT(vector<size_t>&ACGT_count)
{
    size_t ACGT_sum=ACGT_count[0]+ACGT_count[1]+ACGT_count[2]+ACGT_count[3];
    size_t i;
    while (ACGT_sum>10000)
    {
        size_t max_i=5;
        size_t max_val=0;
        for (i=0;i<4;i++)
        {
            if (max_val<ACGT_count[i])
            {
                max_val=ACGT_count[i];
                max_i=i;
            }
        }
        if (max_i<=4) ACGT_count[max_i]--;
        ACGT_sum=ACGT_count[0]+ACGT_count[1]+ACGT_count[2]+ACGT_count[3];
    }
    while (ACGT_sum<10000)
    {
        size_t min_i=5;
        size_t min_val=100000;
        for (i=0;i<4;i++)
        {
            if (min_val>ACGT_count[i])
            {
                min_val=ACGT_count[i];
                min_i=i;
            }
        }
        if (min_i<=4) ACGT_count[min_i]++;
        ACGT_sum=ACGT_count[0]+ACGT_count[1]+ACGT_count[2]+ACGT_count[3];
    }
}

void fasta2freq(const string &inputFastaFile,
    const string &outputFreqFile)
{
    ifstream fp_in;
    ofstream fp_out;
    string line;

    vector<size_t>ACGT_count(4,0);
    
    /* read fasta file */
    size_t i;
    if (inputFastaFile!="-") fp_in.open(inputFastaFile.c_str());
    while ((inputFastaFile!="-")?fp_in.good():cin.good())
    {
        if (inputFastaFile!="-") getline(fp_in,line);
        else getline(cin,line);
        if (line.size()==0||line[0]=='>') continue;
        for (i=0;i<line.size();i++)
        {
            switch(line[i])
            {
                case 'A':
                    ACGT_count[0]++;
                    break;
                case 'C':
                    ACGT_count[1]++;
                    break;
                case 'G':
                    ACGT_count[2]++;
                    break;
                case 'T':
                    ACGT_count[3]++;
                    break;
                case 'U':
                    ACGT_count[3]++;
                    break;
            }
        }
        line.clear();
    }
    if (inputFastaFile!="-") fp_in.close();

    /* write freq file */
    double total=0;
    for (i=0;i<4;i++) total+=ACGT_count[i];
    for (i=0;i<4;i++) ACGT_count[i]=ACGT_count[i]*10000/total;
    adjust_ACGT(ACGT_count);

    if (outputFreqFile=="-") cout<<setiosflags(ios::fixed)
            <<"A "<<setprecision(4)<<0.0001*ACGT_count[0]
            <<" C "<<setprecision(4)<<0.0001*ACGT_count[1]
            <<" G "<<setprecision(4)<<0.0001*ACGT_count[2]
            <<" T "<<setprecision(4)<<0.0001*ACGT_count[3]<<endl;
    else 
    {
        fp_out.open(outputFreqFile.c_str());
        fp_out<<setiosflags(ios::fixed)
            <<"A "<<setprecision(4)<<0.0001*ACGT_count[0]
            <<" C "<<setprecision(4)<<0.0001*ACGT_count[1]
            <<" G "<<setprecision(4)<<0.0001*ACGT_count[2]
            <<" T "<<setprecision(4)<<0.0001*ACGT_count[3]<<endl;
        fp_out.close();
    }
    
    ACGT_count.clear();
    return;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc!=3)
    {
        cerr<<docstring;
        return 0;
    }
    string inputFastaFile=argv[1];
    string outputFreqFile=argv[2];
    fasta2freq(inputFastaFile,outputFreqFile);
    return 0;
}
