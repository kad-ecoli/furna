const char* docstring=""
"cisbp2meme pwn.txt meme.txt\n"
"cisbp2meme pwn.txt meme.txt freq.txt\n";

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
using namespace std;

/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void Split(const string &line, vector<string> &line_vec,
    const char delimiter=' ')
{
    bool within_word = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (line[pos]==delimiter)
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

int main(int argc,char **argv)
{
    string pwn_file ="";
    string out_file ="";
    string freq_file="";

    for (int a=1;a<argc;a++)
    {
        if (pwn_file.size()==0)
            pwn_file=argv[a];
        else if (out_file.size()==0)
            out_file=argv[a];
        else if (freq_file.size()==0)
            freq_file=argv[a];
        else
        {
            cerr<<"ERROR: unknown option "<<argv[a]<<endl;
            return 1;
        }
    }

    if (out_file.size()==0)
    {
        cerr<<docstring;
        return 1;
    }

    string freq_txt="A 0.250 C 0.250 G 0.250 T 0.250";
    
    stringstream buf;
    if (freq_file.size())
    {
        ifstream fp;
        fp.open(freq_file.c_str(),ios::in);
        buf<<fp.rdbuf();
        fp.close();
        freq_txt=buf.str();
        buf.str(string());
    }
    
    string name = pwn_file.substr(pwn_file.find_last_of("/\\") + 1);
    name = name.substr(0, name.find_last_of('.'));


    ifstream fp;
    fp.open(pwn_file.c_str(),ios::in); //ifstream fp(filename,ios::in);
    buf<<fp.rdbuf();
    fp.close();

    vector<string> lines;
    Split(buf.str(),lines,'\n'); 
    buf.str(string());

    if (lines.size()<=1)
    {
        cerr<<"ERROR! "<<pwn_file<<" empty."<<endl;
        return 1;
    }
    buf<<lines.size()-1;
    string nsites=buf.str();
    buf.str(string());

    string result="MEME version 4\n"
        "\n"
        "ALPHABET= ACGT\n"
        "strands: +\n"
        "\n"
        "Background letter frequencies\n"+freq_txt+"\n\n"
        "MOTIF "+name+"\n"
        "letter-probability matrix: alength= 4 w= "+nsites+" nsites= "+nsites+"\n";

    vector<string> line_vec;
    vector<double> row(4,0);
    double total=0;
    size_t i,j;
    for (i=1;i<lines.size();i++)
    {
        Split(lines[i],line_vec,'\t');
        total=0;
        for (j=1;j<line_vec.size();j++)
        {
            row[j-1]=atof(line_vec[j].c_str());
            total+=row[j-1];
        }
        buf<<row[0]/total<<"  "
           <<row[1]/total<<"  "
           <<row[2]/total<<"  "
           <<row[3]/total<<endl;
        for (j=0;j<line_vec.size();j++) line_vec[j].clear(); line_vec.clear();
    }
    result+=buf.str();

    if (out_file=="-")
        cout<<result<<flush;
    else
    {
        ofstream fout;
        fout.open(out_file.c_str(),ios::out);
        fout<<result<<flush;
        fout.close();
    }

    /* clean up */
    result.clear();
    buf.str(string());
    name.clear();
    vector<string>().swap(line_vec);
    vector<double>().swap(row);
    return 0;
}
