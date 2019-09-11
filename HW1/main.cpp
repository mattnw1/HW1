#include <iostream>
#include <fstream> //read and write to files
#include <stdio.h>

using namespace std;


int main(int argc, char const *argv[])
{
  int A = 0;
  int C = 0;
  int T = 0;
  int G = 0;
  int AA = 0;
  int AC = 0;
  int AT = 0;
  int AG = 0;
  int CA = 0;
  int CC = 0;
  int CT = 0;
  int CG = 0;
  int TA = 0;
  int TT = 0;
  int TC = 0;
  int TG = 0;
  int GA = 0;
  int GG = 0;
  int GC = 0;
  int GT = 0;
  ifstream dnafile;
  dnafile.open("DNA_Strings.txt");
  int numOfLines = 0;
  string line;
  int dnaSum, dnaMean, dnaVariance, dnaStdDev;

  if (dnafile.is_open())
  {
    while (getline(dnafile,line))
    {

    //std::cout << "DNA file" << '\n';

      // Count single nucleotides
      for (int i = 0; i < line.size(); ++i)
      {
        // using transform() function and ::toupper in STL
        transform(line.begin(), line.end(), line.begin(), ::toupper); // ensures every DNA seq is uppcase
        to_string(line[i]);

        if(line[i] == 'A')
        {
          ++A;
          //cout << A << endl;
        }

        else if(line[i] == 'C')
        ++C;

        else if(line[i] == 'T')
        ++T;

        else if(line[i] == 'G')
        ++G;
      }

      //cout << A << endl;

      // Count bigram nucleotides
      for (int i = 0; i < line.size(); i+=2)
      {
        // using transform() function and ::toupper in STL
        transform(line.begin(), line.end(), line.begin(), ::toupper); // ensures every DNA seq is uppcase
        to_string(line[i]);

        cout << line[i] << "test" << endl;


      }
    std::cout << line << '\n';
    cout << "Total Number of Gs: " << G <<endl;

    ++numOfLines;


    }
      //compute sum of length

      //compute mean of length

      //compute variance of length

      //compute standard deviation of length

  }
    /* code */
    dnafile.close();
    return 0;
  }
