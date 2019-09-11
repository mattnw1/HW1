#include <iostream>
#include <fstream> //read and write to files
#include <stdio.h>
#include <string>
#include <cmath>



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
  int varianceA;
  int varianceC;
  int varianceT;
  int varianceG;
  double relativeProbA;
  double relativeProbC;
  double relativeProbT;
  double relativeProbG;
  double relativeProbAA;
  double relativeProbAC;
  double relativeProbAT;
  double relativeProbAG;
  double relativeProbCA;
  double relativeProbCC;
  double relativeProbCT;
  double relativeProbCG;
  double relativeProbTA;
  double relativeProbTT;
  double relativeProbTC;
  double relativeProbTG;
  double relativeProbGA;
  double relativeProbGG;
  double relativeProbGC;
  double relativeProbGT;



  ifstream dnafile;
  dnafile.open("DNA_Strings.txt");
  int numOfLines = 0;
  string line;
  int dnaSum, dnaMean;
  int dnaVariance = 0;
  int dnaStdDev;

  if (dnafile.is_open())
  {
    while (getline(dnafile,line))
    {

  //  std::cout << "DNA file" << '\n';

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
          //std::cout << "test" << '\n';
        }

        else if(line[i] == 'C')
        ++C;

        else if(line[i] == 'T')
        ++T;

        else if(line[i] == 'G')
        ++G;
      }

    //  cout << A << endl;

      // Count bigram nucleotides
      for (int i = 0; i < line.size(); i+=2)
      {
        // using transform() function and ::toupper in STL
        transform(line.begin(), line.end(), line.begin(), ::toupper); // ensures every DNA seq is uppcase
        to_string(line[i]);
        string bigram = line.substr(i,2);

        if (bigram == "AA")
          ++AA;

        else if (bigram == "AC")
          ++AC;

        else if (bigram == "AT")
          ++AT;

        else if (bigram == "AG")
          ++AG;

        else if (bigram == "CA")
          ++CA;

        else if (bigram == "CC")
          ++CC;

        else if (bigram == "CT")
          ++CT;

        else if (bigram == "CG")
          ++CG;

        else if (bigram == "TA")
          ++TA;

        else if (bigram == "TT")
          ++TT;

        else if (bigram == "TC")
          ++TC;

        else if (bigram == "TG")
          ++TG;

        else if (bigram == "GA")
          ++GA;

        else if (bigram == "GG")
          ++GG;

        else if (bigram == "GC")
          ++GC;

        else if (bigram == "GT")
          ++GT;
        // cout << bigram << " test"<< endl;
      }
      /*
      cout << AA << " AAs" << endl;
      cout << AC << " ACs" << endl;
      cout << AT << " ATs" << endl;
      cout << AG << " AGs" << endl;
      cout << CA << " CAs" << endl;
      cout << CC << " CCs" << endl;
      cout << CT << " CTs" << endl;
      cout << CG << " CGs" << endl;
      cout << TA << " TAs" << endl;
      cout << TT << " TTs" << endl;
      cout << TC << " TCs" << endl;
      cout << TG << " TGs" << endl;
      cout << GA << " GAs" << endl;
      cout << GG << " GGs" << endl;
      cout << GC << " GCs" << endl;
      cout << GT << " GTs" << endl;
      */
    //std::cout << line << '\n';
    //cout << "Total Number of Gs: " << G <<endl;
    ++numOfLines;
    }
    //compute sum of length
    dnaSum = (A+T+G+C);
    //cout << dnaSum << endl;

    //compute mean of length
    dnaMean = (A+T+G+C)/numOfLines;
    //cout << dnaMean << endl;


    //compute variance of length
    varianceA = A-dnaMean;
    varianceC = C-dnaMean;
    varianceT = T-dnaMean;
    varianceG = G-dnaMean;
    //cout << varianceA << "varAtest" << endl;

    dnaVariance += pow(varianceA, 2);
    //cout << dnaVariance << "test1" <<endl;

    dnaVariance += pow(varianceC, 2);
    //cout << dnaVariance << "test2" <<endl;

    dnaVariance += pow(varianceT, 2);
    //cout << dnaVariance << "test3" <<endl;

    dnaVariance += pow(varianceG, 2);
    //cout << dnaVariance << "test4" <<endl;

    dnaVariance=dnaVariance/4;
    //cout << dnaVariance << "test5" <<endl;

    //compute standard deviation of length
    dnaStdDev = sqrt(dnaVariance);
    //cout << dnaStdDev << "stdDev" <<endl;

    //compute relative probabilities for nucleotides
    relativeProbA = (double)A/dnaSum;
    relativeProbC = (double)C/dnaSum;
    relativeProbG = (double)G/dnaSum;
    relativeProbT = (double)T/dnaSum;
    //cout << relativeProbA << " relativeProbA" << endl;
    //cout << relativeProbC << " relativeProbA" << endl;
    //cout << relativeProbG << " relativeProbA" << endl;
    //cout << relativeProbT << " relativeProbA" << endl;
    //cout << relativeProbA + relativeProbC + relativeProbG + relativeProbT << " totalProb" << endl;

    //compute relative probabilities for bigrams
    relativeProbAA = (double)AA/(dnaSum/2);
    relativeProbAC = (double)AC/(dnaSum/2);
    relativeProbAT = (double)AT/(dnaSum/2);
    relativeProbAG = (double)AG/(dnaSum/2);
    relativeProbCA = (double)CA/(dnaSum/2);
    relativeProbCC = (double)CC/(dnaSum/2);
    relativeProbCT = (double)CT/(dnaSum/2);
    relativeProbCG = (double)CG/(dnaSum/2);
    relativeProbTA = (double)TA/(dnaSum/2);
    relativeProbTT = (double)TT/(dnaSum/2);
    relativeProbTC = (double)TC/(dnaSum/2);
    relativeProbTG = (double)TG/(dnaSum/2);
    relativeProbGA = (double)GA/(dnaSum/2);
    relativeProbGG = (double)GG/(dnaSum/2);
    relativeProbGC = (double)GC/(dnaSum/2);
    relativeProbGT = (double)GT/(dnaSum/2);

    //cout << relativeProbAA << " relativeProbAA" << endl;

  }

    dnafile.close();

    ofstream outputFile;
    outputFile.open("matthewNwerem.txt");
    outputFile << "Temp write to file" << endl;
    return 0;
}
