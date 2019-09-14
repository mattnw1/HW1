#include <iostream>
#include <fstream> //read and write to files
#include <stdio.h>
#include <string>
#include <cmath>
#include <ctime>
#include <time.h>
#include <algorithm>




using namespace std;

int main(int argc, char const *argv[])
{
  cout << "This program is asimple analysis program\nthat will compute basic statistics\nfor a list of DNA strings." << endl << endl;
  cout << "Summary statistics have been printed to a text file, along with 1000 DNA strings\nwhose lengths follow a Gaussian distribution based on previously computed basic statistics." << endl;
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
  double relativeProbA = 0;
  double relativeProbC = 0;
  double relativeProbT = 0;
  double relativeProbG = 0;
  double relativeProbAA = 0;
  double relativeProbAC = 0;
  double relativeProbAT = 0;
  double relativeProbAG = 0;
  double relativeProbCA = 0;
  double relativeProbCC = 0;
  double relativeProbCT = 0;
  double relativeProbCG = 0;
  double relativeProbTA = 0;
  double relativeProbTT = 0;
  double relativeProbTC = 0;
  double relativeProbTG = 0;
  double relativeProbGA = 0;
  double relativeProbGG = 0;
  double relativeProbGC = 0;
  double relativeProbGT = 0;
  const double pi = 3.14159265358979;
  bool anotherList = true;



  ifstream dnafile;
  dnafile.open("DNA_Strings.txt");
  int numOfLines = 0;
  string line;
  int dnaSum = 0;
  int dnaMean = 0;
  int dnaVariance = 0;
  int dnaStdDev = 0;

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
    dnaMean = dnaSum/numOfLines;
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

    dnaVariance=dnaVariance/numOfLines;
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

    //outputFile << "Temp write to file" << endl;
    outputFile << "Name: Matthew Nwerem\nStudentID: 2277158\nInstagram: mattnw\nLinkedIn: linkedin.com/in/matthewnwerem/\nCPSC 350-01\nAssignment 1\n" << endl; //instagram is an ect.
    outputFile << "The sum length of DNA nucleotides is: " << dnaSum << endl;
    outputFile << "The mean length of DNA nucleotides is " << dnaMean << endl;
    outputFile << "The variance length of DNA nucleotides is: " << dnaVariance << endl;
    outputFile << "The standard deviation length of DNA nucleotides is: " << dnaStdDev << endl << endl;
    outputFile << "--Below prints the relative probability of each nucleotide--" << endl << endl;
    outputFile << "A:  " << relativeProbA << endl;
    outputFile << "C:  " << relativeProbC << endl;
    outputFile << "T:  " << relativeProbT << endl;
    outputFile << "G:  " << relativeProbG << endl << endl;
    outputFile << "--Below prints the relative probability of each nucleotide bigram--" << endl << endl;
    outputFile << "AA:  " << relativeProbAA << endl;
    outputFile << "AC:  " << relativeProbAC << endl;
    outputFile << "AT:  " << relativeProbAT << endl;
    outputFile << "AG:  " << relativeProbAG << endl;
    outputFile << "CA:  " << relativeProbCA << endl;
    outputFile << "CC:  " << relativeProbCC << endl;
    outputFile << "CT:  " << relativeProbCT << endl;
    outputFile << "CG:  " << relativeProbCG << endl;
    outputFile << "TA:  " << relativeProbTA << endl;
    outputFile << "TC:  " << relativeProbTC << endl;
    outputFile << "TT:  " << relativeProbTT << endl;
    outputFile << "TG:  " << relativeProbTG << endl;
    outputFile << "GA:  " << relativeProbGA << endl;
    outputFile << "GC:  " << relativeProbGC << endl;
    outputFile << "GT:  " << relativeProbGT << endl;
    outputFile << "GG:  " << relativeProbGG << endl;

      //gaussian distribution
      double randomVariableC;
      double randomVariableD;
      double randomNumber;

      for (int i = 0; i < 1000; i++)
      {
        double randomA = (RAND_MAX - rand())/(double)(RAND_MAX);
        double randomB = (RAND_MAX - rand())/(double)(RAND_MAX);
        string printedLine;
        srand(time(0));
        randomVariableC = (sqrt(-2*log(randomA))*cos(2*pi*randomB));
        randomVariableD = round((dnaStdDev*randomVariableC) + dnaMean);
        for (int j = 0; j < randomVariableD; ++j)
        {
          randomNumber = (RAND_MAX - rand())/(double)(RAND_MAX);
          if(randomNumber < relativeProbA)
          {
            printedLine += "A";
          }
          else if (randomNumber < (relativeProbA + relativeProbC))
          {
            printedLine += "C";
          }
          else if (randomNumber < (relativeProbA + relativeProbC + relativeProbT))
          {
            printedLine += "T";
          }
          else
          {
            printedLine += "G";
          }
        }
        outputFile << printedLine << endl;
      }
      outputFile.close();

      //while (anotherList)
      //{

      cout << "--List Processed--" << endl;
      cout << "Would you like to process another list?\nY = Yes\nAny other character = No" << endl;
      char answerConsole = cin.get();
      if (toupper(answerConsole) == 'Y')
      {
        anotherList = true;
        string newFileName;
        ofstream f;
        cout << "Enter a name for the new file to be created" << endl;
        cin >> newFileName;
        cout << "\n--File Created--" << endl;
        getline(cin, newFileName);
        f.open(newFileName);

        for (int i = 0; i < 1000; i++)
        {
          double randomA = (RAND_MAX - rand())/(double)(RAND_MAX);
          double randomB = (RAND_MAX - rand())/(double)(RAND_MAX);
          string printedLine;
          srand(time(0));
          randomVariableC = (sqrt(-2*log(randomA))*cos(2*pi*randomB));
          randomVariableD = round((dnaStdDev*randomVariableC) + dnaMean);
          for (int j = 0; j < randomVariableD; ++j)
          {
            randomNumber = (RAND_MAX - rand())/(double)(RAND_MAX);
            if(randomNumber < relativeProbA)
            {
              printedLine += "A";
            }
            else if (randomNumber < (relativeProbA + relativeProbC))
            {
              printedLine += "C";
            }
            else if (randomNumber < (relativeProbA + relativeProbC + relativeProbT))
            {
              printedLine += "T";
            }
            else
            {
              printedLine += "G";
            }
          }
          //newFileName << printedLine << endl;
        }
      }
      else
      {
        cout << "Exiting Program..." << endl;
        anotherList = false;
        exit(0);
      }

      //}





    return 0;
}
