#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>


void cell(std::string seq1, std::string seq2)
{
    int m = seq1.length() + 1;
    int n = seq2.length() + 1;
    int i,j;
    //assign penalties:
    int gap = 2;
    //int match = 1;
    //int mismatch = 1;
    //int score = 0;
    int dnaarr[m][n];
    //initialize gap penalty row
    //dnaarr[0][0] = 0;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            dnaarr[i][j] = 0;
        }
    }

    for (i = 1; i < m; i++)
    {
        dnaarr[i][0] = dnaarr[i-1][0] - gap;
    }
    //initialize column with gap penalty
    for (j = 1; j < n; j++)
    {
        dnaarr[0][j] = dnaarr[0][j-1] - gap;
    }
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            std::cout << dnaarr[i][j] << " ";

        }
        std::cout << std::endl;
    }
    // set up penalties
    for(i = 1; i < m; i++)
    {
        for(j = 1; j < n; j++)
        {
            //int add = dnaarr[i][j];
            //int minus = dnaarr[i][j];
            if(seq1[j - 1] == seq2[i - 1])
            {
                int y = std::max(dnaarr[i-1][j],dnaarr[i][j-1]);
                int plus = std::max(dnaarr[i-1][j-1],y);
                dnaarr[i][j] = plus + 1;
            }
            if(seq1[j - 1] != seq2[i - 1])
            {
                int x = std::max(dnaarr[i-1][j],dnaarr[i][j-1]);
                int sub = std::max(dnaarr[i-1][j-1],x);
                dnaarr[i][j] = sub -1;
            }
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            std::cout << dnaarr[i][j] << " ";

        }
        std::cout << std::endl;
    }
    //initialize traceback
    std::string dna1;
    std::string dna2;
    char direct;
    char x;
    char y;
    i=m-1;
    j=n-1;
    while(i>0 && j>0)
    {
        //we start with cell to the left
        direct = 'l';
        int al = dnaarr[i-1][j];
        int gox = i - 1;
        int goy = j;
        // comparison with the other cells
        int diag = dnaarr[i-1][j-1];
        int up = dnaarr[i][j-1];
        if(diag > al)
        {
            // travels diagonal
            direct='d';
            al = diag;
            gox = i - 1;
            goy = j - 1;
        }
        if(up > al)
        {
            // travles up/vertical
            direct='u';
            al = up;
            gox = i;
            goy = j - 1;
        }
        if(direct == 'd')
        {
            x = seq1[gox];
            y = seq2[goy];
        }
        if(direct == 'u')
        {
            x = '-';
            y = seq2[goy];
        }
        if(direct == 'l')
        {
            x = seq1[gox];
            y = '-';
        }
        dna1 = x + dna1;
        dna2 = y + dna2;
        i = gox;
        j = goy;
        std::cout << "coord:" << i << " " << j << std::endl;
    }
    std::cout << dna1 << " ";
    std::cout << std::endl;
    std::cout << dna2 << " ";
    std::cout << std::endl;
    std::cout << std::endl;
}


int main(int argc, char * argv[])
{
    //assign penalties:
    //int match = 1;
    //int mismatch = -1;
    std::string stuff = "dnadata.txt";
    //int i = 0;
    //int j = 0;
    std::string seq1;
    std::string seq2;
    std::cout << "Enter first string of DNA: ";
    std::cin >> seq1;
    std::cout << "Enter second string of DNA: ";
    std::cin >> seq2;
    //string::size_type sz;
    int n = seq1.length();
    int m = seq2.length();
    int dnaarr[m][n];
    cell(seq1,seq2);
    return 0;
}
