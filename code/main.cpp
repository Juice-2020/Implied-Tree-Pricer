#include <fstream>
#include <vector>
#include "Market.h"
#include "TreePricer.h"
#include "EuropeanTrade.h"
#include "AmericanTrade.h"
#include "BarrierTrade.h"
#include <utility>
#include <random>
#include <boost/math/distributions/normal.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
using namespace boost::numeric;
void test1()
{
  std::ifstream fin("simpleMkt.dat");
  if (fin) {
    Market mkt;
    fin >> mkt;
    mkt.Print();
    EuropeanOption callOption1(Call, 90, Date(2016, 1, 1));
    EuropeanOption callOption2(Call, 91, Date(2016, 1, 1));
    EuroCallSpread callSpread(90, 91, Date(2016, 1, 1));

    std::cout << "callOption1 struck at 90:\t"
	      << binomialTree(mkt, callOption1, 100) << std::endl;      
    std::cout << "callOption2 struck at 91:\t"
	      << binomialTree(mkt, callOption2, 100) << std::endl;
    std::cout << "callSpread 90 to 91:\t\t"
	      << binomialTree(mkt, callSpread, 100) << std::endl;
    std::cout << "callOption1 - callOption2:\t" 
	      << binomialTree(mkt, callOption1, 100) - binomialTree(mkt, callOption2, 100) << std::endl;
  }    
  else {
    std::cerr << "no simpleMkt.dat found" << std::endl;
  }
}

// american options
void test2()
{
  std::ifstream fin("simpleMkt.dat");
  if (fin) {
    Market mkt;
    fin >> mkt;
    mkt.Print();
    AmericanOption callOption1(Call, 90, Date(2016, 1, 1));
    AmericanOption callOption2(Call, 91, Date(2016, 1, 1));
    AmerCallSpread callSpread(90, 91, Date(2016, 1, 1));

    std::cout << "callOptionAmer1 struck at 90:\t"
	      << binomialTree(mkt, callOption1, 100) << std::endl;
    std::cout << "callOptionAmer2 struck at 91:\t"
	      << binomialTree(mkt, callOption2, 100) << std::endl;
    std::cout << "callSpreadAmer 90 to 91:\t\t"
	      << binomialTree(mkt, callSpread, 100) << std::endl;
  }
  else {
    std::cerr << "no simpleMkt.dat found" << std::endl;
  }
}

// European vs American 
void testAmer()
{
  std::cout << "\n------------------------ Test American -------------------\n";
  std::ifstream fin("simpleMkt.dat");
  std::ofstream foutA1("amerCall.txt");
  std::ofstream foutE1("euroCall.txt");  
  std::ofstream foutA2("amerPut.txt");
  std::ofstream foutE2("euroPut.txt");  
  if (fin) {
    Market mkt;
    fin >> mkt;
    mkt.Print();
    const int treeTimeSteps = 100;
    for (int i = 0; i <= 200; i++) {
      EuropeanOption euroCall(Call, i, Date(2016, 1, 1));
      AmericanOption amerCall(Call, i, Date(2016, 1, 1));
      foutE1 << i << "\t\t\t"
	     << binomialTree(mkt, euroCall, treeTimeSteps)
	     << std::endl;
      foutA1 << i << "\t\t\t"
	    << binomialTree(mkt, amerCall, treeTimeSteps)
	     << std::endl;	
    }
    for (int i = 0; i <= 200; i++) {
      EuropeanOption euroPut(Put, i, Date(2016, 1, 1));
      AmericanOption amerPut(Put, i, Date(2016, 1, 1));
      foutE2 << i << "\t\t\t"
	     << binomialTree(mkt, euroPut, treeTimeSteps)
	     << std::endl;
      foutA2 << i << "\t\t\t"
	     << binomialTree(mkt, amerPut, treeTimeSteps)
	     << std::endl;
    } 
  }
  else
    std::cerr << "no simpleMkt.dat found" << std::endl;
}

void testBarrier()
{
  std::cout << "\n------------------------ Test Barrier -------------------\n";
  std::ifstream fin("simpleMkt.dat");
  if (fin) {
    Market mkt;
    fin >> mkt;
    mkt.Print();
    EuropeanOption callOption1(Call, 90, Date(2016, 1, 1));
    UpBarrier upBar(120);
    BarrierOption barOption1(upBar, callOption1);
    

    std::cout << "callOption1 struck at 90:\t\t\t"
	      << binomialTree(mkt, callOption1, 100) << std::endl;
    std::cout << "barrierOption1 with up barrier at 110:\t\t"
	      << binomialTree(mkt, barOption1, 100) << std::endl;

    // changing upBar
    // upBar = UpBarrier(110);
    // std::cout << "changing upBar: \t\t\t\t" << binomialTree(mkt, barOption1, 100) << std::endl;
  }
  else {
    std::cerr << "no simpleMkt.dat found" << std::endl;
  }
}

void testBarrier2()
{
  std::cout << "\n------------------------ Test Barrier2 -------------------\n";
  std::ifstream fin("simpleMkt.dat");
  std::ofstream foutE1("barrierTest_euroCall.txt");
  std::ofstream foutB1("barrierTest_bar.txt");
  if (fin) {
    Market mkt;
    fin >> mkt;
    mkt.Print();
    const int treeTimeSteps = 100;
    EuropeanOption euroCall(Call, 80, Date(2016, 1, 1));
    for (int i = 0; i <= 200; i++) {
      BarrierOption barOption1(DownBarrier(i), euroCall); 
      foutE1 << i << "\t\t\t"
	     << binomialTree(mkt, euroCall, treeTimeSteps)
	     << std::endl;
      foutB1 << i << "\t\t\t"
	     << binomialTree(mkt, barOption1, treeTimeSteps)
	     << std::endl;
    }
  }
  else
    std::cerr << "no simpleMkt.dat found" << std::endl;   
}

void testTreePricer()
{
  std::cout << "\n------------------------ Test TreePricer -------------------\n";
  std::ifstream fin("simpleMkt.dat");
  if (fin) {
    Market mkt;
    fin >> mkt;
    mkt.Print();
    CRRBinomialTreePricer crrBTreePricer(1000); // nTimeStep = 100

    EuropeanOption euroCall(Call, 80, Date(2016, 1, 1));

    std::cout << crrBTreePricer.Price(mkt, euroCall) << std::endl;
    std::cout << binomialTree(mkt, euroCall, 100) << std::endl;


    JRRNBinomialTreePricer jrrnBTreePricer(1000); // nTimeStep = 100
    std::cout << jrrnBTreePricer.Price(mkt, euroCall) << std::endl;
    
  }
  else
    std::cerr << "no simpleMkt.dat found" << std::endl; 
}

std::vector<std::string> split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

void ShowArray2(std::vector<double> Nrange, std::vector<double> CRR_price, std::vector<double> Imp_price, double nTenors, int precision = 4)
{
    // print out the results
    //std::cout << std::setw(7) << "| Delta\\T | ";
    std::cout << std::setw(7) << "|   Step_Num   | ";
    for (int i = 0; i < nTenors; i++)
        std::cout << std::fixed << std::setprecision(2) << std::setw(7) << Nrange[i] << " | ";
    std::cout << std::endl;
    std::cout << std::setw(7) << "|  CRR_Pricer  | ";
    for (int i = 0; i < nTenors; i++)
        std::cout << std::fixed << std::setprecision(4) << std::setw(7) << CRR_price[i] << " | ";
    std::cout << std::endl;
    std::cout << std::setw(7) << "| Implied_Tree | ";
    for (int i = 0; i < nTenors; i++)
        std::cout << std::fixed << std::setprecision(4) << std::setw(7) << Imp_price[i] << " | ";
    std::cout << std::endl;
}  


void test_impliedTreePricer()
{
    std::ifstream din;
    din.open("impliedvol.txt");
    if (!din.is_open()) {
        std::cout << "no";
    }
    std::string line;
    ImpVol impvol;

    int i = 0;

    while (getline(din, line))   //整行读取，换行符“\n”区分，遇到文件尾标志eof终止读取
    {
        //std::cout << "原始字符串：" << line << std::endl; //整行输出
        char a = '\t';
        std::vector<std::string> fields;
        fields = split(line, a); //清除掉向量fields中第一个元素的无效字符，并赋值给变量name
        impvol.time[i] = stod(fields[0]); //清除掉向量fields中第二个元素的无效字符，并赋值给变量age
        impvol.strike[i] = stod(fields[1]); //清除掉向量fields中第三个元素的无效字符，并赋值给变量birthday
        impvol.ivol[i] = stod(fields[2]);
        i = i + 1;
    }
    std::ifstream fin("simpleMkt.dat");
    if (fin) {
        Market mkt;
        fin >> mkt;
        mkt.volatility = 0.1444;
        mkt.stockPrice = 1.25;
        mkt.interestRate = 0.01;
        std::vector<double> Nrange(10);
        std::vector<double> CRR_price(10);
        std::vector<double> Imp_price(10);
        cout << "const volatility:\t" << mkt.volatility << endl;
        cout << "stock price:\t\t" << mkt.stockPrice << endl;
        cout << "interest rate:\t\t" << mkt.interestRate << endl;
        cout << "down_out barrier\t" << 1.2 << endl;
        int nTenors = 10;
        std::cout << "\n------------------------------------------------ European Option --------------------------------------------------\n";
        for (int i = 0; i < nTenors; i++) {
            Nrange[i] = 10 + i * 20;
            EuropeanOption euroCall(Call, 1.2, Date(2016, 1, 1));
            Imp_price[i] = impliedTree_pricer(mkt, euroCall, Nrange[i], impvol);
            CRRBinomialTreePricer crrBTreePricer(Nrange[i]);
            //std::cout << crrBTreePricer.Price(mkt, euroCall) << std::endl;
            CRR_price[i] = crrBTreePricer.Price(mkt, euroCall);
            //std::cout << bsformula( C,1.2, 1, 1.25, 0.13, 0.01, 0) << std::endl;
        }
        ShowArray2(Nrange,CRR_price, Imp_price, nTenors);
        std::cout << "\n------------------------------------------------ American Option --------------------------------------------------\n";
        for (int i = 0; i < nTenors; i++) {
            Nrange[i] = 10 + i * 20;
            AmericanOption amCall(Call, 1.2, Date(2016, 1, 1));
            Imp_price[i] = impliedTree_pricer(mkt, amCall, Nrange[i], impvol);
            CRRBinomialTreePricer crrBTreePricer(Nrange[i]);
            //std::cout << crrBTreePricer.Price(mkt, euroCall) << std::endl;
            CRR_price[i] = crrBTreePricer.Price(mkt, amCall);
            //std::cout << bsformula( C,1.2, 1, 1.25, 0.13, 0.01, 0) << std::endl;
        }
        ShowArray2(Nrange, CRR_price, Imp_price, nTenors);
        std::cout << "\n------------------------------------------------ Barrier Option --------------------------------------------------\n";
        for (int i = 0; i < nTenors; i++) {
            Nrange[i] = 10 + i * 20;
            EuropeanOption euroCall(Call, 1.2, Date(2016, 1, 1));
            BarrierOption barOption(DownBarrier(1.2), euroCall);
            Imp_price[i] = impliedTree_pricer(mkt, barOption, Nrange[i], impvol);
            CRRBinomialTreePricer crrBTreePricer(Nrange[i]);
            //std::cout << crrBTreePricer.Price(mkt, euroCall) << std::endl;
            CRR_price[i] = crrBTreePricer.Price(mkt, barOption);
            //std::cout << bsformula( C,1.2, 1, 1.25, 0.13, 0.01, 0) << std::endl;
         }
        ShowArray2(Nrange, CRR_price, Imp_price, nTenors);
    }
}



int main()
{
  //test1();
  //test2();
  //testAmer();
  // testBarrier();
  // testBarrier2();
  test_impliedTreePricer();
  //  testTreePricer();
  return 0;
}
