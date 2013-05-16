#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

//For reference see chapter 10 of the book Beautiful Testing Edited by Tim Riley and Adam Goucher

//REMARK: From here on it is understood that we are working with gaussian distribution functions!

void print(vector<double> a){
  cout << "(";
  for(int i=0; i<a.size()-1; i++)
    cout << a[i] << "," << endl;
  cout << a[a.size()-1] << ")" << endl;
}

//This function is the cumulative distribution function of a gaussian of mean
//mu and standard deviation sigma calculated in y
double nc(double x, double mu=0, double sigma=1){
  if(mu==0 && sigma==1)
    return 0.5*erfc(-x/sqrt(2));
  else{
    double y=(x-mu)/sigma;
    return 0.5*erfc(-y/sqrt(2));
  }
}

void test_nc(){
  cout.precision(20);
  //Expected result for nc(i)
  double v[8]={0.0013498980316300945267,0.022750131948179207200,0.15865525393145705141,0.50000000000000000000,
               0.84134474606854294859,0.97724986805182079280,0.99865010196836990547,0.99996832875816688008};
  //Expected result for nc(i,0,sqrt(0.5))
  double u[8]={0.00001104524849929272068638806,0.002338867490523632918965372,0.07864960352514256532938968,0.5000000000000000000000000,
               0.9213503964748574346706103,0.9976611325094763670810346,0.9999889547515007072793136,0.9999999922913710498599906};
  for(int i=-3; i<5; i++){
    cout << "nc(" << i << ")=" << nc((double)i) << endl;
    cout << " v[" << i+3 << "]=" << v[i+3] << endl;
    cout << "nc(" << i << ",0,sqrt[0.5])=" << nc((double)i,0.,sqrt(0.5)) << endl;
    cout << "             u[" << i+3 << "]=" << u[i+3] << endl;
  }
}

//Function that construct Fn(x_i). Take a large number of samples n.
//For each sample x_i we can compare the actual proportion of samples
//less than x_i to the proportion of samples we would expect to have
//seen. In other words, we will compare the empirical distribution 
//function with the theoretical distribution function.
//The empirical distribution is defined as:
//
//              the number of x_j <= x_i
//  Fn(x_i) = ----------------------------
//                        n
//
//while the theroetical distribution is the cumulative distribution function nc.
//
//The output vector contains Fn(x) with x varying from mu-5*sigma to mu+5*sigma by 0.1*sigma.
//The results are placed in the output vector into ascending order with respect to x.
vector<double> Fn(vector<double> sample, double mu=0, double sigma=1){
  vector<double> out;
  for(int i=0; i<100; i++)
    out.push_back((double)(count_if(sample.begin(),sample.end(),bind2nd(less<double>(), mu+(-5+0.1*i)*sigma)))/sample.size());
  return out;
}

//Function that construct F(x).
//The output vector contains F(x) with x varying from mu-5*sigma to mu+5*sigma by 0.1*sigma.
//The results are placed in the output vector into ascending order with respect to x.
vector<double> F(double mu=0, double sigma=1){
  vector<double> out;
  for(int i=0; i<100; i++)
    out.push_back(nc(mu+(-5+0.1*i)*sigma,mu,sigma));
  return out;
}  

//This function calculates K+ that is defined as sqrt(n)*max[Fn(x)-F(x)]
double Kplus(vector<double> sample, double mu=0, double sigma=1){
  vector<double> teo=F(mu,sigma);
  vector<double> exp=Fn(sample,mu,sigma);
  vector<double> diff(teo.size());
  for(int i=0; i<diff.size(); i++)
    diff[i]=exp[i]-teo[i];
  return (*max_element(diff.begin(),diff.end()))*sqrt(sample.size());
}

//This function calculates K- that is defined as sqrt(n)*max[F(x)-Fn(x)]
double Kminus(vector<double> sample, double mu=0, double sigma=1){
  vector<double> teo=F(mu,sigma);
  vector<double> exp=Fn(sample,mu,sigma);
  vector<double> diff(teo.size());
  for(int i=0; i<diff.size(); i++)
    diff[i]=teo[i]-exp[i];
  return (*max_element(diff.begin(),diff.end()))*sqrt(sample.size());
}

//Functions to count how often the entries of a vector are between 0.07089 and 1.5174
bool is_between (double i) { return (i>0.07089 && i<1.5174); }

double how_often_between(vector<double> sample){
  return (double)(count_if(sample.begin(),sample.end(),is_between))/sample.size();
}

//This is the Kolmogorof_Smirnov test for a set of set of samples. Only Kplus is considered.
double Kolmogorov_Smirnov(vector<vector<double>> samples, double mu=0, double sigma=1){
  vector<double> Kp;
  for(int i=0; i<samples.size(); i++)
    Kp.push_back(Kplus(samples[i],mu,sigma));
  return how_often_between(Kp);
}  

