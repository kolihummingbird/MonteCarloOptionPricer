#include<iostream>
struct CallPayoff 
{
	double K;
	double operator() (double S) const { return std::max(S - K, 0.0);}
};
