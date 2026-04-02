#include<iostream>

double generate TerminalPrice(double S0, double r, double sigma, double T, std::mt19937_64& rng)
{
	std::normal_distribution<double> norm(0.0, 1.0);
	double Z = norm(rng);
	return SO * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);
}