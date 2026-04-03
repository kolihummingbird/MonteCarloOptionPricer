#include<iostream>
#include <random>
#include <cmath>

using namespace std;

std::mt19937 rng(std::random_device{}());
std::normal_distribution<double> normal(0.0, 1.0);

double put_payoff(double S, double K)
{
    return std::max(K - S, 0.0);
}

double call_payoff(double S, double K)
{
    return std::max(S - K, 0.0);
}
double normal_random() 
{
	return normal(rng);
}

double normal_cdf(double x)
{
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

double black_scholes_call(double S, double K, double r, double sigma, double T)
{
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T)
                / (sigma * std::sqrt(T));

    double d2 = d1 - sigma * std::sqrt(T);

    return S * normal_cdf(d1) - K * std::exp(-r * T) * normal_cdf(d2);
}



double simulate_path(double S0, double mu, double sigma, double T, int steps)
{
	double S = S0;
	double dt = T / steps;
	for (int i=0; i < steps; ++i)
	{
		//if (i % 200 == 0) { std::cout << "Intermediate values: " << S << std::endl; }
		double Z = normal_random();
		S *= std::exp((mu - 0.5 * sigma * sigma) * dt
			+ sigma * std::sqrt(dt) * Z);
	}
	return S;
}

double monte_carlo_call(int num_paths, double S, double K, double r, double sigma, double T)
{
	double payoff_sum = 0.0;
	for (int i = 0; i < num_paths; ++i)
	{
		double Z = normal_random();
		double ST = S * std::exp((r - 0.5 * sigma * sigma) * T
			+ sigma * std::sqrt(T) * Z);
		double payoff = std::max(ST - K, 0.0);
		payoff_sum += payoff;
	}
	return std::exp(-r * T) * (payoff_sum / num_paths);
}



int main()
{

	double S = 100.0;
	double K = 100.0;
	double r = 0.05;
	double sigma = 0.2;
	double T = 1.0;
	int steps = 1000;

	//std::cout << "Monte Carlo: " << mc << std::endl;
	//std::cout << "Black-Scholes: " << bs << std::endl;
	//std::cout << "Abs error: " << std::abs(mc - bs) << endl;
	std::cout << "Paths\tMCPrice\tError\t\tError*sqrt(N)" << std::endl;
	for (int i=3;i<=6;++i)
	{
		int paths = std::pow(10,i);
		double mc = monte_carlo_call(paths, S, K, r, sigma, T);
		double bs = black_scholes_call(S, K, r, sigma, T);
		double error = std::abs(mc - bs);

		std::cout << paths <<"\t" << mc << "\t" << error <<"\t"<<error*std::sqrt(paths) << std::endl;
	}

	return 0;
}


/*
double dt = T / steps;

simulate_path(S, r, sigma, T, steps);
std::cout << "Final price: " << S << std::endl;

for (int simulations = 100; simulations <= 10000; simulations *= 10)
{

	//int simulations = 1000;
	double sum = 0.0;

	for (int i = 0; i < simulations; ++i)
	{
		double ST = simulate_path(S, r, sigma, T, steps);
		//sum += simulate_path(S, mu, sigma, T, steps);
		sum += call_payoff(ST, K);
	}

	//double average = sum / simulations;
	double price = sum / simulations;
	price *= std::exp(-0.05 * 1.0);

	//std::cout << "Average price at " << simulations << " simulations : " << average << std::endl;
	std::cout << "Call price (MC): " << price << std::endl;
}
*/
