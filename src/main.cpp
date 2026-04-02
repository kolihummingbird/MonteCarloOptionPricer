#include<iostream>
#include <random>
#include <cmath>

using namespace std;

std::mt19937 rng(std::random_device{}());
std::normal_distribution<double> normal(0.0, 1.0);

double normal_random() 
{
	return normal(rng);
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

int main()
{

	double S = 100;
	double mu = 0.05;
	double sigma = 0.2;
	double T = 1.0;
	int steps = 1000;

	double dt = T / steps;
	simulate_path(S, mu, sigma, T, steps);
	std::cout << "Final price: " << S << std::endl;

	for (int simulations = 100; simulations <= 10000; simulations *= 10)
	{

		//int simulations = 1000;
		double sum = 0.0;
		for (int i = 0; i < simulations; ++i)
		{
			sum += simulate_path(S, mu, sigma, T, steps);
		}
		double average = sum / simulations;

		std::cout << "Average price at " << simulations << " simulations : " << average << std::endl;
	}

	return 0;
}
