#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>

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
double black_scholes_delta(double S, double K, double r,
	double sigma, double T)
{
	double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));

	return normal_cdf(d1);
}

double monte_carlo_call_antithetic(int num_paths, double S,
	double K, double r, double sigma, double T)
{
	double payoff_sum = 0.0;
	for (int i = 0; i < num_paths / 2 ; ++i)
	{
		double Z = normal_random();
		
		double ST1 = S * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);
		double ST2 = S * std::exp((r - 0.5 * sigma * sigma) * T - sigma * std::sqrt(T) * Z);
		
		double payoff1 = std::max(ST1-K,0.0);
		double payoff2 = std::max(ST2 - K, 0.0);
	
		payoff_sum += (payoff1 + payoff2);
	}
	return std::exp(-r * T) * (payoff_sum / num_paths);

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

double monte_carlo_call_stats(int num_paths, double S, 
	double K, double r, double sigma, 
	double T, double& std_error)
{
	double sum = 0.0;
	double sum_sq = 0.0;

	for (int i = 0; i < num_paths; ++i)
	{
		double Z = normal_random();

		double ST = S * std::exp((r - 0.5 * sigma * sigma) * T
			+ sigma * std::sqrt(T) * Z);

		double payoff = std::exp(-r * T) * std::max(ST - K, 0.0);

		sum += payoff;
		sum_sq += payoff * payoff;
	}
	double mean = sum / num_paths;

	double variance = (sum_sq / num_paths) - (mean * mean);
	std_error = std::sqrt(variance / num_paths);
	
	return std::exp(-r * T) * mean;
}


double monte_carlo_delta(int paths, double S, double K, double r, double sigma, double T) 
{
	double eps = 0.01 * S;

	double price_up = monte_carlo_call(paths, S + eps, K, r, sigma, T);
	double price_down = monte_carlo_call(paths, S - eps, K, r, sigma, T);

	return (price_up - price_down) / (2.0 * eps);
	
}

double monte_carlo_delta_crn(int paths, double S, 
	double K, double r, double sigma, double T) //common random numbers
{
	double eps = 0.01 * S;
	double sum_up = 0.0;
	double sum_down = 0.0;

	for (int i = 0; i < paths; ++i)
	{
		double Z = normal_random();
		double ST_up = (S + eps) * std::exp((r-0.5*sigma*sigma)*T +sigma * std::sqrt(T) * Z);
		double ST_down = (S - eps) * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);
		
		sum_up += std::max(ST_up - K, 0.0);
		sum_down += std::max(ST_down - K, 0.0);
	}
	double price_up = std::exp(-r * T) * (sum_up / paths);
	double price_down = std::exp(-r * T) * (sum_down / paths);

	return (price_up - price_down) / (2.0 * eps);
}

double monte_carlo_delta_crn_stats(int paths, double S, 
    double K, double r, double sigma, double T, double& std_error)
{
    double eps = 0.01 * S;
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < paths; ++i)
    {
        double Z = normal_random();

        double ST_up = (S + eps) * std::exp((r - 0.5 * sigma * sigma) 
            * T + sigma * std::sqrt(T) * Z);
        double ST_down = (S - eps) * std::exp((r - 0.5 * sigma * sigma) 
            * T + sigma * std::sqrt(T) * Z);

        double payoff_up = std::max(ST_up - K, 0.0);
        double payoff_down = std::max(ST_down - K, 0.0);

        //finite difference per path
        double delta_i = (payoff_up - payoff_down) / (2.0 * eps);

        sum += delta_i;
        sum_sq += delta_i * delta_i;
    }
    double mean = sum / paths;

    double variance = (sum_sq / paths) - (mean * mean);
    std_error = std::exp(-r * T) * std::sqrt(variance / paths);

    return std::exp(-r * T) * mean;
}

double monte_carlo_delta_pathwise(int paths,
	double S, double K, double r, double sigma, double T)
{
	double sum = 0.0;

	for (int i = 0; i < paths; ++i)
	{
		double Z = normal_random();

		double ST = S * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);

		if (ST > K)
		{
			sum += ST / S;
		}
	}

	return std::exp(-r * T) * (sum / paths);
}

double monte_carlo_delta_pathways_stats(int paths, double S, double K,
	double r, double sigma, double T, double& std_error)
{
	double sum = 0.0;
	double sum_sq = 0.0;

	for (int i = 0; i < paths; ++i)
	{
		double Z = normal_random();

		double ST = S * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);

		double delta_i = 0.0;

		if (ST > K)
		{
			delta_i = ST / S;
		}
		double discounted_delta = std::exp(-r * T) * delta_i;
		sum += discounted_delta;
		sum_sq += discounted_delta * discounted_delta;
	}
	
	double mean = sum / paths;

	double variance = (sum_sq / paths) - (mean * mean);
	std_error = std::sqrt(variance / paths);

	return mean;
}


int main()
{

	double S = 100.0;
	double K = 100.0;
	double r = 0.05;
	double sigma = 0.2;
	double T = 1.0;
	//int steps = 10000; //for discrete
	int paths = 10000;	// for continuous


	double se;
	double price = monte_carlo_call_stats(paths, S, K, r, sigma, T, se);

	double discounted_se = std::exp(-r * T) * se;
	std::cout << "Price: " << price << std::endl;
	std::cout << "95% CI: ["
		<< price - 1.96 * discounted_se << ", "
		<< price + 1.96 * discounted_se << "]"
		<< std::endl;

    double delta_se;
    
	//double delta_mc = monte_carlo_delta(paths, S, K, r, sigma, T);
	
	double delta_bs = black_scholes_delta(S, K, r, sigma, T);
	std::cout << "Delta BS: " << delta_bs << std::endl;
	
	double delta_pw_se;
	double delta_mc_crn = monte_carlo_delta_crn(paths, S, K, r, sigma, T);
	double delta_mc = monte_carlo_delta_crn_stats(paths, S, K, r, sigma, T, delta_se);
	double delta_pw = monte_carlo_delta_pathways_stats(paths, S, K, r, sigma, T, delta_pw_se);

	std::cout << "Method        Delta      StdErr" << std::endl;
	std::cout << "CRN       " << delta_mc << "   " << delta_se << std::endl;
	std::cout << "Pathwise  " << delta_pw << "   " << delta_pw_se << std::endl;

	std::cout << "Delta MC (naive): " << delta_mc << std::endl;
	std::cout << "Delta MC (CRN): " << delta_mc_crn << std::endl;
	std::cout << "Delta Pathwise: " << delta_pw << std::endl;
		
	std::cout << "Delta 95% CI: ["
		<< delta_mc - 1.96 * delta_se << ", "
		<< delta_mc + 1.96 * delta_se << "]"
		<< std::endl;	

	//std::cout << "Monte Carlo: " << mc << std::endl;
	//std::cout << "Black-Scholes: " << bs << std::endl;
	//std::cout << "Abs error: " << std::abs(mc - bs) << endl;
	
	std::cout << std::fixed << std::setprecision(6);
	std::cout << std::setw(10) << "Paths"
			  << std::setw(15) << "MC"
			  << std::setw(15) << "Error"
			  << std::setw(20) << "ScaledErr"
			  << std::setw(15) << "AntiMC"
			  << std::setw(15) << "AntiErr"
			  << std::endl;
	
	double bs = black_scholes_call(S, K, r, sigma, T);
	for (int i=3;i<=7;++i)
	{
		int paths = std::pow(10, i);
		double amc = monte_carlo_call_antithetic(paths, S, K, r, sigma, T);
		double mc = monte_carlo_call(paths, S, K, r, sigma, T);
		double error = std::abs(mc - bs);

		double antierror = std::abs(amc - bs);

		std::cout << std::setw(10) << paths
			 	  << std::setw(15) << mc
			      << std::setw(15) << error
				  << std::setw(20) << error * std::sqrt(paths)
				  << std::setw(15) << amc
				  << std::setw(15) << antierror
				  << std::endl;
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
