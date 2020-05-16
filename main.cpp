#include <iostream>
#include <cmath>

#include "option.h"
#include "montecarlo.h"
#include "random.h"

int main()
{
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "FILE: " << __FILE__ << " DATE: " << __DATE__ << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    // instantiate option
    double K = 50.;
    double S = 40.;
    double T = 1.;
    double r = .02;
    double sigma = .4;
    char type = 'p';
    double price = 12.54;
    VanillaOption<double> myopt(S, K, T, r, sigma, type);

    // instantiate monte carlo parameters
    unsigned int seed = 1023;
    unsigned int iters = 10;
    std::vector<unsigned long> num_paths = {100, 500, 1000, 5000, 10000, 20000, 30000, 40000, 50000};

    // run problem 1
    BoxMuller<double> boxmuller;
    Simulator<double> sim1(myopt, num_paths, iters, boxmuller, seed);
    sim1.run_monte_carlo();
    sim1.print();

    std::cout << "\n=========================================\n\n";

    // run problem 2
    Antithetic<double> antithetic;
    Simulator<double> sim2(myopt, num_paths, iters, antithetic, seed);
    sim2.run_monte_carlo();
    sim2.print();

    return 0;
}