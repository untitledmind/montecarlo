#ifndef MONTECARLO_HH
#define MONTECARLO_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "random.hh"
#include "option.hh"

template <typename Type> class Simulator;

template <typename Type = double>
class Simulator {
    private:
        Type price;
        Type diff;
        unsigned int seed;
        unsigned int iters;
        std::vector<unsigned long> num_paths;
        std::vector<Type> price_results;
        std::vector<Type> stdev_results;
        std::vector<Type> sterr_results;
        std::vector<std::chrono::seconds> timer;
        std::string method;
        Option<Type>* option;
        RandomGenerator<Type>* random;
    public:
        Simulator();
        Simulator(Option<Type>& rhs, RandomGenerator<Type>& _random, std::vector<unsigned long>& _num_paths, unsigned int& _iters, unsigned int& _seed);
        virtual ~Simulator();
        void black_scholes();
        std::vector<Type> monte_carlo(unsigned long& _num_paths) const;
        void run_monte_carlo();
        std::vector<Type> monte_carlo_var_vol(unsigned long& _num_paths) const;
        void run_monte_carlo_var_vol();
        void print() const;
};

template <typename Type>
Simulator<Type>::Simulator() {}

template <typename Type>
Simulator<Type>::Simulator(Option<Type>& _option, RandomGenerator<Type>& _random, std::vector<unsigned long>& _num_paths, unsigned int& _iters, unsigned int& _seed):
    option(&_option), random(&_random), num_paths(_num_paths), iters(_iters), seed(_seed) {
        srand(seed);
    }

// template <typename Type>
// Simulator<Type>::Simulator(Option<Type>& _option, StepwiseVolatility<Type>& _vol_func, RandomGenerator<Type>& _random, std::vector<unsigned long>& _num_paths, unsigned int& _iters, unsigned int& _seed):
//     option(&_option), vol_func(&_vol_func), random(&_random), num_paths(_num_paths), iters(_iters), seed(_seed) {
//         srand(seed);
//     }

template <typename Type>
Simulator<Type>::~Simulator() {}

template <typename Type>
std::vector<Type> Simulator<Type>::monte_carlo(unsigned long& _num_paths) const {
    Type variance = option->get_vol() * option->get_vol() * option->get_expiry();
    Type root_variance = sqrt(variance);
    Type moved_spot = option->get_spot() * exp(option->get_rate() * option->get_expiry() - .5*variance);

    Type factor = 1.001;
    Type base_root_variance = root_variance;

    Type running_sum;
    Type psum1;
    Type psum2;
    Type stdev;
    Type sterr;

    for (unsigned long i=0; i < _num_paths; ++i) {
        root_variance *= factor;
        if (root_variance > 1.1 * base_root_variance)
            root_variance = base_root_variance;

        Type sample = random->generate();       // sample is drawn from pure virtual function based on method
        Type this_spot = moved_spot * exp(root_variance * sample);
        Type payoff = (*option)(this_spot);

        running_sum += payoff;

        psum1 += payoff;
        psum2 += payoff*payoff;
    }
    
    Type mean = running_sum / _num_paths;
    mean *= exp(-option->get_rate() * option->get_expiry());
    stdev = sqrt((exp(-2.*option->get_rate() * option->get_expiry())/(_num_paths-1))*
        (psum2 - (psum1*psum1)/(_num_paths) ));
    sterr = stdev/sqrt(_num_paths);

    // pass statistics of the monte carlo simulation to the caller
    std::vector<Type> result = {mean, stdev, sterr};
    return result;
}


template <typename Type>
void Simulator<Type>::run_monte_carlo() {
    // copy num_paths so monte carlo can alter the copy as it progresses
    std::vector<unsigned long> monte_num_paths = num_paths;
    // srand(seed);

    // double the number of paths if using antithetic method
    if (!random->get_signature().compare(std::string("Antithetic generator"))) {
        std::for_each(monte_num_paths.begin(), monte_num_paths.end(), [](unsigned long& n) {
            n *= 2;
        });
    }
    
    for (int j=0; j < monte_num_paths.size(); ++j) {
        std::cout << "Beginning new trial with " << monte_num_paths[j] << " paths\n";
        auto start = std::chrono::system_clock::now();
        for (unsigned int i=0; i < iters; ++i) {
            
            // time the simulation
            std::vector<Type> results = monte_carlo(monte_num_paths[j]);
            Type mean = results[0];
            Type stdev = results[1];
            Type sterr = results[2];

            std::cout << std::setw(10);
            diff = mean - price;
            std::cout << "Num paths = " << monte_num_paths[j] << '\t' << "Option price = " << mean << "\tdiff = " << diff << "\tstdev = " << stdev << "\tsterr = " << sterr << '\n';
            monte_num_paths[j] *= 2;

            if (i == iters - 1) {               // log only the last estimate for each trial
                price_results.push_back(mean);
                stdev_results.push_back(stdev);
                sterr_results.push_back(sterr);
            }
        }
        auto stop = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        timer.push_back(duration);
        std::cout << '\n';
    }
    method = std::string("Monte Carlo");
}

template <typename Type>
std::vector<Type> Simulator<Type>::monte_carlo_var_vol(unsigned long& _num_paths) const {
    Type mean;
    Type running_sum;
    Type psum1;
    Type psum2;
    Type stdev;
    Type sterr;
    
    // get volatility function data before entering the loop to reduce overhead
    std::vector<Type> dt_vec = option->get_vol_func()->get_dt_vec();
    std::vector<Type> vol_vec = option->get_vol_func()->get_vol_vec();
    int steps = option->get_vol_func()->get_steps();

    for (unsigned long i=0; i < _num_paths; ++i) {
        Type this_spot = option->get_spot();
        for (int i=0; i < steps-1; ++i) {
            Type dt = dt_vec[i];
            Type vol = vol_vec[i];

            Type variance = vol * vol * dt;
            Type root_variance = sqrt(variance);
            // Type factor = 1.001;
            // Type base_root_variance = root_variance;

            // root_variance *= factor;
            // if (root_variance > 1.1 * base_root_variance)
            //     root_variance = base_root_variance;

            Type sample = random->generate();       // sample is drawn from pure virtual function based on method
            this_spot *= exp(option->get_rate() * dt - .5*variance + root_variance * sample);
        }
        Type payoff = (*option)(this_spot);

        running_sum += payoff;
        psum1 += payoff;
        psum2 += payoff*payoff;
    }
    
    mean = running_sum / _num_paths;
    mean *= exp(-option->get_rate() * option->get_expiry());
    stdev = sqrt((exp(-2.*option->get_rate() * option->get_expiry())/(_num_paths-1))*
        (psum2 - (psum1*psum1)/(_num_paths) ));
    sterr = stdev/sqrt(_num_paths);
    
    // pass statistics of the monte carlo simulation to the caller
    std::vector<Type> result = {mean, stdev, sterr};
    return result;
}


template <typename Type>
void Simulator<Type>::run_monte_carlo_var_vol() {
    // copy num_paths so monte carlo can alter the copy as it progresses
    std::vector<unsigned long> monte_num_paths = num_paths;
    // srand(seed);

    // double the number of paths if using antithetic method
    if (!random->get_signature().compare(std::string("Antithetic generator"))) {
        std::for_each(monte_num_paths.begin(), monte_num_paths.end(), [](unsigned long& n) {
            n *= 2;
        });
    }
    
    for (int j=0; j < monte_num_paths.size(); ++j) {
        std::cout << "Beginning new trial with " << monte_num_paths[j] << " paths\n";
        auto start = std::chrono::system_clock::now();
        for (unsigned int i=0; i < iters; ++i) {
            
            // time the simulation
            std::vector<Type> results = monte_carlo_var_vol(monte_num_paths[j]);
            Type mean = results[0];
            Type stdev = results[1];
            Type sterr = results[2];

            std::cout << std::setw(10);
            diff = mean - price;
            std::cout << "Num paths = " << monte_num_paths[j] << '\t' << "Option price = " << mean << "\tdiff = " << diff << "\tstdev = " << stdev << "\tsterr = " << sterr << '\n';
            monte_num_paths[j] *= 2;

            if (i == iters - 1) {               // log only the last estimate for each trial
                price_results.push_back(mean);
                stdev_results.push_back(stdev);
                sterr_results.push_back(sterr);
            }
        }
        auto stop = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        timer.push_back(duration);
        std::cout << '\n';
    }
    method = std::string("Monte Carlo");
}

template <typename Type>
void Simulator<Type>::print() const {
    std::cout << "Method used: " << method << " with " << random->get_signature() << '\n';
    std::cout << "Trials  Time    Result       Std deviation    Standard error\n";
    for (int i=0; i<num_paths.size(); ++i)
        std::cout << num_paths[i] << '\t' << timer[i].count() << " s\t" << price_results[i] << '\t' << stdev_results[i] << '\t' << sterr_results[i] << '\n';
}

#endif