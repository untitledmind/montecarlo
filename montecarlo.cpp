#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>

template <typename Type> class Option;
template <typename Type> class Simulator;
template <typename Type> class RandomGenerator;
template <typename Type> class BoxMuller;       // derived
template <typename Type> class Antithetic;      // derived

template <typename Type = double>
class Option {
    private:
        Type spot;
        Type strike;
        Type expiry;
        Type rate;
        Type vol;
        char its_type;
    public:
        Option() {}
        Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const Type& _vol, char& _its_type):
            spot(_spot), strike(_strike), expiry(_expiry), rate(_rate), vol(_vol), its_type(_its_type) {}
        Option(const Option& rhs): spot(rhs.spot), strike(rhs.strike), expiry(rhs.expiry), rate(rhs.rate), vol(rhs.vol), its_type(rhs.its_type) {}
        virtual ~Option() {}
        Type& operator= (const Option& rhs) {
            if (this == &rhs)
                return *this;
            spot = rhs.spot;
            strike = rhs.strike;
            expiry = rhs.expiry;
            rate = rhs.rate;
            vol = rhs.vol;
            its_type = rhs.its_type;
            return *this;
        }
        Type operator()(Type& spot) const {
            // functor returns payoff based on option type
            switch (its_type) {
                case 'c':
                    return std::max(spot - strike, 0.0);
                case 'p':
                    return std::max(strike - spot, 0.);
                default:
                    throw("Unknown option type");
            };
        }
        Type get_spot() const {return spot;}
        Type get_strike() const {return strike;}
        Type get_expiry() const {return expiry;}
        Type get_rate() const {return rate;}
        Type get_vol() const {return vol;}
        char get_its_type() const {return its_type;}
};

template <typename Type = double>
class RandomGenerator {
    public:
        RandomGenerator() {}
        virtual ~RandomGenerator() {}
        virtual Type generate() = 0;
        virtual std::string get_signature() const = 0;
};

template <typename Type = double>
class BoxMuller: public RandomGenerator<Type> {
    private:
        std::string signature = "Box Muller generator";
    public:
        BoxMuller() {}
        virtual ~BoxMuller() {}
        virtual Type generate() {
            Type x = 0.;
            Type y = 0.;
            Type euclid = 0.;

            do {
                x = 2.*rand() / static_cast<Type>(RAND_MAX) - 1;
                y = 2.*rand() / static_cast<Type>(RAND_MAX) - 1;
                euclid = x*x + y*y;
            }
            while (euclid >= 1.);

            return x * sqrt(-2*log(euclid) / euclid);
        }
        std::string get_signature() const {return signature;}
};

template <typename Type = double>
class Antithetic: public RandomGenerator<Type> {
    private:
        Type epsilon;       // store the last random number generated
        bool odd_even;
        std::string signature = "Antithetic generator";
    public:
        Antithetic() {}
        virtual ~Antithetic() {}
        virtual Type generate() {
            if (odd_even) {
                Type x = 0.;
                Type y = 0.;
                Type euclid = 0.;

                do {
                    x = 2.*rand() / static_cast<Type>(RAND_MAX) - 1;
                    y = 2.*rand() / static_cast<Type>(RAND_MAX) - 1;
                    euclid = x*x + y*y;
                }
                while (euclid >= 1.);
                epsilon = x * sqrt(-2*log(euclid) / euclid);
                odd_even = false;

                return epsilon;
            }
            else {
                odd_even = true;
                return -epsilon;
            }
        }
        std::string get_signature() const {return signature;}
};

template <typename Type = double>
class Simulator {
    private:
        Type price = .593;
        Type diff;
        unsigned int seed;
        unsigned int iters;
        std::vector<unsigned long> num_paths;

        // containers to hold simulation results
        std::vector<Type> price_results;
        std::vector<Type> stdev_results;
        std::vector<Type> sterr_results;
        std::vector<std::chrono::seconds> timer;
        std::string method;
        Option<Type> opt;
        
        RandomGenerator<Type>* random;      // pointer to base random class

    public:
        Simulator() {}
        Simulator(Option<Type>& rhs, std::vector<unsigned long>& _num_paths, unsigned int& _iters, RandomGenerator<Type>& _random, unsigned int& _seed): opt(rhs), num_paths(_num_paths), iters(_iters), random(&_random), seed(_seed) {}
        virtual ~Simulator() {}
        std::vector<Type> monte_carlo(unsigned long& _num_paths) const {
            Type variance = opt.get_vol() * opt.get_vol() * opt.get_expiry();
            Type root_variance = sqrt(variance);
            Type moved_spot = opt.get_spot() * exp(opt.get_rate() * opt.get_expiry() + -.5*variance);

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
                Type payoff = opt(this_spot);

                running_sum += payoff;

                psum1 += payoff;
                psum2 += payoff*payoff;
            }

            Type mean = running_sum / _num_paths;
            mean *= exp(-opt.get_rate() * opt.get_expiry());
            stdev = sqrt((exp(-2.0*opt.get_rate() * opt.get_expiry())/(_num_paths-1))*
                (psum2 - (psum1*psum1)/(_num_paths) ));
            sterr = stdev/sqrt(_num_paths);

            // pass statistics of the monte carlo simulation to the caller
            std::vector<Type> result = {mean, stdev, sterr};
            return result;
        }

        void run_monte_carlo() {
            // copy num_paths so monte carlo can alter the copy as it progresses
            std::vector<unsigned long> monte_num_paths = num_paths;

            // double the number of paths if using antithetic method
            if (!random->get_signature().compare(std::string("Antithetic generator")))
                std::for_each(monte_num_paths.begin(), monte_num_paths.end(), [](unsigned long& n){n *= 2;});
            
            for (int j=0; j < monte_num_paths.size(); ++j) {
                std::cout << "Beginning new trial with " << monte_num_paths[j] << " paths\n";
                for (unsigned int i=0; i < iters; ++i) {
                    srand(seed);
                    
                    // time the simulation
                    auto start = std::chrono::high_resolution_clock::now();

                    std::vector<Type> results = monte_carlo(monte_num_paths[j]);
                    Type mean = results[0];
                    Type stdev = results[1];
                    Type sterr = results[2];

                    auto stop = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

                    std::cout << std::setw(10);
                    diff = mean - price;
                    std::cout << "Num paths = " << monte_num_paths[j] << '\t' << "Option price = " << mean << "\tdiff = " << diff << "\tstdev = " << stdev << "\tsterr = " << sterr << '\n';
                    monte_num_paths[j] *= 2;

                    if (i == iters - 1) {               // log only the last estimate for each trial
                        price_results.push_back(mean);
                        stdev_results.push_back(stdev);
                        sterr_results.push_back(sterr);
                        timer.push_back(duration);
                    }
                }
                std::cout << '\n';
            }
            method = std::string("Monte Carlo");
        }

        void print() const {
            std::cout << "Method used: " << method << " with " << random->get_signature() << '\n';
            std::cout << "Trials  Time    Result       Std deviation    Standard error\n";
            for (int i=0; i<num_paths.size(); ++i) {
                std::cout << num_paths[i] << '\t' << timer[i].count() << " s\t" << price_results[i] << '\t' << stdev_results[i] << '\t' << sterr_results[i] << '\n';
            }
        }
};

int main()
{
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "FILE: " << __FILE__ << " DATE: " << __DATE__ << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    // instantiate option
    double S = 56.;
    double K = 50.;
    double T = .50;
    double r = .08;
    double sigma = .20;
    char type = 'p';
    double price = 5.428;
    Option<double> myopt(S, K, T, r, sigma, type);

    // instantiate monte carlo parameters
    unsigned int seed = 1023;
    unsigned int iters = 10;
    std::vector<unsigned long> num_paths = {100, 1000, 10000, 50000, 100000};

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