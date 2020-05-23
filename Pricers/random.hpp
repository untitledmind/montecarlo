#ifndef RANDOM_HPP
#define RANDOM_HPP

template <typename Type> class RandomGenerator;
template <typename Type> class BoxMuller;
template <typename Type> class Antithetic;

#include <cmath>
#include <string>

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
        const std::string signature = "Box Muller generator";
    public:
        BoxMuller() {}
        virtual ~BoxMuller() {}
        virtual Type generate() {
            Type x = 0.;
            Type y = 0.;
            Type euclid = 0.;

            do {
                x = 2. * rand() / static_cast<Type>(RAND_MAX) - 1;
                y = 2. * rand() / static_cast<Type>(RAND_MAX) - 1;
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
        const std::string signature = "Antithetic generator";
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

#endif