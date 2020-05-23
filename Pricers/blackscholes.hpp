#ifndef BLACKSCHOLES_HPP
#define BLACKSCHOLES_HPP

#include "../Option/option.hpp"
#include <iostream>

template <typename Type = double>
class BlackScholes {
    private:
        Option<Type>* option;
        Type d1;
        Type d2;
        Type NP(Type& x) {
            return (1.0/sqrt(2.0 * 3.1415)* exp(-x*x*0.5));
        }
        Type N(Type& x) {
            Type a1 = 0.319381530;
            Type a2 = -0.356563782;
            Type a3 = 1.781477937;
            Type a4 = -1.821255978;
            Type a5 = 1.330274429;
            Type k;
            k = 1/(1+0.2316419*x);
            if (x >= 0.0)
                return 1 - NP(x)*((a1*k) + (a2* k*k) + (a3*k*k*k) + (a4*k*k*k*k) + (a5*k*k*k*k*k));
            else
                return 1 - N(-x);
        }
    public:
        BlackScholes() {}
        BlackScholes(Option<Type>* option_): option(&option_) {}
        virtual ~BlackScholes() {}
        Type& run() {
            d1 = (log(option->get_spot()/option->get_strike()) + (option->get_rate()-option->get_dividend() + 0.5*(option->get_vol()*option->get_vol())) * option->get_expiry()) / (option->get_vol() * sqrt(option->get_expiry()));
            d2 = d1 - option->get_vol()*sqrt(option->get_expiry());
            switch (option->get_its_type()) {
                case 'c':
                    return option->get_spot()*exp(-option->get_dividend()*option->get_expiry())*N(d1) - option->get_strike()*exp(-option->get_rate()*option->get_expiry())*N(d2);
                case 'p':
                    return option->get_strike()*exp(-option->get_rate()*option->get_expiry())*N(-d2) - option->get_spot()*exp(-option->get_dividend()*option->get_expiry())*N(-d1);
                default:
                    throw("unknown option type")
            };
        }
};

#endif