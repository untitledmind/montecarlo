#ifndef BINTREE_HPP
#define BINTREE_HPP

template <typename Type> class BinomialTree;

#include "../Option/option.hpp"
#include "../QSMatrix_EB/matrix.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

template <typename Type = double>
class BinomialTree {
    private:
        Option<Type>* option;
        Type num_steps;
        QSMatrix<Type> stock;
        QSMatrix<Type> oprices;
        std::vector<Type> stockvec;
        std::vector<Type> opricevec;
        Type value;
        Type dt;
        Type u;
        Type d;
        Type p;
        Type discount;
        Type ex_val;
        Type a;
        int counter;
    public:
        BinomialTree();
        BinomialTree(Option<Type>& _option, const int& _num_steps);
        virtual ~BinomialTree();
        Type get_value();
        Type run_with_matrix();
        Type run_with_matrix_improvement2();
        Type run_with_vector();
        void print() const;
};

template <typename Type>
BinomialTree<Type>::BinomialTree() {}

template <typename Type>
BinomialTree<Type>::BinomialTree(Option<Type>& _option, const int& _num_steps): option(&_option), num_steps(_num_steps) {
    dt = option->get_expiry() / num_steps;
    u = exp(option->get_vol() * sqrt(dt));
    d = 1 / u;
    discount = exp(-option->get_rate() * dt);
    p = (exp(option->get_rate() * dt) - d) / (u - d);
    a = log(option->get_strike()/(option->get_spot() * pow(d, num_steps))) / log(u/d);
}

template <typename Type>
BinomialTree<Type>::~BinomialTree() {}

template <typename Type>
Type BinomialTree<Type>::get_value() {
    return value;
}

template <typename Type>
Type BinomialTree<Type>::run_with_matrix() {
    counter = 0;
    QSMatrix<Type> empty(num_steps+1, num_steps+1, Type(0));
    oprices = empty;
    stock = empty;

    Type strike = option->get_strike();
    Type T = option->get_expiry();
    Type vol = option->get_vol();
    Type r = option->get_rate();
    char its_type = option->get_its_type();

    stock[0][0] = option->get_spot();

    for (int j=1; j <= num_steps; ++j) {
        stock[0][j] = stock[0][j-1] * u;
        ++counter;
    }
    
    for (int j=1; j <= num_steps; ++j)
    for (int i=1; i <= j; ++i) {
        stock[i][j] = stock[i-1][j-1] * d;
        ++counter;
    }

    for (int i=0; i <= num_steps; ++i) {
        oprices[i][num_steps] = (*option)(stock[i][num_steps]);
        ++counter;
    }

    for (int j = num_steps-1; j >= 0; --j)
    for (int i=0; i <= j; ++i) {
        oprices[i][j] = (p * oprices[i][j+1] + (1 - p)*oprices[i+1][j+1]) * discount;
        ex_val = (*option)(stock[i][j]);        // payoff functor defined in option.hpp
        ++counter;
        if (option->get_its_style() == 'a' && (ex_val > oprices[i][j])) {
            oprices[i][j] = ex_val;
            ++counter;
        }
    }
    value = oprices[0][0];
    return value;
}

template <typename Type>
Type BinomialTree<Type>::run_with_matrix_improvement2() {
    counter = 0;
    QSMatrix<Type> empty(num_steps+1, num_steps+1, Type(0));
    oprices = empty;
    stock = empty;

    Type spot = option->get_spot();
    Type strike = option->get_strike();
    Type T = option->get_expiry();
    Type vol = option->get_vol();
    Type r = option->get_rate();
    char its_type = option->get_its_type();

    stock[0][0] = spot;

    for (int j=1; j <= num_steps; ++j) {
        stock[0][j] = stock[0][j-1] * u;
        ++counter;
    }
    
    for (int j=1; j <= num_steps; ++j)
    for (int i=1; i <= j; ++i) {
        stock[i][j] = stock[i-1][j-1] * d;
        ++counter;
    }

    for (int i=0; i <= num_steps; ++i) {
        oprices[i][num_steps] = (*option)(stock[i][num_steps]);
        ++counter;
    }

    for (int j = num_steps-1; j >= 0; --j)
    for (int i=0; i <= fmin(num_steps - a, j); ++i) {
        oprices[i][j] = (p * oprices[i][j+1] + (1 - p)*oprices[i+1][j+1]) * discount;
        ex_val = (*option)(stock[i][j]);        // payoff functor defined in option.hpp
        ++counter;
        if (option->get_its_style() == 'a' && (ex_val > oprices[i][j])) {
            oprices[i][j] = ex_val;
            ++counter;
        }
    }
    value = oprices[0][0];
    return value;
}

template <typename Type>
Type BinomialTree<Type>::run_with_vector() {
    counter = 0;
    stockvec.resize(num_steps + 1);
    opricevec.resize(num_steps + 1);

    Type spot = option->get_spot();
    Type strike = option->get_strike();
    Type T = option->get_expiry();
    Type vol = option->get_vol();
    Type r = option->get_rate();
    char its_type = option->get_its_type();

    stockvec[0] = spot;
    for (int j=1; j <= num_steps; ++j) {
        stock[0][j] = stock[0][j-1] * u;
        ++counter;
    }
    
    for (int j=1; j <= num_steps; ++j)
    for (int i=1; i <= j; ++i) {
        stock[i][j] = stock[i-1][j-1] * d;
        ++counter;
    }

    for (int i=0; i <= num_steps; ++i) {
        oprices[i][num_steps] = (*option)(stock[i][num_steps]);
        ++counter;
    }

    for (int j = num_steps-1; j >= 0; --j)
    for (int i=0; i <= num_steps; ++i) {
        oprices[i][j] = (p * oprices[i][j+1] + (1 - p)*oprices[i+1][j+1]) * discount;
        ex_val = (*option)(stock[i][j]);        // payoff functor defined in option.hpp
        ++counter;
        if (option->get_its_style() == 'a' && (ex_val > oprices[i][j])) {
            oprices[i][j] = (*option)(stock[i][j]);
            ++counter;
        }
    }
    value = opricevec[0];
    return value;
}


template <typename Type>
void BinomialTree<Type>::print() const {
    std::vector<QSMatrix<Type> > trees = {oprices, stock};

    std::cout << "\nPrinting option and stock trees... mind the gap...\n";
    for (QSMatrix<Type> tree : trees) {
        std::cout << "      Step";
        std::cout.setf(ios::fixed);

        for (int j=0; j <= num_steps; ++j) {
            std::cout.width(10);
            std::cout << j;
        }

        std::cout << '\n';
        std::cout << "      Time";
        std::cout.precision(2);
        std::cout.setf(ios::fixed);
        std::cout.setf(ios::showpoint);

        for (int j=0; j <= num_steps; ++j) {
            std::cout.width(10);
            std::cout << j * dt;
        }

        std::cout << '\n';

        for (int i=0; i <= num_steps; ++i) {
            std::cout << '\n';
            for (int k=1; k <= i+1; ++k)
                cout << "          ";
            for (int j=i; j <= num_steps; ++j) {
                std::cout.width(10);
                std::cout << tree[i][j];
            }
        }
        std::cout << '\n';
    }
    std::cout << "\n========================================================\n";
    std::cout << "TRIAL RESULTS\n\n" << (option->get_its_type() == 'c' ? "Call" : "Put");
    std::cout << " option value = " << value << '\n';
    std::cout << "Number of steps = " << counter << '\n';
    std::cout << "u = " << u << '\n';
    std::cout << "d = " << d << '\n';
    std::cout << "K = " << option->get_strike() << '\n';
    std::cout << "a = " << a << '\n';

}

#endif