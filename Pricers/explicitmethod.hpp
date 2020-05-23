#ifndef EXPLICITMETHOD_HPP
#define EXPLICITMETHOD_HPP

template <typename Type> class ExplicitMethod;

#include "../Option/option.hpp"
#include <iostream>
#include <vector>

template <typename Type = double>
class ExplicitMethod {
    private:
        Option<Type>* option;
        Type sigma_sqr;
        Type delta_S;
        std::vector<Type> S_values;
        Type delta_t;
        std::vector<Type> a;
        std::vector<Type> b;
        std::vector<Type> c;
        std::vector<Type> f;
        std::vector<Type> f_next;
        Type r1;
        Type r2;
        Type value;
        const int& S_steps;      //granulatiry of spot price
        const int& t_steps;      //granulatity of time
    public:
        ExplicitMethod();
        ExplicitMethod(Option<Type>& option_, const int& S_steps_, const int& t_steps_);
        virtual ~ExplicitMethod();
        
        Type run();
        void print();
};

template <typename Type>
ExplicitMethod<Type>::ExplicitMethod() {}

template <typename Type>
ExplicitMethod<Type>::ExplicitMethod(Option<Type>& option_, const int& S_steps_, const int& t_steps_): option(&option_), S_steps(S_steps_), t_steps(t_steps_) {
    int M = S_steps;
    int N = t_steps;

    sigma_sqr = option->get_vol() * option->get_vol();
    delta_S = (2.0 * option->get_spot()) / M;
    delta_t = option->get_expiry()/N;
    r1 = 1.0 / (1.0 + option->get_rate()*delta_t);
    r2 = delta_t / (1.0 + option->get_rate()*delta_t);

    S_values.reserve(M + 1);
    a.reserve(M);
    b.reserve(M);
    c.reserve(M);
    f.reserve(M+1);
    f_next.reserve(M+1);
}

template <typename Type>
ExplicitMethod<Type>::~ExplicitMethod() {}

template <typename Type>
Type ExplicitMethod<Type>::run() {
    int M = S_steps;
    int N = t_steps;

    for (unsigned m=0; m <= M; ++m)
        S_values[m] = m * delta_S;
    
    for (unsigned int j=1; j < M; ++j) {
        a[j] = r2 * 0.5*j * (-option->get_rate() + sigma_sqr*j);
        b[j] = r1 * (1.0 - sigma_sqr*j*j*delta_t);
        c[j] = r2 * 0.5*j * (option->get_rate() + sigma_sqr*j);
    }

    for (unsigned m=0; m <= M; ++m)
        f_next[m] = (*option)(S_values[m]);    // payoff at S_values[m]
    
    for (int t = N-1; t >= 0; --t) {
        f[0] = (*option)(S_values[0]);      // modification - functor payoff

        for (unsigned m=1; m < M; ++m)
            f[m] = a[m]*f_next[m-1] + b[m]*f_next[m] + c[m]*f_next[m+1]; // weighted average

        f[M] = (*option)(S_values[M]);      // modification - functor payoff

        for (unsigned m=0; m <= M; ++m)
            f_next[m] = f[m];       // update old price to present price 
    }

    value = f[M/2];
    return value;
}

template <typename Type>
void ExplicitMethod<Type>::print() {
    int M = S_steps;
    int N = t_steps;
    printf("M = %d N = %d\n", M,N);

    std::cout << " ---- S_Values(@ t = T) ----\n";
    for (unsigned m=0; m <= M; ++m) {
        S_values[m] = m * delta_S;
        if (!(m % 10))
            std::cout << '\n';
        printf("%5.2f\t", S_values[m]);
    }

    printf("\ndelta_S = %5.3f delta_t = %5.3f option = %c\n", delta_S, delta_t, option->get_its_type());
    printf("delta_t/(delta_S*delta_S) = %5.3f\n\n",delta_t/(delta_S*delta_S) );

    for (int t = N-1; t >= 0; --t) {
        printf("t = %d  time = %5.3f\n", t, t*delta_t);

        
        for (int i=0; i < M; ++i) {
            if (!(i % 10))
                printf("\nPutPrice[%2d] = ", i);
            printf("%5.2f\t", f[i]);
        }
        printf("\n-----------------\n");
    }

    printf("delta_S*M/2 = %5.3f\n", delta_S*M/2);
    std::cout << "\n Option value = " << value << '\n';
}

#endif