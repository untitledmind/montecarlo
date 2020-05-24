#ifndef OPTION_HPP
#define OPTION_HPP

#include <cmath>
#include <algorithm>
#include "volatility.hpp"

template <typename Type> class Option;
template <typename Type> class VanillaOption;
template <typename Type> class DoubleDigital;
template <typename Type> class SquareRootPayoff;

template <typename Type = double>
class Option {
    private:
        const Type spot;
        const Type strike;
        const Type expiry;
        const Type dividend;
        Type rate;
        Type vol;
        Type price;
        VolatilityFunction<Type>* vol_func;
        const char its_type;
        const char its_style;
        void copy(const Option<Type>& rhs);
    public:
        Option();
        Option(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, const Type& vol_, const char& its_type_, const char& its_style_, const Type& dividend_);
        Option(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, VolatilityFunction<Type>& vol_, const char& its_type_, const char& its_style_, const Type& dividend_);
        Option(const Option& rhs);
        virtual ~Option();
        Type& operator= (const Option& rhs);
        virtual Type operator()(Type& spot) const = 0;

        void set_vol(Type& vol_);
        void set_price(Type& price_);
        const Type& get_spot() const;
        const Type& get_strike() const;
        const Type& get_expiry() const;
        const Type& get_dividend() const;
        const Type& get_rate() const;
        const Type& get_vol() const;
        VolatilityFunction<Type>& get_vol_func() const;
        const char& get_its_type() const;
        const char& get_its_style() const;
};

template <typename Type = double>
class VanillaOption: public Option<Type> {
    public:
        VanillaOption();
        VanillaOption(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, const Type& vol_, const char& its_type_, const char& its_style_, const Type& dividend_ = 0);
        VanillaOption(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, VolatilityFunction<Type>& vol_func_, const char& its_type_, const char& its_style_, const Type& dividend_ = 0);
        VanillaOption(const Option<Type>& rhs);
        virtual ~VanillaOption();
        virtual Type operator() (Type& spot) const;
};

template <typename Type = double>
class DoubleDigital: public Option<Type> {
    private:
        const Type strike2;
    public:
        DoubleDigital();
        DoubleDigital(const Type& spot_, const Type& strike_, const Type& strike2_, const Type& expiry_, const Type& rate_, const Type& vol_, const char& its_type_, const char& its_style_, const Type& dividend_ = 0);
        DoubleDigital(const Type& spot_, const Type& strike_, const Type& strike2_, const Type& expiry_, const Type& rate_, VolatilityFunction<Type>& vol_func_, const char& its_type_, const char& its_style_, const Type& dividend_ = 0);
        DoubleDigital(const Option<Type>& rhs);
        virtual ~DoubleDigital();
        virtual Type operator() (Type& spot) const;
        const Type& get_strike2() const;
};

template <typename Type = double>
class SquareRootPayoff: public Option<Type> {
    public:
        SquareRootPayoff();
        SquareRootPayoff(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, const Type& vol_, const char& its_type_, const char& its_style_, const Type& dividend_ = 0);
        SquareRootPayoff(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, VolatilityFunction<Type>& vol_func_, const char& its_type_, const char& its_style_, const Type& dividend_ = 0);
        SquareRootPayoff(const Option<Type>& rhs);
        virtual ~SquareRootPayoff();
        virtual Type operator() (Type& spot) const;
};

template <typename Type>
void Option<Type>::copy(const Option<Type>& rhs) {
    spot = rhs.spot;
    strike = rhs.strike;
    expiry = rhs.expiry;
    rate = rhs.rate;
    vol = rhs.vol;
    its_type = rhs.its_type;
}

template <typename Type>
Option<Type>::Option() {}

template <typename Type>
Option<Type>::Option(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, const Type& vol_, const char& its_type_, const char& its_style_, const Type& dividend_):
    spot(spot_), strike(strike_), expiry(expiry_), rate(rate_), vol(vol_), its_type(its_type_), its_style(its_style_), dividend(dividend_) {}

template <typename Type>
Option<Type>::Option(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, VolatilityFunction<Type>& _vol_func, const char& its_type_, const char& its_style_, const Type& dividend_):
    spot(spot_), strike(strike_), expiry(expiry_), rate(rate_), vol_func(&_vol_func), its_type(its_type_), its_style(its_style_), dividend(dividend_) {}

template <typename Type>
Option<Type>::Option(const Option& rhs): spot(rhs.spot), strike(rhs.strike), expiry(rhs.expiry), rate(rhs.rate), vol(rhs.vol), its_type(rhs.its_type), its_style(rhs.its_style) {}

template <typename Type>
Option<Type>::~Option() {}

template <typename Type>
Type& Option<Type>::operator= (const Option& rhs) {
    if (this != &rhs)
        copy(rhs);
    return *this;
}

template <typename Type>
void Option<Type>::set_vol(Type& vol_) {
    vol = vol_;
}

template <typename Type>
void Option<Type>::set_price(Type& price_) {
    price = price_;
}

template <typename Type>
const Type& Option<Type>::get_spot() const {
    return spot;
}

template <typename Type>
const Type& Option<Type>::get_strike() const {
    return strike;
}

template <typename Type>
const Type& Option<Type>::get_expiry() const {
    return expiry;
}

template <typename Type>
const Type& Option<Type>::get_dividend() const {
    return dividend;
}

template <typename Type>
const Type& Option<Type>::get_rate() const {
    return rate;
}

template <typename Type>
const Type& Option<Type>::get_vol() const {
    return vol;
}

template <typename Type>
VolatilityFunction<Type>& Option<Type>::get_vol_func() const {
    return *vol_func;
}

template <typename Type>
const char& Option<Type>::get_its_type() const {
    return its_type;
}

template <typename Type>
const char& Option<Type>::get_its_style() const {
    return its_style;
}

template <typename Type>
VanillaOption<Type>::VanillaOption(): Option<Type>::Option() {}

template <typename Type>
VanillaOption<Type>::VanillaOption(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, const Type& vol_, const char& its_type_, const char& its_style_, const Type& dividend_):
    Option<Type>::Option(spot_, strike_, expiry_, rate_, vol_, its_type_, its_style_, dividend_) {}

template <typename Type>
VanillaOption<Type>::VanillaOption(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, VolatilityFunction<Type>& _vol_func, const char& its_type_, const char& its_style_, const Type& dividend_):
    Option<Type>::Option(spot_, strike_, expiry_, rate_, _vol_func, its_type_, its_style_, dividend_) {}

template <typename Type>
VanillaOption<Type>::VanillaOption(const Option<Type>& rhs): Option<Type>::Option(rhs) {}

template <typename Type>
VanillaOption<Type>::~VanillaOption() {}

template <typename Type>
Type VanillaOption<Type>::operator()(Type& spot) const {
    // functor returns payoff based on option type
    switch (this->get_its_type()) {
        case 'c':
            return std::max(spot - this->get_strike(), 0.0);
        case 'p':
            return std::max(this->get_strike() - spot, 0.);
        default:
            throw("Unknown option type");
    };
}

template <typename Type>
DoubleDigital<Type>::DoubleDigital(): Option<Type>::Option() {}

template <typename Type>
DoubleDigital<Type>::DoubleDigital(const Type& spot_, const Type& strike_, const Type& strike2_, const Type& expiry_, const Type& rate_, const Type& vol_, const char& its_type_, const char& its_style_, const Type& dividend_):
    Option<Type>::Option(spot_, strike_, expiry_, rate_, vol_, its_type_, its_style_, dividend_), strike2(strike2_) {}

template <typename Type>
DoubleDigital<Type>::DoubleDigital(const Type& spot_, const Type& strike_, const Type& strike2_, const Type& expiry_, const Type& rate_, VolatilityFunction<Type>& _vol_func, const char& its_type_, const char& its_style_, const Type& dividend_):
    Option<Type>::Option(spot_, strike_, expiry_, rate_, _vol_func, its_type_, its_style_, dividend_), strike2(strike2_) {}

template <typename Type>
DoubleDigital<Type>::DoubleDigital(const Option<Type>& rhs): Option<Type>::Option(rhs) {}

template <typename Type>
DoubleDigital<Type>::~DoubleDigital() {}

template <typename Type>
Type DoubleDigital<Type>::operator()(Type& spot) const {
    // functor returns payoff based on option type
    switch (this->get_its_type()) {
        case 'p': {
            if (spot > this->get_strike() && spot < this->get_strike2())
                return std::max(spot - this->get_strike(), 0.);
            else if (spot > this->get_strike2())
                return this->get_strike2();
            return Type(0);
        }
        case 'c': {
            if (spot < this->get_strike())
                return std::max(this->get_strike() - spot, 0.);
            else if (spot < this->get_strike2())
                return this->get_strike2();
            return Type(0);
        }
        default:
            throw("Unknown option type");
    };
}

template <typename Type>
const Type& DoubleDigital<Type>::get_strike2() const {
    return strike2;
}

template <typename Type>
SquareRootPayoff<Type>::SquareRootPayoff(): Option<Type>::Option() {}

template <typename Type>
SquareRootPayoff<Type>::SquareRootPayoff(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, const Type& vol_, const char& its_type_, const char& its_style_, const Type& dividend_):
    Option<Type>::Option(spot_, strike_, expiry_, rate_, vol_, its_type_, its_style_, dividend_) {}

template <typename Type>
SquareRootPayoff<Type>::SquareRootPayoff(const Type& spot_, const Type& strike_, const Type& expiry_, const Type& rate_, VolatilityFunction<Type>& _vol_func, const char& its_type_, const char& its_style_, const Type& dividend_):
    Option<Type>::Option(spot_, strike_, expiry_, rate_, _vol_func, its_type_, its_style_, dividend_) {}

template <typename Type>
SquareRootPayoff<Type>::SquareRootPayoff(const Option<Type>& rhs): Option<Type>::Option(rhs) {}

template <typename Type>
SquareRootPayoff<Type>::~SquareRootPayoff() {}

template <typename Type>
Type SquareRootPayoff<Type>::operator()(Type& spot) const {
    switch (this->get_its_type()) {
        case 'c':
            return sqrt(std::max(spot - this->get_strike(), 0.));
        case 'p':
            return sqrt(std::max(this->get_strike() - spot, 0.));
        default:
            throw("Unknown option type");
    };
}

#endif