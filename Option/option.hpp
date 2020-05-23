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
        const Type rate;
        Type vol;
        Type price;
        VolatilityFunction<Type>* vol_func;
        const char its_type;
        const char its_style;
        void copy(const Option<Type>& rhs);
    public:
        Option();
        Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, Type& _vol, char& _its_type, char& _its_style);
        Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, VolatilityFunction<Type>& _vol, char& _its_type, char& _its_style);
        Option(const Option& rhs);
        virtual ~Option();
        Type& operator= (const Option& rhs);
        virtual Type operator()(Type& spot) const = 0;

        void set_vol(Type& _vol);
        const Type& get_spot() const;
        const Type& get_strike() const;
        const Type& get_expiry() const;
        const Type& get_rate() const;
        const Type& get_vol() const;
        VolatilityFunction<Type>* get_vol_func();
        const char& get_its_type() const;
        const char& get_its_style() const;
};

template <typename Type = double>
class VanillaOption: public Option<Type> {
    public:
        VanillaOption();
        VanillaOption(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, Type& _vol, char& _its_type, char& _its_style);
        VanillaOption(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, VolatilityFunction<Type>& _vol, char& _its_type, char& _its_style);
        VanillaOption(const Option<Type>& rhs);
        virtual ~VanillaOption();
        virtual Type operator() (Type& spot) const;
};

template <typename Type = double>
class DoubleDigital: public Option<Type> {
    private:
        Type strike2;
    public:
        DoubleDigital();
        DoubleDigital(const Type& _spot, const Type& _strike, const Type& _strike2, const Type& _expiry, const Type& _rate, Type& _vol, char& _its_type, char& _its_style);
        DoubleDigital(const Type& _spot, const Type& _strike, const Type& _strike2, const Type& _expiry, const Type& _rate, VolatilityFunction<Type>& _vol, char& _its_type, char& _its_style);
        DoubleDigital(const Option<Type>& rhs);
        virtual ~DoubleDigital();
        virtual Type operator() (Type& spot) const;
        Type get_strike2() const {return strike2;}
};

template <typename Type = double>
class SquareRootPayoff: public Option<Type> {
    public:
        SquareRootPayoff();
        SquareRootPayoff(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, Type& _vol, char& _its_type, char& _its_style);
        SquareRootPayoff(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, VolatilityFunction<Type>& _vol, char& _its_type, char& _its_style);
        SquareRootPayoff(const Option<Type>& rhs);
        virtual ~SquareRootPayoff();
        virtual Type operator() (Type& spot) const;
};

template <typename Type>
Option<Type>::Option() {}

template <typename Type>
Option<Type>::Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, Type& _vol, char& _its_type, char& _its_style):
    spot(_spot), strike(_strike), expiry(_expiry), rate(_rate), vol(_vol), its_type(_its_type), its_style(_its_style) {}

template <typename Type>
Option<Type>::Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, VolatilityFunction<Type>& _vol_func, char& _its_type, char& _its_style):
    spot(_spot), strike(_strike), expiry(_expiry), rate(_rate), vol_func(&_vol_func), its_type(_its_type), its_style(_its_style) {}

template <typename Type>
Option<Type>::Option(const Option& rhs): spot(rhs.spot), strike(rhs.strike), expiry(rhs.expiry), rate(rhs.rate), vol(rhs.vol), its_type(rhs.its_type), its_style(rhs.its_style) {}

template <typename Type>
Option<Type>::~Option() {}

template <typename Type>
Type& Option<Type>::operator= (const Option& rhs) {
    if (this != &rhs) {
        spot = rhs.spot;
        strike = rhs.strike;
        expiry = rhs.expiry;
        rate = rhs.rate;
        vol = rhs.vol;
        its_type = rhs.its_type;
    }
    return *this;
}

template <typename Type>
void Option<Type>::set_vol(Type& _vol) {
    vol = _vol;
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
const Type& Option<Type>::get_rate() const {
    return rate;
}

template <typename Type>
const Type& Option<Type>::get_vol() const {
    return vol;
}

template <typename Type>
VolatilityFunction<Type>* Option<Type>::get_vol_func() {
    return vol_func;
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
VanillaOption<Type>::VanillaOption(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, Type& _vol, char& _its_type, char& _its_style):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol, _its_type, _its_style) {}

template <typename Type>
VanillaOption<Type>::VanillaOption(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, VolatilityFunction<Type>& _vol_func, char& _its_type, char& _its_style):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol_func, _its_type, _its_style) {}

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
DoubleDigital<Type>::DoubleDigital(const Type& _spot, const Type& _strike, const Type& _strike2, const Type& _expiry, const Type& _rate, Type& _vol, char& _its_type, char& _its_style/*, const Type& _strike2*/):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol, _its_type, _its_style), strike2(_strike2) {}

template <typename Type>
DoubleDigital<Type>::DoubleDigital(const Type& _spot, const Type& _strike, const Type& _strike2, const Type& _expiry, const Type& _rate, VolatilityFunction<Type>& _vol_func, char& _its_type, char& _its_style/*, const Type& _strike2*/):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol_func, _its_type, _its_style), strike2(_strike2) {}

template <typename Type>
DoubleDigital<Type>::DoubleDigital(const Option<Type>& rhs): Option<Type>::Option(rhs) {}

template <typename Type>
DoubleDigital<Type>::~DoubleDigital() {}

template <typename Type>
Type DoubleDigital<Type>::operator()(Type& spot) const {
    // functor returns payoff based on option type
    switch (this->get_its_type()) {
        case 'c':
            if (spot > this->get_strike() && spot < this->get_strike2())
                return std::max(spot - this->get_strike(), 0.);
            else if (spot > this->get_strike2())
                return this->get_strike2();
            return Type(0);
        case 'p':
            if (spot < this->get_strike())
                return std::max(this->get_strike() - spot, 0.);
            else if (spot < this->get_strike2())
                return this->get_strike2();
            return Type(0);
        default:
            throw("Unknown option type");
    };

    // if (spot <= this->get_strike())
    //     return 0.;
    // else if (spot >= this->get_strike2())
    //     return 0.;
    // return 1.;
}

template <typename Type>
SquareRootPayoff<Type>::SquareRootPayoff(): Option<Type>::Option() {}

template <typename Type>
SquareRootPayoff<Type>::SquareRootPayoff(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, Type& _vol, char& _its_type, char& _its_style):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol, _its_type, _its_style) {}

template <typename Type>
SquareRootPayoff<Type>::SquareRootPayoff(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, VolatilityFunction<Type>& _vol_func, char& _its_type, char& _its_style):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol_func, _its_type, _its_style) {}

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