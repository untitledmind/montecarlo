#ifndef OPTION_HPP
#define OPTION_HPP

#include <cmath>
#include <algorithm>

#include "volatility.hpp"

template <typename Type> class Option;
template <typename Type> class VanillaOption;
template <typename Type> class DoubleDigital;

template <typename Type = double>
class Option {
    private:
        const Type spot;
        const Type strike;
        const Type expiry;
        Type rate;
        Type vol;
        const VolatilityFunction<Type>* vol_func;
        const char its_type;
        void copy(const Option<Type>& rhs);
    public:
        Option();
        Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const Type& _vol, char& _its_type);
        Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const VolatilityFunction<Type>& _vol, char& _its_type);
        Option(const Option& rhs);
        virtual ~Option();
        Type& operator= (const Option& rhs);
        virtual Type operator()(Type& spot) const = 0;

        Type set_vol(const Type& _vol);
        Type get_spot() const;
        Type get_strike() const;
        Type get_expiry() const;
        Type get_rate() const;
        Type get_vol() const;
        const VolatilityFunction<Type>* get_vol_func() {return vol_func;}
        char get_its_type() const;
};

template <typename Type = double>
class VanillaOption: public Option<Type> {
    public:
        VanillaOption();
        VanillaOption(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const Type& _vol, char& _its_type);
        VanillaOption(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const VolatilityFunction<Type>& _vol, char& _its_type);
        VanillaOption(const Option<Type>& rhs);
        virtual ~VanillaOption();
        virtual Type operator() (Type& spot) const;
};

template <typename Type = double>
class DoubleDigital: public Option<Type> {
    public:
        DoubleDigital();
        DoubleDigital(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const Type& _vol, char& _its_type);
        DoubleDigital(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const VolatilityFunction<Type>& _vol, char& _its_type);
        DoubleDigital(const Option<Type>& rhs);
        virtual ~DoubleDigital();
        virtual Type operator() (Type& spot) const;
};

template <typename Type>
Option<Type>::Option() {}

template <typename Type>
Option<Type>::Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const Type& _vol, char& _its_type):
    spot(_spot), strike(_strike), expiry(_expiry), rate(_rate), vol(_vol), its_type(_its_type) {}

template <typename Type>
Option<Type>::Option(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const VolatilityFunction<Type>& _vol_func, char& _its_type):
    spot(_spot), strike(_strike), expiry(_expiry), rate(_rate), vol_func(&_vol_func), its_type(_its_type) {}

template <typename Type>
Option<Type>::Option(const Option& rhs): spot(rhs.spot), strike(rhs.strike), expiry(rhs.expiry), rate(rhs.rate), vol(rhs.vol), its_type(rhs.its_type) {}

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
Type Option<Type>::set_vol(const Type& _vol) {
    vol = _vol;
}

template <typename Type>
Type Option<Type>::get_spot() const {
    return spot;
}

template <typename Type>
Type Option<Type>::get_strike() const {
    return strike;
}

template <typename Type>
Type Option<Type>::get_expiry() const {
    return expiry;
}

template <typename Type>
Type Option<Type>::get_rate() const {
    return rate;
}

template <typename Type>
Type Option<Type>::get_vol() const {
    return vol;
}

template <typename Type>
char Option<Type>::get_its_type() const {
    return its_type;
}

template <typename Type>
VanillaOption<Type>::VanillaOption(): Option<Type>::Option() {}

template <typename Type>
VanillaOption<Type>::VanillaOption(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const Type& _vol, char& _its_type):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol, _its_type) {}

template <typename Type>
VanillaOption<Type>::VanillaOption(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const VolatilityFunction<Type>& _vol_func, char& _its_type):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol_func, _its_type) {}

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
DoubleDigital<Type>::DoubleDigital(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const Type& _vol, char& _its_type):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol, _its_type) {}

template <typename Type>
DoubleDigital<Type>::DoubleDigital(const Type& _spot, const Type& _strike, const Type& _expiry, const Type& _rate, const VolatilityFunction<Type>& _vol_func, char& _its_type):
    Option<Type>::Option(_spot, _strike, _expiry, _rate, _vol_func, _its_type) {}

template <typename Type>
DoubleDigital<Type>::DoubleDigital(const Option<Type>& rhs): Option<Type>::Option(rhs) {}

template <typename Type>
DoubleDigital<Type>::~DoubleDigital() {}

template <typename Type>
Type DoubleDigital<Type>::operator()(Type& spot) const {
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

#endif