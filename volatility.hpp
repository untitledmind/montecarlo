#ifndef VOLATILITY_HH
#define VOLATILITY_HH

#include <vector>
#include <utility>
#include <initializer_list>
#include <iostream>
#include <iomanip>

template <typename Type> class VolatilityFunction;
template <typename Type> class StepwiseVolatility;

template <typename Type = double>
class VolatilityFunction {
    public:
        VolatilityFunction() {}
        virtual ~VolatilityFunction() {}
        virtual int get_steps() const = 0;
        virtual std::vector<std::pair<Type,Type> > get_function() const = 0;
        virtual std::vector<std::pair<Type,Type> > get_dt() const = 0;
        virtual std::vector<Type> get_dt_vec() const = 0;
        virtual std::vector<Type> get_vol_vec() const = 0;       
        virtual Type operator() (const Type& t) const = 0;
        virtual std::pair<Type,Type> operator[] (const unsigned int& _step) const = 0;
};

template <typename Type = double>
class StepwiseVolatility: public VolatilityFunction<Type> {
    private:
        int steps;
        std::vector<std::pair<Type,Type> > vol_func;
        std::vector<std::pair<Type,Type> > dt_func;
        std::vector<Type> dt_vec;
        std::vector<Type> vol_vec;
        void copy(const StepwiseVolatility<Type>& rhs);
    public:
        StepwiseVolatility();
        StepwiseVolatility(std::vector<std::pair<Type,Type> >& _vols);
        StepwiseVolatility(std::initializer_list<Type>& _vols);
        StepwiseVolatility(const StepwiseVolatility<Type>& rhs);
        virtual ~StepwiseVolatility();

        int get_steps() const;
        std::vector<std::pair<Type,Type> > get_function() const;
        std::vector<std::pair<Type,Type> > get_dt() const;
        std::vector<Type> get_dt_vec() const;
        std::vector<Type> get_vol_vec() const;

        StepwiseVolatility<Type>& operator= (const StepwiseVolatility<Type>& rhs);
        StepwiseVolatility<Type>& operator= (std::vector<std::pair<Type,Type> >& _vols);
        virtual Type operator() (const Type& t) const;
        virtual std::pair<Type,Type> operator[] (const unsigned int& _step) const;

        template <typename U>
        friend std::ostream& operator<< (std::ostream& os, const StepwiseVolatility<U>& rhs);
};

template <typename Type>
void StepwiseVolatility<Type>::copy(const StepwiseVolatility<Type>& rhs) {
    steps = rhs.steps;
    vol_func = rhs.vol_func;
}

template <typename Type>
StepwiseVolatility<Type>::StepwiseVolatility() {}

template <typename Type>
StepwiseVolatility<Type>::StepwiseVolatility(std::vector<std::pair<Type,Type> >& _vols):
    vol_func(_vols), steps(_vols.size()) {
            // calc change in time and respective volatility
        for (int i=0; i < steps-1; ++i) {
            Type dt = std::get<0>(vol_func[i+1]) - std::get<0>(vol_func[i]);
            Type vol = std::get<1>(vol_func[i]);
            dt_func.push_back( {dt,vol} );
            dt_vec.push_back(dt);
            vol_vec.push_back(vol);
        }
    }

template <typename Type>
StepwiseVolatility<Type>::StepwiseVolatility(std::initializer_list<Type>& _vols):
    vol_func(_vols) {}

template <typename Type>
StepwiseVolatility<Type>::StepwiseVolatility(const StepwiseVolatility<Type>& rhs) {
    if (*this != rhs)
        copy(rhs);
    return *this;
}

template <typename Type>
StepwiseVolatility<Type>::~StepwiseVolatility() {}

template <typename Type>
int StepwiseVolatility<Type>::get_steps() const {
    return steps;
}

template <typename Type>
std::vector<std::pair<Type,Type> > StepwiseVolatility<Type>::get_function() const {
    return vol_func;
}

template <typename Type>
std::vector<std::pair<Type,Type> > StepwiseVolatility<Type>::get_dt() const {
    return dt_func;
}

template <typename Type>
std::vector<Type> StepwiseVolatility<Type>::get_dt_vec() const {
    return dt_vec;
}

template <typename Type>
std::vector<Type> StepwiseVolatility<Type>::get_vol_vec() const {
    return vol_vec;
}


template <typename Type>
StepwiseVolatility<Type>& StepwiseVolatility<Type>::operator= (const StepwiseVolatility<Type>& rhs) {
    if (*this != rhs)
        copy(rhs);
    return *this;
}

template <typename Type>
StepwiseVolatility<Type>& StepwiseVolatility<Type>::operator= (std::vector<std::pair<Type,Type> >& _vols) {
    vol_func = _vols;
    steps = _vols.size();
}

template <typename Type>
Type StepwiseVolatility<Type>::operator() (const Type& t) const {
    for (int i=0; i < steps-1; ++i) {
        std::pair<Type,Type> point;
        if (t < std::get<0>(vol_func[i+1]))
            return std::get<1>(vol_func[i]);
    }
    throw(std::runtime_error("volatility does not exist for that t"));
}

template <typename Type>
std::pair<Type,Type> StepwiseVolatility<Type>::operator[] (const unsigned int& _step) const {
    // used to access dt elements - returns volatility and respective dt (width)
    return dt_func[_step];
}

template <typename U>
std::ostream& operator<< (std::ostream& os, const StepwiseVolatility<U>& rhs) {
    os << "Time" << '\t' << std::setw(6) << "Vol\n";
    for (int i=0; i < rhs.steps; ++i) {
        os << std::get<0>(rhs.vol_func[i]) << '\t' << std::setw(5) << std::get<1>(rhs.vol_func[i]) << '\n';
    }
    return os;
}

#endif
