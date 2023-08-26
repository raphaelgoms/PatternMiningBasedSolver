#include "cec22_functions_wrapper.h"
#include "cec22_test_func.cpp"

#include <cassert>

Cec22FunctionsWrapper::Cec22FunctionsWrapper(int function_id) 
    : ObjectiveFunction<double>() 
{ 
    this->function_id_ = function_id;
}

RealType Cec22FunctionsWrapper::eval(double *x, int vsize) 
{
    assert(this->function_id_ > 0);

    RealType result;
    cec22_test_func(x, &result, vsize, 1, this->function_id_);

    return result;
}

RealType Cec22FunctionsWrapper::eval(vector<double> const &x)
{
    assert(this->function_id_ > 0);

    RealType result;
    cec22_test_func(&vector<double>(x)[0], &result, x.size(), 1, this->function_id_);

    return 0.0;
}


  