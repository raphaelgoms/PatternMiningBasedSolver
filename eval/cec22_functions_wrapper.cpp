#include "cec22_functions_wrapper.h"

#ifdef CEC_2022_BENCHMARK

#include "objective_function.h"
#include "cec22_test_func.c"

#include<iostream>

#include <cassert>

double *OShift,*M,*y,*z,*x_bound, *Rn;
int ini_flag=0,n_flag,func_flag,*SS;
FILE *fpt;


Cec22FunctionsWrapper::Cec22FunctionsWrapper(int function_id) 
    : ObjectiveFunction<double>() 
{ 
    this->function_id_ = function_id;
    //this->getFunctionId();
    
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

    return result;
}

#endif  