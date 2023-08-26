#include "objective_function.h"

template <class T>
string ObjectiveFunction<T>::getFunctionName()
{
    return this->function_name;
}

template <class T>
void ObjectiveFunction<T>::setFunctionName(string function_name)
{
    this->function_name = function_name;
}

template <class T>
RealType ObjectiveFunction<T>::getOptimum()
{
    return this->optimum;
}

template <class T>
void ObjectiveFunction<T>::setOptimum(RealType optimum)
{
    this->optimum = optimum;
}
