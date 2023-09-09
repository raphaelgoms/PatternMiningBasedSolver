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
int ObjectiveFunction<T>::getFunctionId()
{
    return this->function_id_;
    
}

template <class T>
void ObjectiveFunction<T>::setFunctionId(int function_id)
{
    this->function_id_ = function_id;
}

template <>
void ObjectiveFunction<double>::setFunctionId(int function_id)
{
    this->function_id_ = function_id;
}

template <class T>
int ObjectiveFunction<T>::getDimension()
{
    return this->dimension;
}

template <class T>
void ObjectiveFunction<T>::setDimension(int dimension)
{
    this->dimension = dimension;
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
