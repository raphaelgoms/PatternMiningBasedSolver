#ifndef Objective_Function_
#define Objective_Function_

#include<vector>
#include<string>

using namespace std;
typedef double RealType;

template <class T>
class ObjectiveFunction
{
protected:
    string function_name;
    T optimum;

public:
    /// @brief 
    /// @param x vector to be evaluted by this function (represents a solution)
    /// @param vsize size of the solution vector
    /// @return The cost of the function to x vector 
    virtual RealType eval(T *x, int vsize) = 0;

    /// @brief 
    /// @param x vector to be evaluted by this function (represents a solution)
    /// @return The cost of the function to x vector
    virtual RealType eval(vector<T> const &x) = 0;
    
    /// @brief 
    /// @return the function name
    string getFunctionName();

    /// @brief 
    /// @param fname the function name
    void setFunctionName(string function_name);

    /// @brief 
    /// @return the best value to this function 
    RealType getOptimum();

    /// @brief 
    /// @param optimum the best possible value to this function 
    void setOptimum(RealType optimum);
};

#endif