#include "definitions.h"

#ifdef CEC_2014_BENCHMARK
#include "objective_function.h"

class Cec14FunctionsWrapper : public ObjectiveFunction<double> {
    //int function_id_ = -1;

public:
    Cec14FunctionsWrapper(int function_id);

    int getFunctionId();
    
    RealType eval(double *x, int vsize) override;
    RealType eval(vector<double> const &x) override;
};
#endif