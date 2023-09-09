
#include "objective_function.h"

class Cec22FunctionsWrapper : public ObjectiveFunction<double> {
    //int function_id_ = -1;

public:
    Cec22FunctionsWrapper(int function_id);

    int getFunctionId();
    
    RealType eval(double *x, int vsize) override;
    RealType eval(vector<double> const &x) override;
};