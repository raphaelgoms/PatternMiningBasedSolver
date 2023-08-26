#include<iostream>
#include "eval/cec22_functions_wrapper.h"

int main(int argc, char **argv) {

    Cec22FunctionsWrapper f(1);

    std::cout << f.eval({ 0, 0 }) << endl;

    return 0;
}