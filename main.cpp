#include<iostream>
#include "eval/cec22_functions_wrapper.h"
#include "heuristics/de.h"

int g_problem_size = 10;
int g_pop_size = (int)round(g_problem_size * 18);
double g_arc_rate = 2.6;
int g_memory_size = 6;
double g_p_best_rate = 0.11;
int g_function_number;
int g_restart_generation;
double g_min_region;
double g_max_region;
int g_number_of_used_vars = 10;
unsigned int g_max_num_evaluations;
double g_optimum[12] = {300, 400, 600, 800, 900, 1800, 2000, 2200, 2300, 2400, 2600, 2700};

int main(int argc, char **argv) {

    if (g_problem_size == 10)
        g_max_num_evaluations = 200000;
    else
        g_max_num_evaluations = 1000000;

    Cec22FunctionsWrapper *f = new Cec22FunctionsWrapper(1);

    std::cout << f->eval({ 0, 0 }) << endl;

    LSHADE lshade(f);

    lshade.run();

    return 0;
}