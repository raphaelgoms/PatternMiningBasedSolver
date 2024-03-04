/*
  L-SHADE implemented by C++ for Special Session & Competition on Real-Parameter Single Objective Optimization at CEC-2014

  Version: 1.0   Date: 16/Apr/2014
  Written by Ryoji Tanabe (rt.ryoji.tanabe [at] gmail.com)
*/

#ifndef _HEADER_H_
#define _HEADER_H_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <math.h>
#include <map>
#include <set>
#include <tuple>
#include <fstream>
#include <random>
#include <ctime>

#include "../eval/objective_function.h"
#include "../mining/pattern_miner.h" 

//#define M_PI 3.1415926535897932384626433832795029

using namespace std;

typedef  double variable;
typedef variable *Individual;
typedef  double Fitness;

typedef map<int, tuple<double, double>> pattern;

extern int g_function_number;
extern int g_problem_size;
extern unsigned int g_max_num_evaluations;
extern int g_restart_generation;
extern int function_name;

extern int g_pop_size;
extern int g_memory_size;
extern double g_p_best_rate;
extern double g_arc_rate;
extern Fitness g_optimum[12];

extern double g_min_region;
extern double g_max_region;

//void cec22_test_func(double *, double * ,int,int,int);
//void cec21_basic_func(double *, double * ,int,int,int);

void fprintPopulation(vector <Individual> pop, int generation, string basepath);


class searchAlgorithm {
public:
  searchAlgorithm();
  searchAlgorithm(ObjectiveFunction<double> *f);
  virtual Fitness run() = 0;

  int used_gen_count = 0;
  int function_calls = 0;
  double success_rate;
protected:
  ObjectiveFunction<double> *f;

  Fitness evaluateIndividual(Individual individual);
  void evaluatePopulation(const vector<Individual> &pop, vector<Fitness> &fitness);
  void initializeFitnessFunctionParameters();

  void initializeParameters();
  Individual makeNewIndividual();
  void modifySolutionWithParentMedium(Individual child, Individual parent);
  void setBestSolution(const vector<Individual> &pop, const vector<Fitness> &fitness, Individual &bsf_solution, Fitness &bsf_fitness);

  tuple<double, double> lineSearch(Individual current_position, double h, int axis, bool stop_in_first_improve = false);
  bool constructGreedyRandomized(Individual& current_position, set<int> unfixed_positions, double h = 1e-1, double alpha = 0);

  //Return random value with uniform distribution [0, 1)
  inline double randDouble() {
    return (double)rand() / (double) RAND_MAX;
  }

  /*
    Return random value from Cauchy distribution with mean "mu" and variance "gamma"
    http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Cauchy
  */
  inline double cauchy_g(double mu, double gamma) {
    return mu + gamma * tan(M_PI*(randDouble() - 0.5));
  }

  /*
    Return random value from normal distribution with mean "mu" and variance "gamma"
    http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Gauss
  */
  inline double gauss(double mu, double sigma){
    return mu + sigma * sqrt(-2.0 * log(randDouble())) * sin(2.0 * M_PI * randDouble());
  }

  //Recursive quick sort with index array
  template<class VarType>
    void sortIndexWithQuickSort(VarType array[], int first, int last, int index[]) {
    VarType x = array[(first + last) / 2];
    int i = first;
    int j = last;
    VarType temp_var = 0;
    int temp_num = 0;

    while (true) {
      while (array[i] < x) i++;
      while (x < array[j]) j--;
      if (i >= j) break;

      temp_var = array[i];
      array[i] = array[j];
      array[j] = temp_var;

      temp_num = index[i];
      index[i] = index[j];
      index[j] = temp_num;

      i++;
      j--;
    }

    if (first < (i -1)) sortIndexWithQuickSort(array, first, i - 1, index);
    if ((j + 1) < last) sortIndexWithQuickSort(array, j + 1, last, index);
  }

  template<typename T>
  T radomlySelectElement(std::set<T> const& v)
  {
    auto it = v.cbegin();
    int random = rand() % v.size();
    std::advance(it, random);

    return *it;
  }

  int function_number;
  int problem_size;
  variable max_region;
  variable min_region;
  vector<variable> lower_bounds;
  vector<variable> upper_bounds;
  Fitness optimum;
  // acceptable error value
  Fitness epsilon;
  unsigned int max_num_evaluations;
  int pop_size;
};

class LSHADE: public searchAlgorithm {
public:
  LSHADE(ObjectiveFunction<double> *f);
  virtual Fitness run();
  void setSHADEParameters();
  void reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness);
  void operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, int &target, int &p_best_individual, variable &scaling_factor, variable &cross_rate, const vector<Individual> &archive, int &arc_ind_count);

  int arc_size;
  double arc_rate;
  variable p_best_rate;
  int memory_size;
  int reduction_ind_num;
};

class DM_LSHADE: public searchAlgorithm {
public:
  enum PatternUsageStrategy {
    MUTATION_SCHEME,
    INSERT_IN_POP
  };

  enum PatternSelectionStrategy {
    RANDOM,
    QUEUE
  };

  enum SolutionFillingStrategy {
    P_BEST,
    MOMENTUM,
    LINE_SEARCH
  };
    
  enum EliteType {
    BY_GENERATION,
    CROSS_GENERATION
  };

  DM_LSHADE(ObjectiveFunction<double> *f,
    PatternSelectionStrategy pattern_sel_strategy = RANDOM, 
    SolutionFillingStrategy filling_strategy = P_BEST,
    PatternUsageStrategy pattern_usage_strategy = MUTATION_SCHEME,
    EliteType elite_type = BY_GENERATION);

  virtual Fitness run();
  void setSHADEParameters();
  void reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness);
  void operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, int &target, int &p_best_individual, variable &scaling_factor, variable &cross_rate, const vector<Individual> &archive, int &arc_ind_count);
  void operateCurrentToPBest1BinWithArchiveAndXPattern(const vector<Individual>& pop, Individual child, int& target, int& p_best_individual, 
        variable& scaling_factor, variable& cross_rate, const vector<Individual>& archive, int& arc_ind_count, map<int,double> pattern);
  void updateElite(vector<Individual> curr_pop, int* sorted_indexes, double *fitness);

  int arc_size;
  double arc_rate;
  variable p_best_rate;
  int memory_size;
  int reduction_ind_num;

  // new atributes:
  int number_of_patterns = 0;
  string mining_algorithm = "";
  float clusters_count = 0;
  PatternUsageStrategy pattern_usage_strategy = MUTATION_SCHEME;
  PatternSelectionStrategy pattern_sel_strategy = RANDOM;
  SolutionFillingStrategy filling_strategy = P_BEST; 
  EliteType elite_type = BY_GENERATION;
  vector<tuple<Individual, double>> elite;
  int elite_max_size = 30;
};

#endif


