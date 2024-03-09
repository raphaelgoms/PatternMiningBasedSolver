
#include "de.h"


#include <filesystem>
#include <iostream>

using namespace alglib;

using namespace std;

DM_LSHADE::DM_LSHADE(ObjectiveFunction<double> *f, 
  PatternSelectionStrategy pattern_sel_strategy, 
  SolutionFillingStrategy filling_strategy,
  PatternUsageStrategy pattern_usage_strategy,
  EliteType elite_type)
{
  this->f = f;
  this->pattern_sel_strategy = pattern_sel_strategy;
  this->filling_strategy = filling_strategy;
  this->pattern_usage_strategy = pattern_usage_strategy;
  this->elite_type = elite_type;
}

Fitness DM_LSHADE::run()
{
  this->function_calls = 0;
  this->used_gen_count = 0;
  // cout << scientific << setprecision(8);
  initializeParameters();
  setSHADEParameters();

  // cout << pop_size << endl;
  // cout << arc_size << endl;
  // cout << p_best_rate << endl;
  // cout << memory_size << endl;

  bool print = false;
  int generation = 1;

  vector<Individual> pop;
  vector<Fitness> fitness(pop_size, 0);

  vector<Individual> children;
  vector<Fitness> children_fitness(pop_size, 0);

  // initialize population
  for (int i = 0; i < pop_size; i++)
  {
    pop.push_back(makeNewIndividual());
    children.push_back((variable *)malloc(sizeof(variable) * problem_size));
  }

  // evaluate the initial population's fitness values
  evaluatePopulation(pop, fitness);

  Individual bsf_solution = (variable *)malloc(sizeof(variable) * problem_size);
  Fitness bsf_fitness;
  int nfes = 0;

  if ((fitness[0] - optimum) < epsilon)
    fitness[0] = optimum;
  bsf_fitness = fitness[0];
  for (int j = 0; j < problem_size; j++)
    bsf_solution[j] = pop[0][j];
  /////////////////////////////////////////////////////////////////////////
  for (int i = 0; i < pop_size; i++)
  {
    nfes++;

    if ((fitness[i] - optimum) < epsilon)
      fitness[i] = optimum;

    if (fitness[i] < bsf_fitness)
    {
      bsf_fitness = fitness[i];
      for (int j = 0; j < problem_size; j++)
        bsf_solution[j] = pop[i][j];
    }

    // if (nfes % 1000 == 0) {
    //   //      cout << nfes << " " << bsf_fitness - optimum << endl;
    //   cout << bsf_fitness - optimum << endl;
    // }

    if (nfes >= max_num_evaluations)
      break;
  }
  ////////////////////////////////////////////////////////////////////////////

  // for external archive
  int arc_ind_count = 0;
  int random_selected_arc_ind;
  vector<Individual> archive;
  for (int i = 0; i < arc_size; i++)
    archive.push_back((variable *)malloc(sizeof(variable) * problem_size));

  int num_success_params;
  vector<variable> success_sf;
  vector<variable> success_cr;
  vector<variable> dif_fitness;

  // the contents of M_f and M_cr are all initialiezed 0.5
  vector<variable> memory_sf(memory_size, 0.5);
  vector<variable> memory_cr(memory_size, 0.5);

  variable temp_sum_sf;
  variable temp_sum_cr;
  variable sum;
  variable weight;

  // memory index counter
  int memory_pos = 0;

  // for new parameters sampling
  variable mu_sf, mu_cr;
  int random_selected_period;
  variable *pop_sf = (variable *)malloc(sizeof(variable) * pop_size);
  variable *pop_cr = (variable *)malloc(sizeof(variable) * pop_size);

  // for current-to-pbest/1
  int p_best_ind;
  int p_num = alglib::round(pop_size * p_best_rate);
  int *sorted_array = (int *)malloc(sizeof(int) * pop_size);
  Fitness *temp_fit = (Fitness *)malloc(sizeof(Fitness) * pop_size);

  // for linear population size reduction
  int max_pop_size = pop_size;
  int min_pop_size = 4;
  int plan_pop_size;

  int cmp_count = 0, succ_count = 0;
  double mean_gen_cost = 0, best_gen_cost = 1.0e+10;

  std::ofstream success, clusters_stats, gen_costs, conv_speed;

  // success.open("stats/success_rate", std::ofstream::out | std::ofstream::app);
  // gen_costs.open("stats/gen_costs", std::ofstream::out | std::ofstream::app);
  // clusters_stats.open("stats/clusters_stats", std::ofstream::out | std::ofstream::app);
  // conv_speed.open("stats/conv_speed", std::ofstream::out | std::ofstream::app);

  // data mining aux. structures
  map<int, double> _pattern;
  vector<map<int, double>> patterns;
  //vector<tuple<Individual, double>> elite;
  PatternMiner miner;

  int currentPatternIndex = 0;
  // main loop
  while (nfes < max_num_evaluations)
  {
    // cout << nfes << "-" << bsf_fitness - optimum << endl;
    if (bsf_fitness - optimum < 1.0e-8) {
        //cout << "cfo: " << nfes << endl;
        break;
    }

    this->used_gen_count++;

    for (int i = 0; i < pop_size; i++)
      sorted_array[i] = i;
    for (int i = 0; i < pop_size; i++)
      temp_fit[i] = fitness[i];
    sortIndexWithQuickSort(&temp_fit[0], 0, pop_size - 1, sorted_array);

    // update elite set: ========================
    if (elite_type == CROSS_GENERATION) {
      updateElite(pop, sorted_array, temp_fit);
    } else { // BY_GENERATION
      elite.clear();
      for (int i=0; i<p_num; i++) {
          elite.push_back({ pop[sorted_array[i]], fitness[sorted_array[i]] });
      }
    }
    // ==========================================


    // patterns mining: ========================
    vector<vector<double>> points;
    for (int i = 0; i < elite.size(); i++) {
        auto elm = elite[i];
        vector<double> point = vector<double>(get<0>(elm), get<0>(elm) + problem_size);
        points.push_back(point);
    }


    _PatternMiner<double> _pm;
    _pm.mine(points);

    if (mining_algorithm == "xmeans") {
      patterns = miner.extractPatterns(points, lower_bounds, upper_bounds);
    } 
    else {
      int k;
      if (number_of_patterns == 0)
        k = p_num;
      else       
        k = min(number_of_patterns, p_num);
      
      patterns = miner.extractPatterns(points, lower_bounds, upper_bounds, k);
    }
    this->clusters_count += patterns.size();
    //==========================================

    if (pattern_usage_strategy == INSERT_IN_POP) {
      // INSERT IN POP. STRATEGY:
      for (size_t i = 0; i < patterns.size(); i++)
      {
        Individual new_ind = makeNewIndividual();
        
        for (auto var : patterns[i])
          new_ind[var.first] = var.second;

        pop[pop.size() - 1 - i] = new_ind;
      }
    }
    // ========================

    for (int target = 0; target < pop_size; target++)
    {
      // In each generation, CR_i and F_i used by each individual ax_pos are generated by first selecting an index r_i randomly from [1, H]
      random_selected_period = rand() % memory_size;
      mu_sf = memory_sf[random_selected_period];
      mu_cr = memory_cr[random_selected_period];

      // generate CR_i and repair its value
      if (mu_cr == -1)
      {
        pop_cr[target] = 0;
      }
      else
      {
        pop_cr[target] = gauss(mu_cr, 0.1);
        if (pop_cr[target] > 1)
          pop_cr[target] = 1;
        else if (pop_cr[target] < 0)
          pop_cr[target] = 0;
      }

      // generate F_i and repair its value
      do
      {
        pop_sf[target] = cauchy_g(mu_sf, 0.1);
      } while (pop_sf[target] <= 0);

      if (pop_sf[target] > 1)
        pop_sf[target] = 1;


      /**
       * Strategies:
       *  1. pattern and pbest with same id
       *  2. pattern getted from queue
       *  3. pattern and pbest with random(different) ids
      */

      unsigned r = rand(); // LEMBRAR: a distribuição de random numbers influencia
      p_best_ind = sorted_array[r % p_num];
      
      if ( pattern_usage_strategy == INSERT_IN_POP) {
        operateCurrentToPBest1BinWithArchive(pop, &children[target][0], target, p_best_ind, pop_sf[target], pop_cr[target], archive, arc_ind_count);
      } else {
        int chosed_pattern_idx = r % patterns.size();
        switch (pattern_sel_strategy)
        {
        case RANDOM:
          _pattern = patterns[chosed_pattern_idx];
          break;

        case QUEUE:
          _pattern = patterns[currentPatternIndex % patterns.size()];
          currentPatternIndex++;
          break;
        
        default:
          _pattern = patterns[chosed_pattern_idx];
          break;
        }

        operateCurrentToPBest1BinWithArchiveAndXPattern(pop, &children[target][0], target, p_best_ind, pop_sf[target], pop_cr[target], archive, arc_ind_count, _pattern);

      }
    }

    generation++;
    
    // evaluate the children's fitness values
    evaluatePopulation(children, children_fitness);

    /////////////////////////////////////////////////////////////////////////
    // update the bsf-solution and check the current number of fitness evaluations
    //  if the current number of fitness evaluations over the max number of fitness evaluations, the search is terminated
    //  So, this program is unconcerned about L-SHADE algorithm directly
    for (int i = 0; i < pop_size; i++)
    {
      nfes++;
      // if (nfes%5000==0) {
      //     conv_speed << bsf_fitness - optimum << ";";
      // }
      mean_gen_cost += fitness[i];
      // following the rules of CEC 2014 real parameter competition,
      // if the gap between the error values of the best solution found and the optimal solution was 10^{−8} or smaller,
      // the error was treated as 0
      if ((children_fitness[i] - optimum) < epsilon)
        children_fitness[i] = optimum;

      if (children_fitness[i] < bsf_fitness)
      {
        bsf_fitness = children_fitness[i];
        for (int j = 0; j < problem_size; j++)
          bsf_solution[j] = children[i][j];
      }

      if (nfes >= max_num_evaluations)
        break;
    }

    mean_gen_cost /= pop_size;
    ////////////////////////////////////////////////////////////////////////////

    // generation alternation
    for (int i = 0; i < pop_size; i++)
    {
      cmp_count++;
      if (children_fitness[i] == fitness[i])
      {
        fitness[i] = children_fitness[i];
        for (int j = 0; j < problem_size; j++)
          pop[i][j] = children[i][j];
      }
      else if (children_fitness[i] < fitness[i])
      {
        succ_count++;
        dif_fitness.push_back(fabs(fitness[i] - children_fitness[i]));
        fitness[i] = children_fitness[i];
        for (int j = 0; j < problem_size; j++)
          pop[i][j] = children[i][j];

        // successful parameters are preserved in S_F and S_CR
        success_sf.push_back(pop_sf[i]);
        success_cr.push_back(pop_cr[i]);

        // parent vectors ax_pos which were worse than the trial vectors u_i are preserved
        if (arc_size > 1)
        {
          if (arc_ind_count < arc_size)
          {
            for (int j = 0; j < problem_size; j++)
              archive[arc_ind_count][j] = pop[i][j];
            arc_ind_count++;
          }
          // Whenever the size of the archive exceeds, randomly selected elements are deleted to make space for the newly inserted elements
          else
          {
            random_selected_arc_ind = rand() % arc_size;
            for (int j = 0; j < problem_size; j++)
              archive[random_selected_arc_ind][j] = pop[i][j];
          }
        }
      }
    }

    num_success_params = success_sf.size();

    // if numeber of successful parameters > 0, historical memories are updated
    if (num_success_params > 0)
    {
      memory_sf[memory_pos] = 0;
      memory_cr[memory_pos] = 0;
      temp_sum_sf = 0;
      temp_sum_cr = 0;
      sum = 0;

      for (int i = 0; i < num_success_params; i++)
        sum += dif_fitness[i];

      // weighted lehmer mean
      for (int i = 0; i < num_success_params; i++)
      {
        weight = dif_fitness[i] / sum;

        memory_sf[memory_pos] += weight * success_sf[i] * success_sf[i];
        temp_sum_sf += weight * success_sf[i];

        memory_cr[memory_pos] += weight * success_cr[i] * success_cr[i];
        temp_sum_cr += weight * success_cr[i];
      }

      memory_sf[memory_pos] /= temp_sum_sf;

      if (temp_sum_cr == 0 || memory_cr[memory_pos] == -1)
        memory_cr[memory_pos] = -1;
      else
        memory_cr[memory_pos] /= temp_sum_cr;

      // increment the counter
      memory_pos++;
      if (memory_pos >= memory_size)
        memory_pos = 0;

      // clear out the S_F, S_CR and delta fitness
      success_sf.clear();
      success_cr.clear();
      dif_fitness.clear();
    }

    // calculate the population size in the next generation
    plan_pop_size = alglib::round((((min_pop_size - max_pop_size) / (double)max_num_evaluations) * nfes) + max_pop_size);

    if (pop_size > plan_pop_size)
    {
      reduction_ind_num = pop_size - plan_pop_size;
      if (pop_size - reduction_ind_num < min_pop_size)
        reduction_ind_num = pop_size - min_pop_size;

      reducePopulationWithSort(pop, fitness);

      // resize the archive size
      arc_size = pop_size * g_arc_rate;
      if (arc_ind_count > arc_size)
        arc_ind_count = arc_size;

      // resize the number of p-best individuals
      p_num = alglib::round(pop_size * p_best_rate);
      if (p_num <= 1)
        p_num = 2;
    }
  }

  this->clusters_count /= generation;
  this->success_rate = (double)succ_count / cmp_count;
  return bsf_fitness - optimum;
}

void DM_LSHADE::updateElite(vector<Individual> curr_pop, int* sorted_indexes, double* fitness)
{
    int max = elite_max_size > curr_pop.size() ? curr_pop.size() : elite_max_size;
    for (int i = 0; i < max; i++)
    {
        if (elite.size() < elite_max_size) {
            elite.push_back({ curr_pop[sorted_indexes[i]], fitness[i] });
        }
        else
        {
            std::vector<tuple<Individual, double>>::iterator elite_member = elite.end();
            int pos = elite.size()-1;
            while (elite_member != elite.begin() && fitness[i] < get<1>(*(elite_member-1))) {
                elite_member--;
                if (pos > 0) pos--;
            }

            int index = sorted_indexes[i];
            Individual ind = curr_pop[index];
            Fitness fit = fitness[i];

            elite.insert(elite_member, { ind , fit });
        }
    }

    if (elite.size() > elite_max_size)
        elite.resize(elite_max_size);
}


void DM_LSHADE::operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, int &target, int &p_best_individual, variable &scaling_factor, variable &cross_rate, const vector<Individual> &archive, int &arc_ind_count)
{
  int r1, r2;

  do
  {
    r1 = rand() % pop_size;
  } while (r1 == target);
  do
  {
    r2 = rand() % (pop_size + arc_ind_count);
  } while ((r2 == target) || (r2 == r1));

  int random_variable = rand() % problem_size;

  if (r2 >= pop_size)
  {
    r2 -= pop_size;
    for (int i = 0; i < problem_size; i++)
    {
      if ((randDouble() < cross_rate) || (i == random_variable))
      {
        child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - archive[r2][i]);
      }
      else
      {
        child[i] = pop[target][i];
      }
    }
  }
  else
  {
    for (int i = 0; i < problem_size; i++)
    {
      if ((randDouble() < cross_rate) || (i == random_variable))
      {
        child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - pop[r2][i]);
      }
      else
      {
        child[i] = pop[target][i];
      }
    }
  }

  // If the mutant vector violates bounds, the bound handling method is applied
  modifySolutionWithParentMedium(child, pop[target]);
}


void DM_LSHADE::operateCurrentToPBest1BinWithArchiveAndXPattern(const vector<Individual>& pop, Individual child, int& target, int& p_best_individual,
    variable& scaling_factor, variable& cross_rate, const vector<Individual>& archive, int& arc_ind_count, map<int, double> pattern) {
    int r1, r2;

    set<int> unfixed_positions;

    do {
        r1 = rand() % pop_size;
    } while (r1 == target);
    do {
        r2 = rand() % (pop_size + arc_ind_count);
    } while ((r2 == target) || (r2 == r1));

    int random_variable = rand() % problem_size;

    if (r2 >= pop_size) {
        r2 -= pop_size;
        for (int i = 0; i < problem_size; i++) {
            if ((randDouble() < cross_rate) || (i == random_variable)) {
                if (pattern.count(i)) {
                  child[i] = pop[target][i] + scaling_factor * (pattern[i] - pop[target][i]) + scaling_factor * (pop[r1][i] - archive[r2][i]);
                } else {
                  switch (filling_strategy)
                  {   
                  case P_BEST:
                    child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - archive[r2][i]);
                    break;
                  
                  case MOMENTUM:
                    child[i] = pop[target][i];
                    break;
                  
                  case LINE_SEARCH:
                    child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - archive[r2][i]);
                    unfixed_positions.insert(i);
                    break;

                  default:
                    break;
                  }
                }
            }
            else {
                child[i] = pop[target][i];
            }
        }
    }
    else {
        for (int i = 0; i < problem_size; i++) {
            
            if ((randDouble() < cross_rate) || (i == random_variable)) {
                if (pattern.count(i)) {
                  child[i] = pop[target][i] + scaling_factor * (pattern[i] - pop[target][i]) + scaling_factor * (pop[r1][i] - pop[r2][i]);
                } else {
                  switch (filling_strategy)
                  {
                  case P_BEST:
                    child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - pop[r2][i]);
                    break;
                  
                  case MOMENTUM:
                    child[i] = pop[target][i];
                    break;
                  
                  case LINE_SEARCH:
                    child[i] = pop[target][i] + scaling_factor * (pop[p_best_individual][i] - pop[target][i]) + scaling_factor * (pop[r1][i] - pop[r2][i]);
                    unfixed_positions.insert(i);
                    break;

                  default:
                    break;
                  }
                }
            }
            else {
                child[i] = pop[target][i];
            }
        }
    }

    //If the mutant vector violates bounds, the bound handling method is applied
    modifySolutionWithParentMedium(child, pop[target]);

    // Semi-greedy search around child
    if (filling_strategy == LINE_SEARCH) 
      constructGreedyRandomized(child, unfixed_positions);
}


void DM_LSHADE::reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness)
{
  int worst_ind;

  for (int i = 0; i < reduction_ind_num; i++)
  {
    worst_ind = 0;
    for (int j = 1; j < pop_size; j++)
    {
      if (fitness[j] > fitness[worst_ind])
        worst_ind = j;
    }

    pop.erase(pop.begin() + worst_ind);
    fitness.erase(fitness.begin() + worst_ind);
    pop_size--;
  }
}

void DM_LSHADE::setSHADEParameters()
{
  arc_rate = g_arc_rate;
  arc_size = (int)alglib::round(pop_size * arc_rate);
  p_best_rate = g_p_best_rate;
  memory_size = g_memory_size;
}
