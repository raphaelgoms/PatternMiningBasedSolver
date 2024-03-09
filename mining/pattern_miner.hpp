#include "alglib/stdafx.h"
#include "alglib/dataanalysis.h"

#include <vector>
#include <map>
#include <math.h>
#include <random>

#include "clusterizer.hpp"

#define PI 3.1415926535897932384626433832795029

using namespace std;
using namespace alglib;

extern int g_number_of_used_vars;

template <typename T> // Currently T can be a real number or a interval of real numbers
using Pattern = map<int, T>;
using Interval = tuple<double, double>;

template <typename T>
class PatternMiner {
public:
	virtual vector<Pattern<T>> mine(const DataSet& data_set, int nmb_of_patterns=0) {};
};

template <>	
inline vector<Pattern<double>> PatternMiner<double>::mine(const DataSet& data_set, int nmb_of_patterns) {
	Clusterizer clusterizer;
	std::vector<Cluster> clusters;
	
	if (nmb_of_patterns)
		clusters = clusterizer.run(data_set, "kmeans", nmb_of_patterns);	
	else 
		clusters = clusterizer.run(data_set, "xmeans");

	std::vector<Pattern<double>> patterns;
	for (auto cl=clusters.begin(); cl!=clusters.end(); cl++) {
		Pattern<double> p;
		for (size_t k = 0; k < (*cl).centroid.size(); k++) 
			p[k] = (*cl).centroid[k];
		patterns.push_back(p);
	}

    return patterns;
}

template <>
inline vector<Pattern<Interval>> PatternMiner<Interval>::mine(const DataSet& data_seti, int nmb_of_patterns) {
    return vector<Pattern<Interval>>();
}