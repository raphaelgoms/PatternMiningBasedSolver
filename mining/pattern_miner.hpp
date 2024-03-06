#include "alglib/stdafx.h"
#include "alglib/dataanalysis.h"

#include <vector>
#include <map>
#include <math.h>
#include <random>

#define PI 3.1415926535897932384626433832795029

using namespace std;
using namespace alglib;

extern int g_number_of_used_vars;

template <typename T> // Currently T can be a real number or a interval of real numbers
using Pattern = map<int, T>;
using Interval = tuple<double, double>;

template <typename T>
class _PatternMiner {
	vector<Pattern<T>> mine(vector<vector<double>> data_set);
};

class PatternMiner {

	vector<vector<double>> dataSet;
	vector<vector<double>> getClusterPoints(int k, integer_1d_array cidx);

	kmeansreport kMeans_alglib(const vector<vector<double>> &data, int k, vector<double> min, vector<double> max);
	kmeansreport xMeans_alglib(const vector<vector<double>> &data, int kMax, const vector<double> &min, const vector<double> &max);

	int improveStructure(int k, const kmeansreport &clusters, int M, const vector<double> &min, const vector<double> &max);
	double logLikelihood(int R, int Rn, double variance, double M, double K);
 	
	double average(const vector<double> &data);
	double variance(const vector<vector<double>> &data, const vector<double> &centroid);
    double distance(const vector<double> & v1, const vector<double> & v2);
	double distance(double * v1, double *  v2, size_t vsize);

	inline double normalize(const double &value, const double &lower_bound, const double &upper_bound) {
    	return (value - lower_bound) / (upper_bound - lower_bound);
	}
	
public:
	// return the list de clusters as partterns
	vector<map<int, double>> extractPatterns(vector<vector<double>> data_set, const vector<double>& lower_bound, const vector<double>& upper_bound, int k=0);
};

