#include <vector>
#include <map>
#include <math.h>

#define PI 3.1415926535897932384626433832795029

using namespace std;

typedef struct {
    vector<double> points;
    double centroid;
} Cluster;

typedef struct {
    vector<vector<double>> points;
    vector<double> centroid;
	double radius;
} MCluster; // MultiCluster: para dados multidimencionais;

class Xmeans {
	vector<Cluster> kMeans(const vector<double> &data, int k, double min, double max);
	vector<MCluster> kMeans(const vector<vector<double>> &data, int k, vector<double> min, vector<double> max);
	
	int improveStructure(int k, const vector<Cluster> &clusters, int M, double min, double max);
	int improveStructure(int k, const vector<MCluster> &clusters, int M, const vector<double> &min, const vector<double> &max);

	double logLikelihood(int R, int Rn, double variance, double M, double K);
  	double variance(const vector<double> &data);
 	double variance(const vector<vector<double>> &data, const vector<double> &centroid);

  	
	vector<double>getNormalizedStandartDeviations(const  vector<vector<double>> &cluster, 
		const vector<double> &centroid, const vector<double> &lower_bounds, 
		const vector<double> &upper_bounds);

	inline double randDouble() {
    	return (double)rand() / (double) RAND_MAX;
	}

	inline double normalize(const double &value, const double &lower_bound, const double &upper_bound) {
    	return (value - lower_bound) / (upper_bound - lower_bound);
	}

	double max_radius;

public:
	double average(const vector<double> &data);
    double distance(const vector<double> & v1, const vector<double> & v2);
  	

	vector<double> getCentroid(const vector<vector<double>> &cluster);
	double getMaxRadius() { return max_radius; }	
	vector<Cluster> xMeans(const vector<double> &data, int kMax, double min, double max);
	vector<MCluster> xMeans(const vector<vector<double>> &data, int kMax, const vector<double> &min, const vector<double> &max);
	
	map<int, double> extractPatternByVariable(const vector<vector<double>>& elite, const vector<double>& lower_bound, const vector<double>& upper_bound, int k=0);
	map<int, double> extractPatternGX(const vector<vector<double>>& elite, const vector<double>& lower_bound, const vector<double>& upper_bound, int k=0);
	map<int, double> extractPatternRX(const vector<vector<double>>& elite, const vector<double>& lower_bound, const vector<double>& upper_bound, int k=0);
	
	// return the list de clusters as partterns
	vector<map<int, double> > extractPatterns(const vector<vector<double>>& elite, const vector<double>& lower_bound, const vector<double>& upper_bound, int k=0);
};

