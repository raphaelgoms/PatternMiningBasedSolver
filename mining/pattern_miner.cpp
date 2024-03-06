#include "pattern_miner.hpp"
#include <math.h>
#include <iostream>
#include <sstream>

// Recursive quick sort with index array
template <class VarType>
void _sortIndexWithQuickSort(VarType array[], int first, int last, int index[]) { VarType x = array[(first + last) / 2]; int i = first;
	int j = last;
	VarType temp_var = 0;
	int temp_num = 0;

	while (true)
	{
		while (array[i] < x)
			i++;
		while (x < array[j])
			j--;
		if (i >= j)
			break;

		temp_var = array[i];
		array[i] = array[j];
		array[j] = temp_var;

		temp_num = index[i];
		index[i] = index[j];
		index[j] = temp_num;

		i++;
		j--;
	}

	if (first < (i - 1))
		_sortIndexWithQuickSort(array, first, i - 1, index);
	if ((j + 1) < last)
		_sortIndexWithQuickSort(array, j + 1, last, index);
}

template <class T>
void print(vector<vector<T>> data)
{
	for (size_t l = 0; l < data.size(); l++)
	{
		cout << '[';
		for (int m = 0; m < data[l].size(); m++)
		{
			if (m > 0)
				cout << ", ";
			cout << data[l][m];
		}
		cout << ']' << endl;
	};
}

real_2d_array convertToReal2dArray(const vector<vector<double> >& data_set) {
	stringstream sdata;
	sdata << "[";
	for(auto it =data_set.begin();it!=data_set.end();it++)
	{
		if (it != data_set.begin()) {
			sdata << ", ";
		}

		sdata << "[";
		for(int i =0;i<(*it).size();i++)
		{
			if (i>0) {
				sdata << ", ";
			}
			sdata << (*it)[i];
		}
		sdata << "]";
	}
	sdata << "]";

	real_2d_array xy(sdata.str().c_str());
	return xy;
}

vector< map<int, double> > PatternMiner::extractPatterns(vector<vector<double> > data_set, const vector<double>& lower_bound, const vector<double>& upper_bound, int k) {
	
	this->dataSet = data_set;
	vector<map<int, double> > patterns;

	Pattern<double> p;
	p[1] = 0.1;

	Pattern<Interval> p2;
	p2[1] = { 0, 1 };

	kmeansreport rep;
	
	if (k > 0)
		rep = kMeans_alglib(data_set, k, lower_bound, upper_bound);  
	else
		rep = xMeans_alglib(data_set, data_set.size(), lower_bound, upper_bound);

	integer_1d_array cidx = rep.cidx;
    real_2d_array cz = rep.c;

	int problem_size = data_set.front().size();

	int *sorted_array = (int *)malloc(sizeof(int) * problem_size);
  	double *variables_std_dev = (double *)malloc(sizeof(double) * problem_size);

	vector<vector<double>> clst_points;
	for (int i = 0; i < cz.rows(); i++)
	{
		double *centroid = cz[i];
		clst_points = getClusterPoints(i, cidx);

		double min_dev = 1.0e10;
		double avg_dev = 0.0;
		int chosed_var = alglib::randomreal() * problem_size;

		for (size_t j = 0; j < problem_size; j++)
		{	
			variables_std_dev[j] = 0.0;
			for (size_t l = 0; l < clst_points.size(); l++)
			{
				variables_std_dev[j] += pow(clst_points[l][j] - centroid[j], 2.0);
			}

			variables_std_dev[j] = sqrt(variables_std_dev[j]/clst_points.size());
			avg_dev += variables_std_dev[j];

			if (variables_std_dev[j] < min_dev) {
				min_dev = variables_std_dev[j];
				chosed_var = j;
			}

			sorted_array[j] = j;
		} 
		avg_dev /= problem_size;

		// Ordenar variaveis pelo desvio padrão:
		_sortIndexWithQuickSort(&variables_std_dev[0], 0, problem_size - 1, sorted_array);
		
		map<int, double> _pattern;
		int num_of_patterns_to_use = g_number_of_used_vars;
		for (size_t j = 0; j < num_of_patterns_to_use; j++)
		{	
			_pattern.insert(pair<int, double>(sorted_array[j], centroid[sorted_array[j]]));
		}
			
		patterns.push_back(_pattern);
	}

	return patterns;	
}

double PatternMiner::distance(const vector<double> & v1, const vector<double> & v2) {
	double sum = 0;
	for (int i = 0; i < v1.size(); i++)	{
		sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	
	return sqrt(sum);
}

double PatternMiner::distance(double * v1, double *  v2, size_t vsize) {
	double sum = 0;
	for (int i = 0; i < vsize; i++)	{
		sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	
	return sqrt(sum);
}

kmeansreport PatternMiner::kMeans_alglib(const vector<vector<double>> &data, int k, vector<double> min, vector<double> max) {

	clusterizerstate s;
    kmeansreport rep;

	real_2d_array data_set = convertToReal2dArray(data);	
	ae_int_t disttype = 2; // Euclidian distance
	int _k = k;

	do {

		try {
			clusterizercreate(s);
			clusterizersetpoints(s, data_set, disttype);
			//clusterizersetkmeanslimits(s, 5, 0);
			clusterizerrunkmeans(s, _k, rep);
			_k--;
		} catch (alglib::ap_error ap_error) {
			cout << ap_error.msg << endl;
		}

		
	} while (rep.terminationtype == -3);

	return rep;
}

kmeansreport PatternMiner::xMeans_alglib(const vector<vector<double>> &data, int kMax, const vector<double> &min, const vector<double> &max)
{
	int k = 1;
	int k_old = k;
	 
	bool stopSplitting = false;
	int MaxNumberOfIterations = 2;

	int iteration = 0;
	kmeansreport clusters_rep;

	while (!stopSplitting && k < kMax)
	{
		k_old = k;

		clusters_rep = kMeans_alglib(data, k, min, max);

		int add_k = improveStructure(k, clusters_rep, data[0].size(), min, max);

		k += add_k;
		stopSplitting = k_old == k || k >= kMax;
	}

	return kMeans_alglib(data, k_old, min, max);
}

vector<vector<double>> PatternMiner::getClusterPoints(int k, integer_1d_array cidx) {
	vector<vector<double>> cpoints;
	for (size_t i = 0; i < cidx.length(); i++)
	{
		if (cidx[i] == k) 
			cpoints.push_back(dataSet[i]);
	}	
	return cpoints;
 }

int countClustersPoints(int k, integer_1d_array cidx) {
	int count = 0;
	for (size_t i = 0; i < cidx.length(); i++) {
		if (cidx[i] == k) count++;
	}
	return count;
}

int PatternMiner::improveStructure(int k, const kmeansreport &clusters_rep, int M, const vector<double> &min, const vector<double> &max)
{
	vector<double> bic_before_split(k);
	vector<double> bic_after_split(k);

	int clst_n_params = M + 1;
	int add_k = 0;

	int clst_size, subclst_size;
	for (int i = 0; i < k; i++)
	{
		vector<vector<double>> clst_points = getClusterPoints(i, clusters_rep.cidx);
		clst_size = clst_points.size();
		if (clst_size < 2)
		 	continue;

		auto centroid = vector<double>(clusters_rep.c[i], clusters_rep.c[i] + M);
		double clst_variance = variance(clst_points, centroid);

		bic_before_split[i] = logLikelihood(clst_size, clst_size, clst_variance, M, 1) - (clst_n_params / 2.0) * log(clst_points.size());
		
		kmeansreport sub_clusters_rep = kMeans_alglib(clst_points, 2, min, max);

		double log_likelihood = 0.0;
		for (int j = 0; j < 2; j++)
		{
			vector<vector<double> > subclst_points = getClusterPoints(j, sub_clusters_rep.cidx);
			int subclst_size = subclst_points.size();

			if (subclst_size < 1)
			 	continue;

			auto centroid = vector<double>(sub_clusters_rep.c[j], sub_clusters_rep.c[j] + M);
			double subclst_variance = variance(subclst_points, centroid);

			log_likelihood += logLikelihood(clst_size, subclst_size, subclst_variance, M, 2);
		}

		int subclst_n_params = 2 * clst_n_params;
		bic_after_split[i] = log_likelihood - (subclst_n_params / 2.0) * log(clst_size);
		if (bic_before_split[i] < bic_after_split[i])
			add_k += 1;
	}

	return add_k;
}

double PatternMiner::logLikelihood(int R, int Rn, double variance, double M, double K)
{
	double res = Rn * (log(Rn) - log(R) - 0.5 * (log(2 * PI) + M * log(variance) + 1)) + 0.5 * K;
	if (isnan(res)) res = 0.0;
	return res;
}

double PatternMiner::average(const vector<double> &data) {

	if (data.size() == 0) return 0.0;

	double sum = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		sum += data[i];
	}
	return sum / data.size();
}

double PatternMiner::variance(const vector<vector<double> > &data, const vector<double> &centroid) {

	if (data.size() == 0) 
		return 0.0;

	double sum = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		sum += distance(centroid, data[i]);
	}
	return sum / (data.size());
}
