#include "xmeans.h"
#include <math.h>
#include <iostream>

vector<map<int, double> > Xmeans::extractPatterns(const vector<vector<double> >& elite, const vector<double>& lower_bound, const vector<double>& upper_bound, int k=0) {
	
	// Clusteriza elite:
	int elite_size = elite.size();
	int problem_size = elite.front().size();
	vector<MCluster> clusters;
	
	if (k > 0)
		clusters = kMeans(elite, k, lower_bound, upper_bound);
	else {	
		clusters = xMeans(elite, elite_size, lower_bound, upper_bound);
	}

	int idx = -1;
	int maxSize = 0;

	max_radius = 0;
	vector<map<int, double> > patterns;
	//std::cout << elite.size() << " | " << clusters.size() << std::endl;
	for (int i = 0; i < clusters.size(); i++)
	{
		if (clusters[i].points.size()) {
			vector<double> centroid = getCentroid(clusters[i].points);
			map<int, double> pattern;
			for (int i = 0; i < problem_size; i++)
				pattern.insert(pair<int, double>(i, centroid[i]));
			
			patterns.push_back(pattern);
			if (max_radius < clusters[i].radius) {
				max_radius = clusters[i].radius;
			}
		}		
	}

	return patterns;	
}


map<int, double> Xmeans::extractPatternGX(const vector<vector<double> >& elite, const vector<double>& lower_bound, const vector<double>& upper_bound, int k=0) {
	
	// Clusteriza elite:
	int elite_size = elite.size();
	int problem_size = elite.front().size();
	vector<MCluster> clusters;
	
	if (k > 0)
		clusters = kMeans(elite, k, lower_bound, upper_bound);
	else {	
		clusters = xMeans(elite, elite_size, lower_bound, upper_bound);
	}

	int idx = -1;
	int maxSize = 0;
	
	for (int i = 0; i < clusters.size(); i++)
	{
		if (clusters[i].points.size() > maxSize) {
			maxSize = clusters[i].points.size();
			//cout << maxSize << endl;
			idx = i;
		}
	}

	MCluster greatest_cluster = clusters[idx];

	// Extrai padrão
	vector<double> meanVector = getCentroid(greatest_cluster.points);
	//vector<double> stdVector = getNormalizedStandartDeviations(greatest_cluster.points, meanVector, lower_bound, upper_bound);
	
	map<int, double> pattern;

	for (int i = 0; i < problem_size; i++) {
		//std::cout << stdVector[i] << " ";
		//if (stdVector[i] < 0.01)
		pattern[i] = meanVector[i];
	} 
	//std::cout << std::endl;

	return pattern;	
}

map<int, double> Xmeans::extractPatternRX(const vector<vector<double> >& elite, const vector<double>& lower_bound, const vector<double>& upper_bound, int k) {
	
	// Clusteriza elite:
	int elite_size = elite.size();
	int problem_size = elite.front().size();
	vector<MCluster> clusters;
	
	if (k > 0)
		clusters = kMeans(elite, k, lower_bound, upper_bound);
	else {	
		clusters = xMeans(elite, elite_size, lower_bound, upper_bound);
	}

	int idx = -1;
	int maxSize = 0;
	
	idx = (int)round((clusters.size() -1) * randDouble());
	MCluster greatest_cluster = clusters[idx];

	// Extrai padrão
	vector<double> meanVector = getCentroid(greatest_cluster.points);
	vector<double> stdVector = getNormalizedStandartDeviations(greatest_cluster.points, meanVector, lower_bound, upper_bound);
	
	map<int, double> pattern;

	for (int i = 0; i < problem_size; i++) {
		//std::cout << stdVector[i] << " ";
		//if (stdVector[i] < 0.01)
		pattern[i] = meanVector[i];
	} 
	//std::cout << std::endl;

	return pattern;	
}

map<int, double> Xmeans::extractPatternByVariable(const vector<vector<double> >& elite, const vector<double>& lower_bound, const vector<double>& upper_bound, int k)
{
	map<int, double> pattern;
	int g_elite_size = elite.size();
	int g_problem_size = elite.front().size();
	vector<double> data(g_elite_size);
	vector<double> choosed;

	for (int i = 0; i < g_problem_size; i++)
	{
		for (unsigned int j = 0; j < g_elite_size; j++)
		{
			data[j] = elite[j][i] ;
		}

		vector<Cluster> clusters;
	
		if (k > 0)
			clusters = kMeans(data, k, lower_bound[i], upper_bound[i]);
		else {	
			clusters = xMeans(data, g_elite_size, lower_bound[i], upper_bound[i]);
		}

		vector<Cluster> clusters = xMeans(data, g_elite_size, lower_bound[i], upper_bound[i]);
		//cout << "Clusters size: " << clusters.size() << endl;
		int max = -1;
		int gt_cluster;

		for (int l = 0; l < clusters.size(); l++)
		{
			int pt_count = clusters[l].points.size();
			if (max < pt_count)
			{
				max = pt_count;
				gt_cluster = l;
			}
		}
		
		pattern[i] = clusters[gt_cluster].centroid;
	} 

	return pattern;
}

vector<Cluster> Xmeans::kMeans(const vector<double> &data, int k, double min, double max)
{
	int nearest;
	double dist;
	vector<Cluster> clusters(k);

	for (int i = 0; i < k; i++)	{
		Cluster c;
		c.centroid = (max - min) * rand() + min;
		clusters[i] = c;
		clusters[i].points = vector<double>();
	}

	bool converged = false;
	int numberOfIterations = 0;
	int MaxNumberOfIterarions = 100;

	while (!converged && numberOfIterations < MaxNumberOfIterarions)
	{
		for (int i = 0; i < k; i++)
			clusters[i].points = vector<double>();

		for (int i = 0; i < data.size(); i++) {
			double min_dist = INFINITY;
			for (int j = 0; j < k; j++)
			{
				double dist = fabs(data[i] - clusters[j].centroid);
				if (dist < min_dist)
				{
					min_dist = dist;
					nearest = j;
				}
			}

			clusters[nearest].points.push_back(data[i]);
		}

		converged = true;
		for (int i = 0; i < k; i++)
		{
			if (clusters[i].points.size() == 0)
				continue;

			double centroid = average(clusters[i].points);
			if (centroid != clusters[i].centroid)
				converged = false;

			clusters[i].centroid = centroid;
		}

		numberOfIterations++;
	}

	return clusters;
}

double Xmeans::distance(const vector<double> & v1, const vector<double> & v2) {
	double sum = 0;
	for (int i = 0; i < v1.size(); i++)	{
		sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	
	return sqrt(sum);
}

vector<double> Xmeans::getCentroid(const vector<vector<double> > &cluster) {

	int dim = cluster[0].size();
	vector<double> centroid;

	for (int i = 0; i < dim; i++) {

		double sum = 0.0;
		for (int j = 0; j < cluster.size(); j++) {
			sum += cluster[j][i];
		}

		centroid.push_back(sum / cluster.size());
	}

	return centroid;
}

vector<double> Xmeans::getNormalizedStandartDeviations(const vector<vector<double> > &cluster, 
	const vector<double> &centroid, const vector<double> &lower_bounds, 
	const vector<double> &upper_bounds) {

	int dim = centroid.size();
	vector<double> deviations;

	for (int i = 0; i < dim; i++) {

		double acc = 0.0;
		for (int j = 0; j < cluster.size(); j++) {
			double normalized_mean = normalize(centroid[i], lower_bounds[i], upper_bounds[i]);
			double normalized_value = normalize(cluster[j][i], lower_bounds[i], upper_bounds[i]);
			acc += (normalized_value - normalized_mean) *  (normalized_value - normalized_mean);
		}

		deviations.push_back(sqrt(acc / cluster.size()));
	}
	
	return deviations;
}


vector<MCluster> Xmeans::kMeans(const vector<vector<double> > &data, int k, vector<double> min, vector<double> max)
{
	int nearest;
	double dist;
	vector<MCluster> clusters(k);
	int dim = data[0].size();

	for (int i = 0; i < k; i++)	{
		MCluster c;
		c.radius = 0;
		c.centroid = vector<double>();
		for (int j = 0; j < dim; j++) {
			c.centroid.push_back((max[j] - min[j]) * randDouble() + min[j]);
		}
		clusters[i] = c;
		clusters[i].points = vector<vector<double> >();
	}

	bool converged = false;
	int numberOfIterations = 0;
	int MaxNumberOfIterarions = 100;

	while (!converged && numberOfIterations < MaxNumberOfIterarions)
	{
		for (int i = 0; i < k; i++)
			clusters[i].points.clear();

		for (int i = 0; i < data.size(); i++) {
			double min_dist = INFINITY;
			double selected_dist = 0;
			for (int j = 0; j < k; j++)
			{
				double dist = distance(data[i], clusters[j].centroid);
				if (dist < min_dist) {				
					min_dist = dist;
					nearest = j;
					selected_dist = dist;
				}
			}

			if (selected_dist > clusters[nearest].radius) 
				  clusters[nearest].radius = selected_dist;

			clusters[nearest].points.push_back(data[i]);
		}

		converged = true;
		for (int i = 0; i < k; i++)
		{
			
			if (clusters[i].points.size() == 0)
				continue;

			vector<double> centroid = getCentroid(clusters[i].points);
			if (centroid != clusters[i].centroid)
				converged = false;
			
			clusters[i].centroid = centroid;
		}

		numberOfIterations++;
	}

	return clusters;
}

vector<Cluster> Xmeans::xMeans(const vector<double> &data, int kMax, double min, double max)
{
	int k = 1;
	int k_old = k;
	 
	bool stopSplitting = false;
	int MaxNumberOfIterations = 2;

	int iteration = 0;
	vector<Cluster> clusters;

	while (!stopSplitting && k < kMax)
	{
		k_old = k;

		clusters = kMeans(data, k, min, max);
		int add_k = improveStructure(k, clusters, 1, min, max);

		k += add_k;
		stopSplitting = k_old == k || k >= kMax;
	}

	return kMeans(data, k_old, min, max);
}

vector<MCluster> Xmeans::xMeans(const vector<vector<double> > &data, int kMax, const vector<double> &min, const vector<double> &max)
{
	int k = 1;
	int k_old = k;
	 
	bool stopSplitting = false;
	int MaxNumberOfIterations = 2;

	int iteration = 0;
	vector<MCluster> clusters;

	while (!stopSplitting && k < kMax)
	{
		k_old = k;

		clusters = kMeans(data, k, min, max);

		int add_k = improveStructure(k, clusters, 1, min, max);

		k += add_k;
		stopSplitting = k_old == k || k >= kMax;
	}

	return kMeans(data, k_old, min, max);
}

int Xmeans::improveStructure(int k, const vector<MCluster> &clusters, int M, const vector<double> &min, const vector<double> &max)
{
	vector<double> bic_before_split(k);
	vector<double> bic_after_split(k);

	int clst_n_params = M + 1;
	int add_k = 0;

	for (int i = 0; i < k; i++)
	{
		vector<vector<double> > clst_points = clusters[i].points;
		int clst_size = clst_points.size();

		if (clst_size <= 2)
			continue;

		double clst_variance = variance(clst_points, getCentroid(clst_points));

		bic_before_split[i] = logLikelihood(clst_size, clst_size, clst_variance, M, 1) - (clst_n_params / 2.0) * log(clst_size);

		vector<MCluster> sub_clusters = kMeans(clst_points, 2, min, max);

		double log_likelihood = 0.0;
		for (int j = 0; j < 2; j++)
		{
			vector<vector<double> > subclst_points = sub_clusters[j].points;
			int subclst_size = subclst_points.size();

			if (subclst_size <= 2)
				continue;

			double subclst_variance = variance(subclst_points, getCentroid(subclst_points));
			log_likelihood = log_likelihood + logLikelihood(clst_size, subclst_size, subclst_variance, M, 2);
		}
		int subclst_n_params = 2 * clst_n_params;
		bic_after_split[i] = log_likelihood - (subclst_n_params / 2.0) * log(clst_size);
		if (bic_before_split[i] < bic_after_split[i])
			add_k += 1;
	}

	return add_k;
}

int Xmeans::improveStructure(int k, const vector<Cluster> &clusters, int M, double min, double max)
{
	vector<double> bic_before_split(k);
	vector<double> bic_after_split(k);

	int clst_n_params = M + 1;
	int add_k = 0;

	for (int i = 0; i < k; i++)
	{
		vector<double> clst_points = clusters[i].points;
		int clst_size = clst_points.size();

		if (clst_size <= 2)
			continue;

		double clst_variance = variance(clst_points);
		bic_before_split[i] = logLikelihood(clst_size, clst_size, clst_variance, M, 1) - (clst_n_params / 2.0) * log(clst_size);

		vector<Cluster> sub_clusters = kMeans(clst_points, 2, min, max);

		double log_likelihood = 0.0;
		for (int j = 0; j < 2; j++)
		{
			vector<double> subclst_points = sub_clusters[j].points;
			int subclst_size = subclst_points.size();

			double subclst_variance = variance(subclst_points);
			log_likelihood = log_likelihood + logLikelihood(clst_size, subclst_size, subclst_variance, M, 2);
		}

		int subclst_n_params = 2 * clst_n_params;
		bic_after_split[i] = log_likelihood - (subclst_n_params / 2.0) * log(clst_size);
		if (bic_after_split[i] > bic_before_split[i])
			add_k += 1;
	}

	return add_k;
}

double Xmeans::logLikelihood(int R, int Rn, double variance, double M, double K)
{
	double res = Rn * (log(Rn) - log(R) - 0.5 * (log(2 * PI) + M * log(variance) + 1)) + 0.5 * K;
	if (isnan(res)) res = 0.0;
	return res;
}

double Xmeans::average(const vector<double> &data) {

	if (data.size() == 0) return 0.0;

	double sum = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		sum += data[i];
	}
	return sum / data.size();
}

double Xmeans::variance(const vector<double> &data) {

	if (data.size() == 0) return 0.0;

	double mean = average(data);

	double sum = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		sum += (data[i] - mean) * (data[i] - mean);
	}
	return sum / data.size();
}

double Xmeans::variance(const vector<vector<double> > &data, const vector<double> &centroid) {

	if (data.size() == 0) 
		return 0.0;

	double sum = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		sum += distance(centroid, data[i]);
	}
	return sum / (data.size());
}
