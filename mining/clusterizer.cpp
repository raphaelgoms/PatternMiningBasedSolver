#include "clusterizer.hpp"
#include <sstream>
#include <iostream>

std::vector<std::vector<double>> Clusterizer::getClusterPoints(const DataSet& ds, int k, integer_1d_array cidx) {
	std::vector<std::vector<double>> cpoints;
    for (size_t i = 0; i < cidx.length(); i++) {
		if (cidx[i] == k) 
			cpoints.push_back(ds[i]);
	}	
	return cpoints; 
}

kmeansreport Clusterizer::xmeans_alglib(const DataSet &ds, int k_max)
{
    return kmeansreport();
}

kmeansreport Clusterizer::kmeans_alglib(const DataSet &ds, int k)
{
    clusterizerstate s;
    kmeansreport rep;

	real_2d_array ar = dataSetToReal2dArray(ds);	
	ae_int_t disttype = 2; // Euclidian distance
	int _k = k;

	do {
		try {
			clusterizercreate(s);
			clusterizersetpoints(s, ar, disttype);
			//clusterizersetkmeanslimits(s, 5, 0);
			clusterizerrunkmeans(s, _k, rep);
			_k--;
		} catch (alglib::ap_error ap_error) {
			std::cout << ap_error.msg << std::endl;
		}		
	} while (rep.terminationtype == -3);

	return rep;
}

std::vector<Cluster> Clusterizer::kmeansreportToClusters(const DataSet& ds, kmeansreport ks)
{
    std::vector<Cluster> clusters;
    for (int i = 0; i < ks.c.rows(); i++) {
        Cluster cl {
            i, 
            std::vector<double>(ks.c[i], ks.c[i] + ds.front().size()), 
            getClusterPoints(ds, cl.id, ks.cidx)
        };

        clusters.push_back(cl);
    }
    
    return clusters;
}

real_2d_array Clusterizer::dataSetToReal2dArray(const DataSet &ds)
{
    std::stringstream sdata;
	sdata << "[";
	for(auto it =ds.begin();it!=ds.end();it++) {
		if (it != ds.begin()) 
			sdata << ", ";

		sdata << "[";
		for(int i =0;i<(*it).size();i++) {
			if (i>0) 
                sdata << ", ";
			sdata << (*it)[i];
		}
		sdata << "]";
	}
	sdata << "]";

	return real_2d_array(sdata.str().c_str());    
}

Clusterizer::Clusterizer()
{
}

Clusterizer::Clusterizer(DataSet ds)
{
    this->ds = ds;
}

void Clusterizer::setDataSet(DataSet ds)
{
    this->ds = ds;
}

std::vector<Cluster> Clusterizer::run(std::string algorithm, int k)
{
    return std::vector<Cluster>();
}

std::vector<Cluster> Clusterizer::run(DataSet ds, std::string algorithm, int k)
{
    if (algorithm == "kmeans") 
        return kmeansreportToClusters(ds, kmeans_alglib(ds, k));
    
    return std::vector<Cluster>();
}
