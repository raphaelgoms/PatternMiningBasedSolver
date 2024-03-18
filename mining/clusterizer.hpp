#include <vector>
#include <string>

#include "alglib/stdafx.h"
#include "alglib/dataanalysis.h"

using namespace alglib;
using DataSet = std::vector<std::vector<double>>;

struct Cluster
{
    int id;
    std::vector<double> centroid;
    std::vector<std::vector<double>> members;
};

class Clusterizer { // ???? Transform this in a interface ????
private:
    DataSet ds;

    // functions wrappers to algorithms from alglib:
    kmeansreport xmeans_alglib(const DataSet& ds, int k_max=-1);
    kmeansreport kmeans_alglib(const DataSet& ds, int k=3);

    // auxiliar methods
    std::vector<Cluster> kmeansreportToClusters(const DataSet& ds, kmeansreport ks);
    real_2d_array dataSetToReal2dArray(const DataSet& ds);
    std::vector<std::vector<double>> getClusterPoints(const DataSet &ds, int k, integer_1d_array cidx);
    
public:
    Clusterizer();
    Clusterizer(DataSet ds);

    std::vector<Cluster> run(std::string algorithm="kmeans", int k=3);
    std::vector<Cluster> run(DataSet ds, std::string algorithm="kmeans", int k=3);

    void setDataSet(DataSet ds);
};