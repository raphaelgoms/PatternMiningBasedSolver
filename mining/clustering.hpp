#include <vector>
#include <string>

using DataSet = std::vector<std::vector<double>>;

typedef struct Cluster
{
    int id;
    std::vector<double> centroid;
    std::vector<std::vector<double>> members;
};

class Cluterizer {
private:
    DataSet ds;    

public:
    Clusterizer();
    Clusterizer(DataSet ds);

    std::vector<Cluster> run(std::string algorithm="kmeans", int k=3);
    std::vector<Cluster> run(DataSet ds, std::string algorithm="kmeans", int k=3);

    void setDataSet(DataSet ds);
};