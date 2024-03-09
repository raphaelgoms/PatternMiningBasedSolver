#include <set>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <map>
#include <any>
#include <nlohmann/json.hpp>

using namespace std;

class ProgramArgsParser
{
private:
    int args_count;
    char **arg_vector;
    set<string> args_names;
    map<string, string> arg_types;
    map<string, any> arg_values;
    nlohmann::json args_json;
    
public:
    ProgramArgsParser(/* args */);
    ~ProgramArgsParser();

    
    void parseArgs(int argc, char **argv, bool printargs=false);
    any getArgValue(string arg);
};


