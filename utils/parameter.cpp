#include "parameter.h"

#include <iostream>
#include <fstream>
#include <string>


using namespace std;

ProgramArgs::ProgramArgs(/* args */)
{
    ifstream file;
    file.open("settings.json");
    file >> args_json;

    for (auto&p:args_json["parameters"]) {
        if (p.contains("values")) {
            for (auto&v:p["values"]) {
                args_names.insert(v);
                arg_types[v] = p["type"];
            } 
        } else { 
            args_names.insert(p["name"]);
            arg_types[p["name"]] = p["type"];
        }
    }
}

ProgramArgs::~ProgramArgs()
{
}

void ProgramArgs::parseArgs(int argc, char **argv, bool printargs)
{
    if (printargs) {
        cout << "--------------------------------" << endl;
        cout << "\tProgram Arguments" << endl;
        cout << "--------------------------------" << endl;
    }

    int i = 0;
    while (++i < argc)
    {
        for (string argname: args_names) {

            string __argname = "--" + argname;
            if (!strcmp(__argname.c_str(), argv[i]))
            {
                if (arg_types[argname] == "integer") {
                    arg_values[argname] = atoi(argv[++i]);
                    if (printargs) std::cout << argname << ": " << std::any_cast<int>(arg_values[argname]) << endl;
                }
                else if (arg_types[argname] == "real") {
                    float arg = atof(argv[++i]);
                    arg_values[argname] = arg;
                    if (printargs) std::cout << argname << ": " << std::any_cast<float>(arg_values[argname]) << endl;
                }
                else if (arg_types[argname] == "AlgorithmCode") {
                    arg_values[argname] = argname;
                    if (printargs) std::cout << arg_types[argname] << ": " << std::any_cast<string>(arg_values[argname]) << endl;
                }
                else if (arg_types[argname] == "EliteTypeCode") {
                    arg_values[argname] = argname;
                    if (printargs) std::cout << arg_types[argname] << ": " << std::any_cast<string>(arg_values[argname]) << endl;
                }
                else if (arg_types[argname] == "PatternUsageStrategy") {
                    arg_values[argname] = argname;
                    if (printargs) std::cout << arg_types[argname] << ": " << std::any_cast<string>(arg_values[argname]) << endl;
                }
                else if (arg_types[argname] == "PatternSelectionStrategy") {
                    arg_values[argname] = argname;
                    if (printargs) std::cout << arg_types[argname] << ": " << std::any_cast<string>(arg_values[argname]) << endl;
                }
                else if (arg_types[argname] == "SolutonFillingStrategy") {
                    arg_values[argname] = argname;
                    if (printargs) std::cout << arg_types[argname] << ": " << std::any_cast<string>(arg_values[argname]) << endl;
                }
            }
        }
    }

     if (printargs) {
        cout << "--------------------------------" << endl;
    }

}

auto ProgramArgs::getArgs()
{
}
