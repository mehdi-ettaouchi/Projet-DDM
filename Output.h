#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "Problem.h"



class Output
{
    private:
        Problem* P_;
    public:
        Output(Problem*);
        ~Output(){};
        void Save_sol(std::string st);
        void splot_solution(std::string);
};


#endif