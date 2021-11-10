#ifndef _IO_CPP

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "Output.h"
#include <stdio.h>

using namespace std;

/*-----------------------------------------------------
Classe de sauvegarde et de construction du fichier à charger dans gnuplot directement pour plus d'accessibilité.
*/

Output::Output(Problem *Pb) : P_(Pb)
{
  std::cout << "Classe IO initié" << std::endl;
}

//Sauvegarde de solution
void Output::Save_sol(std::string st)
{
  int Nx = P_->get_Nx();
  int Ny = P_->get_Ny();
  std::vector<double> sol(Nx * Ny);
  sol = P_->get_sol();
  double x, y, dx = P_->get_dx(), dy = P_->get_dy();
  if (sol.size() > Nx * Ny) //in seq !=
  {
    cout << "solution de mauvaise taille " << sol.size() << endl;
  }
  else
  {
    ofstream myfile;
    myfile.open(st);
    for (int i = 0; i < Ny; i++)
    {
      for (int j = 0; j < Nx; j++)
      {
        x = j * dx;
        y = i * dy;
        myfile << x << " " << y << " " << sol[j + i * Nx] << endl;
      }
    }
    myfile.close();
  }
}

void Output::splot_solution(std::string sol_file_name)
{
  //Pour ne pas refaire sur gnuplot un tracé plusieurs fois, surtout en testant.
  //suffit de load "cmd_file_line" dans gnuplot
  string command_filename = "cmd_file_line";
  ofstream command_unit;
  command_unit.open(command_filename.c_str());
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Sol --->'\n";
  command_unit << "set title 'Solution'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "splot '" << sol_file_name
               << "' using 1:2:3 w p \n";
}

#define _IO_CPP
#endif
