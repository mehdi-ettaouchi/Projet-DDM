#ifndef _BC_CPP

#include "Problem.h"
#include "BC.h"
#include <math.h>
#include <cmath>

/*-----------------------------------------------------
Classe des differentes fonctions numériques utilisées. fonctions de bords, terme source et conditions initiales.
*/

#define _USE_MATH_DEFINES

BC::BC(int Nx, int Ny, double Lx, double Ly) : Lx_(Lx), Ly_(Ly), Nx_(Nx), Ny_(Ny)
{
  std::cout << "La lecture des fonctions relatives au problème" << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
}

double BC::Source_term(const double x, const double y, const double t, const int cas) const
{
  if (cas == 4)
  {
    return 2 * (y - y * y + x - x * x);
  }
  else if (cas == 5)
  {
    return sin(x) + cos(y);
  }
  else if (cas == 6)
  {
    return exp(-(x - Lx_ / 2) * (x - Lx_ / 2)) * exp(-(y - Ly_ / 2) * (y - Ly_ / 2)) * cos(M_PI * t / 2);
  }
  else
  {
    std::cout << "choix indisponible " << cas << std::endl;
  }
}

double BC::Initial_condition(const double x, const double y, const double t) const
{
  return 0.;
}

double BC::Dirichlet_Function0(const double x, const double y, const double t, int cas) const
{
  if (cas == 4)
  {
    return 0.;
  }
  else if (cas == 5)
  {
    return sin(x) + cos(y);
  }
  else if (cas == 6)
  {
    return 0.;
  }
  else
  {
    std::cout << "fonctions de dirichlet non sécifié" << std::endl;
  }
}

double BC::Dirichlet_Function1(const double x, const double y, const double t, int cas) const
{
  if (cas == 4)
  {
    return 0.;
  }
  else if (cas == 5)
  {
    return sin(x) + cos(y);
  }
  else if (cas == 6)
  {
    return 1.;
  }
  else
  {
    std::cout << "fonctions de dirichlet non sécifié" << std::endl;
  }
  return 0.;
}

double BC::Exact_solution(const double x, const double y, const double t) const
{
  return 0.;
}

#define _BC_CPP
#endif
