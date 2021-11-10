#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include <iomanip>
#include <string>
#include "Problem.h"
#include "GradConj.h"
#include "BC.h"
#include "Output.h"
#include "Readfile.h"
#include "mpi.h"

using namespace std;

//macros pour la présentation et le debuggage

#define bloc std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
#define SHOW(a) std::cout << #a << std::endl;

//Fonction de réprensetation de vecteurs de la librairie vector.h

void print_vector(std::vector<double> x)
{
  int n = x.size();
  cout << "le vecteur de taille " << n << endl;

  for (int i = 0; i < n; i++)
  {
    cout << x[i] << " ";
  }
  cout << endl;
  cout << "----------------------------------" << endl;
}

//Fonctions de représentation de matrice (vecteur de vecteurs) de forme abrégée (liste de vecteurs) ou verbose (matrice carrée)

void print_matrix(std::vector<std::vector<double>> A)
{
  cout << "la diagonale de la matrice:" << endl;
  cout << "[ ";
  for (int i = 0; i < A[0].size(); i++)
  {
    cout << A[0][i] << " ";
  }
  cout << " ]" << endl;
  cout << "la sur diagonale de la matrice:" << endl;
  cout << "[ ";
  for (int i = 0; i < A[1].size(); i++)
  {
    cout << A[1][i] << " ";
  }
  cout << " ]" << endl;
  cout << "la sur-sur diagonale de la matrice:" << endl;
  cout << "[ ";
  for (int i = 0; i < A[2].size(); i++)
  {
    cout << A[2][i] << " ";
  }
  cout << " ]" << endl;
  cout << endl;
}
void print_matrix_verbose(std::vector<std::vector<double>> A)
{
  cout << "version verbose de la matrice creuse" << endl;
  std::cout << "-------------------------------------------------" << std::endl;
  if (A.size() != 3)
  {
    cout << "mauvais format de matrice" << endl;
  }
  else
  {
    //Retrouvant Nx et Ny
    int Ny = A[0].size() / (A[0].size() - A[2].size());
    int Nx = A[0].size() / Ny;
    cout << Nx << " " << Ny << endl;
    for (int i = 0; i < Nx * Ny; i++)
    {
      for (int j = 0; j < Nx * Ny; j++)
      {
        if (i == j)
        {
          cout << A[0][i] << " ";
        }
        else if ((j == i + 1) || (j == i - 1))
        {
          cout << A[1][i] << " ";
        }
        else if (j == Nx + i)
        {
          cout << A[2][i];
        }
        else
        {
          cout << "0 ";
        }
      }
      cout << endl;
    }
  }
}

//fonctions d'équilibrage de charges
vector<int> charge(int n, int Np, int me)
{
  int limite = n - Np * (n / Np);
  vector<int> res(2);

  if (me < limite)
  {
    res[0] = me * (n / Np + 1);
    res[1] = res[0] + n / Np;
  }
  else
  {
    res[0] = limite * (n / Np + 1) + (me - limite) * (n / Np);
    res[1] = res[0] + (n / Np) - 1;
  }
  return res;
}

int main(int argc, char **argv)
{

  // Démarrage du chrono
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

  if (argc < 2)
  {
    cout << "Ooops, forgot to give the data_file name" << endl;
    exit(0);
  }

  //Récuperation le nom du fichier d'entrée
  const string data_file_name = argv[1];
  bloc
          cout
      << data_file_name << endl;
  Readfile *Rf = new Readfile(data_file_name);
  Rf->Read_data_file();
  // ------------------------------------------------------------

  //Récupération des données du problème
  double Lx = Rf->Get_Lx(), Ly = Rf->Get_Ly(), D = Rf->Get_D(), deltat = Rf->Get_dt(), tf = Rf->Get_tfinal(), alpha = Rf->Get_alpha(), beta = Rf->Get_beta();
  int Nx = Rf->Get_Nx(), Ny = Rf->Get_Ny(), Nt = 4, n = Rf->Get_n(), cas = Rf->Get_cas();
  // Vérification de la suffisance des conditions de bords
  if(alpha==0 && beta==0){
    std::cout<< " conditions de bords insuffisantes"<< endl;
    return 0; 
  }
  //Nt ou delta t jouent le meme role
  int nb_iter = 0;
  // Initialisation
  MPI_Init(&argc, &argv);
  int rank, Np;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &Np);
  Nx /= Np;
  // Conditons de recouvrement initiales
  std::vector<double> stencil1(Ny - 2, 0.0), stencil2(Ny - 2, 0.0);
  ofstream SolFile, SpeedupFile;
  std::vector<double> y(Nx * Ny, 2.);
  SpeedupFile.open("Speedup_Pour Np=" + std::to_string(Np));
  SpeedupFile << "# n  speedup" << endl;
  int nbegin,nend;
  if (alpha!=0){
    nbegin=3;
    nend=Nx-1;
  }
  else {
    nbegin=1;
    nend=Nx;
  }
  for (n = nbegin; n < nend; n++)
  {
    double t1 = MPI_Wtime();
    double domXsize = Lx / (Np - n * (Np - 1) / (Nx + 1));
    double deltax = domXsize / (Nx + 1.0);
    double deltay = Ly / (Ny + 1);
    double recouv = n * domXsize / (Nx + 1);
    bool Converged = false;
    // Déclaration d'une classe gc pour pouvoir utiliser les ops sur matrices/vecteurs
    std::vector<std::vector<double>> B(3);
    std::vector<double> g(Nx * Ny, 4.);
    GradConj mc(B, g, Nx, Ny);
    // Réinitialisation du nombre d'itérations
    nb_iter=0;
    // Boucle de convergence de la méthode de recouvrement
    while (!Converged)
    {
      //Initialisation des différentes classes nécessaires à la résolution dans la classe problem
      BC bc = BC(Nx, Ny, domXsize, Ly);
      Problem P = Problem(&bc, Nx, Ny, Nt, domXsize, Ly, deltat, rank, Np, n, recouv, Lx, alpha, beta, stencil1, stencil2);
      //Résolution du problème cas
      P.Solve_problem(cas, tf);
      y = P.get_sol();
      // Exception du  cas à 1 proc
      if (Np == 1)
      {
        Converged = true;
        break;
      }
      // construction des stencils
      std::vector<double> stencilL = P.getLeftStencil();
      std::vector<double> stencilR = P.getRightStencil();
      // Communication entre procs
      if (rank != 0 && rank != Np - 1)
      {
        MPI_Status status;
        MPI_Send(&stencilL[0], Ny - 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        MPI_Send(&stencilR[0], Ny - 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&stencilL[0], Ny - 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&stencilR[0], Ny - 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
      }
      else
      {
        if (rank == 0)
        {
          MPI_Status status;
          MPI_Send(&stencilR[0], Ny - 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
          MPI_Recv(&stencilR[0], Ny - 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
        }
        if (rank == Np - 1)
        {
          MPI_Status status;
          MPI_Send(&stencilL[0], Ny - 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
          MPI_Recv(&stencilL[0], Ny - 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
      }
      // Rajout des termes dépendent de la solution du meme proc
      if (alpha != 0)
      {
        for (int i = 1; i < Ny - 1; i++)
        {
          stencilL[i - 1] +=y[i*Nx+2]-(2.*beta*deltax/alpha)*y[i*Nx+1];
          stencilR[i-1]+=y[i*Nx+Nx-3]-(2.*beta*deltax/alpha)*y[i*Nx+Nx-2];
        }
      }
      // Actualiser le nombre d'itérations
      nb_iter++;
      //vérification de la condition de l'erreur
      double err = 0;
      if (rank != 0 && rank != Np - 1)
      {
        err += max(mc.normMax(mc.sum(stencil1, stencilL, -1)), mc.normMax(mc.sum(stencil2, stencilR, -1)));
      }
      else
      {
        if (rank == 0)
        {
          err += mc.normMax(mc.sum(stencil2, stencilR, -1));
        }
        else
        {
          err += mc.normMax(mc.sum(stencil1, stencilL, -1));
        }
      }
      MPI_Allreduce(&err, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      if (err < pow(10, -4))
      {
      std:
        cout << " le nombre d'itérations pour converger est de " << nb_iter << endl;
        Converged = true;
      }
      else
      {
        std::cout << "on est à l'iteration " << nb_iter << "avec une erreur de " << err << endl;
        stencil1 = stencilL;
        stencil2 = stencilR;
      }
    }

    // Sauvegarde de la solution
    SolFile.open("Sol" + std::to_string(rank) + ".dat");
    SolFile << "# x      y    u(x,y)" << endl;
    for (int i = 0; i < Ny; i++)
    {

      for (int k = rank * (n / Np); k < Nx - (Np - rank - 1) * (n / Np); k++)
      {
        SolFile << rank * deltax * (Nx - n) + k * deltax << " " << i * deltay << " " << y[i * Nx + k] << endl;
      }
    }
    SolFile.close();

    // Fin du chrono
    double t2 = MPI_Wtime();
    double s = t2 - t1;
    MPI_Allreduce(&s, &s, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // Sauvegarder le speedup
    if (rank == 0)
    {

      SpeedupFile << n << " " << 6.92 / s << endl;
    }
  }
  SpeedupFile.close();
  // Affichage du résultat
  // std::cout << "-------------------------------------------------" << std::endl;
  // cout << "Cela a pris " << (t2 - t1) << " seconds"
  //      << "pour le proc " << rank << endl;
  MPI_Finalize();
  return 0;
}
