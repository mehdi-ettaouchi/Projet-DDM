#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include <mpi.h>
#include <vector>

#include "Problem.h"
#include "GradConj.h"
#include "BC.h"
#include "Output.h"
#include "Readfile.h"

/*-----------------------------------------------------
Une implémentation parallèledu gradient conjugué et de la résolution de l'équation. Les commentaires dans ce code pourront porter des précisions éventuelles au rapport.

//Sommaire de la parralèlisation//

On parallèlise sur Nx*Ny, on ne porte donc aucune contrainte sur Ny et le nombre de proc (Ny peut etre indivisible par Nproc), on ditribue ainsi la charge totlae qui est Nx*Ny au procs et effectue les opérations séparables (somme, norme) indpendemment. A chaque itération pour le produit, on aura besoin (selon les cas) d'un nombre d'éléments différents du vecteur x. Ceci peut etre fait des procs imméditement voisins ou plus loin. Le code est adapté à faire cela en "mémorisant" la charge relative de chaque proc et en la comparant à la posiition de l'indice qui balaye le bloc actuel. En arrivant à lapartie où il faut calculer le pas à chaque itération une réduction (en somme de carrée) du pas dans tous les procs All_Reduce. A la convergence, chaque proc écrit sa propre solution dans un fichier à part.
*/

using namespace std;


//Macros de présentation et de debuggage
#define bloc std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
#define SHOW(a) std::cout << #a << std::endl;
#define debug(i) std::cout << "debug" << i << std::endl;
#define siz(a) std::cout << a.size() << std::endl;

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
  cout << "la sous diagonale de la matrice:" << endl;
  cout << "[ ";
  for (int i = 0; i < A[3].size(); i++)
  {
    cout << A[3][i] << " ";
  }
  cout << " ]" << endl;
  cout << "la sous-sous diagonale de la matrice:" << endl;
  cout << "[ ";
  for (int i = 0; i < A[4].size(); i++)
  {
    cout << A[4][i] << " ";
  }
  cout << " ]" << endl;
  cout << endl;
}

void print_matrix_verbose(std::vector<std::vector<double>> A)
{
  cout << "version verbose de la matrice creuse" << endl;
  std::cout << "-------------------------------------------------" << std::endl;
  if (A.size() != 5)
  {
    cout << "mauvais format de matrice" << endl;
  }
  else
  {
    //retrouvant Nx et Ny
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

//Fonction de distribution optimale de charge entre les procs
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
  int i(0);

  MPI_Init(&argc, &argv);

  double t1(0.), t2(0.);
  t1 = MPI_Wtime();

  //récuperation du nom du fichier d'entrée 
  const string data_file_name = argv[1];
  bloc
  cout<< data_file_name << endl;
  Readfile *Rf = new Readfile(data_file_name);
  Rf->Read_data_file();
  // ------------------------------------------------------------
  //Récupération des données du problème
  double Lx = Rf->Get_Lx(), Ly = Rf->Get_Ly(), D = Rf->Get_D(), deltat = Rf->Get_dt(), tf = Rf->Get_tfinal(), t(0.);
  int Nx = Rf->Get_Nx(), Ny = Rf->Get_Ny(), Nt = 4; //Nt ou deltat jouent le meme role et s'adaptent entre elles si l'une est renseignée au lieu de l'autre.


  //initialisation du parallélisme et des différentes variables dont on aura besoin (tag, indices, pas d'espace, variables de procs..)
  int me, Np, tag, input, begin, end;
  double deltax = Lx / (Nx + 1), deltay = Ly / (Ny + 1);
  double betax = -deltat / pow(deltax, 2), betay = -deltat / pow(deltay, 2);
  double alpha = 1 - 2 * (betax + betay);
  tag = 100;

  std::vector<int> v(2);
  //v Contiendra les charges


  BC bc = BC(Nx, Ny, Lx, Ly);
  Problem P = Problem(&bc, Nx, Ny, Nt, Lx, Ly, deltat);

  MPI_Comm_size(MPI_COMM_WORLD, &Np);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  v = charge(Nx * Ny, Np, me);

  int size = v[1] - v[0] + 1;
  printf("size : %d", size);

  std::vector<std::vector<double>> C(5, vector<double>(size));
  std::vector<double> x(size, 0.), w(size, 2.), f(size, 0.), h(size, 0.);

  //on souhaite faire une somme de deux vecteurs, suffit d'appeler les fonctions somme sur la partition du vecteur selon charge, de meme pour norme et produit scalaire etc

  //connaitre le rang sur la matrice NxNy x NxNy
  int rang = 0;
  for (int k = 0; k < me; k++)
  {
    rang += charge(Nx * Ny, Np, k)[1] - charge(Nx * Ny, Np, k)[0] + 1;
  }

  // construction du second membre --------------------------
  int reste, quotient, reste0, quotient0;
  reste = rang % Nx;
  quotient = (rang - reste) / Nx;
  reste0 = reste;
  quotient0 = quotient;

  int cas = Rf->Get_cas();
  std::cout << "le cas utilisé est" << cas << std::endl;
  //Que le source terme 
  for (int iter = 0; iter < size; iter++)
  {
    f[iter] = bc.Source_term(reste * deltax, quotient * deltay, t, cas);

    //Cette boucle sera utiliser pas mal de fois pour remettre les indices dans l'ordre
    reste += 1;
    if (reste > Nx - 1)
    {
      reste = 0;
      quotient++;
    }
  }
  

  //La construction des conditions au bords:
  int rangp = rang + size;
  int rs1 = rang % Nx;
  int q1 = (rang - rs1) / Nx;
  int rs2 = rangp % Nx;
  int q2 = (rangp - rs2) / Nx;
  int box = q2 - q1;
  /*   cout<<"this is box "<<box<<endl;
    cout<<"this is rangp "<<rang<<endl;
    cout<<"this is rangp "<<rangp<<endl;
    cout<<"this is q1 and q2 "<<q1<<" "<<q2<<endl;
    cout<<"this is rs1 and rs2 "<<rs1<<" "<<rs2<<endl; */
  int j;
  for (int i = q1; i <= q2; i++)
  {
    for (int j = 0; j < Nx; j++)
    {
      if (i == q1 && j >= rs1 || i == q2 && j <= rs2 || i != q1 && i != q2)
      {
        if (i == 0 && j == 0)
        {
          f[i * Nx + j - rang] += D * bc.Dirichlet_Function0(j * deltax, i * deltay, t, cas) / (deltax * deltax) + D * bc.Dirichlet_Function1(j * deltax, i * deltay, t, cas) / (deltay * deltay);
        }
        if (i == 0 && j != 0 && j != Nx - 1)
        {
          f[i * Nx + j - rang] += D * bc.Dirichlet_Function1(j * deltax, i * deltay, t, cas) / (deltay * deltay);
        }
        if (i == 0 && j == Nx - 1)
        {
          f[i * Nx + j - rang] += D * bc.Dirichlet_Function0(j * deltax, i * deltay, t, cas) / (deltax * deltax) + D * bc.Dirichlet_Function1(j * deltax, i * deltay, t, cas) / (deltay * deltay);
        }
        if (i == Ny - 1 && j == 0)
        {
          f[i * Nx + j - rang] += D * bc.Dirichlet_Function0(j * deltax, i * deltay, t, cas) / (deltax * deltax) + D * bc.Dirichlet_Function1(j * deltax, i * deltay, t, cas) / (deltay * deltay);
        }
        if (i == Ny - 1 && j == Nx - 1)
        {
          f[i * Nx + j - rang] += D * bc.Dirichlet_Function0(j * deltax, i * deltay, t, cas) / (deltax * deltax) + D * bc.Dirichlet_Function1(j * deltax, i * deltay, t, cas) / (deltay * deltay);
        }
        if (i == Ny - 1 && j != 0 && j != Nx - 1)
        {
          f[i * Nx + j - rang] += D * bc.Dirichlet_Function1(j * deltax, i * deltay, t, cas) / (deltay * deltay);
        }
        if (j == 0 && i != 0 && i != Ny - 1)
        {
          f[i * Nx + j - rang] += D * bc.Dirichlet_Function0(j * deltax, i * deltay, t, cas) / (deltax * deltax);
        }
        if (j == Nx - 1 && i != 0 && i != Ny - 1)
        {
          f[i * Nx + j - rang] += D * bc.Dirichlet_Function0(j * deltax, i * deltay, t, cas) / (deltax * deltax);
        }
      }
    }
  }

  //construction de la matrice, cette fois liste (vecteur) à 5 vecteurs.
  for (int i = 0; i < size; i++)
  {
    if (rang + i - Nx < 0)
    {
      C[0][i] = 0;
    }
    else
    {
      C[0][i] = betay;
    }
    if (rang + i + Nx < Nx * Ny)
    {
      C[4][i] = betay;
    }
    else
    {
      C[4][i] = 0;
    }
    if (rang + i - 1 < 0 || ((rang + i) % Nx) == 0)
    {
      C[1][i] = 0;
    }
    else
    {
      C[1][i] = betax;
    }
    if (rang + i + 1 > Nx * Ny || ((rang + i) % Nx) == Nx - 1)
    {
      C[3][i] = 0;
    }
    else
    {
      C[3][i] = betax;
    }
    C[2][i] = alpha;
  }

  //
  //ébut de la zone integrée gradient conjuguée-----------------
  int n = size;
  int k_ = Nx * Ny;
  std::vector<std::vector<double>> A(C);
  std::vector<double> r1(n), b1(n), p1(n), temp(n);
  double alpha1;
  double gamma;
  double alphan, alphad, gamman;
  std::vector<double> rSuivant(n);
  std::vector<double> xSuivant(n);
  std::vector<double> z1(n);

  double beta;
  int nb_iterat_;
  MPI_Status Status;
  int q, r;
  bool a, b, c, d;
  q = Nx / size;
  r = Nx - q * size;

  for (double it_t = 0; it_t <= tf; it_t += deltat)
  {
    j = 0;
    nb_iterat_ = 0;
    b1 = GradConj::sum(x, GradConj::prod_scal(f, deltat), 1);
    r1 = b1;
    p1 = r1; // calcul du residu
    beta = pow(GradConj::norm(r1), 2);
    MPI_Allreduce(&beta, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    beta = sqrt(beta);
    while (j < k_)
    {
      //Les vecteurs dénominés tool serveront de vecteurs temp pour inverser les vecteurs envoyés ou rvus qui sont en effet stacker en pile dans les variables d'entrée.
      vector<double> z(Nx, 0.), y(Nx, 0.), prod(size, 0.), ztool(size), ztool1(r);
      x = p1;

      //me va envoyer ses éléments aux procs qui en a besoin
      for (int k = 1; k < q + 2; k++)
      {
        if (me + k < Np)
        {
          if (k == q + 1)
          {
            MPI_Send(&x[size - r], r, MPI_DOUBLE, me + k, 0, MPI_COMM_WORLD);
            MPI_Recv(&y[(k - 1) * size], r, MPI_DOUBLE, me + k, 0, MPI_COMM_WORLD, &Status);
          }
          else
          {
            MPI_Send(&x[0], size, MPI_DOUBLE, me + k, Np, MPI_COMM_WORLD);

            MPI_Recv(&y[(k - 1) * size], size, MPI_DOUBLE, me + k, 0, MPI_COMM_WORLD, &Status);
          }
        }
        if (me - k > -1)
        {
          if (k == q + 1)
          {
            MPI_Send(&x[0], r, MPI_DOUBLE, me - k, 0, MPI_COMM_WORLD);

            MPI_Recv(&z[(k - 1) * size], r, MPI_DOUBLE, me - k, 0, MPI_COMM_WORLD, &Status);

            for (int itr = 0; itr < r; itr++)
            {
              if (r - itr < itr)
              {
                z[(k - 1) * size + itr] = ztool1[(k - 1) * size + r - itr - 1];
              }
              else
              {
                ztool1[(k - 1) * size + itr] = z[(k - 1) * size + itr];
                z[(k - 1) * size + itr] = z[(k - 1) * size + r - itr - 1];
              }
            }
          }
          else
          {
            MPI_Send(&x[0], size, MPI_DOUBLE, me - k, Np, MPI_COMM_WORLD);

            MPI_Recv(&z[(k - 1) * size], size, MPI_DOUBLE, me - k, 0, MPI_COMM_WORLD, &Status);

            for (int itr = 0; itr < size; itr++)
            {
              if (size - itr < itr)
              {
                z[(k - 1) * size + itr] = ztool[k * size - itr - 1];
              }
              else
              {
                ztool[(k - 1) * size + itr] = z[(k - 1) * size + itr];
                z[(k - 1) * size + itr] = z[k * size - itr - 1];
              }
            }
          }
        }
      }

      //Produit matrice vecteur adapté à la matrice creuse
      for (int i = 0; i < size; i++)
      {
        a = (i - Nx > -1);
        b = (i + Nx < size);
        c = (i > 0);
        d = (i + 1 < size);

        prod[i] += C[2][i] * x[i];
        if (a)
        {
          prod[i] += C[0][i] * x[i - Nx];
        }
        else
        {
          prod[i] += C[0][i] * z[Nx - i - 1];
        }
        if (c)
        {
          prod[i] += C[1][i] * x[i - 1];
        }
        else
        {
          prod[i] += C[1][i] * z[0];
        }
        if (b)
        {
          prod[i] += C[4][i] * x[i + Nx];
        }
        else
        {
          prod[i] += C[4][i] * y[Nx + i - size];
        }
        if (d)
        {
          prod[i] += C[3][i] * x[i + 1];
        }
        else
        {
          prod[i] += C[3][i] * y[0];
        }
      }
      z1 = prod;
      if (nb_iterat_ > 0)
      {
        x = xSuivant;
      }
      else
      {
        for (int i = 0; i < size; i++)
        {
          x[i] = 0.;
        }
      }

      alphan = pow(beta, 2);
      alphad = GradConj::dot_product(z1, p1);

      //Réduction du pas entre les procs, une réduction a été aussi au début de chaque itération
      MPI_Allreduce(&alphad, &alphad, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      alpha1 = alphan / alphad;

      xSuivant = GradConj::sum(x, GradConj::prod_scal(p1, alpha1), 1);
      rSuivant = GradConj::sum(r1, GradConj::prod_scal(z1, alpha1), -1);

      gamman = GradConj::dot_product(rSuivant, rSuivant);

      MPI_Allreduce(&gamman, &gamman, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


      gamma = gamman / alphan;
      p1 = GradConj::sum(rSuivant, GradConj::prod_scal(p1, gamma), 1);
      x = xSuivant;
      r1 = rSuivant;
      beta = sqrt(gamman);

      nb_iterat_ = nb_iterat_ + 1;
      if (beta < pow(10, -10))
      {
        break;
      }
      j++;
    }
  }


  //Plot des solutions de chaque proc
  Output io = Output(&P);

  ofstream myfile;
  string st = "solution_from_proc" + to_string(me) + ".txt";
  myfile.open(st);
  double x1, y1;

  for (int j = quotient0; j < quotient + 1; j++)
  {
    if (j > quotient0)
    {
      reste0 = 0;
    }
    for (int i = reste0; i < Nx; i++)
    {

      x1 = i * deltax;
      y1 = j * deltay;
      bool a = (f[i + j * Nx - rang] > -1);
      myfile << x1 << " " << y1 << " " << a * x[i + j * Nx - rang] << endl;
      if (j == quotient && i > reste)
      {
        break;
      }
    }
  }

  myfile.close();

  t2 = MPI_Wtime();

  bloc

  cout<< "The thread " << me << " computed the solution in " << t2 - t1 << endl;

  MPI_Finalize();

  return 0;
}
