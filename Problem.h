#ifndef PBDEF
#define PBDEF


#include <fstream>
#include <iostream>
#include <vector>

#include "GradConj.h"
#include "BC.h"


class Problem
{
private:
    const BC* functions_;
    int Nx_,Ny_,Nt_,n_;
    double Lx_,Ly_,t_,Lx0_;
    double deltax_=Lx_/(Nx_+1.);
    double deltay_=Ly_/(Ny_+1.);
    double deltat_=t_/(Nt_+1.);
    double alpha_,beta_;
    const double D_=1.;
    int procID_,r_,Np_;
    std::vector<double> F_,Bd_;
    std::vector<std::vector<double>> A_;
    std::vector<double> sol_;
    std::vector<double> stencil1_,stencil2_;



public:
    //constructor
    Problem(BC* functions,int , int , int , double , double, double,int,int,int,double,double,double,double,std::vector<double>,std::vector<double>);
    //destructor
    ~Problem(){};
    //problem related
    std::vector<std::vector<double>> Construct_Matrix();
    void Construct_F(int,double, std::vector<double>&);
    void Construct_Bd(int,double);
    std::vector<double> Constructe_Bd(int,double);
    void Solve_problem(int cas, double tf);
    //get solution and parameters
    std::vector<double> get_sol();
    std::vector<double> getLeftStencil();
    std::vector<double> getRightStencil();
    double get_Lx() const {return Lx_;} ;
    double get_Ly() const {return Ly_;} ;
    double get_Nx() const {return Nx_;} ;
    double get_Ny() const {return Ny_;} ;
    double get_dx() const {return deltax_;} ;
    double get_dy() const {return deltay_;} ;
    double get_n() const {return n_;};
    double get_alpha() const {return alpha_;};
    double get_beta() const {return beta_;};
    std::vector<double> get_stencil1() const {return stencil1_ ;};

};


#endif
