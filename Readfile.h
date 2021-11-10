#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>

class Readfile
{
private:
  std::string file_name_;
  double tfinal_, dt_, Lx_, Ly_, D_,n_,alpha_,beta_;
  int Nx_, Ny_, cas_;

  std::string results_;

  bool print_info_;
  bool if_tfinal_;
  bool if_dt_;
  bool if_results_;
  bool if_D_;
  bool if_n_;
  bool if_Lx_;
  bool if_Ly_;
  bool if_Nx_;
  bool if_Ny_;
  bool if_cas_;
  bool if_alpha_;
  bool if_beta_;

public:
  //constructor
  Readfile(std::string file_name);

  //ones to define
  void Read_data_file();

  void Assembel_sol_file(int);

  //header defined
  void Adapt_dt(double dt) { dt_ = dt; };
  const bool Print_info() const { return print_info_; };
  const double Get_tfinal() const { return tfinal_; };
  const double Get_dt() const { return dt_; };
  const double Get_D() const { return D_; };
  const double Get_Nx() const { return Nx_; };
  const double Get_Ny() const { return Ny_; };
  const double Get_Lx() const { return Lx_; };
  const double Get_Ly() const { return Ly_; };
  const double Get_n() const { return n_;};
  const int Get_cas() const { return cas_; };
  const int Get_beta() const { return beta_; };
  const int Get_alpha() const { return alpha_; };
  const std::string Get_results() const { return results_; };
};

#define _DATA_FILE_H
#endif