#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include <functional>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <vector>
#include "preconditioner.hpp"

std::vector<double>
cgsolve(const CooMatrix<double>&   A,
	const std::vector<double>& b, const int& k) {
    
  assert((NbCol(A)==NbRow(A)) &&
	 (b.size()==NbCol(A)) );
 
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    p   = r;
  auto   Ap   = A*p;    
  double r2   = std::pow(Norm(r),2);  
  double eps2 = 1e-6;//(1e-8)*r2;
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;
    
  //std::size_t niter = 0;    
  int niter = 0;
  //while( (r2>eps2 && niter++<1000) && niter<=k){
  //while(r2>eps2 && niter<=k){
  while( niter++<1000 && niter<=k){
      
    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = r2/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;    
    r2    = std::pow(Norm(r),2);
    beta  = r2/(alpha*pAp);
    p     = beta*p+r;

    //niter++;
    // if((niter%50)==0){
    //   std::cout << std::left << std::setw(7) << niter << "\t";
    //   std::cout << Norm(r) << std::endl;
    // }
    
  }
    
  return x;
    
}


std::vector<double>
cgsolve_modifie(const CooMatrix<double>&   A,
	const std::vector<double>& b, const int& k, const std::vector<double>& u) {
    
  assert((NbCol(A)==NbRow(A)) &&
	 (b.size()==NbCol(A)) );
 
  auto    x   = u;//std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    p   = r;
  auto   Ap   = A*p;    
  double r2   = std::pow(Norm(r),2);  
  double eps2 = 1e-6;//(1e-8)*r2;
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;
    
  int niter = 0;    
  //while( (r2>eps2 && niter++<1000) && niter<=k){
  //while(r2>eps2 && niter<=k){
  while( niter++<1000 && niter<=k){
      
    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = r2/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;    
    r2    = std::pow(Norm(r),2);
    beta  = r2/(alpha*pAp);
    p     = beta*p+r;

    //niter++;
    // if((niter%50)==0){
    //   std::cout << std::left << std::setw(7) << niter << "\t";
    //   std::cout << Norm(r) << std::endl;
    // }
    
  }
    
  return x;
    
}

std::vector<double> PCGSolver(const CooMatrix<double>& A, const std::vector<double>& b, const int& k) {


  auto P = CholeskyPrec(A);
  auto x = std::vector<double>(b.size(),0.);
  auto r = b - A*x;
  auto r_new = r;
  auto z = P*r;
  auto z_new = z;
  auto p = z;
  auto Ap = std::vector<double>(b.size(),0.);

  double r2 = std::pow(Norm(r),2);
  double eps2 = 1e-6;//(1e-8)*r2;
  eps2 *= std::abs((b|b));      
  double alpha,beta,pAp;
    
  int niter = 0;
  while ((r2>eps2 && niter++<1000) & niter<k) {
    Ap = A*p;
    pAp = std::real((p|Ap));
    alpha = (r|z)/pAp;
    x += alpha*p;
    r_new -= alpha*Ap;
    z_new = P*r_new;
    beta = (r_new|z_new)/(r|z);
    p = z_new + beta*p;
    z = z_new;
    r = r_new;
    // r2 = std::pow(Norm(r),2);
  }

  return x;
}




#endif
