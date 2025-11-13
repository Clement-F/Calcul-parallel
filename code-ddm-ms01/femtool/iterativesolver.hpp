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
	const std::vector<double>& b, const int& k_temp=1000) {
    
  std::size_t k =  std::size_t(k_temp);
    
  assert((NbCol(A)==NbRow(A)) &&
	 (b.size()==NbCol(A)) );
 
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    p   = r;
  auto   Ap   = A*p;    
  double r2   = std::pow(Norm(r),2);  
  double eps2 = (1e-8)*r2;
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;
    
  std::size_t niter = 0;    
  while(( r2>eps2 && niter++<1000) && niter<k){
      
    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = r2/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;    
    r2    = std::pow(Norm(r),2);
    beta  = r2/(alpha*pAp);
    p     = beta*p+r;
      
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
    }
      
  }
    
  return x;
    
}

// std::vector<double>
// PCGSolver(const CooMatrix<double>&   A, const CooMatrix<double>&   Q, const std::vector<double>& b, 
//   const int& k=1000, const double& epsi=10e-6)
// {
//   auto    x   = std::vector<double>(b.size(),0.);
//   auto    r   = b - A*x;
//   auto    z   = Q*r;
//   auto    p   = z;
//   auto   Ap   = A*p;   

//   auto r_prec = r;  auto z_prec = z;
//   double alpha, beta;
//   double r2   = std::pow(Norm(r),2);  

//   std::size_t niter = 0;  
//   while((r2>epsi && niter++<1000) && niter<k)
//   {
//     r_prec = r; z_prec = z;
//     Ap   = A*p;  
//     alpha = (r|z)/(p|Ap);
//     x     = x + alpha*p;
//     r     = r - alpha*Ap;
//     z     = Q*r;
//     beta  = (r|z)/(r_prec|z_prec);
//     p     = z + beta*p; 
//     r2    = std::pow(Norm(r),2);
    
//     if((niter%50)==0){
//       std::cout << std::left << std::setw(7) << niter << "\t";
//       std::cout << Norm(r) << std::endl;
//     }
//   }
//   return x;
    
// }



std::vector<double>
PCGSolver(const CooMatrix<double>&   A, const std::vector<double>& b, 
  const int& k_temp=1000, const double& epsi=10e-6)
{

  std::size_t k =  std::size_t(k_temp);

  CholeskyPrec Q(A);
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b - A*x;
  auto    z   = Q*r;
  auto    p   = z;
  auto   Ap   = A*p;   

  auto r_prec = r;  auto z_prec = z;
  double alpha, beta;
  double r2   = std::pow(Norm(r),2);  
  double eps2 = epsi*r2;
  eps2       *= std::abs((b|b));      

  std::size_t niter = 0;  
  while((r2>eps2 && niter++<1000) && niter<k)
  {
    r_prec = r; z_prec = z;
    Ap   = A*p;  
    alpha = (r|z)/(p|Ap);
    x     = x + alpha*p;
    r     = r - alpha*Ap;
    z     = Q*r;
    beta  = (r|z)/(r_prec|z_prec);
    p     = z + beta*p; 
    r2    = std::pow(Norm(r),2);
    
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
    }
  }
  return x;
    
}


#endif
