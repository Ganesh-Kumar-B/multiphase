#pragma once

#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<fstream>
#include<algorithm>
#include <sstream>
#include<string>
#include "lbmD3Q19.h"
#include "GRID_3D.h"
#define PI 3.14159265



template<typename T, typename T1>
void collide(Grid_N_C_3D<T> &grid,Grid_N_C_3D<T> &rho,Grid_N_C_3D<T> &pnid,Grid_N_C_3D<T> &fnid,Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho,
            lbmD3Q19<T1> &lb,double beta,double tau, double TbyTc, double kappa, int t ){

    Grid_N_C_3D<T>  laplacian_pnidplusfnidbyrho            (grid.n_x,grid.n_y,grid.n_z,1,1);
    
    double feq_Node[19] = {0}, ux = 0, uy = 0, uz = 0;

    double rho_critical = 1.0, T_critical = lb.theta0/TbyTc ; 
    double b = 1.0/(3.0*rho_critical), a = b*T_critical*27.0/8.0;
    kappa = kappa*a;   //!scale it properly with the delxs

    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){
            
                double Rho = 0;
                get_moments(grid, lb,  ux, uy, uz,Rho, i, j ,k); 
                rho.Node(i,j,k) = Rho;
            }
        }
    }
    Periodic(rho);

    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){

                //> laplacian of Rho  
                laplacian_rho.Node(i,j,k)  = 0.0;
                double del_t = 1.0;
                double Coeff = (2.0/(del_t*del_t*lb.theta0));

                for(int dv = 0; dv< grid.d_v; dv++)
                    laplacian_rho.Node(i,j,k) += lb.W[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;

                laplacian_rho.Node(i,j,k) = Coeff * ( laplacian_rho.Node(i,j,k) - rho.Node(i,j,k));

                //> gradient of Rho
                double grad_rhox = 0.0,grad_rhoy = 0.0,grad_rhoz = 0.0;
                del_t = 1.0;
                Coeff = (1.0/(del_t*lb.theta0));

                for(int dv = 0; dv< grid.d_v; dv++){
                    grad_rhox += lb.W[dv]*lb.Cx[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
                    grad_rhoy += lb.W[dv]*lb.Cy[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
                    grad_rhoz += lb.W[dv]*lb.Cz[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
                }
                double grad_rho = Coeff* (grad_rhox + grad_rhoy + grad_rhoz);


                //> pnid
                pnid.Node(i,j,k) = (rho.Node(i,j,k)*rho.Node(i,j,k)*b *lb.theta0 )/(1.0 - rho.Node(i,j,k)*b)  - 
                                a * rho.Node(i,j,k)*rho.Node(i,j,k) ;
                
                //> Fnid
                fnid.Node(i,j,k) =  - a * rho.Node(i,j,k)*rho.Node(i,j,k)   -   rho.Node(i,j,k)*lb.theta0 *log(1.0 - rho.Node(i,j,k)*b)
                                    -kappa * 0.5*grad_rho*grad_rho
                                    ;
                
                //> munid
                munid.Node(i,j,k) = -lb.theta0*log(1.0 - rho.Node(i,j,k)*b) ;
                munid.Node(i,j,k) += rho.Node(i,j,k)*b*lb.theta0/(1.0 - rho.Node(i,j,k)*b);
                munid.Node(i,j,k) -= 2.0*rho.Node(i,j,k)*a 
                                    -kappa*laplacian_rho.Node(i,j,k)
                                    ;
            }
        }    
    }
    Periodic(laplacian_rho);
    Periodic(pnid);
    Periodic(fnid);
    Periodic(munid); 

    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){

      
                laplacian_pnidplusfnidbyrho.Node(i,j,k)  = 0.0;
                double del_t = 1.0;
                double Coeff = (2.0/(del_t*del_t*lb.theta0));

                for(int dv = 0; dv< grid.d_v; dv++)
                    laplacian_pnidplusfnidbyrho.Node(i,j,k) += lb.W[dv]*((pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Node(i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])) / rho.Node(i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])) ;

                laplacian_pnidplusfnidbyrho.Node(i,j,k) = Coeff * ( laplacian_pnidplusfnidbyrho.Node(i,j,k) - ((pnid.Node(i,j,k) + fnid.Node(i,j,k)) / rho.Node(i,j,k)));
            }
      
        }
    }

    Periodic(laplacian_pnidplusfnidbyrho);
    




    //< Collision
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){
                
                double Rho = 0.0;

                //> CHEMICAL POTENTIAL FORMULATION 
                double grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;
                double Fx = 0, Fy = 0, Fz = 0.0;
                double del_t = 1.0;
                double Coeff_grad = (1.0/(del_t*lb.theta0));

                for(int dv = 0; dv< grid.d_v; dv++){
                    grad_mux += lb.W[dv]*lb.Cx[dv]*munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
                    grad_muy += lb.W[dv]*lb.Cy[dv]*munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
                    grad_muz += lb.W[dv]*lb.Cz[dv]*munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
                }


                Fx = - Coeff_grad*(grad_mux);
                Fy = - Coeff_grad*(grad_muy);  
                Fz = - Coeff_grad*(grad_muz);  


                //> MECHANICAL FORMULATION


                // // > direct              
                double Fx  = 0.0, Fy = 0.0, Fz = 0.0;
                double del_t = 1.0;
                double Coeff_grad = (1.0/(del_t*lb.theta0));


                for(int dv = 0; dv < grid.d_v; dv++){

                    Fx += Coeff_grad*(               
                                        lb.W[dv]*lb.Cx[dv]*( (pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Node( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
                                        )
                                        - kappa*Coeff_grad*( lb.W[dv]*lb.Cx[dv]*( laplacian_rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
                                        ;

                    Fy += Coeff_grad*(               
                                        lb.W[dv]*lb.Cy[dv]*( (pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Node( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
                                        )
                                        - kappa*Coeff_grad*( lb.W[dv]*lb.Cy[dv]*( laplacian_rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
                                        ;

                     Fz += Coeff_grad*(               
                                        lb.W[dv]*lb.Cz[dv]*( (pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Node( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
                                        )
                                        - kappa*Coeff_grad*( lb.W[dv]*lb.Cz[dv]*( laplacian_rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
                                        ;

                }
                
                
                // //> 4th order corrections
                double coeff_4th_grad = (lb.theta0*del_t*del_t)/2.0;

                for(int dv = 0; dv < grid.d_v; dv++){
                    Fx -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cx[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
                    Fy -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cy[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
                    Fz -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cz[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
                }

                Fx = -1.0*Fx;
                Fy = -1.0*Fy;
                Fz = -1.0*Fz;
                
                

                get_moments(grid, lb,  ux, uy, uz,Rho, i, j, k, Fx, Fy , Fz);            //for the node
                get_equi(feq_Node,lb, ux, uy,uz, Rho);

                //> normal
                for (int dv = 0; dv< 19; dv++){
                    grid.Node(i,j,k,dv) =  grid.Node(i,j,k,dv) + 2.0* beta*(feq_Node[dv] - grid.Node(i,j,k,dv))
                                        + 2.0 *beta * tau*lb.thetaInverse * rho.Node(i,j)* lb.W[dv] * (Fx * lb.Cx[dv] + Fy * lb.Cy[dv] + Fz * lb.Cz[dv])
                                        ;
                }       

                //> with entropic 
                //<this works properly
                // double x_i[9];
                // for(int dv = 0; dv< grid.d_v; dv++)
                //     x_i[dv] = feq_Node[dv]/grid.Node(i,j,dv) - 1.0;
                
                // double alpha = 0;
                // for(int dv = 0; dv< grid.d_v; dv++){
                //     alpha = 2.0;
                //     if( std::fabs(x_i[dv]) > 0.0001){
                //         calculateAlpha(lb,x_i,&grid.Node(i,j,0),beta,alpha);
                //         break;
                //     }
                // }

                // for (int dv = 0; dv< 9; dv++){
                //     grid.Node(i,j,dv) =  grid.Node(i,j,dv) + alpha* beta*(feq_Node[dv] - grid.Node(i,j,dv))
                //                         + (1 - 0.5*alpha*beta)*lb.thetaInverse * rho.Node(i,j)* lb.W[dv] * (Fx * lb.Cx[dv] + Fy * lb.Cy[dv] + Fz * lb.Cz[dv] );
                // }

                //>entopic 3rd type 
                // double feq_zeroforce[9] ={};

                // double x_i[9];
                // for(int dv = 0; dv< grid.d_v; dv++)
                //     x_i[dv] = feq_Node[dv]/grid.Node(i,j,dv) - 1.0;
                
                // double alpha = 0;
                // for(int dv = 0; dv< grid.d_v; dv++){
                //     alpha = 2.0;
                //     if( std::fabs(x_i[dv]) > 0.0001){
                //         calculateAlpha(lb,x_i,&grid.Node(i,j,0),beta,alpha);
                //         break;
                //     }
                // }

                // get_moments(grid, lb,  ux, uy,Rho, i, j, 0, 0 );            //for the node
                // get_equi(feq_zeroforce,lb, ux, uy, Rho);

                // for (int dv = 0; dv< 9; dv++){
                //     grid.Node(i,j,dv) =  grid.Node(i,j,dv) + alpha* beta*(feq_zeroforce[dv] - grid.Node(i,j,dv))
                //                         +feq_Node[dv] -  feq_zeroforce[dv] ;
                // }
            
            }

        }
    }
}
    


template<typename T>
void calculateAlpha(lbmD3Q19<T> &lbModel,T* x_i,T* f_i,T beta,T& alpha)
{
  double a(0.0), b(0.0), c(0.0),oneBySix(1.0/6.0);
  double ximin(0.0), ximax(0.0);

  alignas(32) T xSq [19];
  alignas(32) T fxSq[19];

  for(int dv = 0;dv<19;dv++)
  { 
    ximin = std::min(ximin, x_i[dv])  ;
    ximax = std::max(ximax, x_i[dv])  ;
  }

  for(int dv = 0;dv<19;dv++)
  {
    xSq[dv]  = x_i[dv]*x_i[dv] ;
    fxSq[dv] = f_i[dv]*x_i[dv]*x_i[dv] ;

    if(x_i[dv]<0.0)
      a += fxSq[dv]*x_i[dv]*0.5 ;

    b += fxSq[dv]*0.5 ;
    c += fxSq[dv]/(1.0 + 0.5*x_i[dv]) ;
  }

  T alphaMax = -1.0/(beta*ximin);

  T k;
  if(a<0 && b>0 && c>0)
    k = (b-sqrt(b*b - 4.0*a*c))/(2.0*a);
  else
    k = 1.5;

  a = 0.0;b=0.0;c=0.0;
//   T kBeta   = k*beta;
  T beta2   = beta*beta;
  T fourByK = 4.0/k;
  T hBeta = 0.0;

  for(int dv = 0; dv < 19; dv++)
  {
    if(x_i[dv]<0.0)
    {
      a += fxSq[dv]*x_i[dv]*beta2*( 1.0/6.0 - hBeta*x_i[dv]/12.0 + hBeta*hBeta*x_i[dv]*x_i[dv]/20 - hBeta*hBeta*hBeta*x_i[dv]*x_i[dv]*x_i[dv]/5.0 );
      b += fxSq[dv]*0.5;
    }

    if(x_i[dv]>0.0)
      b += f_i[dv]*( (x_i[dv]*x_i[dv]*0.5) - beta2*(x_i[dv]*x_i[dv]*x_i[dv]/15.0)* ( (4.0/(fourByK+x_i[dv])) + (2.0/(fourByK+2.0*x_i[dv])) + (4.0/(fourByK+3.0*x_i[dv])) ));

    c += f_i[dv]*(60.0*x_i[dv]*x_i[dv] + 60.0*x_i[dv]*x_i[dv]*x_i[dv] + 11.0*x_i[dv]*x_i[dv]*x_i[dv]*x_i[dv])/( 60.0 + 90.0*x_i[dv] + 36.0*x_i[dv]*x_i[dv] + 3.0*x_i[dv]*x_i[dv]*x_i[dv]);
  }

  T  h;
  if(a<0 && b>0 && c>0)
    h = (b-std::sqrt(b*b - 4.0*a*c))/(2.0*a);
  else
    h = 2.1;

  a = 0.0;
  hBeta = h*beta;

  for(int dv = 0; dv < 19; dv++)
  {
    if(x_i[dv]<0.0)
    {
      a += f_i[dv]*x_i[dv]*x_i[dv]*x_i[dv]*beta2*( 1.0/6.0 - hBeta*x_i[dv]/12.0 + hBeta*hBeta*x_i[dv]*x_i[dv]/20 - hBeta*hBeta*hBeta*x_i[dv]*x_i[dv]*x_i[dv]/5.0 );
    }
  }

  if(a<0 && b>0 && c>0)
    alpha = 2.0*c/(b+sqrt(b*b - 4.0*a*c));

  if(alpha > alphaMax)
  {
    if ( alphaMax > 1.0)
      alpha = 0.5*(1.0+alphaMax) ;
    else
      alpha = 0.95*alphaMax;
  }
}



template<typename T>
void get_equi(double *feq , lbmD3Q19<T> &lb, double ux, double uy, double uz, double rho){


    double u2 = ux*ux + uy*uy + uz*uz;
    double a1=0;
    double first,second, third,feq0=0;
    for (int dv = 0; dv< 19; dv++){

        feq0 = rho*lb.W[dv];

        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv] + uz*lb.Cz[dv])*lb.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb.thetaInverse;
        feq[dv] = feq0*(1+ first + second + third);    

        
    }
         
}





template<typename T,typename T1>
void get_moments(Grid_N_C_3D<T> &grid, lbmD3Q19<T1> &lb,double &Ux, double &Uy, double &Uz,double &Rho,  int X, int Y, int Z, double Fx = 0, double Fy = 0, double Fz = 0){ ///node or cell 0-Node 1- cell
   Ux  = 0.0;
   Uy  = 0.0;
   Uz  = 0.0;
   Rho = 0.0;


    for(int dv = 0; dv <grid.d_v; dv++){
        Ux  += grid.Node(X,Y,Z,dv)*lb.Cx[dv];
        Uy  += grid.Node(X,Y,Z,dv)*lb.Cy[dv];
        Uz  += grid.Node(X,Y,Z,dv)*lb.Cz[dv];
        Rho += grid.Node(X,Y,Z,dv);
    }  

    Ux = Ux/Rho + 0.5*Fx;
    Uy = Uy/Rho + 0.5*Fy;
    Uz = Uz/Rho + 0.5*Fz;  

}

//period no of waves
template<typename T, typename T1>
void initialization(Grid_N_C_3D<T> &grid,lbmD3Q19<T1> &lb,double Rho_mean ,double amplitude, double period){


    double Feq_node[19] = {0},Rho = 0.0;
    double x,y,z
           ;    ///distance between nodes 
    
    double ux_node =0, uy_node = 0, uz_node = 0;
    
    double k_w = 2*M_PI* period;
    
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){

                x = ((double)i)/ grid.n_x;
                y = ((double)j)/ grid.n_y;
                z = ((double)k)/ grid.n_z;

                Rho = Rho_mean + amplitude * sin(k_w* x)* sin(k_w*y)*sin(k_w*z);

                get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho);

                for (int dv = 0; dv<grid.d_v; dv++){
                    grid.Node(i,j,k,dv) = Feq_node[dv]
                   
                    ;
                }
            }
        }
    }
    
    
        
}
        



template<typename T>
void printMass(Grid_N_C_3D<T> &grid){    
    double a = 0;
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){

                for (int dv = 0; dv< grid.d_v; dv++){
                    a += grid.Node(i,j,k,dv) ;
                }
            }
        }
    }
std::cout<<"Total mass =  "<<a<<std::endl;

}


;