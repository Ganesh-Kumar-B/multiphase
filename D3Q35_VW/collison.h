#pragma once

#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<fstream>
#include<algorithm>
#include <sstream>
#include<string>
#include "lbmD3Q35.h"
#include "GRID_3D.h"
#include "multiphase.h"
#define PI 3.14159265



template<typename T, typename T1>
void collide(Grid_N_C_3D<T> &grid,
            lbmD3Q35<T1> &lb,double beta,double tau, double TbyTc, double kappa, int t ){

    Grid_N_C_3D<T>  laplacian_pnidplusfnidbyrho     (grid.n_x,grid.n_y,grid.n_z,2,1);
    Grid_N_C_3D<T>  rho                             (grid.n_x,grid.n_y,grid.n_z,2,1);   
    Grid_N_C_3D<T>  pnid                            (grid.n_x,grid.n_y,grid.n_z,2,1);   
    Grid_N_C_3D<T>  fnid                            (grid.n_x,grid.n_y,grid.n_z,2,1);       
    Grid_N_C_3D<T>  munid                           (grid.n_x,grid.n_y,grid.n_z,2,1);   
    Grid_N_C_3D<T>  laplacian_rho                   (grid.n_x,grid.n_y,grid.n_z,2,1);   
    Grid_N_C_3D<T>  laplacian_fnid                  (grid.n_x,grid.n_y,grid.n_z,2,1);   
    Grid_N_C_3D<T>  gradient_rho                    (grid.n_x,grid.n_y,grid.n_z,2,3);   //3 components

    double feq_Node[35] = {0},feq_Cell[35] = {0}, ux = 0, uy = 0, uz = 0;



    Multiphase_terms(grid,rho,pnid, fnid, munid,laplacian_rho,laplacian_fnid,gradient_rho,lb,TbyTc,kappa );



    //< Collision
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){
                
                double Rho = 0.0;
                double Fx = 0, Fy = 0, Fz = 0;

                // Multiphase_Force_Node(grid,rho,pnid, fnid, munid,laplacian_rho,lb,Fx,Fy,Fz,i,j,k );
                Multiphase_Force_eta_Node(grid,rho,pnid, fnid, munid,laplacian_rho,gradient_rho,lb,Fx,Fy,Fz,i,j,k,kappa );
                


                get_moments_Node(grid, lb,  ux, uy, uz,Rho, i, j, k, Fx, Fy , Fz);            //for the node
                get_equi(feq_Node,lb, ux, uy,uz, Rho);


                // //> normal
                // for (int dv = 0; dv< grid.d_v; dv++){
                //     grid.Node(i,j,k,dv) =  grid.Node(i,j,k,dv) + 2.0* beta*(feq_Node[dv] - grid.Node(i,j,k,dv))
                //                         + 2.0 *beta * tau*lb.thetaInverse * rho.Node(i,j,k)* lb.W[dv] * (Fx * lb.Cx[dv] + Fy * lb.Cy[dv] + Fz * lb.Cz[dv])
                //                         ;
                // }
                //> with entropic 
                //<this works properly
                //$for node
                double x_i[35];
                for(int dv = 0; dv< grid.d_v; dv++)
                    x_i[dv] = feq_Node[dv]/grid.Node(i,j,dv) - 1.0;
                
                double alpha = 0;
                for(int dv = 0; dv< grid.d_v; dv++){
                    alpha = 2.0;
                    if( std::fabs(x_i[dv]) > 0.0001){
                        calculateAlpha(lb,x_i,&grid.Node(i,j,0),beta,alpha);
                        break;
                    }
                }

                for (int dv = 0; dv< 35; dv++){
                    grid.Node(i,j,k,dv) =  grid.Node(i,j,k,dv) + alpha* beta*(feq_Node[dv] - grid.Node(i,j,k,dv))
                                        + (1 - 0.5*alpha*beta)*lb.thetaInverse * rho.Node(i,j,k)* lb.W[dv] * (Fx * lb.Cx[dv] + Fy * lb.Cy[dv] + Fz * lb.Cz[dv] );
                }       




                //< CELLS    
                Rho = 0.0;

                
                // Multiphase_Force_Cell(grid,rho,pnid, fnid, munid,laplacian_rho,lb,Fx,Fy,Fz,i,j,k );
                Multiphase_Force_eta_Cell(grid,rho,pnid, fnid, munid,laplacian_rho,gradient_rho,lb,Fx,Fy,Fz,i,j,k,kappa );

                get_moments_Cell(grid, lb,  ux, uy, uz,Rho, i, j, k, Fx, Fy , Fz);            //for the node
                get_equi(feq_Cell,lb, ux, uy,uz, Rho);


                //> normal
                // for (int dv = 0; dv< grid.d_v; dv++){
                //     grid.Cell(i,j,k,dv) =  grid.Cell(i,j,k,dv) + 2.0* beta*(feq_Cell[dv] - grid.Cell(i,j,k,dv))
                //                         + 2.0 *beta * tau*lb.thetaInverse * rho.Cell(i,j,k)* lb.W[dv] * (Fx * lb.Cx[dv] + Fy * lb.Cy[dv] + Fz * lb.Cz[dv])
                //                         ;
                // }       



                //$ for Cell
                
                for(int dv = 0; dv< grid.d_v; dv++)
                    x_i[dv] = feq_Cell[dv]/grid.Cell(i,j,dv) - 1.0;
                
                alpha = 0;
                for(int dv = 0; dv< grid.d_v; dv++){
                    alpha = 2.0;
                    if( std::fabs(x_i[dv]) > 0.0001){
                        calculateAlpha(lb,x_i,&grid.Cell(i,j,0),beta,alpha);
                        break;
                    }
                }

                for (int dv = 0; dv< 35; dv++){
                    grid.Cell(i,j,k,dv) =  grid.Cell(i,j,k,dv) + alpha* beta*(feq_Cell[dv] - grid.Cell(i,j,k,dv))
                                        + (1 - 0.5*alpha*beta)*lb.thetaInverse * rho.Cell(i,j,k)* lb.W[dv] * (Fx * lb.Cx[dv] + Fy * lb.Cy[dv] + Fz * lb.Cz[dv] );
                }



            }
        }
    }
}
    


template<typename T>
void calculateAlpha(lbmD3Q35<T> &lbModel,T* x_i,T* f_i,T beta,T& alpha)
{
  double a(0.0), b(0.0), c(0.0),oneBySix(1.0/6.0);
  double ximin(0.0), ximax(0.0);

  alignas(32) T xSq [35];
  alignas(32) T fxSq[35];

  for(int dv = 0;dv<35;dv++)
  { 
    ximin = std::min(ximin, x_i[dv])  ;
    ximax = std::max(ximax, x_i[dv])  ;
  }

  for(int dv = 0;dv<35;dv++)
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

  for(int dv = 0; dv < 35; dv++)
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

  for(int dv = 0; dv < 35; dv++)
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
void get_equi(double *feq , lbmD3Q35<T> &lb, double ux, double uy, double uz, double rho){


    double u2 = ux*ux + uy*uy + uz*uz;
    double a1=0;
    double first,second, third,feq0=0;
    for (int dv = 0; dv< 35; dv++){

        feq0 = rho*lb.W[dv];

        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv] + uz*lb.Cz[dv])*lb.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb.thetaInverse;
        feq[dv] = feq0*(1+ first + second + third);    

        
    }
         
}





template<typename T,typename T1>
void get_moments_Node(Grid_N_C_3D<T> &grid, lbmD3Q35<T1> &lb,double &Ux, double &Uy, double &Uz,double &Rho,  int X, int Y, int Z, double Fx = 0, double Fy = 0, double Fz = 0){ ///node or cell 0-Node 1- cell
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

template<typename T,typename T1>
void get_moments_Cell(Grid_N_C_3D<T> &grid, lbmD3Q35<T1> &lb,double &Ux, double &Uy, double &Uz,double &Rho,  int X, int Y, int Z, double Fx = 0, double Fy = 0, double Fz = 0){ ///node or cell 0-Cell 1- cell
    Ux  = 0.0;
    Uy  = 0.0;
    Uz  = 0.0;
    Rho = 0.0;



    for(int dv = 0; dv <grid.d_v; dv++){
        Ux  += grid.Cell(X,Y,Z,dv)*lb.Cx[dv];
        Uy  += grid.Cell(X,Y,Z,dv)*lb.Cy[dv];
        Uz  += grid.Cell(X,Y,Z,dv)*lb.Cz[dv];
        Rho += grid.Cell(X,Y,Z,dv);
    }  

    Ux = Ux/Rho + 0.5*Fx;
    Uy = Uy/Rho + 0.5*Fy;
    Uz = Uz/Rho + 0.5*Fz;  

}



//period no of waves
template<typename T, typename T1>
void initialization(Grid_N_C_3D<T> &grid,lbmD3Q35<T1> &lb,double Rho_mean ,double amplitude, double period){


    double Feq_node[35] = {0},Feq_cell[35] = {0},Rho = 0.0;
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
                   
                    grid.Node(i,j,k,dv) = Feq_node[dv];
                    
                }
                
                x = (((double)i) + 0.5)/ grid.n_x;
                y = (((double)j) + 0.5)/ grid.n_y;
                z = (((double)k) + 0.5)/ grid.n_z;

                Rho = Rho_mean + amplitude * sin(k_w* x)* sin(k_w*y)*sin(k_w*z);

                get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho);

                for (int dv = 0; dv<grid.d_v; dv++){
                    
                    grid.Cell(i,j,k,dv) = Feq_node[dv];

                    

                }


                //>shock initialization
                // if(i < grid.n_x/2){
                //     Rho = 1.2;

                //     get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho);

                //     for (int dv = 0; dv<grid.d_v; dv++){
                        
                //         grid.Node(i,j,k,dv) = Feq_node[dv];
                //         grid.Cell(i,j,k,dv) = Feq_node[dv];
                        
                        
                //     }
                    
                // }
                // else{

                //     Rho = 0.8;

                //     get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho);

                //     for (int dv = 0; dv<grid.d_v; dv++){
                        
                //         grid.Node(i,j,k,dv) = Feq_node[dv];
                //         grid.Cell(i,j,k,dv) = Feq_node[dv];

                //     }

                // }

                //>single initialization
                // grid.Node(10,10,10,26) = 1.0;


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
                    a += grid.Cell(i,j,k,dv) ;

                }
            }
        }
    }
std::cout<<"Total mass =  "<<a<<std::endl;

}


;




               
    




