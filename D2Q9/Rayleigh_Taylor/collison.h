#pragma once

#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<fstream>
#include<algorithm>
#include <sstream>
#include<string>
#include "lbmD2Q9.h"
#include "GRID_2D.h"
#include "multiphase.h"
#define PI 3.14159265



template<typename T, typename T1>
void collide(Grid_N_C_2D<T> &grid,
            lbmD2Q9<T1> &lb9, real beta,real tau, real TbyTc, real kappa, int t,Grid_N_C_2D<T> &Force, real g ){

    Grid_N_C_2D<T>  laplacian_pnidplusfnidbyrho     (grid.n_x,grid.n_y,1,1);
    Grid_N_C_2D<T>  rho                             (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  pnid                            (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  fnid                            (grid.n_x,grid.n_y,1,1);       
    Grid_N_C_2D<T>  munid                           (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  laplacian_rho                   (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  laplacian_fnid                  (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  gradient_rho                    (grid.n_x,grid.n_y,1,2);   //2 components

    real feq_Node[9] = {0},



    ux = 0, uy = 0;

    real eta =  0;   //   0 ----> fourth order   1-----> second order 



    real rho_critical = 1.0, T_critical = lb9.theta0/TbyTc ; 
    real b = 0.521772/(rho_critical), a = b*T_critical/0.377332;
    kappa = kappa*a;

    Multiphase_terms(grid,Force,rho,pnid, fnid, munid,laplacian_rho,laplacian_fnid,gradient_rho,lb9,TbyTc,kappa, a, b );

    //       //  first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            
        
            real Rho = 0.0;

            Multiphase_Force_Node(grid,rho,pnid, fnid, munid,laplacian_rho,lb9,Force,i,j, kappa, a, b, g );   

        

            
            get_moments_Node(grid, lb9,  ux, uy,Rho, i, j, Force );            //for the node
            get_equi(feq_Node ,lb9, ux, uy, Rho);

        
            // // //> normal
            for (int dv = 0; dv< grid.d_v; dv++){
                grid.Node(i,j,dv) =  grid.Node(i,j,dv) + 2.0* beta*(feq_Node[dv] - grid.Node(i,j,dv))
                                    + 2.0 *beta * tau*lb9.thetaInverse *  feq_Node[dv] * (Force.Node(i,j,0) * (lb9.Cx[dv] - ux) + Force.Node(i,j,1) * (lb9.Cy[dv] - uy) );
                                    ;
            }


        }
    }


}
    







template<typename T, typename T1>
void initialization(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,lbmD2Q9<T1> &lb9, real phi_l ,real phi_h, real rho_l , real rho_h, real a, real b){

    real Feq_node[9] = {0},Geq_node[9] = {0},rho = 1.0;
    real x,y;

    real phi = 0;
    
    real ux = 0,  uy =  0;
    
    real u1 = 0.0, u2 = 0.0;

    real p = 0;
    
    for(int i = gridf.nbx; i <= gridf.nex ; i++){
        for(int j = gridf.nby;j <= gridf.ney ; j++){
            
            x = ((real)i)/ gridf.n_x ;
            y = ((real)j)/ gridf.n_x ;


            phi = phi_h;

            if(y>1){
                phi = phi_l;
            }
            
            

            get_equi_f(Feq_node   ,lb9, ux, uy, phi);   //need phi, u 

            for (int dv = 0; dv< gridf.d_v; dv++){
                gridf.Node(i,j,dv) = Feq_node[dv];
            }

            rho = rho_l + ((phi - phi_l)/(phi_h - phi_l)) *(rho_h - rho_l);

            real eta = rho*b/4.0;

            p = rho *lb9.theta0*(1 + eta + eta*eta - eta*eta*eta)/pow(1 - eta, 3)  - a *rho*rho ;

            get_equi_g (Geq_node   ,lb9, ux, uy, rho,p);                               // need p and rho and u

            for(int dv = 0; dv< gridg.d_v; dv++){
                gridg.Node(i,j,dv) = Geq_node[dv];
            }




        }
    }
}



template<typename T, typename T1>
void initialization_equilibrium_profile_y(Grid_N_C_2D<T> &grid,lbmD2Q9<T1> &lb,real rho_liq,real rho_gas ){

	real Feq_node[9] = {0},Rho = 0.0;
    real x,y
           ;    ///distance between nodes 
    
    real  x_0 = 0.0;
    real  y_0 = 0.0;

    real phi_l = +1.0;
    real phi_h = -1.0;

    real phi;
    real ux_node = 0.0, uy_node = 0.0;
    
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){

            x = ((real)i)/ grid.n_x - x_0;
            y = ((real)j)/ grid.n_x - y_0;
            z = ((real)k)/ grid.n_x - z_0;

            phi = (tanh((y - 1 - 0.005*cos(2.0*M_PI*x))/(sqrt(2) * (1.0/grid.n_x))));

            Rho = rho_gas + (phi - phi_l)/(phi_h - phi_l) *(rho_liq - rho_gas);
            

            get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho);

            for (int dv = 0; dv<grid.d_v; dv++)
                grid.Node(i,j,k,dv) = Feq_node[dv];


        }
    }
    
}



template<typename T>
void get_equi(real *feq , lbmD2Q9<T> &lb, real ux, real uy, real rho){


    real u2 = ux*ux + uy*uy ;
    real a1 = 0;
    real first,second, third,feq0=0;
    for (int dv = 0; dv< 9; dv++){

        feq0 = rho*lb.W[dv];

        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv] )*lb.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb.thetaInverse;
        feq[dv] = feq0*(1+ first + second + third);    

        
    }
         
}





template<typename T,typename T1>
void get_moments_Node(Grid_N_C_2D<T> &grid, lbmD2Q9<T1> &lb,real &Ux, real &Uy, real &Rho,  int X, int Y, Grid_N_C_2D<T> &Force){ ///node or cell 0-Node 1- cell
    Ux  = 0.0;
    Uy  = 0.0;
    Rho = 0.0;


    for(int dv = 0; dv <grid.d_v; dv++){
        Ux  += grid.Node(X,Y,dv)*lb.Cx[dv];
        Uy  += grid.Node(X,Y,dv)*lb.Cy[dv];
        Rho += grid.Node(X,Y,dv);
    }  

    Ux = Ux/Rho + 0.5*Force.Node(X,Y,0);
    Uy = Uy/Rho + 0.5*Force.Node(X,Y,1);


}











template<typename T>
void printMass(Grid_N_C_2D<T> &gridf){    
    real a = 0;
    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            for (int dv = 0; dv< 9; dv++){
                a += gridf.Node(i,j,dv) ;
            }
        }
    }
    std::cout<<"   "<<a<<std::endl;

}
