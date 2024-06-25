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
#include "advection.h"
#include "GRID_2D.h"
#include "multiphase.h"
#define PI 3.14159265


template<typename T, typename T1>
void collide(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,Grid_N_C_2D<T> &rho, Grid_N_C_2D<T> &phi,Grid_N_C_2D<T> &mu,Grid_N_C_2D<T> &laplacian_phi,
            lbmD2Q9<T1> &lb,real beta,real tau,real tauphi, real TbyTc, int t,Grid_N_C_2D<T> &Force, real kappa,real gamma_s, real A ){

        
    real    feq_Node[9] = {0}, 
            geq_Node[9] = {0},
            ux = 0, uy = 0;

    real eta =0;   //   0 ----> fourth order   1-----> second order 
    
    Multiphase_terms(gridf, gridg,lb,rho,phi , mu,laplacian_phi,kappa, A );

    //       //  first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            
            real Rho = 0.0;

            Multiphase_Force_Node(gridf,phi , mu,lb,Force,i,j );             										
            
            get_moments_Node_f(gridf, lb,  ux, uy, Rho, i, j, Force);            //for the node
            get_equi_f(feq_Node,lb, ux, uy, Rho, phi.Node(i,j),mu.Node(i,j));

            // // //> normal
            for (int dv = 0; dv< gridf.d_v; dv++){
                gridf.Node(i,j,dv) =  gridf.Node(i,j,dv) + (1.0/tau )*(feq_Node[dv] - gridf.Node(i,j,dv))
                                    +(1.0 - (1.0 / (2.0* tau)  ) )*lb.thetaInverse * rho.Node(i,j)* lb.W[dv] * (Force.Node(i,j,0) * lb.Cx[dv] + Force.Node(i,j,1) * lb.Cy[dv] + Force.Node(i,j,1) * lb.Cz[dv])
                                    ;
            }


        }
    }



    for(int i = 0 + gridg.noghost; i < gridg.n_x_node - (gridg.noghost) ; i++){
        for(int j = 0 + gridg.noghost;j < gridg.n_y_node - (gridg.noghost) ; j++){    
            
            real Rho = 0.0;
            
            get_moments_Node_f(gridf, lb,  ux, uy,Rho, i, j,Force);

            get_equi_g(feq_Node,lb, ux, uy, phi.Node(i,j),gamma_s,mu.Node(i,j));

            // //> normal
            for (int dv = 0; dv< gridg.d_v; dv++){
                gridg.Node(i,j,dv) =  gridg.Node(i,j,dv) + (1.0/tauphi)*(feq_Node[dv] - gridg.Node(i,j, dv))
                                    ;
            }

        }
    }


}
    







template<typename T, typename T1>
void initialization(Grid_N_C_2D<T> &gridf,lbmD2Q9<T1> &lb9, real Rho){

    real Feq_node[9] = {0};

    real phi = 0;
    
    real ux = 0,  uy =  0;
    

    
    for(int i = gridf.nbx; i <= gridf.nex ; i++){
        for(int j = gridf.nby;j <= gridf.ney ; j++){
            

            get_equi(Feq_node,lb9,ux,uy,Rho);

            for (int dv = 0; dv<gridf.d_v; dv++)
                gridf.Node(i,j,dv) = Feq_node[dv];

        }
    }
}





template<typename T>
void get_equi_f(real *feq , lbmD2Q9<T> &lb, real ux, real uy, real &rho, real &phi, real &mu){

    real u2 = ux*ux + uy*uy ;
    real a1=0;
    real first,second, third, chem_pot , feq0=0;
    real sum = 0;
    for (int dv = 1; dv< 9; dv++){

        feq0 = rho*lb.W[dv];

        chem_pot = (phi* mu) / (rho* lb.theta0) ;

        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv] )*lb.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb.thetaInverse;
        feq[dv] = feq0*(1+ chem_pot +first + second + third);    

        sum += feq[dv];
    }

    feq[0] = rho - sum;
         
}


template<typename T>
void get_equi_g(real *feq , lbmD2Q9<T> &lb, real ux, real uy,  real &phi,  real &gamma, real &mu){


    real u2 = ux*ux + uy*uy ;
    real first,second, third, chem_pot ,feq0=0;
    real sum = 0;

    for (int dv = 1; dv< 9; dv++){

        feq0 = lb.W[dv];

        chem_pot = (gamma* mu) / (lb.theta0) ;

        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv] )*lb.thetaInverse;
        second = 0.5*(first * first) ;
        third = - 0.5*u2*lb.thetaInverse;
        feq[dv] = feq0*( chem_pot + phi*  first +phi * second + phi * third);    

        sum += feq[dv];
    }
    
    
    feq[0] = phi - sum;

}






template<typename T,typename T1>
void get_moments_Node_f(Grid_N_C_2D<T> &grid, lbmD2Q9<T1> &lb,real &Ux, real &Uy, real &Rho,  int X, int Y, Grid_N_C_2D<T> &Force){ ///node or cell 0-Node 1- cell
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




template<typename T,typename T1>
void get_moments_Node_g(Grid_N_C_2D<T> &grid, lbmD2Q9<T1> &lb,real &phi,  int X, int Y){ ///node or cell 0-Node 1- cell
    
    phi = 0.0;

    for(int dv = 0; dv <grid.d_v; dv++){
        
        phi += grid.Node(X,Y,dv);
    }  

}


template<typename T>
void printMass(Grid_N_C_2D<T> &grid){    
    real a = 0;
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){

            for (int dv = 0; dv< grid.d_v; dv++){
                a += grid.Node(i,j,dv) ;

            }
            
        }
    }
std::cout<<"Total mass =  "<<a<<std::endl;

}






template<typename T, typename T1>
void initialization_2D_droplet(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,Grid_N_C_2D<T> &rho,
Grid_N_C_2D<T> &phi,Grid_N_C_2D<T> &mu,Grid_N_C_2D<T> &laplacian_phi,lbmD2Q9<T1> &lb,real Rho_mean, real kappa, real gamma_s, real A ){

	real Feq_node[9] = {0};
    real x,y;     ///distance between nodes 
    

    real  x_0 = 0.5;
    real  y_0 = 0.5;

    real ux_node = 0.0, uy_node = 0.0;
    

    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){
        for(int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){

            x = ((real)i)/ gridf.n_x - x_0;
            y = ((real)j)/ gridf.n_y - y_0;
            
            phi.Node(i,j) = -1.0;

            if( x * x  + y * y < 0.25*0.25 ){

                phi.Node(i,j) = 1.0;

            }

        

            
        }
    }

    Periodic_left_Right(phi); 
    Periodic_top_bottom(phi); 


    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){
        for(int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){


            //> laplacian of phi  
            
            real del_t = 1.0;
            real Coeff = (2.0/(del_t*del_t*lb.theta0));
            
            laplacian_phi.Node(i,j)  = 0.0;

            for(int dv = 0; dv< 9; dv++)
                laplacian_phi.Node(i,j) += lb.W[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] ) ;

            

            laplacian_phi.Node(i,j) = Coeff * ( laplacian_phi.Node(i,j) - phi.Node(i,j));

            
            
            //>-------------------------------------------<\\


            //$------------------Node

            real    mu = - A * phi.Node(i,j) +  A * phi.Node(i,j) *  phi.Node(i,j) * phi.Node(i,j)   ;
                    mu -= kappa*laplacian_phi.Node(i,j);


            get_equi_f(Feq_node,lb,ux_node,uy_node,Rho_mean, phi.Node(i,j),mu);

            for (int dv = 0; dv<gridf.d_v; dv++)
                gridf.Node(i,j,dv) = Feq_node[dv];

            get_equi_g(Feq_node,lb,ux_node,uy_node ,phi.Node(i,j), gamma_s, mu);

            for (int dv = 0; dv<gridf.d_v; dv++)
                gridg.Node(i,j,dv) = Feq_node[dv];
            

                
        }
    }
}




