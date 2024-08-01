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
            lbmD2Q9<T1> &lb9, real beta,real tau, real TbyTc,real a,real b ,real kappa, real &sigma, int t,Grid_N_C_2D<T> &Force, Grid_N_C_2D<T> &P_tensor ,real g, real dx , real dt ){

    Grid_N_C_2D<T>  laplacian_pnidplusfnidbyrho     (grid.n_x,grid.n_y,1,1);
    Grid_N_C_2D<T>  rho                             (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  pnid                            (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  fnid                            (grid.n_x,grid.n_y,1,1);       
    Grid_N_C_2D<T>  munid                           (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  laplacian_rho                   (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  laplacian_munid                 (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  laplacian_fnid                  (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  gradient_rho                    (grid.n_x,grid.n_y,1,2);   //2 components



    real feq_Node[9] = {0},

    ux = 0, uy = 0;

    real eta =  0;   //   0 ----> fourth order   1-----> second order 



    sigma = 0;

    Multiphase_terms(grid,Force,rho,pnid, fnid, munid,laplacian_rho,laplacian_munid,laplacian_fnid,gradient_rho,P_tensor,lb9,TbyTc,kappa,sigma, a, b, dx, dt );

    //       //  first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            
        
            real Rho = 0.0;


            Multiphase_Force_Node(grid,rho,pnid, fnid, munid,laplacian_rho,lb9,Force,i,j, kappa, a, b, g,dx,dt );   //this gives force density
            // Multiphase_Force_Node_4th_order(grid,rho,pnid, fnid, munid,laplacian_rho,laplacian_munid,lb9,Force,i,j, kappa, a, b, g,dx,dt );   //this gives force density

            
            // Multiphase_Force_P(grid,rho,pnid, fnid, munid,laplacian_rho,P_tensor,lb9,Force,i,j, kappa, a, b, g,dx,dt );   //this gives force density



            get_moments_Node(grid, lb9,  ux, uy,Rho, i, j, Force,dx,dt );            //for the node
            get_equi(feq_Node ,lb9, ux, uy, Rho);


            // // //> normal
            for (int dv = 0; dv< grid.d_v; dv++){
                grid.Node(i,j,dv) =  grid.Node(i,j,dv) + 2.0* beta*(feq_Node[dv] - grid.Node(i,j,dv))
                                    + (1.0 - beta)*dt*lb9.thetaInverse *  feq_Node[dv] * (Force.Node(i,j,0) * (lb9.Cx[dv] - ux) + Force.Node(i,j,1) * (lb9.Cy[dv] - uy) );
                                    ;
            }



    


            // // // //> normal //from the he paper
            // for (int dv = 0; dv< grid.d_v; dv++){
            //     grid.Node(i,j,dv) =  grid.Node(i,j,dv) + (1.0/tau)*(feq_Node[dv] - grid.Node(i,j,dv))
            //                         - (1.0 - 0.5/tau)*lb9.thetaInverse *  feq_Node[dv] * (Force.Node(i,j,0) * (lb9.Cx[dv] - ux) + Force.Node(i,j,1) * (lb9.Cy[dv] - uy) );
            //                         ;
            // }


        }
    }


}
    









template<typename T, typename T1>
void initialization(Grid_N_C_2D<T> &gridf,lbmD2Q9<T1> &lb9, real Rho, real rho_liq, real rho_gas){

    real Feq_node[9] = {0};

    real phi = 0;
    
    real ux = 0,  uy =  0;
    
    real x,y;

    real x_0 = 0.0;
    real y_0 = 0.0;

    
    for(int i = gridf.nbx; i <= gridf.nex ; i++){
        for(int j = gridf.nby;j <= gridf.ney ; j++){
            
            x = ((real)i)/ gridf.n_x - x_0;
            y = ((real)j)/ gridf.n_x - y_0;


            Rho = (rho_liq + rho_gas)* 0.5 + (rho_liq - rho_gas) *0.5* tanh(x);

            get_equi(Feq_node,lb9,ux,uy,Rho);

            for (int dv = 0; dv<gridf.d_v; dv++)
                gridf.Node(i,j,dv) = Feq_node[dv];

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



            phi = (tanh((y - 2.0 - 0.00*cos(2.0*M_PI*x))/(sqrt(2) * (1.0/grid.n_x))));
            


            Rho = rho_gas + (phi - phi_l)/(phi_h - phi_l) *(rho_liq - rho_gas);
            


            get_equi(Feq_node,lb,ux_node,uy_node,Rho);

            for (int dv = 0; dv<grid.d_v; dv++)
                grid.Node(i,j,dv) = Feq_node[dv];


        }
    }
}





template<typename T, typename T1>
void initialization_ellipse(Grid_N_C_2D<T> &grid,lbmD2Q9<T1> &lb,real rho_liq,real rho_gas ){

	real Feq_node[9] = {0},Rho = 0.0;
    real x,y
           ;    ///distance between nodes 
    
    real  x_0 = 0.0;
    real  y_0 = 0.0;

    real phi_l = -1.0;
    real phi_h = +1.0;

    real phi;
    real ux_node = 0.0, uy_node = 0.0;
    
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){

            
            x = ((real)i)/ grid.n_x - x_0;
            y = ((real)j)/ grid.n_x - y_0;

            if(y > 2.5){
                phi = tanh( (0.15 - sqrt( (x -0.5)*(x - 0.5 ) + 1.0*(y - 3.0)*(y - 3.0) )  )/
                            (sqrt(2.0) * (1.2/ grid.n_x) )  
                            );

                Rho = rho_gas + (phi - phi_l)/(phi_h - phi_l) *(rho_liq - rho_gas);

                

                get_equi(Feq_node,lb,ux_node,uy_node,Rho);

                for (int dv = 0; dv<grid.d_v; dv++)
                    grid.Node(i,j,dv) = Feq_node[dv];
            }


        }
    }
}



template<typename T, typename T1>
void initialization_circle(Grid_N_C_2D<T> &grid,lbmD2Q9<T1> &lb,real rho_liq,real rho_gas, real R ){

	real Feq_node[9] = {0},Rho = 0.0;
    real x,y
           ;    ///distance between nodes 
    
    real  x_0 = 0.0;
    real  y_0 = 0.0;

    real phi_l = -1.0;
    real phi_h = +1.0;

    real phi;
    real ux_node = 0.0, uy_node = 0.0;
    
    // for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
    //     for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){

    //         x = ((real)i)/ grid.n_x - x_0;
    //         y = ((real)j)/ grid.n_x - y_0;
            
    //         phi = tanh( (R*R - sqrt( (x -0.5)*(x - 0.5 ) + 0.5*(y - 0.5)*(y - 0.5) )  ) / 
    //                     (sqrt(2.0) * (1.0/ grid.n_x) )  
    //                     );

    //         Rho = rho_gas + (phi - phi_l)/(phi_h - phi_l) *(rho_liq - rho_gas);



    //         get_equi(Feq_node,lb,ux_node,uy_node,Rho);

    //         for (int dv = 0; dv<grid.d_v; dv++)
    //             grid.Node(i,j,dv) = Feq_node[dv];


    //     }
    // }



    x_0 = 0.5, y_0 = 0.5;
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            
            
            
            x = ((real)i)/ grid.n_x - x_0;
            y = ((real)j)/ grid.n_x - y_0;
            

            Rho = rho_gas;
            if(x*x + y*y < R*R){
                Rho = rho_liq;
            }
            

            get_equi(Feq_node,lb,ux_node,uy_node,Rho);

            for (int dv = 0; dv<grid.d_v; dv++)
                grid.Node(i,j,dv) = Feq_node[dv];


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
void get_moments_Node(Grid_N_C_2D<T> &grid, lbmD2Q9<T1> &lb,real &Ux, real &Uy, real &Rho,  int X, int Y, Grid_N_C_2D<T> &Force,real dx, real dt){ ///node or cell 0-Node 1- cell
    Ux  = 0.0;
    Uy  = 0.0;
    Rho = 0.0;


    for(int dv = 0; dv <grid.d_v; dv++){
        Ux  += grid.Node(X,Y,dv)*lb.Cx[dv];
        Uy  += grid.Node(X,Y,dv)*lb.Cy[dv];
        Rho += grid.Node(X,Y,dv);
    }  

    Ux = Ux/Rho + 0.5*Force.Node(X,Y,0)*dt;
    Uy = Uy/Rho + 0.5*Force.Node(X,Y,1)*dt;


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
