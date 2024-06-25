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
#include "collison.h"

template<typename T, typename T1>
void Multiphase_terms(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg, lbmD2Q9<T1> &lb,  Grid_N_C_2D<T> &rho, Grid_N_C_2D<T> &phi, Grid_N_C_2D<T> &munid, 
            Grid_N_C_2D<T> &laplacian_phi, real kappa, real A
            ){


    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            
                get_moments_Node_g(gridg, lb, phi.Node(i,j), i, j ); 


        }
    }

    Periodic_left_Right(phi); 
    Periodic_top_bottom(phi); 



    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){


                //> laplacian of phi  
                
                real del_t = 1.0;
                real Coeff = (2.0/(del_t*del_t*lb.theta0));
                
                laplacian_phi.Node(i,j)  = 0.0;

                for(int dv = 0; dv< 9; dv++)
                    laplacian_phi.Node(i,j) += lb.W[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv]  ) ;


                laplacian_phi.Node(i,j) = Coeff * ( laplacian_phi.Node(i,j) - phi.Node(i,j));

                //>---------------------





                //> munid
                munid.Node(i,j)   = - A * phi.Node(i,j) +  A * phi.Node(i,j) *  phi.Node(i,j) * phi.Node(i,j)   ;
                
                munid.Node(i,j) -= kappa*laplacian_phi.Node(i,j);


            
        }    
    }

    Periodic_left_Right(munid); 
    Periodic_top_bottom(munid); 


}



template<typename T, typename T1>
void Multiphase_Force_Node(Grid_N_C_2D<T> &grid, Grid_N_C_2D<T> &phi,   Grid_N_C_2D<T> &mu, 
            lbmD2Q9<T1> &lb, Grid_N_C_2D<T> &Force, int i, int j){



    Force.Node(i,j,0) = 0.0;
    Force.Node(i,j,1) = 0.0;


    real del_t = 1.0;
    real Coeff_grad = (1.0/(del_t*lb.theta0));

    for(int dv = 0; dv< 9; dv++){
         Force.Node(i,j,0) += lb.W[dv]*lb.Cx[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] ) ;
         Force.Node(i,j,1) += lb.W[dv]*lb.Cy[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] ) ;
    }


    Force.Node(i,j,0) = mu.Node(i,j)*Coeff_grad*( Force.Node(i,j,0));
    Force.Node(i,j,1) = mu.Node(i,j)*Coeff_grad*( Force.Node(i,j,1));  

}


;