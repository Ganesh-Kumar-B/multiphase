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

enum coodinates{X,Y,Z};

template<typename T, typename T1>
void Multiphase_terms(Grid_N_C_3D<T> &gridf,Grid_N_C_3D<T> &gridg, lbmD3Q35<T1> &lb,  Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &phi, Grid_N_C_3D<T> &munid, 
            Grid_N_C_3D<T> &laplacian_phi, real kappa, real A
            ){


    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            for(int k = 0 + gridf.noghost;k < gridf.n_z_node - (gridf.noghost) ; k++){
            
                get_moments_Node_g(gridg, lb, phi.Node(i,j,k), i, j ,k); 

                get_moments_Cell_g(gridg, lb, phi.Cell(i,j,k), i, j ,k); 

            }
        }
    }


    Periodic(phi);



    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            for(int k = 0 + gridf.noghost;k < gridf.n_z_node - (gridf.noghost) ; k++){


                //> laplacian of phi  
                
                real del_t = 1.0;
                real Coeff = (2.0/(del_t*del_t*lb.theta0));
                
                laplacian_phi.Node(i,j,k)  = 0.0;

                for(int dv = 0; dv< 27; dv++)
                    laplacian_phi.Node(i,j,k) += lb.W[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] , k + (int)lb.Cz[dv] ) ;

                for(int dv = 27; dv< 35;dv++)
                    laplacian_phi.Node(i,j,k) += lb.W[dv]*phi.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;

                laplacian_phi.Node(i,j,k) = Coeff * ( laplacian_phi.Node(i,j,k) - phi.Node(i,j,k));

                laplacian_phi.Cell(i,j,k)  = 0.0;
                
                for(int dv = 0; dv< 27; dv++)
                    laplacian_phi.Cell(i,j,k) += lb.W[dv]*phi.Cell( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] , k + (int)lb.Cz[dv] ) ;

                for(int dv = 27; dv< 35;dv++)
                    laplacian_phi.Cell(i,j,k) += lb.W[dv]*phi.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;

                laplacian_phi.Cell(i,j,k) = Coeff * ( laplacian_phi.Cell(i,j,k) - phi.Cell(i,j,k));
                //>---------------------






                //> munid
                munid.Node(i,j,k)   = - A * phi.Node(i,j,k) +  A * phi.Node(i,j,k) *  phi.Node(i,j,k) * phi.Node(i,j,k)   ;
                
                munid.Node(i,j,k) -= kappa*laplacian_phi.Node(i,j,k);

    

                munid.Cell(i,j,k) = - A * phi.Cell(i,j,k) +  A * phi.Cell(i,j,k) *  phi.Cell(i,j,k) * phi.Cell(i,j,k)   ;

                munid.Cell(i,j,k) -= kappa*laplacian_phi.Cell(i,j,k);

            }
        }    
    }

    Periodic(munid); 

}






template<typename T, typename T1>
void Multiphase_Force_Node(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &phi,   Grid_N_C_3D<T> &mu, 
            lbmD3Q35<T1> &lb, Grid_N_C_3D<T> &Force, int i, int j,int k){



    Force.Node(i,j,k,0) = 0.0;
    Force.Node(i,j,k,1) = 0.0;
    Force.Node(i,j,k,2) = 0.0;


    real del_t = 1.0;
    real Coeff_grad = (1.0/(del_t*lb.theta0));

    for(int dv = 0; dv< 27; dv++){
         Force.Node(i,j,k,0) += lb.W[dv]*lb.Cx[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv]  , k + (int)lb.Cz[dv]) ;
         Force.Node(i,j,k,1) += lb.W[dv]*lb.Cy[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv]  , k + (int)lb.Cz[dv]) ;
         Force.Node(i,j,k,2) += lb.W[dv]*lb.Cz[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv]  , k + (int)lb.Cz[dv]) ;
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        Force.Node(i,j,k,0) += lb.W[dv]*lb.Cx[dv]*phi.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv] , k + (int)lb.CzF[dv]) ;
        Force.Node(i,j,k,1) += lb.W[dv]*lb.Cy[dv]*phi.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv] , k + (int)lb.CzF[dv]) ;
        Force.Node(i,j,k,2) += lb.W[dv]*lb.Cz[dv]*phi.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv] , k + (int)lb.CzF[dv]) ;
    }

    Force.Node(i,j,k,0) = mu.Node(i,j,k)*Coeff_grad*( Force.Node(i,j,k,0));
    Force.Node(i,j,k,1) = mu.Node(i,j,k)*Coeff_grad*( Force.Node(i,j,k,1));  
    Force.Node(i,j,k,2) = mu.Node(i,j,k)*Coeff_grad*( Force.Node(i,j,k,2));  

}



template<typename T, typename T1>
void Multiphase_Force_Cell(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &phi,   Grid_N_C_3D<T> &mu, 
            lbmD3Q35<T1> &lb, Grid_N_C_3D<T> &Force, int i, int j,int k){

    Force.Cell(i,j,k,0) = 0.0;
    Force.Cell(i,j,k,1) = 0.0;
    Force.Cell(i,j,k,2) = 0.0;


    real del_t = 1.0;
    real Coeff_grad = (1.0/(del_t*lb.theta0));
    


    del_t = 1.0;
    Coeff_grad = (1.0/(del_t*lb.theta0));

    for(int dv = 0; dv< 27; dv++){
        Force.Cell(i,j,k,0) += lb.W[dv]*lb.Cx[dv]*phi.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        Force.Cell(i,j,k,1) += lb.W[dv]*lb.Cy[dv]*phi.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        Force.Cell(i,j,k,2) += lb.W[dv]*lb.Cz[dv]*phi.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
    }

    for(int dv = 27; dv< 35; dv++){
        Force.Cell(i,j,k,0) += lb.W[dv]*lb.Cx[dv]*phi.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
        Force.Cell(i,j,k,1) += lb.W[dv]*lb.Cy[dv]*phi.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
        Force.Cell(i,j,k,2) += lb.W[dv]*lb.Cz[dv]*phi.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
    }


    Force.Cell(i,j,k,0) = mu.Cell(i,j,k)* Coeff_grad*(Force.Cell(i,j,k,0));
    Force.Cell(i,j,k,1) = mu.Cell(i,j,k)* Coeff_grad*(Force.Cell(i,j,k,1));  
    Force.Cell(i,j,k,2) = mu.Cell(i,j,k)* Coeff_grad*(Force.Cell(i,j,k,2));  

    
    }





