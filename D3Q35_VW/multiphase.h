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

template<typename T, typename T1>
void Multiphase_terms(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &pnid, Grid_N_C_3D<T> &fnid, Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho,
            lbmD3Q35<T1> &lb, double TbyTc, double kappa){

    double ux = 0, uy = 0, uz = 0;

    double rho_critical = 1.0, T_critical = lb.theta0/TbyTc ; 
    double b = 1.0/(3.0*rho_critical), a = b*T_critical*27.0/8.0;
    kappa = kappa*a;


    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){
            
                get_moments_Node(grid, lb,  ux, uy, uz, rho.Node(i,j,k), i, j ,k); 

                get_moments_Cell(grid, lb,  ux, uy, uz, rho.Cell(i,j,k), i, j ,k); 

            }
        }
    }


    Periodic(rho);



    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){



                //> laplacian of Rho  
                
                double del_t = 1.0;
                double Coeff = (2.0/(del_t*del_t*lb.theta0));
                
                laplacian_rho.Node(i,j,k)  = 0.0;

                for(int dv = 0; dv< 27; dv++)
                    laplacian_rho.Node(i,j,k) += lb.W[dv]*rho.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] , k + (int)lb.Cz[dv] ) ;

                for(int dv = 27; dv< 35;dv++)
                    laplacian_rho.Node(i,j,k) += lb.W[dv]*rho.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;

                laplacian_rho.Node(i,j,k) = Coeff * ( laplacian_rho.Node(i,j,k) - rho.Node(i,j,k));


                laplacian_rho.Cell(i,j,k)  = 0.0;
                
                for(int dv = 0; dv< 27; dv++)
                    laplacian_rho.Cell(i,j,k) += lb.W[dv]*rho.Cell( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] , k + (int)lb.Cz[dv] ) ;

                for(int dv = 27; dv< 35;dv++)
                    laplacian_rho.Cell(i,j,k) += lb.W[dv]*rho.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;

                laplacian_rho.Cell(i,j,k) = Coeff * ( laplacian_rho.Cell(i,j,k) - rho.Cell(i,j,k));





                //> pnid
                pnid.Node(i,j,k) = (rho.Node(i,j,k)*rho.Node(i,j,k)*b *lb.theta0 )/(1.0 - rho.Node(i,j,k)*b)  - 
                                a * rho.Node(i,j,k)*rho.Node(i,j,k) ;

                pnid.Cell(i,j,k) = (rho.Cell(i,j,k)*rho.Cell(i,j,k)*b *lb.theta0 )/(1.0 - rho.Cell(i,j,k)*b)  - 
                                a * rho.Cell(i,j,k)*rho.Cell(i,j,k) ;

                
                //> Fnid
                fnid.Node(i,j,k) =  - a * rho.Node(i,j,k)*rho.Node(i,j,k)   -   rho.Node(i,j,k)*lb.theta0 *log(1.0 - rho.Node(i,j,k)*b)
                                    ;

                fnid.Cell(i,j,k) =  - a * rho.Cell(i,j,k)*rho.Cell(i,j,k)   -   rho.Cell(i,j,k)*lb.theta0 *log(1.0 - rho.Cell(i,j,k)*b)
                                    
                                    ;
                
                
                //> munid
                munid.Node(i,j,k) = -lb.theta0*log(1.0 - rho.Node(i,j,k)*b) ;
                munid.Node(i,j,k) += rho.Node(i,j,k)*b*lb.theta0/(1.0 - rho.Node(i,j,k)*b);
                munid.Node(i,j,k) -= 2.0*rho.Node(i,j,k)*a; 
                munid.Node(i,j,k) -= kappa*laplacian_rho.Node(i,j,k);


                munid.Cell(i,j,k) = -lb.theta0*log(1.0 - rho.Cell(i,j,k)*b) ;
                munid.Cell(i,j,k) += rho.Cell(i,j,k)*b*lb.theta0/(1.0 - rho.Cell(i,j,k)*b);
                munid.Cell(i,j,k) -= 2.0*rho.Cell(i,j,k)*a ;
                munid.Cell(i,j,k) -= kappa*laplacian_rho.Cell(i,j,k);





            }
        }    
    }
    Periodic(laplacian_rho);
    Periodic(pnid);
    Periodic(fnid);
    Periodic(munid); 
}






template<typename T, typename T1>
void Multiphase_Force_Node(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &pnid, Grid_N_C_3D<T> &fnid, Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho,
            lbmD3Q35<T1> &lb, double &Fx, double &Fy, double &Fz, int i, int j,int k){

    //> CHEMICAL POTENTIAL FORMULATION 
    double grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;
    Fx = 0, Fy = 0, Fz = 0.0;
    double del_t = 1.0;
    double Coeff_grad = (1.0/(del_t*lb.theta0));

    for(int dv = 0; dv< 27; dv++){
        grad_mux += lb.W[dv]*lb.Cx[dv]*munid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        grad_muy += lb.W[dv]*lb.Cy[dv]*munid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        grad_muz += lb.W[dv]*lb.Cz[dv]*munid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        grad_mux += lb.W[dv]*lb.Cx[dv]*munid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;
        grad_muy += lb.W[dv]*lb.Cy[dv]*munid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;
        grad_muz += lb.W[dv]*lb.Cz[dv]*munid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;
    }

    Fx = - Coeff_grad*(grad_mux);
    Fy = - Coeff_grad*(grad_muy);  
    Fz = - Coeff_grad*(grad_muz);  


    // // //> MECHANICAL FORMULATION
    // // > direct              
    // double Fx  = 0.0, Fy = 0.0, Fz = 0.0;
    // double del_t = 1.0;
    // double Coeff_grad = (1.0/(del_t*lb.theta0));


    // for(int dv = 0; dv < grid.d_v; dv++){

    //     Fx += Coeff_grad*(               
    //                         lb.W[dv]*lb.Cx[dv]*( (pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Node( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
    //                         )
    //                         - kappa*Coeff_grad*( lb.W[dv]*lb.Cx[dv]*( laplacian_rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
    //                         ;

    //     Fy += Coeff_grad*(               
    //                         lb.W[dv]*lb.Cy[dv]*( (pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Node( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
    //                         )
    //                         - kappa*Coeff_grad*( lb.W[dv]*lb.Cy[dv]*( laplacian_rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
    //                         ;

    //      Fz += Coeff_grad*(               
    //                         lb.W[dv]*lb.Cz[dv]*( (pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Node( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
    //                         )
    //                         - kappa*Coeff_grad*( lb.W[dv]*lb.Cz[dv]*( laplacian_rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
    //                         ;

    // }
    
    
    // // // 4th order corrections
    // // double coeff_4th_grad = (lb.theta0*del_t*del_t)/2.0;

    // // for(int dv = 0; dv < grid.d_v; dv++){
    // //     Fx -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cx[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
    // //     Fy -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cy[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
    // //     Fz -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cz[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
    // // }

    // Fx = -1.0*Fx;
    // Fy = -1.0*Fy;
    // Fz = -1.0*Fz;
    //>--------------------------


}

template<typename T, typename T1>
void Multiphase_Force_Cell(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &pnid, Grid_N_C_3D<T> &fnid, Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho,
            lbmD3Q35<T1> &lb, double &Fx, double &Fy, double &Fz, int i, int j,int k){

    //> CHEMICAL POTENTIAL FORMULATION 
    double grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;
    Fx = 0, Fy = 0, Fz = 0.0;
    double del_t = 1.0;
    double Coeff_grad = (1.0/(del_t*lb.theta0));
    
    
    //> CHEMICAL POTENTIAL FORMULATION 
                 grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;
                 Fx = 0, Fy = 0, Fz = 0.0;
                 del_t = 1.0;
                 Coeff_grad = (1.0/(del_t*lb.theta0));

                for(int dv = 0; dv< 27; dv++){
                    grad_mux += lb.W[dv]*lb.Cx[dv]*munid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
                    grad_muy += lb.W[dv]*lb.Cy[dv]*munid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
                    grad_muz += lb.W[dv]*lb.Cz[dv]*munid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
                }

                for(int dv = 27; dv< 35; dv++){
                    grad_mux += lb.W[dv]*lb.Cx[dv]*munid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
                    grad_muy += lb.W[dv]*lb.Cy[dv]*munid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
                    grad_muz += lb.W[dv]*lb.Cz[dv]*munid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
                }


                Fx = - Coeff_grad*(grad_mux);
                Fy = - Coeff_grad*(grad_muy);  
                Fz = - Coeff_grad*(grad_muz);  

                //> MECHANICAL FORMULATION
                // // > direct              
                // Fx  = 0.0, Fy = 0.0, Fz = 0.0;
                // del_t = 1.0;
                // Coeff_grad = (1.0/(del_t*lb.theta0));


                // for(int dv = 0; dv < grid.d_v; dv++){

                //     Fx += Coeff_grad*(               
                //                         lb.W[dv]*lb.Cx[dv]*( (pnid.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Cell( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
                //                         )
                //                         - kappa*Coeff_grad*( lb.W[dv]*lb.Cx[dv]*( laplacian_rho.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
                //                         ;

                //     Fy += Coeff_grad*(               
                //                         lb.W[dv]*lb.Cy[dv]*( (pnid.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Cell( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
                //                         )
                //                         - kappa*Coeff_grad*( lb.W[dv]*lb.Cy[dv]*( laplacian_rho.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
                //                         ;

                //      Fz += Coeff_grad*(               
                //                         lb.W[dv]*lb.Cz[dv]*( (pnid.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]))/ rho.Cell( i+ lb.Cx[dv] ,j + lb.Cy[dv], k + lb.Cz[dv])  )
                //                         )
                //                         - kappa*Coeff_grad*( lb.W[dv]*lb.Cz[dv]*( laplacian_rho.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])))
                //                         ;

                // }
                
                
                // // // //> 4th order corrections
                // // //  coeff_4th_grad = (lb.theta0*del_t*del_t)/2.0;

                // // // for(int dv = 0; dv < grid.d_v; dv++){
                // // //     Fx -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cx[dv]*laplacian_pnidplusfnidbyrho.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
                // // //     Fy -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cy[dv]*laplacian_pnidplusfnidbyrho.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
                // // //     Fz -= coeff_4th_grad* Coeff_grad*(lb.W[dv]*lb.Cz[dv]*laplacian_pnidplusfnidbyrho.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv],k + lb.Cz[dv])) ;
                // // // }

                // Fx = -1.0*Fx;
                // Fy = -1.0*Fy;
                // Fz = -1.0*Fz;
                //> ------------------------------

    
    
    
    
    
    
    
    
    
    }