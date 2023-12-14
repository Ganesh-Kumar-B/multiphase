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
#define PI 3.14159265



template<typename T, typename T1>
void collide(Grid_N_C_2D<T> &grid,Grid_N_C_2D<T> &rho,Grid_N_C_2D<T> &pnid,Grid_N_C_2D<T> &fnid,Grid_N_C_2D<T> &munid, Grid_N_C_2D<T> &laplacian_rho,
            lbmD2Q9<T1> &lb,double beta,double tau, double TbyTc, double kappa, int t ){

    double feq_Node[9] = {0}, ux = 0, uy = 0, G[9]={0};

    double rho_critical = 1.0, T_critical = lb.theta0/TbyTc ; 
    double b = 1.0/(3.0*rho_critical), a = b*T_critical*27.0/8.0;
    kappa = kappa*a;

    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            
            double Rho = 0;
            get_moments(grid, lb,  ux, uy,Rho, i, j); 
            rho.Node(i,j) = Rho;

        }
    }
    Periodic(rho);

   for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
     
            //> laplacian of Rho  
            laplacian_rho.Node(i,j)  = 0.0;
            double del_t = 1.0;
            double Coeff = (2.0/(del_t*del_t*lb.theta0));

            for(int dv = 0; dv< grid.d_v; dv++)
                laplacian_rho.Node(i,j) += lb.W[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;

            laplacian_rho.Node(i,j) = Coeff * ( laplacian_rho.Node(i,j) - rho.Node(i,j));

            //> gradient of Rho
            double grad_rhox = 0.0,grad_rhoy = 0.0;
            del_t = 1.0;
            Coeff = (1.0/(del_t*lb.theta0));

            for(int dv = 0; dv< grid.d_v; dv++){
                grad_rhox += lb.W[dv]*lb.Cx[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;
                grad_rhoy += lb.W[dv]*lb.Cy[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;
            }
            double grad_rho = Coeff* (grad_rhox + grad_rhoy);


            //> pnid
            pnid.Node(i,j) = (rho.Node(i,j)*rho.Node(i,j)*b *lb.theta0 )/(1.0 - rho.Node(i,j)*b)  - 
                            a * rho.Node(i,j)*rho.Node(i,j) ;
            
            //> Fnid
            fnid.Node(i,j) =  - a * rho.Node(i,j)*rho.Node(i,j)   -   rho.Node(i,j)*lb.theta0 *log(1.0 - rho.Node(i,j)*b)
                                -kappa * 0.5*grad_rho*grad_rho
                                ;
            
            //> munid
            munid.Node(i,j) = -lb.theta0*log(1.0 - rho.Node(i,j)*b) ;
            munid.Node(i,j) += rho.Node(i,j)*b*lb.theta0/(1.0 - rho.Node(i,j)*b);
            munid.Node(i,j) -= 2.0*rho.Node(i,j)*a 
                                -kappa*laplacian_rho.Node(i,j)
                                ;
        }    
    }

    Periodic(laplacian_rho);
    Periodic(pnid);
    Periodic(fnid);
    Periodic(munid); 





    //< Collision
    /// first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            
            double Rho = 0.0;

            //> CHEMICAL POTENTIAL FORMULATION 
            double grad_mux = 0.0, grad_muy = 0.0;
            double Fx = 0, Fy = 0;
            double del_t = 1.0;
            double Coeff_grad = (1.0/(del_t*lb.theta0));

            for(int dv = 0; dv< grid.d_v; dv++){
                grad_mux += lb.W[dv]*lb.Cx[dv]*munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;
                grad_muy += lb.W[dv]*lb.Cy[dv]*munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;

            }

            Fx = - Coeff_grad*(grad_mux);
            Fy = - Coeff_grad*(grad_muy);    
            

            //> MECHANICAL FORMULATION
            //> SECOND ORDER 
            // double Fx  = 0, Fy = 0;
            // double del_t = 1.0;
            // double Coeff_grad = (1.0/(del_t*lb.theta0));

            // for(int dv = 0; dv < grid.d_v; dv++){

            //     Fx += Coeff_grad*(  (1.0/rho.Node(i,j))               
            //                         *( lb.W[dv]*lb.Cx[dv]*( pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) + munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) )) 
            //                          +
            //                         (pnid.Node(i,j) + munid.Node(i,j))
            //                         *( lb.W[dv]*lb.Cx[dv]*(1.0/rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv])))
            //                         )
            //                         - kappa*Coeff_grad*( lb.W[dv]*lb.Cx[dv]*( laplacian_rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv])))
            //                         ;

            //     std::cout<<(  (1.0/rho.Node(i,j))               
            //                         *( lb.W[dv]*lb.Cx[dv]*( pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) + munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) )) 
            //                          +
            //                         (pnid.Node(i,j) + munid.Node(i,j))
            //                         *( lb.W[dv]*lb.Cx[dv]*(1.0/rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv])))
            //                         )<<std::endl;

            //     Fy += Coeff_grad*(  (1.0/rho.Node(i,j))               *( lb.W[dv]*lb.Cy[dv]*( pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) + munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ))  +
            //                         (pnid.Node(i,j) + munid.Node(i,j))*( lb.W[dv]*lb.Cy[dv]*(1.0/rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv])))
            //                         )
            //                         - kappa*Coeff_grad*( lb.W[dv]*lb.Cy[dv]*( laplacian_rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv])))
            //                         ;                
                                
                                
            // }

            Fx = 0;
            Fy = 0;
            // std::cout<<Fx<<std::endl;

            get_moments(grid, lb,  ux, uy,Rho, i, j, Fx, Fy );            //for the node
            get_equi(feq_Node,lb, ux, uy, Rho);

            for (int dv = 0; dv< 9; dv++){
                 
                grid.Node(i,j,dv) =  grid.Node(i,j,dv) + 2.0* beta*(feq_Node[dv] - grid.Node(i,j,dv))
                                    + 2.0 *beta * tau*lb.thetaInverse * rho.Node(i,j)* lb.W[dv] * (Fx * lb.Cx[dv] + Fy * lb.Cy[dv]) ;
         
            } 
        }
    }
}
    





template<typename T>
void get_equi(double feq[9], lbmD2Q9<T> &lbD2Q9, double ux, double uy, double rho){


    double u2 = ux*ux + uy*uy;
     double a1=0;
    double first,second, third,feq0=0;
    for (int dv = 0; dv< 9; dv++){

        feq0 = rho*lbD2Q9.W[dv];

        first  = (ux*lbD2Q9.Cx[dv] + uy*lbD2Q9.Cy[dv])*lbD2Q9.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lbD2Q9.thetaInverse;
        feq[dv] = feq0*(1+ first + second + third);    


        
    }
         
}





template<typename T,typename T1>
void get_moments(Grid_N_C_2D<T> &lbgrid, lbmD2Q9<T1> &lbD2Q9,double &Ux, double &Uy,double &Rho,  int X, int Y, double Fx = 0  , double Fy = 0){ ///node or cell 0-Node 1- cell
   Ux  = 0.0;
   Uy  = 0.0;
   Rho = 0.0;


    for(int dv = 0; dv <9; dv++){
        Ux  += lbgrid.Node(X,Y,dv)*lbD2Q9.Cx[dv];
        Uy  += lbgrid.Node(X,Y,dv)*lbD2Q9.Cy[dv];
        Rho += lbgrid.Node(X,Y,dv);
    }  

    Ux = Ux/Rho + 0.5*Fx;
    Uy = Uy/Rho + 0.5*Fy;  

}

//period no of waves
template<typename T, typename T1>
void initialization(Grid_N_C_2D<T> &lbgrid,lbmD2Q9<T1> &lbD2Q9,double Rho_mean ,double amplitude, double period){


    double Feq_node[9] = {0},Rho;
    double x,y
           ;    ///distance between nodes 
    
    double ux_node =0, uy_node = 0;
    
    double k = 2*M_PI* period;
    
    for(int i = 0 + lbgrid.noghost; i < lbgrid.n_x_node - (lbgrid.noghost); i++){
        for(int j = 0 + lbgrid.noghost; j < lbgrid.n_y_node - (lbgrid.noghost); j++){
            
            x = ((double)i)/ lbgrid.n_x;
            y = ((double)j)/ lbgrid.n_y;

            Rho = Rho_mean + amplitude * sin(k* x)* sin(k*y);

            get_equi(Feq_node,lbD2Q9,ux_node,uy_node,Rho);

            for (int dv = 0; dv<9; dv++){
                lbgrid.Node(i,j,dv) = Feq_node[dv]; 
            }
        }
    }
    
    
        
}
        



template<typename T>
void printMass(Grid_N_C_2D<T> &lbgrid){    
    double a = 0;
    for(int i = 0 + lbgrid.noghost; i < lbgrid.n_x_node - (lbgrid.noghost) ; i++){
        for(int j = 0 + lbgrid.noghost;j < lbgrid.n_y_node - (lbgrid.noghost) ; j++){
            for (int dv = 0; dv< 9; dv++){
                a += lbgrid.Node(i,j,dv) ;
            }
        }}
std::cout<<"   "<<a<<std::endl;

}


;