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
void collide(Grid_N_C_2D<T> &grid,Grid_N_C_2D<T> &rho,Grid_N_C_2D<T> &pnid,Grid_N_C_2D<T> &fnid,Grid_N_C_2D<T> &munid,
            lbmD2Q9<T1> &lb,double beta,double tau, double TbyTc, double kappa ){

    double feq_Node[9] = {0},feq_Cell[9]={0}, ux = 0, uy = 0, G[9]={0};

    double rho_critical = 1.0, T_critical = lb.theta0/TbyTc ; 
    double b = 1.0/(3.0*rho_critical), a = b*T_critical*27.0/8.0;


    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            
            double Rho = 0;
            get_moments(grid, lb,  ux, uy,Rho, i, j); 
            rho.Node(i,j) = Rho;

        }
    }

    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
     
            //> laplacian of Rho  
            double laplacian_rho  = 0.0;
            double del_t = 1.0;
            double Coeff = (2.0/(del_t*del_t*lb.theta0));

            for(int dv = 0; dv< grid.d_v; dv++)
                laplacian_rho += lb.W[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;

            laplacian_rho = Coeff * ( laplacian_rho - rho.Node(i,j));

            //> gradient of Rho
            // double grad_rhox = 0.0,grad_rhoy = 0.0;
            // del_t = 1.0;
            // Coeff = (1.0/(del_t*lb.theta0));

            // for(int dv = 0; dv< grid.d_v; dv++){
            //     grad_rhox += lb.W[dv]*lb.Cx[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;
            //     grad_rhoy += lb.W[dv]*lb.Cy[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;
            // }
            // double grad_rho = Coeff* (grad_rhox + grad_rhoy);


            // //> pnid
            // pnid.Node(i,j) = (rho*rho*b *lb.theta0 )/(1.0 - rho*b)  - 
            //                 a * rho*rho ;
            
            // //> Fnid
            // fnid.Node(i,j) =  - a * rho*rho   -   rho*lb.theta0 *log(1.0 - rho*b)
            //              // - myVDW.kappa * 0.5*pow((myLattice[iX+1].rho -myLattice[iX-1].rho)/(2.0),2)
            //                     ;
            
            //> munid
            munid.Node(i,j) = -lb.theta0*log(1.0 - rho.Node(i,j)*b) ;
            munid.Node(i,j) += rho.Node(i,j)*b*lb.theta0/(1.0 - rho.Node(i,j)*b);
            munid.Node(i,j) -= 2.0*rho.Node(i,j)*a 
                                -kappa*laplacian_rho
                                ;


        }    
    }





/// first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){

            double Rho = 0.0;
            double grad_mux = 0.0, grad_muy = 0.0;
            double del_t = 1.0;
            double Coeff_grad = (1.0/(del_t*lb.theta0));

            for(int dv = 0; dv< grid.d_v; dv++){
                grad_mux += lb.W[dv]*lb.Cx[dv]*munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;
                grad_muy += lb.W[dv]*lb.Cy[dv]*munid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv]) ;
            }

            grad_mux = Coeff_grad*(grad_mux);
            grad_muy = Coeff_grad*(grad_muy);    

            

            get_moments(grid, lb,  ux, uy,Rho, i, j);            //for the node
            get_equi(feq_Node,lb, ux, uy, Rho);

            for (int dv = 0; dv< 9; dv++){
                 
                grid.Node(i,j,dv) =  grid.Node(i,j,dv) + 2.0* beta*(feq_Node[dv] - grid.Node(i,j,dv))
                                    + 2.0 *beta * tau*lb.thetaInverse * rho.Node(i,j)* lb.W[dv] * (grad_mux * lb.Cx[dv] + grad_muy * lb.Cy[dv]) ;
         
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
void get_moments(Grid_N_C_2D<T> &lbgrid, lbmD2Q9<T1> &lbD2Q9,double &Ux, double &Uy,double &Rho,  int X, int Y){ ///node or cell 0-Node 1- cell
   Ux  = 0.0;
   Uy  = 0.0;
   Rho = 0.0;


    for(int dv = 0; dv <9; dv++){
        Ux  += lbgrid.Node(X,Y,dv)*lbD2Q9.Cx[dv];
        Uy  += lbgrid.Node(X,Y,dv)*lbD2Q9.Cy[dv];
        Rho += lbgrid.Node(X,Y,dv);
    }  

    Ux = Ux/Rho;
    Uy = Uy/Rho;  

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