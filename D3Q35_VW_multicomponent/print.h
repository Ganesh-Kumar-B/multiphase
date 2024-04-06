#pragma once

#include "lbmD3Q35.h"
#include "GRID_3D.h"
#include <fstream>
#include <iostream>
#include<stdio.h>
#include <sys/stat.h>
#include <sys/types.h>






template<typename T, typename T1>
void printdata(lbmD3Q35<T1> &lbModel,  Grid_N_C_3D<T> &gridLB,  int step, real u0)
{

    real u_inv,nx_inv, ny_inv;
    u_inv = 1/ u0;
    
    nx_inv = 1/(real)gridLB.n_x;
    ny_inv = 1/(real)gridLB.n_y;


    std::vector<int> line;bool isPresent;
    T u1,u2,u3,u4,um, rho1,rho2,del=0.05;
    std::ofstream file;
    char fileName[250];
    sprintf(fileName,"./Result/velocity_%d.txt",step) ;
    file.open(fileName) ;
    file<<"x,y,z,ux,uy,uz,rho"<<std::endl;

    for(int i = 0 + gridLB.noghost; i < gridLB.n_x_node - (gridLB.noghost); i++){ 
	    for (int j = 0 + gridLB.noghost; j < gridLB.n_y_node - (gridLB.noghost); j++){
            for(int k = 0 + gridLB.noghost; k < gridLB.n_z_node - (gridLB.noghost); k++){ 


                get_moments_Node(gridLB,lbModel,u1, u2,u3, rho1, i,j,k);
                
                file<<i<<"," <<j  <<","<<k  <<","<<u1<<","<<u2<<","<<u3<<","<<rho1<<std::endl;
                
                get_moments_Cell(gridLB,lbModel,u1, u2,u3, rho1, i,j,k);
                
                file<<i+0.5<<"," <<j+0.5  <<","<<k+0.5<<","<<u1<<","<<u2<<","<<u3<<","<<rho1<<std::endl;
                

            }
        }
    }
}

template<typename T, typename T1>
void print_vtk(lbmD3Q35<T1> &lb,  Grid_N_C_3D<T> &gridf ,Grid_N_C_3D<T> &gridg,  int step, real u0, real theta, Grid_N_C_3D<T> &Force)
{
   
    T u1,u2,u3,u4,um, rho1,rho2,del=0.05;

    std::ofstream file;
    char fileName[250];
    char foldername[250];
    sprintf(foldername,"Result_%.2f",theta);
    mkdir(foldername,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    sprintf(fileName,"./Result_%.2f/velocity_%d.vtk", theta,step) ;
    file.open(fileName);

    file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_POINTS"<<std::endl;
    
    file<<"DIMENSIONS "<<1*gridf.n_x<<" "<<1*gridf.n_y<<" "<<1*gridf.n_z<<std::endl;
    
    file<<"ORIGIN "<<0<<" "<<0<<" "<<0<<std::endl;
    file<<"SPACING "<<1<<" "<<1<<" "<<1<<std::endl;

    file<<"POINT_DATA "<<1*gridf.n_x*1*gridf.n_y*1*gridf.n_z<<std::endl;
    file<<"SCALARS density double\nLOOKUP_TABLE default"<<std::endl;

    for(int k = 0 + gridf.noghost; k < gridf.n_z_node - (gridf.noghost); k++){
	    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
            for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

                get_moments_Node_f(gridf,lb,u1, u2,u3, rho1, i,j,k,Force);
                file<<rho1<<std::endl;

            }
        }
    }

    file<<"SCALARS phi double\nLOOKUP_TABLE default"<<std::endl;


    for(int k = 0 + gridf.noghost; k < gridf.n_z_node - (gridf.noghost); k++){
	    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
            for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 
                
                real phi = 0;

                get_moments_Node_g(gridg,lb,phi,i,j,k);

                file<< phi <<std::endl;

            }
        }
    }


    file<<"VECTORS velocity double"<<std::endl;

    for(int k = 0 + gridf.noghost; k < gridf.n_z_node - (gridf.noghost); k++){
	    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
            for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

                get_moments_Node_f(gridf,lb,u1, u2,u3, rho1, i,j,k,Force);

                file<<u1<<" "<<u2<<" "<<u3<<std::endl;
                
            }
        }
    }


}
