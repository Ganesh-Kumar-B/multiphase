#pragma once

#include "lbmD3Q35.h"
#include "GRID_3D.h"
#include <fstream>
#include <iostream>
#include<stdio.h>







template<typename T, typename T1>
void printdata(lbmD3Q35<T1> &lbModel,  Grid_N_C_3D<T> &gridLB,  int step, double u0)
{

    double u_inv,nx_inv, ny_inv;
    u_inv = 1/ u0;
    
    nx_inv = 1/(double)gridLB.n_x;
    ny_inv = 1/(double)gridLB.n_y;


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
void print_vtk(lbmD3Q35<T1> &lb,  Grid_N_C_3D<T> &grid,  int step, double u0)
{


    T u1,u2,u3,u4,um, rho1,rho2,del=0.05;

    std::ofstream file;
    char fileName[250];
    sprintf(fileName,"./Result/velocity_%d.vtk",step) ;
    file.open(fileName) ;
    file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_GRID"<<std::endl;
    
    file<<"DIMENSIONS "<<2*grid.n_x<<" "<<1*grid.n_y<<" "<<1*grid.n_z<<std::endl;
    
    // file<<"ORIGIN "<<0<<" "<<0<<" "<<0<<std::endl;
    // file<<"SPACING "<<1<<" "<<1<<" "<<1<<std::endl;


    file<<"POINTS "<<2*grid.n_x*1*grid.n_y*1*grid.n_z<<" double"<<std::endl;

    
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){ 
	    for (int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){ 


                file<<i<<" "<<j     <<" "<<k     <<std::endl;
                file<<i +0.5<<" "<<j +0.5<<" "<<k +0.5<<std::endl;


            }
        }
    }
        

    file<<"POINT_DATA "<<2*grid.n_x*1*grid.n_y*1*grid.n_z<<std::endl;
    file<<"SCALARS density double 1\nLOOKUP_TABLE default"<<1<<std::endl;

    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){ 
	    for (int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){ 

                get_moments_Node(grid,lb,u1, u2,u3, rho1, i,j,k);

                file<<rho1<<std::endl;

                get_moments_Cell(grid,lb,u1, u2,u3, rho1, i,j,k);

                file<<rho1<<std::endl;
            }
        }
    }

    // file<<"VECTORS velocity double"<<std::endl;

    // for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){ 
	//     for (int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
    //         for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){ 

    //             get_moments_Node(grid,lb,u1, u2,u3, rho1, i,j,k);

    //             file<<u1<<" "<<u2<<" "<<u3<<std::endl;


    //         }
    //     }
    // }



}
