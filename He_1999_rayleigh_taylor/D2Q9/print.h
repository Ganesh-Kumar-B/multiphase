#ifndef PRINT
#define PRINT

#include "lbmD2Q9.h"
#include "GRID_2D.h"
#include <fstream>
#include <iostream>

#include <sys/stat.h>
#include <sys/types.h>




template<typename T, typename T1>
void print_vtk(lbmD2Q9<T1> &lb,  Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,  int step, double u0,double kappa , double theta)
{
    T u1,u2,u3,u4,um, rho1,rho2,del=0.05;

    std::ofstream file;
    char fileName[250];
    char foldername[250];
    sprintf(foldername,"kappare1000_%0.6f_%.2f",kappa,theta);
    mkdir(foldername,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    sprintf(fileName,"./kappare1000_%0.6f_%.2f/velocity_%d.vtk",kappa, theta,step) ;
    file.open(fileName);
    // file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_GRID"<<std::endl;
    file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_POINTS"<<std::endl;
    
    file<<"DIMENSIONS "<<1*gridf.n_x<<" "<<1*gridf.n_y<<std::endl;
    
    file<<"ORIGIN "<<0<<" "<<0<<" "<<0<<std::endl;
    file<<"SPACING "<<1<<" "<<1<<" "<<1<<std::endl;


    // file<<"POINTS "<<1*gridf.n_x*1*gridf.n_y*1*gridf.n_z<<" double"<<std::endl;

    
    // for(int k = 0 + gridf.noghost; k < gridf.n_z_node - (gridf.noghost); k++){
	//     for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
    //         for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

    //             file<<i<<" "<<j     <<" "<<k     <<std::endl;
    //             file<<i +0.5<<" "<<j +0.5<<" "<<k +0.5<<std::endl;


    //         }
    //     }
    // }
        

    file<<"POINT_DATA "<<1*gridf.n_x*1*gridf.n_y<<std::endl;
    file<<"SCALARS density double 1\nLOOKUP_TABLE default"<<1<<std::endl;

    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

            get_moments(gridf,lb,u1, u2, rho1, i,j);
            file<<rho1<<std::endl;

        }
    }
    

    file<<"VECTORS velocity double"<<std::endl;

    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

            get_moments(gridf,lb,u1, u2, rho1, i,j);

            file<<u1/u0<<" "<<u2/u0<<std::endl;

        }
    }
    



}








#endif