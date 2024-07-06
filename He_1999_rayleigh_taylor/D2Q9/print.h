#ifndef PRINT
#define PRINT

#include "lbmD2Q9.h"
#include "GRID_2D.h"
#include <fstream>
#include <iostream>

#include <sys/stat.h>
#include <sys/types.h>




template<typename T, typename T1>
void print_vtk(lbmD2Q9<T1> &lb9,  Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,Grid_N_C_2D<T> &grad_psi_rho,
                double step, double u0,double kappa,
                double phi_l ,double phi_h, double rho_l , double rho_h,Grid_N_C_2D<T> &Force,const std::string &name )
{
    T u1,u2,p,um, rho1,rho2,del=0.05;


    double phi = 0.0, rho = 0.0;

    std::ofstream file;
    char fileName[250];
    char foldername[250];

    

    sprintf(foldername,"%s_%.2f_",name.c_str(),kappa);
    mkdir(foldername,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    sprintf(fileName,"./%s_%.2f_/velocity_%.6f.vtk",name.c_str(),kappa,0.01*step) ;
    file.open(fileName);



    // file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_GRID"<<std::endl;
    file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_POINTS"<<std::endl;
    
    file<<"DIMENSIONS "<<1*gridf.n_x<<" "<<1*gridf.n_y<<" "<<1<<std::endl;
    
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
        

    // file<<"POINT_DATA "<<1*gridf.n_x*1*gridf.n_y<<std::endl;
    // file<<"SCALARS density double 1\nLOOKUP_TABLE default"<<1<<std::endl;

    // for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
    //     for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

    //         get_moments(gridf,lb9,u1, u2, rho1, i,j);
    //         file<<rho1<<std::endl;

    //     }
    // }

    file<<"POINT_DATA "<<1*gridf.n_x*1*gridf.n_y<<std::endl;
    file<<"SCALARS phi double 1\nLOOKUP_TABLE default"<<std::endl;

    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 
            
            
            get_phi(gridf,lb9,phi,i,j);
            file<<phi<<std::endl;

        }
    }
    

    file<<"SCALARS rho double 1\nLOOKUP_TABLE default"<<std::endl;

    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

            
            get_phi(gridf,lb9,phi,i,j);

            rho = rho_l + ((phi - phi_l)/(phi_h - phi_l)) *(rho_h - rho_l);

            file<<rho<<std::endl;

        }
    }


    file<<"SCALARS pressure double 1\nLOOKUP_TABLE default"<<std::endl;

    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 
            
            get_phi(gridf,lb9,phi,i,j);

            rho = rho_l + ((phi - phi_l)/(phi_h - phi_l)) *(rho_h - rho_l);

            get_P_and_u(gridg,grad_psi_rho, lb9,  u1, u2,p,rho, i, j, Force);

            file<<p<<std::endl;

        }
    }



    file<<"VECTORS velocity double"<<std::endl;

    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

            get_phi(gridf,lb9,phi,i,j);

            rho = rho_l + ((phi - phi_l)/(phi_h - phi_l)) *(rho_h - rho_l);

            get_P_and_u(gridg,grad_psi_rho, lb9,  u1, u2,p,rho, i, j, Force);

            file<<u1<<" "<<u2<<std::endl;
        }
    }
    



}








#endif