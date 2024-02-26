#pragma once

#include"lbmD3Q35.h"
#include"GRID_3D.h"
#include"collison.h"
#include<iomanip>



template<typename T, typename T1>
void collide(Grid_N_C_3D<T> &grid,  lbmD3Q35<T1> &lb, real u0){

    double feq[35] ={0};

            int out_Node_last[] =       { dV_ZERO_P2_ZERO ,
                                        dV_P1_P1_ZERO   ,dV_M1_P1_ZERO  ,dV_ZERO_P1_P1  ,dV_ZERO_P1_M1  ,
                                        dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2    ,
                                        dV_PH1_PH1_PH1  ,dV_MH1_PH1_PH1 ,dV_PH1_PH1_MH1 , dV_MH1_PH1_MH1
                                        };

            int out_Node_secondLast[] = {dV_ZERO_P2_ZERO,
                                        dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2,  
                                        };
                                        
            int cell_last[]           = {dV_ZERO_P2_ZERO ,
                                        dV_P1_P1_ZERO   ,dV_M1_P1_ZERO  ,dV_ZERO_P1_P1  ,dV_ZERO_P1_M1  ,
                                        dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2    ,};

            int cell_secondlast[] =     {dV_ZERO_P2_ZERO,
                                        dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2,  
                                        }
                                        ;
    get_equi(feq,lb,u0,0.0,0.0,1.0);


    for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
        for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){ 

            //top wall
            double G_num_cell = 0, G_den_cell = 0, G_num_node = 0, G_den_node = 0,G_node = 0,G_cell = 0;





        }
    }    





}   