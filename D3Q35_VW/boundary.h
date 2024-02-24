#pragma once

#include"lbmD3Q35.h"
#include"GRID_3D.h"
#include<iomanip>



template<typename T, typename T1>
void collide(Grid_N_C_3D<T> &grid,  lbmD3Q35<T1> &lb, real u0){

    for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
        for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 




        }
    }    





}