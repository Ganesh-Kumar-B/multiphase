#ifndef ADVECTION
#define ADVECTION
//till second order
// // enum velocityDir{ 
// // dV_ZERO_ZERO,
// // dV_P1_ZERO, dV_ZERO_P1, dV_M1_ZERO, dV_ZERO_M1,
// // dV_P3_ZERO, dV_ZERO_P3, dV_M3_ZERO, dV_ZERO_M3,
// // dV_PH1_PH1, dV_MH1_PH1, dV_MH1_MH1, dV_PH1_MH1
// // }; 

#include<iomanip>
#include "GRID_2D.h"
#include "lbmD2Q9.h"




template<typename T>
void advection_D2Q9(Grid_N_C_2D<T> &gridf){


    for(int i = 0 + gridf.noghost ; i < gridf.n_x_node - (gridf.noghost);i++){
        for(int j = 0 + gridf.noghost ; j < gridf.n_y_node - (gridf.noghost);j++){ 
           
            
            gridf.Node(i,j,dV_ZERO_M1) = gridf.Node(i  ,j+1, dV_ZERO_M1 ); 
            gridf.Node(i,j,dV_M1_ZERO) = gridf.Node(i+1,j , dV_M1_ZERO);


            gridf.Node(i,j,dV_M1_M1) = gridf.Node(i+1,j+1,dV_M1_M1);
            gridf.Node(i,j,dV_M1_P1) = gridf.Node(i+1,j-1,dV_M1_P1);

        }
    }



    for(int i = gridf.n_x_node - (gridf.noghost +1); i> gridf.noghost-1; i--){                                                                                        //  #pragma omp parallel for shared(b)
        for(int j = gridf.n_y_node - (gridf.noghost +1); j> gridf.noghost-1; j--){

            gridf.Node(i,j,dV_ZERO_P1) = gridf.Node(i  ,j-1, dV_ZERO_P1 );
            gridf.Node(i,j,dV_P1_ZERO) = gridf.Node(i-1,j  , dV_P1_ZERO );

            gridf.Node(i,j,dV_P1_P1) = gridf.Node(i-1,j-1,dV_P1_P1);
            gridf.Node(i,j,dV_P1_M1) = gridf.Node(i-1,j+1,dV_P1_M1);

        }
    }


}




template<typename T>
void BB_top(Grid_N_C_2D<T> &gridf, lbmD2Q9<T> &lb9,real u0){

    int topN_last[]       =   { dV_ZERO_P1,  dV_P1_P1,   dV_M1_P1 };

    for(int i = 0 + gridf.nbx ; i <= gridf.nex;i++){
        
        int j = gridf.ney;
        for(int l = 0; l < 3; l++) 
                gridf.Node(i + (int)lb9.Cx[topN_last[l]], j +  (int)lb9.Cy[topN_last[l]], oppdV[topN_last[l]]) = gridf.Node(i,j, topN_last[l]);

    }
}

template<typename T>
void BB_bottom(Grid_N_C_2D<T> &gridf, lbmD2Q9<T> &lb9,real u0){

    int bottomN_last[]       =   { dV_ZERO_M1, dV_P1_M1,dV_M1_M1};

    for(int i = 0 + gridf.nbx ; i <= gridf.nex;i++){
        
        int j = gridf.nby;
        for(int l = 0; l < 3; l++) 
            gridf.Node(i + (int)lb9.Cx[bottomN_last[l]], j +  (int)lb9.Cy[bottomN_last[l]], oppdV[bottomN_last[l]]) = gridf.Node(i,j, bottomN_last[l]);

    }

}


template<typename T>
void BB_left(Grid_N_C_2D<T> &gridf, lbmD2Q9<T> &lb9,real u0){

    int leftN_last[]       =   {  dV_M1_ZERO,dV_M1_M1,dV_M1_P1};

    for(int j = 0 + gridf.nby ; j <= gridf.ney;   j++){
        int i = gridf.nbx;

        for(int l = 0; l < 3; l++) 
            gridf.Node(i + (int)lb9.Cx[leftN_last[l]], j +  (int)lb9.Cy[leftN_last[l]], oppdV[leftN_last[l]]) = gridf.Node(i,j, leftN_last[l]);


    }
}


template<typename T>
void BB_right(Grid_N_C_2D<T> &gridf, lbmD2Q9<T> &lb9,real u0){

    int rightN_last[]       =   {  dV_P1_ZERO,dV_P1_M1,dV_P1_P1 };

    for(int j = 0 + gridf.nby ; j <= gridf.ney;   j++){
        int i = gridf.nex;

        for(int l = 0; l < 3; l++) 
            gridf.Node(i + (int)lb9.Cx[rightN_last[l]], j +  (int)lb9.Cy[rightN_last[l]], oppdV[rightN_last[l]]) = gridf.Node(i,j, rightN_last[l]);


    }

}





template<typename T>
void Periodic_left_Right(Grid_N_C_2D<T> &gridf){

    for(int j = 0 ; j < gridf.n_y_node;j++ ){
        for(int dv = 0; dv<gridf.d_v; dv++){

            gridf.Node(0    ,j      ,dv)          = gridf.Node((gridf.n_x_node-gridf.noghost) -1, j, dv);
            // gridf.Node(gridf.n_x_node - gridf.noghost, j , dv) = gridf.Node(1,j,dv);

            gridf.Node(gridf.nex+1, j,dv) = gridf.Node(1, j,dv);
        }

    }

}





template<typename T>
void Periodic_top_bottom(Grid_N_C_2D<T> &gridf){

    for(int i = 0 ; i < gridf.n_x_node;i++ ){
        for(int dv = 0; dv<gridf.d_v; dv++){

            gridf.Node(i,0,dv) = gridf.Node(i,(gridf.n_y_node-gridf.noghost) -1, dv);
            gridf.Node(i,(gridf.n_y_node - gridf.noghost), dv) = gridf.Node(i,1,dv);

        }
    }
}



template<typename T>
void Grad_zero_left_Right(Grid_N_C_2D<T> &gridf){

    for(int j = 0 ; j < gridf.n_y_node;j++ ){
        for(int dv = 0; dv<gridf.d_v; dv++){

            gridf.Node(0,j,dv)                  = gridf.Node(gridf.nbx, j, dv);
            gridf.Node(gridf.nex +1 , j , dv)   = gridf.Node(gridf.nex, j, dv);

        }
    }
}



template<typename T>
void Grad_zero_top_bottom(Grid_N_C_2D<T> &gridf){

    for(int i = 0 ; i < gridf.n_x_node;i++ ){
        for(int dv = 0; dv<gridf.d_v; dv++){

            gridf.Node(i,0           , dv)  = gridf.Node(i,gridf.nby, dv);
            gridf.Node(i,2           , dv)  = gridf.Node(i,gridf.nby, dv);



            gridf.Node(i,gridf.ney +1, dv)  = gridf.Node(i,gridf.ney, dv);
            gridf.Node(i,gridf.ney -1, dv)  = gridf.Node(i,gridf.ney, dv);


        }
    }
}







template<typename T>
void Diffuse_BB_(Grid_N_C_2D<T> &grid,lbmD2Q9<T> &lbd2q9,real u0){    ///based on prasianakis based on the incoming populations
real nx; real ny,ux ,uy = 0,rho = 1.0,G = 0; //top
int j;u0= u0;
real feq[9] = {0};
for(int i = 1; i< grid.n_x_node-1 ; i++){


    //top
    nx = 0;  ny = -1;
    ux = u0; uy = 0;
    j = grid.n_y_node - grid.noghost;
    get_equi(feq,lbd2q9,ux,uy,rho);

    for(int dv = 0; dv<9; dv++){

    G = (grid.Node(i +0, j -1, dV_ZERO_P1) + grid.Node(i+1, j-1, dV_M1_P1) + grid.Node(i -1, j -1,dV_P1_P1))/
        (feq[dV_ZERO_M1] + feq[dV_M1_M1] + feq[dV_P1_M1]);

    if(lbd2q9.Cx[dv]*nx + lbd2q9.Cy[dv]*ny >0){        
        grid.Node(i,j,dv) = feq[dv] *G;

    }}}

//bottom
for(int i = 1; i< grid.n_x_node-1 ; i++){

    for(int dv = 0; dv< 9; dv++){
    nx = 0; ny = 1;
    ux = 0; uy = 0;
    j = 0;
    get_equi(feq, lbd2q9, ux, uy, rho);

    G = (grid.Node(i + 0, j +1, dV_ZERO_M1) + grid.Node(i+1,j+1,dV_M1_M1) + grid.Node(i-1, j+1, dV_P1_M1))/
        (feq[dV_ZERO_P1]  +  feq[dV_P1_P1]  + feq[dV_M1_P1]);

    if(lbd2q9.Cx[dv]*nx + lbd2q9.Cy[dv]*ny >0){
        
        grid.Node(i, j , dv) = feq[dv] *G;
        // std::cout<<G<<" "<<i<<" "<<dv<<std::endl;

    }

}}



// for(int dv= 0; dv<9; dv++){

//     std::cout<<grid.Node(3,grid.n_y_node-1,dv)<<std::endl;
// }

}







#endif