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




template<typename T>
void advection_D2Q9(Grid_N_C_2D<T> &grid){


    for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
    for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 
           
            
            grid.Node(i,j,dV_ZERO_M1) = grid.Node(i  ,j+1, dV_ZERO_M1 ); 
            grid.Node(i,j,dV_M1_ZERO) = grid.Node(i+1,j , dV_M1_ZERO);


            grid.Node(i,j,dV_M1_M1) = grid.Node(i+1,j+1,dV_M1_M1);
            grid.Node(i,j,dV_M1_P1) = grid.Node(i+1,j-1,dV_M1_P1);


}}

        for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){                                                                                        //  #pragma omp parallel for shared(b)
        for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){

            grid.Node(i,j,dV_ZERO_P1) = grid.Node(i  ,j-1, dV_ZERO_P1 );
            grid.Node(i,j,dV_P1_ZERO) = grid.Node(i-1,j  , dV_P1_ZERO );

            grid.Node(i,j,dV_P1_P1) = grid.Node(i-1,j-1,dV_P1_P1);
            grid.Node(i,j,dV_P1_M1) = grid.Node(i-1,j+1,dV_P1_M1);

}}
}

template<typename T>
void cavity_walls(Grid_N_C_2D<T> &grid,double u0){

//bottom wall
for(int i = grid.noghost; i< grid.n_x_node - grid.noghost ; i++){
int j = 0;

grid.Node(i,j,dV_ZERO_P1) = grid.Node(i,j+1,dV_ZERO_M1);
grid.Node(i,j,dV_P1_P1)   = grid.Node(i+1,j+1, dV_M1_M1);
grid.Node(i,j,dV_M1_P1)   = grid.Node(i-1,j+1, dV_P1_M1);
}

//top
for(int i = grid.noghost; i< grid.n_x_node - grid.noghost ; i++){
int j = grid.n_y_node - grid.noghost;
grid.Node(i,j,dV_ZERO_M1) = grid.Node(i,j-1, dV_ZERO_P1) ;
grid.Node(i,j,dV_M1_M1)   = grid.Node(i-1, j-1, dV_P1_P1) - (1.0/6.0)*u0;
grid.Node(i,j,dV_P1_M1)   = grid.Node(i+1, j -1, dV_M1_P1)+ (1.0/6.0)*u0;

}


//left
for(int j = grid.noghost; j< grid.n_y_node - grid.noghost ; j++){
int i = 0;

grid.Node(i,j,dV_P1_ZERO) = grid.Node(i +1, j , dV_M1_ZERO);
grid.Node(i,j,dV_P1_P1)   = grid.Node(i +1,j +1, dV_M1_M1);
grid.Node(i,j,dV_P1_M1)   = grid.Node(i+1, j -1, dV_M1_P1);


}


//right
for(int j = grid.noghost; j< grid.n_y_node - grid.noghost ; j++){
int i = grid.n_x_node - grid.noghost;

grid.Node(i,j,dV_M1_ZERO) = grid.Node(i -1, j , dV_P1_ZERO);
grid.Node(i,j,dV_M1_P1)   = grid.Node(i -1, j+1, dV_P1_M1);
grid.Node(i,j,dV_M1_M1)   = grid.Node(i-1, j-1, dV_P1_P1);



}


}






template<typename T>
void Periodic_left_Right(Grid_N_C_2D<T> &grid){

    for(int j = 0 ; j < grid.n_y_node;j++ ){
        for(int dv = 0; dv<9; dv++){

            grid.Node(0,j,dv) = grid.Node((grid.n_x_node-grid.noghost) -1, j, dv);
            grid.Node((grid.n_x_node - grid.noghost), j , dv) = grid.Node(1,j,dv);

        }

    }


}



template<typename T>
void Periodic_top_bottom(Grid_N_C_2D<T> &grid){

      for(int i = 0 ; i < grid.n_x_node;i++ ){
        for(int dv = 0; dv<9; dv++){

            grid.Node(i,0,dv) = grid.Node(i,(grid.n_y_node-grid.noghost) -1, dv);
            grid.Node(i,(grid.n_y_node - grid.noghost), dv) = grid.Node(i,1,dv);

        }

      }}





template<typename T>
void Diffuse_BB_(Grid_N_C_2D<T> &grid,lbmD2Q9<T> &lbd2q9,double u0){    ///based on prasianakis based on the incoming populations
double nx; double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
int j;u0= u0;
double feq[9] = {0};
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





template<typename T>
void Diffuse_B1(Grid_N_C_2D<T> &grid,lbmD2Q9<T> &lbd2q9,double u0){   //looping over the outgoing populations

    double nx; double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
int j,ic= 0 ,jc = 0;u0= u0;
double feq[9] = {0};



for(int i = 2; i< grid.n_x_node-2;i++){
    j = 0;
    for(int dv = 0; dv<9;dv ++){
        grid.Node(i,j,dv) = 0;
        
    }

    j = grid.n_y_node - 1;
    for(int dv =0 ; dv<9; dv++){
        grid.Node(i,j,dv) = 0;
    }
}



    /////===============top============///////
    nx = 0;  ny = -1;
    ux = u0; uy = 0;
    j = grid.n_y_node - grid.noghost -1;
    get_equi(feq, lbd2q9, ux, uy, rho);
     G = 1/(feq[dV_ZERO_M1] + feq[dV_M1_M1] + feq[dV_P1_M1]);

    int dv = dV_ZERO_P1;
    for(int i = 2; i< grid.n_x_node -2;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    dv = dV_M1_P1;
    for(int i = 3; i< grid.n_x_node -1;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    dv = dV_P1_P1;
    for(int i = 1; i< grid.n_x_node -3;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    

    ///////////////////////=========bottom============////////////////
    nx = 0;  ny = 1;
    ux = 0.0; uy = 0;
    j = 1;
    get_equi(feq, lbd2q9, ux, uy, rho);
     G = 1/(feq[dV_ZERO_P1] + feq[dV_M1_P1] + feq[dV_P1_P1]);

    dv = dV_ZERO_M1;
    for(int i = 2; i< grid.n_x_node -2;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    dv = dV_M1_M1;
    for(int i = 3; i< grid.n_x_node -1;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    dv = dV_P1_M1;
    for(int i = 1; i< grid.n_x_node -3;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }       


    
    ///////=============top
    j = grid.n_y_node - grid.noghost;
    for(int i = 1; i< grid.n_x_node -grid.noghost-1;i++){



        grid.Node(i,j,dV_ZERO_M1)   = (-1.0/0.6)*(grid.Node(i,j-1,dV_ZERO_P1) - grid.Node(i,j,dV_ZERO_M1)) + grid.Node(i,j-1,dV_ZERO_P1);
        grid.Node(i,j,dV_P1_M1)     = (-1.0/0.6)*(grid.Node(i,j-1,dV_P1_P1) - grid.Node(i,j,dV_P1_M1)) + grid.Node(i,j-1,dV_P1_P1);
        grid.Node(i,j,dV_M1_M1)     = (-1.0/0.6)*(grid.Node(i,j-1,dV_M1_P1) - grid.Node(i,j,dV_M1_M1)) + grid.Node(i,j-1,dV_M1_P1);

        // for(int dv=0; dv<grid.d_v;dv++){
        //     grid.Node(i,j,dv) = 0.8*(grid.Node((int)lbd2q9.Cx[oppdV[dv]],(int)lbd2q9.Cy[oppdV[dv],oppdV[dv]])-
        //                              grid.Node(i,j,dv)) + grid.Node((int)lbd2q9.Cx[oppdV[dv]],(int)lbd2q9.Cy[oppdV[dv],oppdV[dv]]);
        // }
    
    }
    /////================bottom
    j = 0;
    for(int i = 1; i< grid.n_x_node -grid.noghost-1;i++){
        grid.Node(i,j,dV_ZERO_P1)   = (-1.0/0.6)*(grid.Node(i,j+1,dV_ZERO_M1)  - grid.Node(i,j,dV_ZERO_P1))+ grid.Node(i,j+1,dV_ZERO_M1);
        grid.Node(i,j,dV_P1_P1)     = (-1.0/0.6)*(grid.Node(i,j+1,dV_P1_M1)    - grid.Node(i,j,dV_P1_P1))  + grid.Node(i,j+1,dV_P1_M1);
        grid.Node(i,j,dV_M1_P1)     = (-1.0/0.6)*(grid.Node(i,j+1,dV_M1_M1)    - grid.Node(i,j,dV_M1_P1))  + grid.Node(i,j+1,dV_M1_M1);
    }

}



//note /////////===========================================================1ST ORDER============================================//
//note /////////===========================================================1ST ORDER============================================//
//note /////////===========================================================1ST ORDER============================================//
//note /////////===========================================================1ST ORDER============================================//
//note /////////===========================================================1ST ORDER============================================//

//1st order
template<typename T>
void Diff_non_align_general(Grid_N_C_2D<T> &grid,lbmD2Q9<T> &lbd2q9,double u0,double del){ //1st order

    double nx; double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
int j,ic= 0 ,jc = 0;u0= u0;
double feq[9] = {0};


for(int i = 2; i< grid.n_x_node-2;i++){
    j = 0;
    for(int dv = 0; dv<9;dv ++){
        grid.Node(i,j,dv) = 0;        
    }

    j = grid.n_y_node - 1;
    for(int dv =0 ; dv<9; dv++){
        grid.Node(i,j,dv) = 0;
    }
}



    /////===============top============///////
    nx = 0;  ny = -1;
    ux = u0; uy = 0;
    j = grid.n_y_node - grid.noghost -1;
    get_equi(feq, lbd2q9, ux, uy, rho);
     G = 1/(feq[dV_ZERO_M1] + feq[dV_M1_M1] + feq[dV_P1_M1]);

    int dv = dV_ZERO_P1;
    for(int i = 2; i< grid.n_x_node -2;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    dv = dV_M1_P1;
    for(int i = 3; i< grid.n_x_node -1;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    dv = dV_P1_P1;
    for(int i = 1; i< grid.n_x_node -3;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }



    ///////////////////////=========bottom============////////////////
    nx = 0;  ny = 1;
    ux = 0.0; uy = 0;
    j = 1;
    get_equi(feq, lbd2q9, ux, uy, rho);
     G = 1/(feq[dV_ZERO_P1] + feq[dV_M1_P1] + feq[dV_P1_P1]);

    dv = dV_ZERO_M1;
    for(int i = 2; i< grid.n_x_node -2;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    dv = dV_M1_M1;
    for(int i = 3; i< grid.n_x_node -1;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }

    dv = dV_P1_M1;
    for(int i = 1; i< grid.n_x_node -3;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G;
                }
            }
    }       
std::cout<<std::setprecision(10);
std::cout<<grid.Node(50,grid.n_y_node - grid.noghost    ,dV_M1_M1)<<" x 0.8 "<<std::endl;
std::cout<<grid.Node(50,grid.n_y_node - grid.noghost-1  ,dV_M1_P1)<<"x = 0"<<std::endl;

    

    ///////=============top
    j = grid.n_y_node - grid.noghost;
    for(int i = 1; i< grid.n_x_node -grid.noghost-1;i++){


        grid.Node(i,j,dV_ZERO_M1)   = (-1.0/(2*del))*(grid.Node(i,j-1,dV_ZERO_P1) - grid.Node(i,j,dV_ZERO_M1)) + grid.Node(i,j-1,dV_ZERO_P1);
        grid.Node(i,j,dV_P1_M1)     = (-1.0/(2*del))*(grid.Node(i,j-1,dV_P1_P1) - grid.Node(i,j,dV_P1_M1)) + grid.Node(i,j-1,dV_P1_P1);
        grid.Node(i,j,dV_M1_M1)     = (-1.0/(2*del))*(grid.Node(i,j-1,dV_M1_P1) - grid.Node(i,j,dV_M1_M1)) + grid.Node(i,j-1,dV_M1_P1);

        // for(int dv=0; dv<grid.d_v;dv++){
        //     grid.Node(i,j,dv) = 0.8*(grid.Node((int)lbd2q9.Cx[oppdV[dv]],(int)lbd2q9.Cy[oppdV[dv],oppdV[dv]])-
        //                              grid.Node(i,j,dv)) + grid.Node((int)lbd2q9.Cx[oppdV[dv]],(int)lbd2q9.Cy[oppdV[dv],oppdV[dv]]);
        // }
    
    }


    for(int i = 3; i< grid.n_x_node -1;i++){
    
      j = grid.n_y_node - grid.noghost;  
        std::cout<< grid.Node(i,j,4)<<" "<<i<<" "<<4 <<std::endl;
        std::cout<< grid.Node(i,j,7)<<" "<<i<<" "<<7 <<std::endl;
        std::cout<< grid.Node(i,j,8)<<" "<<i<<" "<<8 <<std::endl;
    }

    /////================bottom
    j = 0;
    for(int i = 1; i< grid.n_x_node -grid.noghost-1;i++){
        grid.Node(i,j,dV_ZERO_P1)   = (-1.0/(2*del))*(grid.Node(i,j+1,dV_ZERO_M1)  - grid.Node(i,j,dV_ZERO_P1))+ grid.Node(i,j+1,dV_ZERO_M1);
        grid.Node(i,j,dV_P1_P1)     = (-1.0/(2*del))*(grid.Node(i,j+1,dV_P1_M1)    - grid.Node(i,j,dV_P1_P1))  + grid.Node(i,j+1,dV_P1_M1);
        grid.Node(i,j,dV_M1_P1)     = (-1.0/(2*del))*(grid.Node(i,j+1,dV_M1_M1)    - grid.Node(i,j,dV_M1_P1))  + grid.Node(i,j+1,dV_M1_M1);
    }

    std::cout<<grid.Node(50,grid.n_y_node - grid.noghost,dV_M1_M1)<<" x "<<std::endl;
} 



//note //////===============================================2ND ORDER================================///////////////////
//note //////===============================================2ND ORDER================================///////////////////
//note //////===============================================2ND ORDER================================///////////////////
//note //////===============================================2ND ORDER================================///////////////////
template<typename T>
void Diff_2nd_order_non_align(Grid_N_C_2D<T> &grid,lbmD2Q9<T> &lbd2q9,double u0,double del){ 

    double nx,y1,w1,w2; 
    double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
    int j,ic= 0 ,jc = 0;u0= u0;
    double feq[9] = {0};

    y1 = grid.Node(50,0,dV_P1_P1);
    // y1 = grid.Node(50,0,dV_M1_M1);

    w1 = ((1)/(2*del-1))*((1-2*del)/(2*del-1-2*del)); //first one
    w2 = ((1-(2*del-1))/(2*del - (2*del-1)))*((1)/(2*del));  ////for calculating 0.8
    // std::cout<<std::setprecision(15);
    

    // // std::cout<< "y1 = x=-0.2: "<<grid.Node(50,0,dV_ZERO_P1)<<std::endl;                                                                     ////print
    // std::cout<< "y1 = x=-0.2: "<<y1<<std::endl;                                                                     ////print


    for(int i = 2; i< grid.n_x_node-2;i++){
    j = 0;
    for(int dv = 0; dv<9;dv ++) {
        grid.Node(i,j,dv) = grid.Node(i,j,dv)*((1)/(2*del-1))*((1-2*del)/(2*del-1-2*del));
    }

    j = grid.n_y_node -1;
    for(int dv =0 ; dv<9; dv++){
        grid.Node(i,j,dv) = grid.Node(i,j,dv)*((1)/(2*del-1))*((1-2*del)/(2*del-1-2*del));
    }
    }

    /////===============top============///////
    nx = 0;  ny = -1;
    ux = u0; uy = 0;
    j = grid.n_y_node - grid.noghost -1;
    get_equi(feq, lbd2q9, ux, uy, rho);
     G = 1/(feq[dV_ZERO_M1] + feq[dV_M1_M1] + feq[dV_P1_M1]);

    int dv = dV_ZERO_P1;
    for(int i = 2; i< grid.n_x_node -2;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G*((1-(2*del-1))/(2*del - (2*del-1)))*((1)/(2*del));
                }
            }
    }

    dv = dV_M1_P1;
    for(int i = 3; i< grid.n_x_node -1;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G*((1-(2*del-1))/(2*del - (2*del-1)))*((1)/(2*del));
                }
            }
    }

    dv = dV_P1_P1;
    for(int i = 1; i< grid.n_x_node -3;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G*((1-(2*del-1))/(2*del - (2*del-1)))*((1)/(2*del));
                }
            }
    }


    ///////////////////////=========bottom============////////////////
    nx = 0;  ny = 1;
    ux = 0.0; uy = 0;
    j = 1;
    get_equi(feq, lbd2q9, ux, uy, rho);
     G = 1/(feq[dV_ZERO_P1] + feq[dV_M1_P1] + feq[dV_P1_P1]);

    dv = dV_ZERO_M1;
    for(int i = 2; i< grid.n_x_node -2;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G*((1-(2*del-1))/(2*del - (2*del-1)))*((1)/(2*del));
                }
            }
        }



    dv = dV_M1_M1;
    for(int i = 3; i< grid.n_x_node -1;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G*((1-(2*del-1))/(2*del - (2*del-1)))*((1)/(2*del));
                }
            }
    }

    
    dv = dV_P1_M1;
    for(int i = 1; i< grid.n_x_node -3;i++){
        ic = i + (int)lbd2q9.Cx[dv]; jc = j + (int)lbd2q9.Cy[dv]; 
        for(int dvc = 0; dvc <9; dvc++){
                if(lbd2q9.Cx[dvc]*nx + lbd2q9.Cy[dvc]*ny > 0){
                grid.Node(ic,jc,dvc) += feq[dvc]*grid.Node(i,j,dv) *G*((1-(2*del-1))/(2*del - (2*del-1)))*((1)/(2*del));
                }
            }
    }       

        // std::cout<< "y2 = x= 0.8: "<<(grid.Node(50,0,dV_P1_P1) - (w1*y1))/w2<<std::endl;                                                              //pr
        // // std::cout<< "y2 = x= 0.8: "<<(grid.Node(50,grid.n_y_node - grid.noghost,dV_M1_M1) - (w1*y1))/w2<<std::endl;                                                              //pr
 


        // std::cout<< "y2 = x = 0: "<<grid.Node(50,1,dV_P1_M1) <<std::endl;                                                                             //pr
        // // std::cout<< "y2 = x = 0: "<<grid.Node(50,grid.n_y_node - grid.noghost-1,dV_M1_P1) <<std::endl;                                                                             //pr



    ///////=============top
    j = grid.n_y_node - grid.noghost;
    for(int i = 1; i< grid.n_x_node -grid.noghost-1;i++){

        grid.Node(i,j,dV_ZERO_M1) = grid.Node(i,j,dV_ZERO_M1)+ grid.Node(i,j-1,dV_ZERO_P1)*((1-(2*del-1))/(0 - (2*del-1)))*((1-2*del)/(0-2*del));
        grid.Node(i,j,dV_P1_M1)   = grid.Node(i,j,dV_P1_M1)  + grid.Node(i,j-1,dV_P1_P1)  *((1-(2*del-1))/(0 - (2*del-1)))*((1-2*del)/(0-2*del));
        grid.Node(i,j,dV_M1_M1)   = grid.Node(i,j,dV_M1_M1)  + grid.Node(i,j-1,dV_M1_P1)  *((1-(2*del-1))/(0 - (2*del-1)))*((1-2*del)/(0-2*del));

        // for(int dv=0; dv<grid.d_v;dv++){
        //     grid.Node(i,j,dv) = 0.8*(grid.Node((int)lbd2q9.Cx[oppdV[dv]],(int)lbd2q9.Cy[oppdV[dv],oppdV[dv]])-
        //                              grid.Node(i,j,dv)) + grid.Node((int)lbd2q9.Cx[oppdV[dv]],(int)lbd2q9.Cy[oppdV[dv],oppdV[dv]]);
        // }
    
    }


  

    /////================bottom
    j = 0;
    for(int i = 1; i< grid.n_x_node -grid.noghost-1;i++){
        grid.Node(i,j,dV_ZERO_P1)   = grid.Node(i,j,dV_ZERO_P1) + grid.Node(i,j+1,dV_ZERO_M1)*((1-(2*del-1))/(0 - (2*del-1)))*((1-2*del)/(0-2*del));
        grid.Node(i,j,dV_P1_P1)     = grid.Node(i,j,dV_P1_P1)   + grid.Node(i,j+1,dV_P1_M1)  *((1-(2*del-1))/(0 - (2*del-1)))*((1-2*del)/(0-2*del));
        grid.Node(i,j,dV_M1_P1)     = grid.Node(i,j,dV_M1_P1)   + grid.Node(i,j+1,dV_M1_M1)  *((1-(2*del-1))/(0 - (2*del-1)))*((1-2*del)/(0-2*del));
    }
// 
    // std::cout<<grid.Node(50,0,dV_P1_P1)<<"x "<<std::endl;                                                                                               //pr
    // std::cout<<grid.Node(50,grid.n_y_node - grid.noghost,dV_M1_M1)<<"x "<<std::endl; 
} 



template<typename T>
void G_calc(Grid_N_C_2D<T> &grid,lbmD2Q9<T> &lbd2q9){ //indirectly contain the local normal of aerofoil

    double feq[9] = {0};
    double ux = 0; double uy = 0; double rho = 1.0;
    get_equi(feq,lbd2q9, ux, uy, rho);

for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
    for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost); j++){ 
        if(grid.marker_node(i,j) == 0){ //enter for only the fluid node
        for(int dv = 0; dv<9; dv++){


            if(grid.marker_node(i + (int)lbd2q9.Cx[dv], j+ (int) lbd2q9.Cy[dv]) == 1){ // if the landing site is solid

                grid.g_inv(i + (int)lbd2q9.Cx[dv], j+ (int) lbd2q9.Cy[dv]) += feq[dv];

            }


        }

        }
}}

///////////==========printing the population addition======================/
// for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
//     for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost); j++){ 

//         if(grid.marker_node(i,j) == 0){ //enter for only the fluid node
//         for(int dv = 0; dv<9; dv++){


//             if(grid.marker_node(i + (int)lbd2q9.Cx[dv], j+ (int) lbd2q9.Cy[dv]) == 1){ // if the landing site is solid

//                 std::cout<<i + (int)lbd2q9.Cx[dv]<<","<<j+ (int) lbd2q9.Cy[dv]<<","<<grid.g_inv(i + (int)lbd2q9.Cx[dv], j+ (int) lbd2q9.Cy[dv])<<","<<dv<<std::endl;

//             }
//         }}
//     }}



}


template<typename T>
void Diffuse_B_Solid(Grid_N_C_2D<T> &grid,lbmD2Q9<T> &lbd2q9,double u0){ 

int i_b,j_b;

double feq[9] = {0};
get_equi(feq,lbd2q9,u0,0.0,1.0);

for(int i = 270 ; i < 350;i++){
    for(int j = 190 ; j <220; j++){
        if(grid.marker_node(i,j) ==1){

            for(int dv = 0; dv<9; dv++){
                grid.Node(i,j,dv) = 0.0;
            }
        }
    }}



for(int i = 270 ; i < 350;i++){
    for(int j = 190 ; j <220; j++){ 
        if(grid.marker_node(i,j) == 0){ //enter for only the fluid node
        for(int dv = 0; dv<9; dv++){

            i_b = i + (int)lbd2q9.Cx[dv]; j_b = j+ (int) lbd2q9.Cy[dv];

            if(grid.marker_node(i_b, j_b ) == 1){ // if the landing site is solid

                for(int dvc = 0; dvc<9; dvc++){
                    

                    if(grid.marker_node(i_b+(int)lbd2q9.Cx[dvc],j_b+ (int) lbd2q9.Cy[dvc]) ==0){

                        grid.Node(i_b,j_b,dvc) +=   feq[dvc] * grid.Node(i,j,dv) /grid.g_inv(i_b,j_b);


                   }
            } }  } }
}}

}


template<typename T, typename T1>
void slip_wall(Grid_N_C_2D<T> &lbgrid,lbmD2Q9<T1> &lbD3Q9){







}





template<typename T, typename T1>
void slip_wall(Grid_N_C_2D<T> &lbgrid,lbmD2Q9<T1> &lbD3Q9, double ux, double uy,double rho){

    int j ;
    for(int i = 1; i< lbgrid.n_x_node - lbgrid.noghost -1 ;i++){
        
        j = lbgrid.n_y_node - lbgrid.noghost;
        
        lbgrid.Node(i,j,dV_ZERO_M1) = lbgrid.Node(i,j -1, dV_ZERO_P1);
        lbgrid.Node(i,j,dV_P1_M1)   = lbgrid.Node(i,j-1 , dV_P1_P1);
        lbgrid.Node(i,j,dV_M1_M1)   = lbgrid.Node(i,j -1, dV_M1_P1);
        
        // lbgrid.Node(i,j,dV_P1_M1)   = lbgrid.Node(i+1,j-1 , dV_M1_P1);
        // lbgrid.Node(i,j,dV_M1_M1)   = lbgrid.Node(i-1,j-1 , dV_P1_P1);

        j = 0;

        lbgrid.Node(i,j,dV_ZERO_P1) = lbgrid.Node(i,j+1, dV_ZERO_M1);
        lbgrid.Node(i,j,dV_P1_P1)   = lbgrid.Node(i,j+1, dV_P1_M1 );
        lbgrid.Node(i,j,dV_M1_P1)   = lbgrid.Node(i,j+1, dV_M1_M1);
    }

}


template<typename T, typename T1>
void inlet(Grid_N_C_2D<T> &lbgrid,lbmD2Q9<T1> &lbD3Q9, double ux, double uy,double rho){

double feq_Node[9] = {0};

get_equi(feq_Node,lbD3Q9,ux,uy,rho);


    for(int i = lbgrid.noghost; i <= lbgrid.noghost  ;i++){
        for(int j = lbgrid.noghost; j<lbgrid.n_y_node - lbgrid.noghost;j++){

            for(int dv = 0; dv <9; dv++){
                lbgrid.Node(i,j,dv) = feq_Node[dv];
            }
                        

            
        }
        }
}






template<typename T, typename T1>
void outlet(Grid_N_C_2D<T> &lbgrid,lbmD2Q9<T1> &lbD3Q9, double ux, double uy,double rho){


double feq_Node[9] = {0};


get_equi(feq_Node,lbD3Q9,ux,uy,rho);
// get_equi(feq_Cell,lbD3Q9,ux,uy,rho);

    for(int i = lbgrid.n_x_node - lbgrid.noghost -4; i < lbgrid.n_x_node - lbgrid.noghost  ;i++){
        for(int j = lbgrid.noghost; j<lbgrid.n_y_node - lbgrid.noghost;j++){

            
            // for(int dv = 0; dv<9; dv++){
            //     lbgrid.Node(i,j,dv) = feq_Node[dv];
            // }

            for(int dv = 0; dv< 9; dv++){
                lbgrid.Node(i,j,dv) = lbgrid.Node(i-1, j,dv);

            }


        }
        }
}






/////////////==========printing the population addition======================/
// for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
//     for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost); j++){ 

//         if(grid.marker_node(i,j) == 0){ //enter for only the fluid node
//         for(int dv = 0; dv<9; dv++){


//             if(grid.marker_node(i + (int)lbd2q9.Cx[dv], j+ (int) lbd2q9.Cy[dv]) == 1){ // if the landing site is solid

//                 std::cout<<i + (int)lbd2q9.Cx[dv]<<","<<j+ (int) lbd2q9.Cy[dv]<<","<<grid.g_inv(i + (int)lbd2q9.Cx[dv], j+ (int) lbd2q9.Cy[dv])<<","<<dv<<std::endl;

//             }
//         }}
//     }}















// template<typename T>
// void advection(Grid_N_C_2D<T> &grid){
    
//         double a,b;
//         for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
//             a = grid.Cell(i,grid.noghost -1,dV_MH1_PH1);
//             // std::cout<<a<<" "<<i<<std::endl;
//         for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost); j++){ 
            
//             grid.Node(i,j,dV_ZERO_M1) = grid.Node(i  ,j+1, dV_ZERO_M1 ); //for SC1 -Node
//             grid.Node(i,j,dV_M1_ZERO) = grid.Node(i+1,j  , dV_M1_ZERO );

//             grid.Cell(i,j,dV_ZERO_M1) = grid.Cell(i  ,j+1, dV_ZERO_M1 ); //for SC1 -cell
//             grid.Cell(i,j,dV_M1_ZERO) = grid.Cell(i+1,j  , dV_M1_ZERO );

//             ///half velocities
//             //std::cout<<grid.Node(i,j,dV_MH1_MH1)<<"before "<<i<<j <<std::endl;
            
//             grid.Node(i,j,dV_MH1_MH1) = grid.Cell(i  ,j  , dV_MH1_MH1);
//             grid.Node(i,j,dV_MH1_PH1) = a;     
//             a = grid.Cell(i,j,dV_MH1_PH1);       

//             grid.Cell(i,j,dV_MH1_MH1) = grid.Node(i+1,j+1, dV_MH1_MH1);                        
//             grid.Cell(i,j,dV_MH1_PH1) = grid.Node(i+1,j  , dV_MH1_PH1);

//             //3 velocities

//             grid.Node(i,j,dV_ZERO_M3) = grid.Node(i  ,j+3, dV_ZERO_M3);
//             grid.Node(i,j,dV_M3_ZERO) = grid.Node(i+3,j  , dV_M3_ZERO);

//             grid.Cell(i,j,dV_ZERO_M3) = grid.Cell(i  ,j+3, dV_ZERO_M3);
//             grid.Cell(i,j,dV_M3_ZERO) = grid.Cell(i+3,j  , dV_M3_ZERO);
//             //std::cout<<grid.Node(i,j,dV_MH1_MH1) <<" "<<i<<j<<std::endl;
//         } 
//     }

//         for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){
//             b = grid.Node(i,grid.n_y_node - grid.noghost,dV_PH1_MH1);
//         for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){

//             //SC1

//             grid.Node(i,j, dV_ZERO_P1) = grid.Node(i  ,j-1, dV_ZERO_P1);
//             grid.Node(i,j, dV_P1_ZERO) = grid.Node(i-1,j  , dV_P1_ZERO);

//             grid.Cell(i,j,dV_ZERO_P1) = grid.Cell(i  ,j-1, dV_ZERO_P1 );
//             grid.Cell(i,j,dV_P1_ZERO) = grid.Cell(i-1,j  , dV_P1_ZERO );
      
//             ///half velocities
//            // std::cout<<grid.Node(i,j,dV_PH1_PH1)<<"before "<<i<<j <<std::endl;


//             grid.Cell(i,j,dV_PH1_MH1) = b;
//             grid.Cell(i,j,dV_PH1_PH1) = grid.Node(i  ,j  , dV_PH1_PH1);
//             b=grid.Node(i,j,dV_PH1_MH1);

//             grid.Node(i,j,dV_PH1_MH1) = grid.Cell(i-1,j  , dV_PH1_MH1);
            
//             grid.Node(i,j,dV_PH1_PH1) = grid.Cell(i-1,j-1, dV_PH1_PH1);



//             //std::cout<<grid.Cell(i,j,dV_PH1_PH1) <<" "<<i<<j<<std::endl;

//             //SC3
//             grid.Node(i,j, dV_ZERO_P3) = grid.Node(i  ,j-3, dV_ZERO_P3);
//             grid.Node(i,j, dV_P3_ZERO) = grid.Node(i-3,j  , dV_P3_ZERO);


//             grid.Cell(i,j,dV_ZERO_P3) = grid.Cell(i  ,j-3, dV_ZERO_P3 );
//             grid.Cell(i,j,dV_P3_ZERO) = grid.Cell(i-3,j  , dV_P3_ZERO );

        
//         }}

// }

// template<typename T>
// void advection1(Grid_N_C_2D<T> &grid){


//         #pragma omp parallel 
//         {


//         #pragma omp for schedule(static)
//         for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
//         for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 

//             grid.Node(i,j,dV_ZERO_M1) = grid.Node(i  ,j+1, dV_ZERO_M1 ); //for SC1 -Node
//             grid.Node(i,j,dV_M1_ZERO) = grid.Node(i+1,j  , dV_M1_ZERO );

//             grid.Node(i,j,dV_ZERO_M3) = grid.Node(i  ,j+3, dV_ZERO_M3);
//             grid.Node(i,j,dV_M3_ZERO) = grid.Node(i+3,j  , dV_M3_ZERO);

//             grid.Node(i,j  ,dV_MH1_MH1) = grid.Cell(i  ,j   , dV_MH1_MH1);
//             grid.Node(i,j  ,dV_MH1_PH1) = grid.Cell(i  ,j -1, dV_MH1_PH1);

//         }
//         }


//         #pragma omp for schedule(static)
//         for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
//         for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 

//             grid.Cell(i,j,dV_ZERO_M1) = grid.Cell(i  ,j+1, dV_ZERO_M1 ); //for SC1 -cell
//             grid.Cell(i,j,dV_M1_ZERO) = grid.Cell(i+1,j  , dV_M1_ZERO );

//             grid.Cell(i,j,dV_ZERO_M3) = grid.Cell(i  ,j+3, dV_ZERO_M3);
//             grid.Cell(i,j,dV_M3_ZERO) = grid.Cell(i+3,j  , dV_M3_ZERO);

//             grid.Cell(i,j,dV_MH1_MH1) = grid.Node(i+1,j+1, dV_MH1_MH1);                        
//             grid.Cell(i,j,dV_MH1_PH1) = grid.Node(i+1,j  , dV_MH1_PH1);
            

//         }
//         }

//         #pragma omp for schedule(static)
//         for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){                                                                                        //  #pragma omp parallel for shared(b)
//         for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){
            
//             grid.Cell(i,j,dV_ZERO_P1) = grid.Cell(i  ,j-1, dV_ZERO_P1 );
//             grid.Cell(i,j,dV_P1_ZERO) = grid.Cell(i-1,j  , dV_P1_ZERO );

//             grid.Cell(i,j,dV_ZERO_P3) = grid.Cell(i  ,j-3, dV_ZERO_P3 );
//             grid.Cell(i,j,dV_P3_ZERO) = grid.Cell(i-3,j  , dV_P3_ZERO );

//             grid.Cell(i,j,dV_PH1_MH1) = grid.Node(i,j+1,dV_PH1_MH1);
//             grid.Cell(i,j,dV_PH1_PH1) = grid.Node(i,j  ,dV_PH1_PH1);

//         }}


//         #pragma omp for schedule(static)
//         for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){                                                                                        //  #pragma omp parallel for shared(b)
//         for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){
            
//             grid.Node(i,j, dV_ZERO_P1) = grid.Node(i  ,j-1, dV_ZERO_P1);
//             grid.Node(i,j, dV_P1_ZERO) = grid.Node(i-1,j  , dV_P1_ZERO);

//             grid.Node(i,j, dV_ZERO_P3) = grid.Node(i  ,j-3, dV_ZERO_P3);
//             grid.Node(i,j, dV_P3_ZERO) = grid.Node(i-3,j  , dV_P3_ZERO);

//             grid.Node(i,j,dV_PH1_MH1) = grid.Cell(i-1,j  ,dV_PH1_MH1);
//             grid.Node(i,j,dV_PH1_PH1) = grid.Cell(i-1,j-1,dV_PH1_PH1);

//         }}


//         }

// }



// template<typename T>
// void periodic(Grid_N_C_2D<T> &grid){

  
//     // //---------------------------------top bottom---------------------------------------//can it be applied separately?????
//     // for(int i = 0; i<grid.n_x_node; i++){
//     //     ///-------top to bottom----///
//     //     //SC1
//     //     grid.Node(i,grid.noghost -1,dV_ZERO_P1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1); ///copying to first adjacent ghost node from last node physical domain
//     //     grid.Cell(i,grid.noghost -1,dV_ZERO_P1) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1); 

//     //     //SC3
//     //     grid.Node(i,0,dV_ZERO_P3) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P3);         ///*****not generalized
//     //     grid.Cell(i,0,dV_ZERO_P3) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P3);

//     //     grid.Node(i,1,dV_ZERO_P3) = grid.Node(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P3);
//     //     grid.Cell(i,1,dV_ZERO_P3) = grid.Cell(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P3);

//     //     grid.Node(i,2,dV_ZERO_P3) = grid.Node(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_P3);
//     //     grid.Cell(i,2,dV_ZERO_P3) = grid.Cell(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_P3);

//     //     //SC-1/2
//     //     grid.Cell(i,grid.noghost -1,dV_PH1_PH1) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_PH1_PH1); //top cell to bottom
//     //     grid.Cell(i,grid.noghost -1,dV_MH1_PH1) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_MH1_PH1);

//     //     ///-------bottom to top------///
//     //     //SC -1
//     //     grid.Node(i,grid.n_y_node - grid.noghost,dV_ZERO_M1) = grid.Node(i, grid.noghost,dV_ZERO_M1);
//     //     grid.Cell(i,grid.n_y_node - grid.noghost,dV_ZERO_M1) = grid.Cell(i, grid.noghost,dV_ZERO_M1);

//     //     //SC3
//     //     grid.Node(i,grid.n_y_node - (grid.noghost),dV_ZERO_M3)   = grid.Node(i,grid.noghost,dV_ZERO_M3);
//     //     grid.Cell(i,grid.n_y_node - (grid.noghost),dV_ZERO_M3)   = grid.Cell(i,grid.noghost,dV_ZERO_M3);

//     //     grid.Node(i,grid.n_y_node - (grid.noghost-1),dV_ZERO_M3) = grid.Node(i,grid.noghost+1,dV_ZERO_M3);
//     //     grid.Cell(i,grid.n_y_node - (grid.noghost-1),dV_ZERO_M3) = grid.Cell(i,grid.noghost+1,dV_ZERO_M3);
        
//     //     grid.Node(i,grid.n_y_node - (grid.noghost-2),dV_ZERO_M3) = grid.Node(i,grid.noghost+2,dV_ZERO_M3);
//     //     grid.Cell(i,grid.n_y_node - (grid.noghost-2),dV_ZERO_M3) = grid.Cell(i,grid.noghost+2,dV_ZERO_M3);


//     //     //SC-1/2
//     //     grid.Node(i,grid.n_y_node - (grid.noghost),dV_PH1_MH1) = grid.Node(i,grid.noghost,dV_PH1_MH1);
//     //     grid.Node(i,grid.n_y_node - (grid.noghost),dV_MH1_MH1) = grid.Node(i,grid.noghost,dV_MH1_MH1);

//     // }

//     //--------------------------------------------------left right-----------------------------------------------//

//     for(int j =0; j<grid.n_y_node  ;j++){
//         //--------right to left-------//
//         grid.Node(grid.noghost -1,j,dV_P1_ZERO) = grid.Node( grid.n_x_node - (grid.noghost +1),j,dV_P1_ZERO); ///copying to first adjacent ghost node from last node physical domain
//         grid.Cell(grid.noghost -1,j,dV_P1_ZERO) = grid.Cell( grid.n_x_node - (grid.noghost +1),j,dV_P1_ZERO);

//         //SC3
//         grid.Node(0,j,dV_P3_ZERO) = grid.Node(grid.n_x_cell - (grid.noghost +3),j,dV_P3_ZERO);         
//         grid.Cell(0,j,dV_P3_ZERO) = grid.Cell(grid.n_x_cell - (grid.noghost +3),j,dV_P3_ZERO);

//         grid.Node(1,j,dV_P3_ZERO) = grid.Node(grid.n_x_cell - (grid.noghost +2),j,dV_P3_ZERO);
//         grid.Cell(1,j,dV_P3_ZERO) = grid.Cell(grid.n_x_cell - (grid.noghost +2),j,dV_P3_ZERO);

//         grid.Node(2,j,dV_P3_ZERO) = grid.Node(grid.n_x_cell - (grid.noghost +1),j,dV_P3_ZERO);
//         grid.Cell(2,j,dV_P3_ZERO) = grid.Cell(grid.n_x_cell - (grid.noghost +1),j,dV_P3_ZERO);

//         //BCC-1/2
//         grid.Cell(grid.noghost -1,j,dV_PH1_PH1) = grid.Cell(grid.n_x_cell - (grid.noghost +1),j,dV_PH1_PH1); 
//         grid.Cell(grid.noghost -1,j,dV_PH1_MH1) = grid.Cell(grid.n_x_cell - (grid.noghost +1),j,dV_PH1_MH1);

//         ///-------left to right-----///
//         //SC -1
//         grid.Node(grid.n_x_cell - grid.noghost,j,dV_M1_ZERO) = grid.Node(grid.noghost,j,dV_M1_ZERO);
//         grid.Cell(grid.n_x_cell - grid.noghost,j,dV_M1_ZERO) = grid.Cell(grid.noghost,j,dV_M1_ZERO);

//         //SC3
//         grid.Node(grid.n_x_cell - (grid.noghost)  ,j,dV_M3_ZERO)   = grid.Node(grid.noghost  ,j,dV_M3_ZERO);
//         grid.Cell(grid.n_x_cell - (grid.noghost)  ,j,dV_M3_ZERO)   = grid.Cell(grid.noghost  ,j,dV_M3_ZERO);

//         grid.Node(grid.n_x_cell - (grid.noghost-1),j,dV_M3_ZERO)   = grid.Node(grid.noghost+1,j,dV_M3_ZERO);
//         grid.Cell(grid.n_x_cell - (grid.noghost-1),j,dV_M3_ZERO)   = grid.Cell(grid.noghost+1,j,dV_M3_ZERO);

//         grid.Node(grid.n_x_cell - (grid.noghost-2),j,dV_M3_ZERO)   = grid.Node(grid.noghost+2,j,dV_M3_ZERO);
//         grid.Cell(grid.n_x_cell - (grid.noghost-2),j,dV_M3_ZERO)   = grid.Cell(grid.noghost+2,j,dV_M3_ZERO);


//         //BCC-1/2
//         grid.Node(grid.n_x_cell - (grid.noghost),j,dV_MH1_PH1) = grid.Node(grid.noghost,j,dV_MH1_PH1);
//         grid.Node(grid.n_x_cell - (grid.noghost),j,dV_MH1_MH1) = grid.Node(grid.noghost,j,dV_MH1_MH1);

//     }

// }





// template<typename T>
// void bounce_back_prep(Grid_N_C_2D<T> &grid){
// //________________________________________________bounceback_________________________________________________________//
   
//     for(int i =1; i<grid.n_x_node-1 ; i++){
        
//         //------SC1
//         //top________(grid.ny_node - grid.noghost) --->> 1st ghost node
//         grid.Node(i, grid.n_y_node - (grid.noghost),dV_ZERO_P1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1);
//         //taking the last cell as the boundary
//         grid.Cell(i, grid.n_y_node - (grid.noghost),dV_ZERO_M1) = grid.Cell(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P1);

//         //bottom
//         grid.Node(i,grid.noghost -1,dV_ZERO_M1) = grid.Node(i,grid.noghost,dV_ZERO_M1);
//         grid.Cell(i,grid.noghost -1,dV_ZERO_M1) = grid.Cell(i,grid.noghost,dV_ZERO_M1);


//         //----BCC1/2
//         //top
//         grid.Node(i, grid.n_y_node - (grid.noghost),dV_PH1_PH1) = grid.Cell(i-1, grid.n_y_node - (grid.noghost+1),dV_PH1_PH1);
//         grid.Node(i, grid.n_y_node - (grid.noghost),dV_MH1_PH1) = grid.Cell(i  , grid.n_y_node - (grid.noghost+1),dV_MH1_PH1);
//         //bottom
//         grid.Cell(i,grid.noghost -1,dV_MH1_MH1) = grid.Node(i+1,grid.noghost,dV_MH1_MH1);
//         grid.Cell(i,grid.noghost -1,dV_PH1_MH1) = grid.Node(i,grid.noghost,dV_PH1_MH1);

//         //---SC3
//         //top
//         grid.Node(i,  grid.n_y_node - (grid.noghost)  ,dV_ZERO_P3) =  grid.Node(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_P3);
//         grid.Node(i,  grid.n_y_node - (grid.noghost-1),dV_ZERO_P3) =  grid.Node(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P3);
//         grid.Node(i,  grid.n_y_node - (grid.noghost-2),dV_ZERO_P3) =  grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P3);
                   
//         grid.Cell(i,  grid.n_y_cell - (grid.noghost)  ,dV_ZERO_P3) =  grid.Cell(i, grid.n_y_cell - (grid.noghost +3),dV_ZERO_P3);
//         grid.Cell(i,  grid.n_y_cell - (grid.noghost-1),dV_ZERO_P3) =  grid.Cell(i, grid.n_y_cell - (grid.noghost +2),dV_ZERO_P3);
//         grid.Cell(i,  grid.n_y_cell - (grid.noghost-2),dV_ZERO_P3) =  grid.Cell(i, grid.n_y_cell - (grid.noghost +1),dV_ZERO_P3);

//         //bottom
//         grid.Node(i, grid.noghost - 3,dV_ZERO_M3) =  grid.Node(i, grid.noghost   ,dV_ZERO_M3);
//         grid.Node(i, grid.noghost - 2,dV_ZERO_M3) =  grid.Node(i, grid.noghost +1,dV_ZERO_M3);
//         grid.Node(i, grid.noghost - 1,dV_ZERO_M3) =  grid.Node(i, grid.noghost +2,dV_ZERO_M3);

//         grid.Cell(i, grid.noghost - 3,dV_ZERO_M3) =  grid.Cell(i, grid.noghost   ,dV_ZERO_M3);
//         grid.Cell(i, grid.noghost - 2,dV_ZERO_M3) =  grid.Cell(i, grid.noghost +1,dV_ZERO_M3);
//         grid.Cell(i, grid.noghost - 1,dV_ZERO_M3) =  grid.Cell(i, grid.noghost +2,dV_ZERO_M3);        
 
//     }
    

// }



// template<typename T>
// void bounce_back_wall(Grid_N_C_2D<T> &grid){




//     for(int i =1; i<grid.n_x_node-1 ; i++){

//     //SC1
//     //top
//     grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_M1) = grid.Node(i, grid.n_y_node - (grid.noghost),dV_ZERO_P1);
//     grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_M1) = grid.Cell(i, grid.n_y_node - (grid.noghost),dV_ZERO_P1);
//     //bottom
//     grid.Node(i,grid.noghost,dV_ZERO_P1) = grid.Node(i,grid.noghost -1,dV_ZERO_M1);
//     grid.Cell(i,grid.noghost,dV_ZERO_P1) = grid.Cell(i,grid.noghost -1,dV_ZERO_M1);
//     //BCC1/2
//     //top
//     grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_MH1_MH1) = grid.Node(i+1, grid.n_y_node - (grid.noghost),dV_PH1_PH1);
//     grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_PH1_MH1) = grid.Node(i  , grid.n_y_node - (grid.noghost),dV_MH1_PH1);
//     //bottom
//     grid.Node(i,grid.noghost,dV_PH1_PH1) = grid.Cell(i-1,grid.noghost -1,dV_MH1_MH1);
//     grid.Node(i,grid.noghost,dV_MH1_PH1) = grid.Cell(i,grid.noghost -1,dV_PH1_MH1);

//     ///SC3
//     //top
//     grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost)  ,dV_ZERO_P3);
//     grid.Node(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost-1),dV_ZERO_P3);
//     grid.Node(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost-2),dV_ZERO_P3);

//     grid.Cell(i, grid.n_y_cell - (grid.noghost +1),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost)  ,dV_ZERO_P3);
//     grid.Cell(i, grid.n_y_cell - (grid.noghost +2),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost-1),dV_ZERO_P3);
//     grid.Cell(i, grid.n_y_cell - (grid.noghost +3),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost-2),dV_ZERO_P3);

//     //bottom
//     grid.Node(i, grid.noghost  ,dV_ZERO_P3) = grid.Node(i, grid.noghost - 1,dV_ZERO_M3);
//     grid.Node(i, grid.noghost+1,dV_ZERO_P3) = grid.Node(i, grid.noghost - 2,dV_ZERO_M3);
//     grid.Node(i, grid.noghost+2,dV_ZERO_P3) = grid.Node(i, grid.noghost - 3,dV_ZERO_M3);
    
//     grid.Cell(i, grid.noghost  ,dV_ZERO_P3) = grid.Cell(i, grid.noghost - 1,dV_ZERO_M3);
//     grid.Cell(i, grid.noghost+1,dV_ZERO_P3) = grid.Cell(i, grid.noghost - 2,dV_ZERO_M3);
//     grid.Cell(i, grid.noghost+2,dV_ZERO_P3) = grid.Cell(i, grid.noghost - 3,dV_ZERO_M3);

//     }
// }


// //////////////--------------------------------------------/working bounce back--------------------------------------------------------//////////
// template<typename T>
// void bounce_back(Grid_N_C_2D<T> &grid){
    

//     for ( int i = 1 ;i< grid.n_x_node -1; i++){

//     //TOP
//     //SC1
//      grid.Node(i, grid.n_y_node - (grid.noghost),dV_ZERO_M1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1);
//      grid.Cell(i, grid.n_y_cell - (grid.noghost),dV_ZERO_M1) = grid.Cell(i, grid.n_y_cell - (grid.noghost+2),dV_ZERO_P1);


//     //SC3
//     grid.Node(i, grid.n_y_node - (grid.noghost  ),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P3);
//     grid.Node(i, grid.n_y_node - (grid.noghost-1),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P3);
//     grid.Node(i, grid.n_y_node - (grid.noghost-2),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_P3);

//     grid.Cell(i, grid.n_y_cell - (grid.noghost  ),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost +2),dV_ZERO_P3);
//     grid.Cell(i, grid.n_y_cell - (grid.noghost-1),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost +3),dV_ZERO_P3);
//     grid.Cell(i, grid.n_y_cell - (grid.noghost-2),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost +4),dV_ZERO_P3);


//     ///SC1/2
// // //noslip bounce back
//     grid.Node(i + 1, grid.n_y_node - (grid.noghost), dV_MH1_MH1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_PH1_PH1);     //////top right and top left corners are bouncing back in to the domain
//     grid.Node(i - 1, grid.n_y_node - (grid.noghost), dV_PH1_MH1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_MH1_PH1);

//         //slip bounceback
//     // grid.Node(i , grid.n_y_node - (grid.noghost), dV_PH1_MH1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_PH1_PH1);     
//     // grid.Node(i , grid.n_y_node - (grid.noghost), dV_MH1_MH1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_MH1_PH1);



//     //Bottom
//     //SC1
   
//     grid.Node(i, grid.noghost - 1,dV_ZERO_P1) = grid.Node(i,grid.noghost +1, dV_ZERO_M1);
//     grid.Cell(i, grid.noghost - 1,dV_ZERO_P1) = grid.Cell(i,grid.noghost   , dV_ZERO_M1);



//     //SC3
//     grid.Node(i, grid.noghost - 1,dV_ZERO_P3) =  grid.Node(i, grid.noghost+1, dV_ZERO_M3);
//     grid.Node(i, grid.noghost - 2,dV_ZERO_P3) =  grid.Node(i, grid.noghost+2, dV_ZERO_M3);
//     grid.Node(i, grid.noghost - 3,dV_ZERO_P3) =  grid.Node(i, grid.noghost+3, dV_ZERO_M3);

//     grid.Cell(i, grid.noghost - 1,dV_ZERO_P3) =  grid.Cell(i, grid.noghost   , dV_ZERO_M3);
//     grid.Cell(i, grid.noghost - 2,dV_ZERO_P3) =  grid.Cell(i, grid.noghost +1, dV_ZERO_M3);
//     grid.Cell(i, grid.noghost - 3,dV_ZERO_P3) =  grid.Cell(i, grid.noghost +2, dV_ZERO_M3);

//     //SC1/2
//     // //noslip bounceback
//     grid.Cell(i + 1, grid.noghost-1, dV_MH1_PH1) =  grid.Cell(i, grid.noghost, dV_PH1_MH1);            //bottom right and left also bounces back                    
//     grid.Cell(i - 1, grid.noghost-1, dV_PH1_PH1) =  grid.Cell(i, grid.noghost, dV_MH1_MH1);

//     // slip bounce back
//     // grid.Cell(i , grid.noghost-1, dV_PH1_PH1) =  grid.Cell(i, grid.noghost, dV_PH1_MH1);
//     // grid.Cell(i , grid.noghost-1, dV_MH1_PH1) =  grid.Cell(i, grid.noghost, dV_MH1_MH1);

// }

// }







// template<typename T, typename T1>
// void inlet(Grid_N_C_2D<T> &lbgrid,lbmD2Q13<T1> &lbD3Q13, double ux, double uy,double rho){


// double feq_Node[13] = {0},feq_Cell[13] = {0};



// get_equi(feq_Node,lbD3Q13,ux,uy,rho);
// get_equi(feq_Cell,lbD3Q13,ux,uy,rho);

//     for(int i = lbgrid.noghost; i <= lbgrid.noghost +2 ;i++){
//         for(int j = lbgrid.noghost; j<lbgrid.n_y_node - lbgrid.noghost;j++){

//             for(int dv = 0; dv <13; dv++){
//                 lbgrid.Node(i,j,dv) = feq_Node[dv];
//             }
                        
//             for(int dv = 0; dv <13; dv++){
//                 lbgrid.Cell(i,j,dv) = feq_Cell[dv];
//             }
            
//         }
//         }
// }


// template<typename T, typename T1>
// void outlet(Grid_N_C_2D<T> &lbgrid,lbmD2Q13<T1> &lbD3Q13, double ux, double uy,double rho){


// double feq_Node[13] = {0},feq_Cell[13] = {0};


// get_equi(feq_Node,lbD3Q13,ux,uy,rho);
// get_equi(feq_Cell,lbD3Q13,ux,uy,rho);

//     for(int i = lbgrid.n_x_node - lbgrid.noghost -3; i < lbgrid.n_x_node - lbgrid.noghost  ;i++){
//         for(int j = lbgrid.noghost; j<lbgrid.n_y_node - lbgrid.noghost;j++){

//             for(int dv = 0; dv <13; dv++){
//                 lbgrid.Node(i,j,dv) = feq_Node[dv];
//             }
                        
//             for(int dv = 0; dv <13; dv++){
//                 lbgrid.Cell(i,j,dv) = feq_Cell[dv];
//             }
            


//             // for(int dv = 0; dv< 13; dv++){
//             //     lbgrid.Node(i,j,dv) = lbgrid.Node(i-1, j,dv);

//             // }
//             // for(int dv = 0; dv<13; dv++){
//             //    lbgrid.Cell(i,j,dv) = lbgrid.Cell(i-1, j,dv);
//             // }


//         }
//         }
// }



;



#endif