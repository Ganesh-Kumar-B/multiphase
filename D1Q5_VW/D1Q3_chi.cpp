/*
  A Simple 1D D1Q5 code for MultiPhase flow
 */
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
/*****Max Size of Array *******/
#define MAX_SIZE  2000 
#define N_DV 5
 
enum direction
{
    ZERO, /* Alias for 0 */
    DX, /* Alias for 1 */
    DMX, /* Alias for -1 */
    D3X, /* Alias for 3 */
    DM3X, /* Alias for -3 */
};

typedef  struct nodeD1Q3
{
    double rho;
    double rhoOld;
    double pNid;
    double FNid;
    double chi;
    double muA;
    double tauM;
    double beta2;
    double betaM;
    double oneM2beta;
    double vel;
    //new foHelmholtz
    double fHelmholtz;
    //define a forcing term
    double Force;
    //define pressure
    double surface;
    double f[N_DV]; 
    double fEq[N_DV]; 
    double f0[N_DV]; 
    double Fh[N_DV]; 
    double df[N_DV];
} nodeD1Q3;

typedef  nodeD1Q3 latticeArr[MAX_SIZE]  ;
 
typedef  struct latticeD1Q3
{
    double weight[N_DV];
    double dvD1Q3[N_DV];
    int sgn[N_DV];
    double T0;
    double T0Inv;
    double c;
    double c2;
} latticeD1Q3;


typedef  struct nonIdealParam  
{
    double a;
    double b;
    double alpha1;
    double kappa;
    double rho0;
}nonIdealParam;



void getLatticeD1Q3(double c, latticeD1Q3 *myD1Q3);
void printRho(latticeArr lattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW,int nX, double beta,int iX_begin, int iX_end, int time,double c);
void getFeqPQuad(double fEq[N_DV], latticeD1Q3 myD1Q3,  double rho,  double vel);
void createBoundaryPeriodic(latticeArr lattice, int nX );
void collideWorking (latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW, int nX, double beta,double tau, int time,double c);
void advect(latticeArr lattice, int nX);
void calculateAlpha( latticeArr myLattice, double& alpha, int i, int time);
void initializePerturbPeriodic(latticeArr lattice, latticeD1Q3 myD1Q3,  int nX, double vel, double rhoMean, double ampDisturb, double periodDisturb);


main()
{
  
    latticeArr myLattice;
	latticeD1Q3 myD1Q3;
	double kappabar, TbyTc,rho0byrhoc;
	double c, rho,u, beta,  T_critical, tau, dt, dX, rho_critical;
	nonIdealParam myVDW;
	int iX,nX,time, finalTime;
	c = 1;
	getLatticeD1Q3(c, &myD1Q3);



    nX = 300;
	beta = 0.7;
	TbyTc = 0.96;
	rho0byrhoc = 0.95;
	kappabar = 0.0001;//0.0625;
	dX = 1.0/(nX-1.0);
	dt =  dX/c;
	
	tau = (1.0-beta)/beta *dt*0.5;   
	printf("\n beta=%lf kn =%lf \n",beta, tau);
	
    T_critical = myD1Q3.T0/TbyTc ;
    rho_critical = 1.0;

    myVDW.b = 1.0/(3.0*rho_critical);
    myVDW.a = myVDW.b*T_critical*27.0/8.0  ;

	myVDW.kappa = kappabar*myVDW.a*dX*dX;
	myVDW.rho0 = rho0byrhoc*rho_critical;	

    finalTime = 990000    ;

	initializePerturbPeriodic(  myLattice,   myD1Q3,    nX, 0.0, myVDW.rho0, 0.001,2 );//0.001,2
	createBoundaryPeriodic(  myLattice,  nX );
    
	for(time =0; time<= finalTime ; time++)
    {      
       for( int iX = nX+2  ; iX >=3 ; iX--) 
       {
            myLattice[iX].rho = myLattice[iX].f[ZERO] + myLattice[iX].f[DX] + myLattice[iX].f[DMX] + myLattice[iX].f[D3X] + myLattice[iX].f[DM3X];
            //      //     myLattice[iX].vel = (myD1Q3.dvD1Q3[DX]*(myLattice[iX].f[DX]-myLattice[iX].f[DMX]))/myLattice[iX].rho ; 
            //   double Force_i    = -(myLattice[iX+1].pNid - myLattice[iX-1].pNid)/(dX*2.0);
            //   myLattice[iX].vel = (myLattice[iX].vel + 0.5*dt*Force_i)/myLattice[iX].rho ; //the velocity is actually not needed this thing is for computing the kappa term gradient
       }
       /* Set Periodicity for computing gradients */
	   myLattice[0     ].vel  = myLattice[nX     ].vel ;
       myLattice[1     ].vel  = myLattice[nX + 1 ].vel ;
       myLattice[2     ].vel  = myLattice[nX + 2 ].vel ;
       myLattice[nX + 3].vel  = myLattice[3      ].vel ;
       myLattice[nX + 4].vel  = myLattice[4      ].vel ;
       myLattice[nX + 5].vel  = myLattice[5      ].vel ;        

       myLattice[0].rho = myLattice[nX].rho ;
       myLattice[1].rho = myLattice[nX + 1].rho;
       myLattice[2].rho = myLattice[nX + 2].rho;
       myLattice[nX + 3].rho = myLattice[3].rho;
       myLattice[nX + 4].rho = myLattice[4].rho;
       myLattice[nX + 5].rho = myLattice[5].rho;
        
          
        
	   collideWorking (  myLattice,   myD1Q3, myVDW,   nX,  beta,tau,time,c);  
       
	   createBoundaryPeriodic ( myLattice,  nX );
	   advect(   myLattice,       nX);

       if(time % 9900== 0) {
         printRho(myLattice,myD1Q3,myVDW,nX,beta,3,nX+2,time,c);
       }
/*  ______________________________________________________________   */	  
	}
    printRho(myLattice,myD1Q3,myVDW,nX,beta,3,nX+2,finalTime,c);
}


void   getLatticeD1Q3(double c, latticeD1Q3 *myD1Q3)
{
        double tmp;
		tmp =sqrt(10.0);
		myD1Q3->T0 = c*c/( (5.0/3.0) + (sqrt(10) /3.0)  );
		printf("\nTemp = %lf\n", myD1Q3->T0);
		myD1Q3->T0Inv = 1.0/myD1Q3->T0;
		myD1Q3->c = c; 
	   	myD1Q3->c2 = c*c;

		myD1Q3->weight[ZERO] =  (1.0/720.0)*( 64.0*(4.0+ sqrt(10)))    ;
		myD1Q3->weight[DX]   =  myD1Q3->weight[DMX]  = (1.0/720.0)*( 27.0*(8.0 - sqrt(10)))     ;
        myD1Q3->weight[D3X]  =  myD1Q3->weight[DM3X] = (1.0/720.0)*( 1.0*(16.0 - 5*sqrt(10)))    ;



		myD1Q3->dvD1Q3[ZERO] = 0.0;
		myD1Q3->dvD1Q3[DX]   = c;
		myD1Q3->dvD1Q3[DMX]  = -1.0 * c;
		myD1Q3->dvD1Q3[D3X]   = 3.0*c;
		myD1Q3->dvD1Q3[DM3X]  = -3.0 * c;
		myD1Q3->sgn[ZERO] = 0;
		myD1Q3->sgn[DX]   = 1;
		myD1Q3->sgn[DMX]  = -1 ;
		myD1Q3->sgn[D3X]   = 1;
		myD1Q3->sgn[DM3X]  = -1;
		printf("\n w0=%lf w1= %lf  \n", myD1Q3->weight[ZERO], myD1Q3->weight[DX]);
		printf("\n rho=%.16lf  T0=%.16lf\n", myD1Q3->weight[ZERO]+2.0*(myD1Q3->weight[DX] + myD1Q3->weight[DM3X]),myD1Q3->T0 );
		return;
}
		

		
		                    /*                 Periodic boundary conditions                      */
/*
 *  Periodicity:
 *
 *         NX-2.5  NX-1.5  NX-0.5   0.5     1.5    2.5    NX-2.5    NX-1.5  NX-0.5   0.5      1.5       2.5
 *          |-------|-------|--------|------|------|--------|--------|-------|--------|--------|---------|  
 * array:   0       1       2        3      4      5       NX       NX+1    NX+2     NX+3     NX+4      NX+5
 *                                           
 *
 *   *********************************************************************************/       
	
void createBoundaryPeriodic(latticeArr myLattice, int nX )
{
    int iX,  dv; 
    /* Periodic BC */
    for (dv = 0; dv < N_DV; dv++) 
    {
        myLattice[0].f[dv] = myLattice[nX     ].f[dv];
        myLattice[1].f[dv] = myLattice[nX + 1 ].f[dv];
        myLattice[2].f[dv] = myLattice[nX + 2 ].f[dv];

        myLattice[nX + 3].f[dv] = myLattice[3].f[dv];
        myLattice[nX + 4].f[dv] = myLattice[4].f[dv];
        myLattice[nX + 5].f[dv] = myLattice[5].f[dv];   
    }
    return;
}
	
void advect(latticeArr myLattice, int nX)
{
    int iX, end;
    end =nX+2;
    for(iX =3; iX<=end;iX++)
    {
        myLattice[iX] .f[DMX]   =myLattice[iX+1] .f[DMX]; 
        myLattice[iX] .f[DM3X]  =myLattice[iX+3] .f[DM3X]; 

    }
    for(iX =end; iX>=3;iX--)
    {
        myLattice[iX] .f[DX]    =myLattice[iX-1] .f[DX]; 
        myLattice[iX] .f[D3X]   =myLattice[iX-3] .f[D3X]; 

    }
    return;
} 

void initializePerturbPeriodic(latticeArr myLattice, latticeD1Q3 myD1Q3,  int nX, double vel, double rhoMean, double ampDisturb, double periodDisturb)
{
    int iX;
    double coord, rhoIn, kx;
    FILE *fpt;
    fpt = fopen("iC.dat","w");
    /*Initialise the lattice */
    kx = 2*M_PI* periodDisturb;
    for(iX =3; iX<=nX+2;iX++)
    {
        coord= (iX-2.5)/nX;
        rhoIn =rhoMean + ampDisturb*sin(kx*coord);
        myLattice[iX].rho =rhoIn;
        /*Inititalise f using equilibrium*/
        getFeqPQuad(myLattice[iX].f,myD1Q3, rhoIn, vel);
        fprintf(fpt,"%.10lf  %.10lf \n", coord, rhoIn);
    }
    createBoundaryPeriodic(  myLattice,   nX );
    fclose(fpt);
 }

void printRho(latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW, int nX, double beta,int iX_begin, int iX_end, int time, double c)
{
    FILE *fpt;
    int iX;
    char *fileName;
    double rho, vel;
    fileName = ( char  * ) malloc( (size_t) (45 * sizeof( char ) ) );
  
    if (!fileName)
    { 
        fprintf(stderr, "Memory allocation failure");
        exit(0);
    }
    sprintf(fileName,"./outputs/rho%d.dat", time);
    
    fpt=fopen(fileName,"w");

    fprintf(fpt,  " # iX \t rho \t vel \t P \t mu \n");
    
    for(iX =iX_begin; iX<iX_end;iX++)
    {
        double P = 0.0, mu = 0.0, dx=0.0, temp;
        dx = 1.0/(nX-1.0);
    
        P = myLattice[iX].rho*myD1Q3.T0/(1.0-myLattice[iX].rho*myVDW.b);
        P += -myVDW.a*myLattice[iX].rho*myLattice[iX].rho;
        P += -myVDW.kappa*myLattice[iX].rho*(myLattice[iX+1].rho+myLattice[iX-1].rho-2.0*myLattice[iX].rho)/(dx*dx);
        P += -0.5*myVDW.kappa*(myLattice[iX+1].rho - myLattice[iX-1].rho)*(myLattice[iX+1].rho - myLattice[iX-1].rho)/(4.0*dx*dx);

        //correction
        //P += -(0.5/beta - 0.25)*temp*temp/myLattice[iX].rho; //term hurts
        //P += 0.25*myVDW.kappa*(myLattice[iX+1].P+myLattice[iX-1].P-2.0*myLattice[iX].P)/(dx*dx); //term may be fine
        //P += -1.0/12.0*myVDW.kappa*(myLattice[iX+1].rho+myLattice[iX-1].rho-2.0*myLattice[iX].rho)/(dx*dx); // this term is good

        mu = myD1Q3.T0*log(myLattice[iX].rho) ;
        mu += -1.5*myD1Q3.T0*log(2.0*M_PI*myD1Q3.T0);
        mu += -2.0*myVDW.a*myLattice[iX].rho -myD1Q3.T0*log(1.0-myLattice[iX].rho*myVDW.b);
        mu += myLattice[iX].rho*myD1Q3.T0*myVDW.b/(1.0-myLattice[iX].rho*myVDW.b);
        mu += - myVDW.kappa*(myLattice[iX-1].rho + myLattice[iX+1].rho - 2.0*myLattice[iX].rho)/(dx*dx);

        fprintf(fpt,  "%d,%.7lf,%.7lf,%.7lf,%.7lf,%.7lf\n",iX, myLattice[iX].rho, myLattice[iX].vel, P, mu,myLattice[iX].chi);
    }
    fclose(fpt);
}

void collideWorking(latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW,  int nX, double beta,double tau, int time, double c)
{
    int iX,dv;
    double  fact,quad, quad1, tmp;
    double fEq[N_DV],fc[N_DV], dRho,dPnid;
    double  vel,alpha, lapRho, dmuA,dmuR, mass,sum, df, dt,dx;
    double rhoReduced ,tmp2,fact2,tauM,betaM,rhoRedby4,g;
    mass = 0.0;
    dx = 1.0/(nX-1.0);
    dt = dx / c;

//    	//define pressure
//     for( iX = nX+2 ; iX >=3 ; iX--)
//     {
// //         myLattice[iX].pNid = myVDW.b*myLattice[iX].rho*myLattice[iX].rho*myD1Q3.T0/(1-myVDW.b*myLattice[iX].rho) - myLattice[iX].rho*myLattice[iX].rho*myVDW.a ;
//        myLattice[iX].pNid = myLattice[iX].chi*myLattice[iX].rho*myD1Q3.T0 ;
//         //surface energy
//     } 



    //$pressure from VW EOS
        for( iX = nX+2 ; iX >=3 ; iX--)
    {
        myLattice[iX].pNid = (myLattice[iX].rho*myLattice[iX].rho* myVDW.b *myD1Q3.T0  )/(1.0 - myLattice[iX].rho*myVDW.b)  - 
                            myVDW.a * myLattice[iX].rho*myLattice[iX].rho ;
                  // pnid =  (rho**2 b theta) / (1 - rho * b) - a rho**2 
    }



    //$ F_nid  free energy
    for( iX = nX+2 ; iX >=3 ; iX--)
    {   
        myLattice[iX].FNid =  - myVDW.a * myLattice[iX].rho*myLattice[iX].rho   -   myLattice[iX].rho*myD1Q3.T0 *log(1.0 - myLattice[iX].rho*myVDW.b)
                            // - myVDW.kappa * 0.5*pow((myLattice[iX+1].rho -myLattice[iX-1].rho)/(dx*2.0),2) //: this will not come 
        ;
    }



    myLattice[0     ].FNid = myLattice[nX       ].FNid ;
    myLattice[1     ].FNid = myLattice[nX + 1   ].FNid ;
    myLattice[2     ].FNid = myLattice[nX + 2   ].FNid ;
    myLattice[nX + 3].FNid = myLattice[3        ].FNid ;
    myLattice[nX + 4].FNid = myLattice[4        ].FNid ;
    myLattice[nX + 5].FNid = myLattice[5        ].FNid ;    
    


    
    myLattice[0     ].pNid = myLattice[nX       ].pNid ;
    myLattice[1     ].pNid = myLattice[nX + 1   ].pNid ;
    myLattice[2     ].pNid = myLattice[nX + 2   ].pNid ;
    myLattice[nX + 3].pNid = myLattice[3        ].pNid ;
    myLattice[nX + 4].pNid = myLattice[4        ].pNid ;
    myLattice[nX + 5].pNid = myLattice[5        ].pNid ;    

   	//calculate chemical potential
   	for( iX = nX+2; iX >=3; iX--)
    {
        //from continuous derivative
        myLattice[iX].muA  = -myD1Q3.T0*log(1.0 - myLattice[iX].rho*myVDW.b) ;
        myLattice[iX].muA += myLattice[iX].rho*myVDW.b*myD1Q3.T0/(1.0 - myLattice[iX].rho*myVDW.b);
        myLattice[iX].muA -= 2.0*myLattice[iX].rho*myVDW.a ;

    }

	myLattice[0     ].muA = myLattice[nX    ].muA    ;
    myLattice[1     ].muA = myLattice[nX + 1].muA ;
    myLattice[2     ].muA = myLattice[nX + 2].muA ;
    myLattice[nX + 3].muA = myLattice[3     ].muA ;
    myLattice[nX + 4].muA = myLattice[4     ].muA ;
    myLattice[nX + 5].muA = myLattice[5     ].muA ;

    for( iX = nX+2; iX >=3; iX--)
    {

        myLattice[iX].muA = myLattice[iX].muA
                            + (dx*dx/6.0)*(myLattice[iX -1].muA + myLattice[iX +1].muA - 2.0*myLattice[iX].muA)/(dx*dx)                               
                        ;    


        double kappa1 = myVDW.kappa; 

        // // # correction from chemical potential
        // kappa1 =   myVDW.kappa 
        //             -   (dx*dx/6.0)    *(   
        //                                 (myLattice[iX].rho * myVDW.b * myVDW.b * myD1Q3.T0)/(pow(1.0 - myLattice[iX].rho * myVDW.b ,2)) 
        //                             +   (myVDW.b*myD1Q3.T0)/( 1.0 - myLattice[iX].rho * myVDW.b  )   - 2.0*myVDW.a
        //                                 )
        //  ;

        myLattice[iX].muA -= kappa1*(myLattice[iX-1].rho + myLattice[iX+1].rho - 2.0*myLattice[iX].rho)/(dx*dx) ;

    }

    myLattice[0     ].muA = myLattice[nX    ].muA    ;
    myLattice[1     ].muA = myLattice[nX + 1].muA ;
    myLattice[2     ].muA = myLattice[nX + 2].muA ;
    myLattice[nX + 3].muA = myLattice[3     ].muA ;
    myLattice[nX + 4].muA = myLattice[4     ].muA ;
    myLattice[nX + 5].muA = myLattice[5     ].muA ;



    //  calcualte the surface terms contribution from kappa
    for( iX = nX+2; iX >=3; iX--)
    {   

        // std::cout<<"before "<<myVDW.kappa<<std::endl;

        double kappa1 = myVDW.kappa ;

        // #correction from pressure
        // double kappa1 =   myVDW.kappa 
        //             -   (dx*dx/6.0)    *(   
        //                                 (myLattice[iX].rho*myLattice[iX].rho * myVDW.b * myVDW.b * myD1Q3.T0)/(pow(1.0 - myLattice[iX].rho * myVDW.b ,2)) 
        //                                 +   (2.0*myVDW.b*myLattice[iX].rho*myD1Q3.T0)/( 1.0 - myLattice[iX].rho * myVDW.b  )   - 2.0*myVDW.a*myLattice[iX].rho
        //                                 )
        // ;

        
        //  std::cout<<"after  "<<kappa1<<std::endl;

        myLattice[iX].surface = kappa1*  (myLattice[iX-1].rho + myLattice[iX+1].rho - 2.0*myLattice[iX].rho)/(dx*dx) ;
    }

    myLattice[0].surface = myLattice[nX].surface ;
    myLattice[1].surface = myLattice[nX + 1 ].surface ;
    myLattice[2].surface = myLattice[nX + 2 ].surface ;
    myLattice[nX + 3].surface = myLattice[3].surface ;
    myLattice[nX + 4].surface = myLattice[4].surface ;
    myLattice[nX + 5].surface = myLattice[5].surface ;
    
    
    

    //calculate force
    for( iX = nX+2; iX >=3 ; iX--)
    {
        mass += myLattice[iX].rho;
        // // a = 1 represents 2nd order a = 4/3 represents 4th order 
        double a = 1.0; double b = 1.0 - a;


        myLattice[iX].Force =   myLattice[iX].rho*(myLattice[iX+1].muA -myLattice[iX-1].muA)/(dx*2.0) *a + 
                                myLattice[iX].rho*(myLattice[iX+2].muA -myLattice[iX-2].muA)/(dx*4.0) *b 
        ;



        //#presure formulation

        // myLattice[iX].Force = (myLattice[iX+1].pNid - myLattice[iX-1].pNid)/(dx*2.0) - ((myLattice[iX + 1 ].surface - myLattice[iX -1 ].surface) )  / (dx*2.0)
        // ;



        //#second order verification of pnid fnid with the this is just to check whether the pnid and fnid are entered correctly since munid = (pnid + fnid)/(rho)
        // double del_muA_second_order = (1.0/(2.0 * dx))*((myLattice[iX+1].pNid  + myLattice[iX+1].FNid)/(myLattice[iX +1].rho) - (myLattice[iX-1].pNid  + myLattice[iX-1].FNid)/(myLattice[iX -1].rho)
        //                                                 -(myVDW.kappa* (myLattice[iX + 1 ].surface - myLattice[iX -1 ].surface) )
        //                                                 )  ;


        // myLattice[iX].Force = myLattice[iX].rho*(del_muA_second_order);



        

        //  //#fourth order to second order 
        // double del_muA_fourth_order =   (myLattice[iX  ].pNid  + myLattice[iX  ].FNid)*
        //                                 (   (a/(2.0*dx)) *( (1.0 /myLattice[iX +1].rho ) - (1.0 /myLattice[iX -1].rho )) +   
        //                                     (b/(4.0*dx)) *( (1.0 /myLattice[iX +2].rho ) - (1.0 /myLattice[iX -2].rho ))  ) +
                                       
        //                                 (1.0 /myLattice[iX   ].rho )*
        //                                 (   (a/(2.0*dx)) *( ( myLattice[iX+1].pNid  + myLattice[iX+1].FNid ) - (myLattice[iX-1].pNid  + myLattice[iX-1].FNid)) +   
        //                                     (b/(4.0*dx)) *( ( myLattice[iX+2].pNid  + myLattice[iX+2].FNid ) - (myLattice[iX-2].pNid  + myLattice[iX-2].FNid))  ) 

        //                 - (1.0/ (2.0 *dx))*(myVDW.kappa* (myLattice[iX + 1 ].surface - myLattice[iX -1 ].surface) )
        //                 ;
 
        // myLattice[iX].Force = myLattice[iX].rho*(del_muA_fourth_order); //total force

        // # -------------------------------------

        myLattice[iX].Force = -myLattice[iX].Force /myLattice[iX].rho; //force density

        myLattice[iX].vel = (myD1Q3.dvD1Q3[DX]*(myLattice[iX].f[DX]-myLattice[iX].f[DMX]) + 3*myLattice[iX].f[D3X] - 3*myLattice[iX].f[DM3X]  ) ;

        myLattice[iX].vel = myLattice[iX].vel/myLattice[iX].rho +  0.5*dt*(myLattice[iX].Force);
        
        getFeqPQuad( myLattice[iX].fEq ,   myD1Q3,    myLattice[iX].rho,   myLattice[iX].vel); 
    } 

    for( iX = nX+2  ; iX >=3 ; iX--)
    {
	    alpha = 2.0;
      	    // calculateAlpha( myLattice, alpha, iX, time);
	
	    double chi = 0.0;

        rhoReduced = myVDW.b*myLattice[iX].rho;
        chi = rhoReduced/(1.0-rhoReduced) - myVDW.a*myLattice[iX].rho/myD1Q3.T0;
        tauM = tau;//*myVDW.rho0/myLattice[iX].rho *(1 + rhoReduced*(5.0/8.0 + rhoReduced*(0.2869+ rhoReduced*(0.1103+ 0.0386*rhoReduced))));
//             tauM = tauM/(1.0+chi);
	    betaM = dt/(2.0*tauM + dt);
      
        for(dv = 0; dv < N_DV; dv++)
        {
            myLattice[iX].f[dv]  = myLattice[iX] .f[dv] + 2.0*betaM*(myLattice[iX].fEq[dv] - myLattice[iX] .f[dv]); 
            myLattice[iX].f[dv] += alpha*betaM*tauM*myD1Q3.T0Inv*myLattice[iX].rho*myD1Q3.weight[dv]*myD1Q3.dvD1Q3[dv]*myLattice[iX].Force;
        }      
    }
      
    if(time % 9900 == 0)
    printf("\nAt time =%d mass = %.16lf ", time, mass);



}

void getFeqPQuad(double fEq[N_DV], latticeD1Q3 myD1Q3,  double rho,  double vel)
{
   
    double  fact, fact3,quad, quad1, tmp,tmp3, cubic;
    double  ratio, ratio3, quad3, cubic1, vel1;
    ratio = myD1Q3.T0Inv*myD1Q3.c2;
    ratio3 =  myD1Q3.T0Inv*myD1Q3.c2 *3*3;
    quad = myD1Q3.T0Inv *0.5*  vel*vel ;
  
    fEq[ZERO]= myD1Q3.weight[ZERO] *rho*  (1.0  -  quad) ;

    fact    =  myD1Q3.T0Inv * vel*myD1Q3.dvD1Q3[DX];
    fact3   =  myD1Q3.T0Inv * vel*myD1Q3.dvD1Q3[D3X];

    quad1  = quad *(ratio  -1.0); 
    quad3  = quad *(ratio3 -1.0);
    tmp = rho* myD1Q3.weight[DX];
    tmp3 = rho* myD1Q3.weight[D3X];

    fEq[DX]= tmp*(1.0 + fact + quad1   ); 
    fEq[DMX]= tmp*(1.0 - fact + quad1 ); 
  
    fEq[D3X]= tmp3*(1.0 + fact3 + quad3   ); 
    fEq[DM3X]= tmp3*(1.0 - fact3 + quad3   ); 
    return;
}


void calculateAlpha( latticeArr myLattice, double& alpha, int i,int time)
{
      alpha = 2.0; 
      
      double a(0.), b(0.), c(0.), ximax;
      double x_i[N_DV];
      
      for(int dv = 0; dv<N_DV; dv++)
      {
	    x_i[dv] = myLattice[i].fEq[dv]/myLattice[i].f[dv]-1.0;
      }
	
      ximax = 0.0;
      for(int dv = 0; dv<N_DV; dv++)
      {
	    ximax = std::max(fabs(ximax),x_i[dv]);
      }
      
      if(ximax>0.01)
      {
	  for(int dv = 0; dv<N_DV; dv++)
	  {
            if(x_i[dv]<0.0)
            {
            a += myLattice[i].f[dv]*x_i[dv]*x_i[dv]*x_i[dv]*0.5;
            c += myLattice[i].f[dv]*x_i[dv]*x_i[dv];
            }
            
            else
            {
            c += myLattice[i].f[dv]*x_i[dv]*x_i[dv]/(1.0+x_i[dv]);
            }
            
            b += myLattice[i].f[dv]*x_i[dv]*x_i[dv]*0.5;
	  }
	  alpha = (b-sqrt(b*b - 4.0*a*c))/(2.0*a);       
      }
      
//       if( (time==3260) && (i==46)) 
//       {
// 	printf("\n f_i  %.16lf  %.16lf %.16lf   ", myLattice[i].f[DX],myLattice[i].f[ZERO],myLattice[i].f[DMX]);
// 	printf("\n fEq  %.16lf  %.16lf %.16lf   ", myLattice[i].fEq[DX],myLattice[i].fEq[ZERO],myLattice[i].fEq[DMX]);
// 	printf("\n x_i  %.16lf  %.16lf %.16lf   ", x_i[DX],x_i[ZERO],x_i[DMX]);
// 	printf("\n f*x  %.16lf  %.16lf %.16lf   ", myLattice[i].f[DX]*x_i[DX],myLattice[i].f[ZERO]*x_i[ZERO],myLattice[i].f[DMX]*x_i[DMX]);
// 	printf("\n abc  %.16lf  %.16lf %.16lf  alpha %.16lf ", a, b ,c, alpha);
// 	printf("\n %16lf diff in f and feq",myLattice[i].f[ZERO] + myLattice[i].f[DX] + myLattice[i].f[DMX] - myLattice[i].fEq[ZERO] - myLattice[i].fEq[DX] - myLattice[i].fEq[DMX]);
// // 	std::cout << "\n"<< b*b - 4.0*a*c << "\n";
//       }
	for(int dv = 0; dv<N_DV; dv++)
	{
	  if( (myLattice[i].f[dv]<0.0) && (a>=0.0) )
	  {
	    alpha = 2.0;
	    printf(" forcing alpha=2 here %d \t %d \t %.16lf \t %.16lf \n",i,dv, myLattice[i].f[dv],a);
	    break;
	  }
// 	  if (time==3427) std::cout << myLattice[i] .f[dv]<< "  " << myLattice[i] .fEq[dv]<< "  " << alpha << "  "<< dv << "  "<<i << "\n";
	}
     
}






// template<typename dataType1>
// void calculateAlpha(lbmRD3Q41<dataType1> &lbModel,dataType1* x_i,dataType1* f_i,dataType1 beta,int index,dataType1& alpha)
// {
//   double a(0.0), b(0.0), c(0.0),oneBySix(1.0/6.0);
//   double ximin(0.0), ximax(0.0);

//   alignas(32) dataType1 xSq[41];
//   alignas(32) dataType1 fxSq[41];

//   for(int dv = 0;dv<41;dv++)
//   {
//     ximin = std::min(ximin, x_i[dv])  ;
//     ximax = std::max(ximax, x_i[dv])  ;
//   }

//   for(int dv = 0;dv<41;dv++)
//   {
//     xSq[dv]  = x_i[dv]*x_i[dv] ;
//     fxSq[dv] = f_i[dv]*x_i[dv]*x_i[dv] ;

//     if(x_i[dv]<0.0)
//       a += fxSq[dv]*x_i[dv]*0.5 ;

//     b += fxSq[dv]*0.5 ;
//     c += fxSq[dv]/(1.0 + 0.5*x_i[dv]) ;
//   }

//   dataType1 alphaMax = -1.0/(beta*ximin);

//   dataType1 k;
//   if(a<0 && b>0 && c>0)
//     k = (b-sqrt(b*b - 4.0*a*c))/(2.0*a);
//   else
//     k = 1.5;

//   a = 0.0;b=0.0;c=0.0;
// //   dataType1 kBeta   = k*beta;
//   dataType1 beta2   = beta*beta;
//   dataType1 fourByK = 4.0/k;
//   dataType1 hBeta = 0.0;

//   for(int dv = 0; dv < 41; dv++)
//   {
//     if(x_i[dv]<0.0)
//     {
//       a += fxSq[dv]*x_i[dv]*beta2*( 1.0/6.0 - hBeta*x_i[dv]/12.0 + hBeta*hBeta*x_i[dv]*x_i[dv]/20 - hBeta*hBeta*hBeta*x_i[dv]*x_i[dv]*x_i[dv]/5.0 );
//       b += fxSq[dv]*0.5;
//     }

//     if(x_i[dv]>0.0)
//       b += f_i[dv]*( (x_i[dv]*x_i[dv]*0.5) - beta2*(x_i[dv]*x_i[dv]*x_i[dv]/15.0)* ( (4.0/(fourByK+x_i[dv])) + (2.0/(fourByK+2.0*x_i[dv])) + (4.0/(fourByK+3.0*x_i[dv])) ));

//     c += f_i[dv]*(60.0*x_i[dv]*x_i[dv] + 60.0*x_i[dv]*x_i[dv]*x_i[dv] + 11.0*x_i[dv]*x_i[dv]*x_i[dv]*x_i[dv])/( 60.0 + 90.0*x_i[dv] + 36.0*x_i[dv]*x_i[dv] + 3.0*x_i[dv]*x_i[dv]*x_i[dv]);
//   }

//   dataType1  h;
//   if(a<0 && b>0 && c>0)
//     h = (b-std::sqrt(b*b - 4.0*a*c))/(2.0*a);
//   else
//     h = 2.1;

//   a = 0.0;
//   hBeta = h*beta;

//   for(int dv = 0; dv < 41; dv++)
//   {
//     if(x_i[dv]<0.0)
//     {
//       a += f_i[dv]*x_i[dv]*x_i[dv]*x_i[dv]*beta2*( 1.0/6.0 - hBeta*x_i[dv]/12.0 + hBeta*hBeta*x_i[dv]*x_i[dv]/20 - hBeta*hBeta*hBeta*x_i[dv]*x_i[dv]*x_i[dv]/5.0 );
//     }
//   }

//   if(a<0 && b>0 && c>0)
//     alpha = 2.0*c/(b+sqrt(b*b - 4.0*a*c));

//   if(alpha > alphaMax)
//   {
//     if ( alphaMax > 1.0)
//       alpha = 0.5*(1.0+alphaMax) ;
//     else
//       alpha = 0.95*alphaMax;
//   }
// }





// //just some tests
//     getFeqPQuad(myLattice[5].fEq ,   myD1Q3,   1.0,   0.0025231234);
//     std::cout<<myLattice[5].fEq[1] - myLattice[5].fEq[2] + 3*myLattice[5].fEq[3] -3* myLattice[5].fEq[4]<<std::endl; 