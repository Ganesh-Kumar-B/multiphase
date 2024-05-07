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
#define MAX_SIZE  5000 
#define N_DV 5
 
enum direction
{
    ZERO, /* Alias for 0 */
    DX, /* Alias for 1 */
    DMX, /* Alias for -1 */
    D3X, /* Alias for 3 */
    DM3X, /* Alias for -3 */
};

typedef  struct nodeD1Q5
{
    double rho;
    double rhoOld;
    double pNid;
    double FNid;
    double chi;
    double muR;
    double muA;
    double dmuA;
    double dmuR;
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
} nodeD1Q5;

typedef  nodeD1Q5 latticeArr[MAX_SIZE]  ;
 
typedef  struct latticeD1Q5
{
    double weight[N_DV];
    double dvD1Q5[N_DV];
    int sgn[N_DV];
    double T0;
    double T0Inv;
    double c;
    double c2;
} latticeD1Q5;


typedef  struct nonIdealParam  
{
    double a;
    double b;
    double alpha1;
    double kappa;
    double rho0;
} nonIdealParam ;

void getLatticeD1Q5(double c, latticeD1Q5 *myD1Q5);
void printRho(latticeArr lattice, latticeD1Q5 myD1Q5, nonIdealParam myVDW,int nX, double beta,int iX_begin, int iX_end, int time,double c);
void getFeqPQuad(double fEq[N_DV], latticeD1Q5 myD1Q5,  double rho,  double vel);

void createBoundaryPeriodic(latticeArr lattice, int nX );
void collideWorking (latticeArr myLattice, latticeD1Q5 myD1Q5, nonIdealParam myVDW, int nX,double dx, double beta,double tau, int time,double c);
void advect(latticeArr lattice, int nX);
void calculateAlpha( latticeArr myLattice, double& alpha, int i, int time);
void calculateAlpha1( latticeArr myLattice, double& alpha, int i, int time, double beta);

void initializePerturbPeriodic(latticeArr lattice, latticeD1Q5 myD1Q5,  int nX, double vel, double rhoMean, double ampDisturb, double periodDisturb);


main()
{
    

    latticeArr myLattice;
	latticeD1Q5 myD1Q5;
	double kappabar, TbyTc,rho0byrhoc;
	double c, rho,u, beta,  T_critical, tau, dt, dX, rho_critical;
	nonIdealParam myVDW;
	int iX,nX,time, finalTime;

    myVDW.b = 4.0;
    myVDW.a = 1.0 ;

    T_critical      = (0.377332*myVDW.a)/myVDW.b;
	TbyTc           = 0.90   ;
    double T_actual = TbyTc* T_critical;

    rho_critical    = 0.521772/myVDW.b;
	rho0byrhoc      = 1.0;
	myVDW.rho0      = rho0byrhoc*rho_critical;	



    c = sqrt( T_actual *  ((5.0/3.0) + (sqrt(10.0) /3.0)));

    std::cout<<" c      :   "<<c<<std::endl;
	
    getLatticeD1Q5(c, &myD1Q5);
    nX = 500;

	kappabar = 1.0;//0.0625;
            
	dX = pow(myVDW.b,1.0/3.0);
	dt = dX/c;
	
	beta = 0.5; 
	tau = (1.0-beta)/beta *dt*0.5;   
    
    std::cout<<" TbyTc  :   "<<TbyTc<<std::endl;
    std::cout<<"dx      :   "<<dX<<std::endl;
    std::cout<<"dt      :   "<<dt<<std::endl;
    std::cout<<"c       :   "<<myD1Q5.c<<std::endl;
    std::cout<<"tau     :   "<<tau<<std::endl;


	myVDW.kappa = kappabar*myVDW.a*dX*dX;
    std::cout<<"kappa "<< myVDW.kappa << std::endl;


    finalTime = 100000    ;

	initializePerturbPeriodic(  myLattice,   myD1Q5,    nX, 0.0, myVDW.rho0, 0.001, 2 );//0.001,2

	createBoundaryPeriodic(  myLattice,  nX );
    
	for(time =0; time<= finalTime ; time++)
    {      
        for( int iX = nX+2  ; iX >=3 ; iX--) 
        {
            myLattice[iX].rho = myLattice[iX].f[ZERO] + myLattice[iX].f[DX] + myLattice[iX].f[DMX] + myLattice[iX].f[D3X] + myLattice[iX].f[DM3X];
            
        }
        /* Set Periodicity for computing gradients */
        myLattice[0     ].vel   = myLattice[nX     ].vel ;
        myLattice[1     ].vel   = myLattice[nX + 1 ].vel ;
        myLattice[2     ].vel   = myLattice[nX + 2 ].vel ;
        myLattice[nX + 3].vel   = myLattice[3      ].vel ;
        myLattice[nX + 4].vel   = myLattice[4      ].vel ;
        myLattice[nX + 5].vel   = myLattice[5      ].vel ;        

        myLattice[0].rho        = myLattice[nX].rho ;
        myLattice[1].rho        = myLattice[nX + 1].rho;
        myLattice[2].rho        = myLattice[nX + 2].rho;
        myLattice[nX + 3].rho   = myLattice[3].rho;
        myLattice[nX + 4].rho   = myLattice[4].rho;
        myLattice[nX + 5].rho   = myLattice[5].rho;


        collideWorking (  myLattice,   myD1Q5, myVDW,   nX,dX,  beta,tau,time,c);  
        
        createBoundaryPeriodic ( myLattice,  nX );
        advect(   myLattice,       nX);

        if(time % 1000== 0) {
            printRho(myLattice,myD1Q5,myVDW,nX,beta,3,nX+2,time,c);
        }
/*  ______________________________________________________________   */	  
	}
    printRho(myLattice,myD1Q5,myVDW,nX,beta,3,nX+2,finalTime,c);
}


void   getLatticeD1Q5(double c, latticeD1Q5 *myD1Q5)
{
        
		myD1Q5->T0 = c*c/( (5.0/3.0) + (sqrt(10) /3.0)  );
		printf("\nTemp = %lf\n", myD1Q5->T0);
		myD1Q5->T0Inv = 1.0/myD1Q5->T0;
		myD1Q5->c = c; 
	   	myD1Q5->c2 = c*c;

		myD1Q5->weight[ZERO] =  (1.0/720.0)*( 64.0*(4.0+ sqrt(10)))    ;
		myD1Q5->weight[DX]   =  myD1Q5->weight[DMX]  = (1.0/720.0)*( 27.0*(8.0 - sqrt(10)))     ;
        myD1Q5->weight[D3X]  =  myD1Q5->weight[DM3X] = (1.0/720.0)*( 1.0*(16.0 - 5*sqrt(10)))    ;



		myD1Q5->dvD1Q5[ZERO ] = 0.0;
		myD1Q5->dvD1Q5[DX   ]   = c;
		myD1Q5->dvD1Q5[DMX  ]  = -1.0 * c;
		myD1Q5->dvD1Q5[D3X  ]   = 3.0*c;
		myD1Q5->dvD1Q5[DM3X]  = -3.0 * c;

		myD1Q5->sgn[ZERO] = 0;
		myD1Q5->sgn[DX]   = 1;
		myD1Q5->sgn[DMX]  = -1 ;
		myD1Q5->sgn[D3X]   = 1;
		myD1Q5->sgn[DM3X]  = -1;
		printf("\n w0=%lf w1= %lf  \n", myD1Q5->weight[ZERO], myD1Q5->weight[DX]);
		printf("\n rho=%.16lf  T0=%.16lf\n", myD1Q5->weight[ZERO]+2.0*(myD1Q5->weight[DX] + myD1Q5->weight[DM3X]),myD1Q5->T0 );
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

void initializePerturbPeriodic(latticeArr myLattice, latticeD1Q5 myD1Q5,  int nX, double vel, double rhoMean, double ampDisturb, double periodDisturb)
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
        myLattice[iX].rho = rhoIn;
        /*Inititalise f using equilibrium*/
        getFeqPQuad(myLattice[iX].f,myD1Q5, rhoIn, vel);
        fprintf(fpt,"%.10lf  %.10lf \n", coord, rhoIn);
    }
    createBoundaryPeriodic(  myLattice,   nX );
    fclose(fpt);
 }

void printRho(latticeArr myLattice, latticeD1Q5 myD1Q5, nonIdealParam myVDW, int nX, double beta,int iX_begin, int iX_end, int time, double c)
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

    fprintf(fpt,  "iX \t rho \t vel \t mu \n");
    
    for(iX =iX_begin; iX<iX_end;iX++)
    {
        double P = 0.0, mu = 0.0, dx=0.0, temp;


        fprintf(fpt,  "%d,%.7lf,%.7lf,%.7lf\n",iX, myLattice[iX].rho/myVDW.rho0, myLattice[iX].vel, myLattice[iX].muA);
    }
    fclose(fpt);
}

void collideWorking(latticeArr myLattice, latticeD1Q5 myD1Q5, nonIdealParam myVDW,  int nX,double dx, double beta,double tau, int time, double c)
{
    int iX,dv;
    double  fact,quad, quad1, tmp;
    double fEq[N_DV],fc[N_DV], dRho,dPnid;
    double  vel,alpha, lapRho, dmuA,dmuR, mass,sum, df, dt;
    double rhoReduced ,tmp2,fact2,tauM,betaM,rhoRedby4,g;



    mass = 0.0;
    dt = dx / c;



   	//calculate chemical potential
   	for( iX = nX+2; iX >=3; iX--)
    {
        double eta_EOS = myLattice[iX].rho * myVDW.b / 4.0;

        myLattice[iX].muA = myD1Q5.T0*(3.0*pow(eta_EOS,3) - 9.0 *eta_EOS*eta_EOS + 8.0 *eta_EOS )/(pow(1 - eta_EOS,3));
        myLattice[iX].muA -= 2.0*myLattice[iX].rho*myVDW.a ;

        myLattice[iX].muA -= myVDW.kappa*(myLattice[iX-1].rho + myLattice[iX+1].rho - 2.0*myLattice[iX].rho)/(dx*dx) ;

    }
	myLattice[0     ].muA = myLattice[nX    ].muA    ;
    myLattice[1     ].muA = myLattice[nX +1 ].muA ;
    myLattice[2     ].muA = myLattice[nX +2 ].muA ;
    myLattice[nX + 3].muA = myLattice[3     ].muA ;
    myLattice[nX + 4].muA = myLattice[4     ].muA ;
    myLattice[nX + 5].muA = myLattice[5     ].muA ;



    //calculate force
    for( iX = nX+2; iX >=3 ; iX--)
    {
        mass += myLattice[iX].rho;
		
        myLattice[iX].Force = myLattice[iX].rho*(myLattice[iX+1].muA -myLattice[iX-1].muA)/(dx*2.0);
        
        myLattice[iX].Force = -myLattice[iX].Force /myLattice[iX].rho; //force density

        myLattice[iX].vel = (myD1Q5.dvD1Q5[DX]*(myLattice[iX].f[DX]-myLattice[iX].f[DMX]) + 3*myLattice[iX].f[D3X] - 3*myLattice[iX].f[DM3X]  ) ;

        myLattice[iX].vel = myLattice[iX].vel/myLattice[iX].rho +  0.5*dt*(myLattice[iX].Force);
        
        getFeqPQuad( myLattice[iX].fEq ,   myD1Q5,    myLattice[iX].rho,   myLattice[iX].vel); 
    } 
    for( iX = nX+2  ; iX >=3 ; iX--)
    {
	    alpha = 2.0;
       
        for(dv = 0; dv < N_DV; dv++)
        {
            myLattice[iX].f[dv]  = myLattice[iX] .f[dv] + alpha*beta*(myLattice[iX].fEq[dv] - myLattice[iX] .f[dv]); 
            myLattice[iX].f[dv] += alpha*beta*tau*myD1Q5.T0Inv*myLattice[iX].rho*myD1Q5.weight[dv]*myD1Q5.dvD1Q5[dv]*myLattice[iX].Force;
        }      
    }
      
    if(time % 1000 == 0)
    printf("\nAt time =%d mass = %.16lf ", time, mass);
}

void getFeqPQuad(double fEq[N_DV], latticeD1Q5 myD1Q5,  double rho,  double vel)
{
   
    double u2 = vel*vel ;

    double first,second, third,feq0=0;

    for (int dv = 0; dv< 5; dv++){

        feq0 = rho*myD1Q5.weight[dv];

        first  = (vel*myD1Q5.dvD1Q5[dv])*myD1Q5.T0Inv;
        second = 0.5*(first * first);
        third = -0.5*u2*myD1Q5.T0Inv ;

        fEq[dv] = feq0*(1+ first + second + third);    

    }

    return;
}





