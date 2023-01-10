//
//  main.cpp
//  HartreeFock
//
//  Created by Junlan Lu on 2018/11/18.
//  Run on a core i7 with XCode 9.4.1 on MacOS 10.14
//

#include <iostream>
#include <cmath>
#include "bigarrayt.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
using namespace std;
/* Defines constant RMAX to be the maximum radius of the lattice, ABhor to be bohr radius, etc.
 */
double RMAX=5, e2 =14.409, h1,h2m=3.81795,ABhor = 0.5299,zstar=2.0-5.0/16.0;

/*Implements the trapezoid integration technique to integrate an integrand dat from 0 to RMAX
 */
double integrate_trapezoid(arrayt<double> dat){
    int i,n=dat.n1();
    double dr = RMAX/n;
    double I = 0.5*dr*(dat(0)+dat(n-1));
    for(i = 1;i<n; i++){
        I = I + dat(i)*dr;
    }
    return I;
}

/* normalizes the R wavefunction
 */
void normalize(arrayt<double> &R){
    int i,n=R.n1();
    double I;
    arrayt<double> R2(n);
    for(i=0;i<n;i++){
        R2(i) = R(i)*R(i);
    }
    I = 1/sqrt(integrate_trapezoid(R2));
    for(i=0;i<n;i++){
        R(i) = R(i)*I;
    }
}

/* solves for the density given the radial wavefunction
 */
void calculateRho(arrayt<double> &Rho, arrayt<double> &R){
    int i,n = R.n1();
    double r,dr = RMAX/n;
    Rho(0)=0;
    for(i=1;i<n;i++){
        r = i*dr;
        Rho(i)=R(i)*R(i)/(2*M_PI*r*r);
    }
}

/* solves for the potential term given the density by solving Poisson's equation
 */
void calculatePhi(arrayt<double> &Rho,arrayt<double> &Phi){
    int i, n=Rho.n1();
    double r,dr = RMAX/n, dr2=dr*dr, m;
    arrayt<double> phi(n), integrand(n);
    for(i=0;i<n;i++){
        r = i*dr;
        integrand(i) = r*Rho(i);
    }
    phi(0) = 0;
    phi(1) = 4*M_PI*e2*dr*integrate_trapezoid(integrand);
    for(i=2;i<n;i++){
        r = i*dr;
        phi(i) = (dr2/12)*(-4*M_PI*e2*r)*(Rho(i)+10*Rho(i-1)+Rho(i-2))-phi(i-2)+2*phi(i-1);
    }
    
    Phi(0)=0;
    for(i=1;i<n;i++){
        r = i*dr;
        Phi(i) = phi(i)/r;
    }
    
    m=(phi(n-1)-phi(n-11))/(10*dr);
    for(i=1;i<n;i++){
        r = i*dr;
        Phi(i) = Phi(i)-m;
    }
}


/* calculates the total energy of the 2-electron system, given the density, potential,
  and radial wave function
 */
double calculateE(arrayt<double> &Phi, arrayt<double> &Rho, arrayt<double> &R, double z){
    int i, n = R.n1();
    double int1=0,int2=0,int3=0,E,r,dr = RMAX/n,DR,KEN,VEN,VEE;
    arrayt<double> integrand(n),integrand1(n-1), integrand2(n-1),integrand3(n-1);
    for(i=0;i<n-1;i++){
        r=i*dr;
        DR=(R(i+1)-R(i))/dr;
        integrand1(i)= DR*DR;
        integrand2(i)= Rho(i)*r;
        integrand3(i)= Phi(i)*Rho(i)*r*r;
    }
    int1=integrate_trapezoid(integrand1);
    KEN=2*h2m*int1;
    int2=integrate_trapezoid(integrand2);
    VEN=-4*M_PI*z*e2*int2;
    int3=integrate_trapezoid(integrand3);
    VEE=M_PI*int3;
    E = KEN+VEE+VEN;
    return E;
}

/* given an eigenvalue epsilon, calculate R(RMAX) by using the numerov method
 */
double calculateRMax(arrayt<double> &Phi, arrayt<double> &R, double epsilon, double z){
    int n=Phi.n1(),i;
    double Rmax, r, dr = RMAX/n, const1=dr*dr/(12*h2m), term1, term2, term3;
    R(0) = 0;
    R(1) = 1;
    for(i=2;i<n;i++){
        r = i*dr;
        term1 = 2*(1-5*const1*(epsilon+z*e2/(r-dr)-Phi(i-1)/2))*R(i-1);
        term2 = (1+const1*(epsilon+z*e2/(r-2*dr)-Phi(i-2)/2))*R(i-2);
        if(i==2) term2 = 0;
        term3 = (1+const1*(epsilon+z*e2/r-Phi(i)/2));
        R(i) = (term1 - term2)/ term3;
    }
    normalize(R);
    Rmax=R(n-1);
    return Rmax;
}

/*bisect function takes in a specified function, a bracket x1-x2, and a tolerance value. It returns
 -1 if a root is not bracketed and the calculated root within the tolerance specifation if
 the root is found. Uses the bisection method.
 */
double bisect(double(*f)(arrayt<double> &, arrayt<double> &, double, double), double x1, double x2, double tol, arrayt<double> &Phi, arrayt<double> &R, double z){
    double f1, f2, f3 = -1, x3;
    f1 = f(Phi, R, x1, z);
    f2 = f(Phi, R, x2, z);
    if((f1 * f2 > 0)){
        return 1;
        
    }
    do {
        x3 = 0.5 * (x1 + x2);
        f3 = f(Phi, R, x3, z);
        if(f1 * f3 > 0){
            x1 = x3;
            f1 = f3;
        }
        else if (f2 * f3 > 0){
            x2 = x3;
            f2 = f3;
        }
    } while (abs(f3) > tol);
    
    return x3;
}

/* implements the shooting method to solve for R(r)
 */
double shoot(arrayt<double> &Phi, arrayt<double> &R, double z){
    int count=0, maxcount = 1e5;
    double epsilon = 1, trial1, trial2, increment = 1, tol = 1e-5;
    trial1 = -200;
    trial2 = trial1+increment;
    while(epsilon > 0){
        epsilon = bisect(calculateRMax, trial1, trial2, tol, Phi, R, z);
        trial1 = trial1 + increment;
        trial2 = trial2 + increment;
        count++;
        if(count > maxcount){
            cout << "Cannot find any eigenvalues" << endl;
            break;
        }
    }
    calculateRMax(Phi, R, epsilon, z);
    return epsilon;
}

/* test function to test the CalculateE() method
 */
void testCalculateE(){
    int i,n=1e4;
    double const1=2*sqrt(zstar/ABhor)*zstar/ABhor,dr=RMAX/n,r,E, z=2.0;
    arrayt<double> R(n),Phi(n),Rho(n);
    for(i=0;i<n;i++){
        r=i*dr;
        R(i) = const1*r*exp(-zstar*r/ABhor);
    }
    normalize(R);
    calculateRho(Rho, R);
    calculatePhi(Rho, Phi);
    E=calculateE(Phi,Rho,R,z);

    cout<<"Calculated ground state energy of Helium: "<< E << endl;
    cout<<"Analytical value of ground state energy of Helium: "<<-e2*(z*z-5*z/8+25.0/256)/ABhor<<endl;
    
    ofstream fp1;
    fp1.open("testCalculateE.dat");
    if( fp1.fail() ) {
        cout << "cannot open file" << endl;
    }
    for(i=0;i<n;i++){
        r=i*dr;
        fp1 << setw(15)<< r << setw(15) << R(i)*R(i) << setw(15) << Rho(i) << setw(15) << Phi(i) << endl;
    }
    fp1.close();
}

/* test function to test the CalculateRMAX() method for an arbitrary epsilon
 */
void testCalculateRMax(){
    int i,n = 1e4;
    double dr = RMAX/n,r,z=2.0;
    arrayt<double> Phi(n), R(n), Rho(n);
    for(i=0; i < n ;i++){
        Phi(i)=0;
        R(i)=2;
    }
    calculateRMax(Phi, R, -20,z);
    ofstream fp1;
    fp1.open("testRMax.dat");
    if( fp1.fail() ) {
        cout << "cannot open file" << endl;
    }
    for(i=0;i<n;i++){
        r=i*dr;
        fp1 << setw(15)<< r << setw(15) << R(i)*R(i) << endl;
    }
    fp1.close();
}

/* test function to test the shoot() method for the 1 electron case
 */
void testShoot(){
    int i,n = 1e4;
    double epsilon, E, dr = RMAX/n,r;
    arrayt<double> Phi(n), R(n), Rho(n);
    for(i=0; i < n ;i++){
        Phi(i)=0;
    }
    epsilon = shoot(Phi, R, 2);
    
    calculateRho(Rho, R);
    E = calculateE(Phi, Rho, R, 2);
    cout << "Energy eigenvalue for z = 2: " <<epsilon << endl;
    ofstream fp1;
    fp1.open("testShoot.dat");
    if( fp1.fail() ) {
        cout << "cannot open file" << endl;
    }
    for(i=0;i<n;i++){
        r = i*dr;
        fp1 << setw(15) << r<< setw(15) << R(i)*R(i) << setw(15) << Rho(i) << setw(15) << Phi(i) << endl;
    }
    fp1.close();
}
/* print the arrays into an output .dat file
 */
void print(arrayt<double> &R, arrayt<double> &Rho, arrayt<double> &Phi, int count){
    int i, n=R.n1();
    double dr=RMAX/n, r;
    ofstream fp;
    string str ="sol2_"+to_string(count)+".dat";
    fp.open(str);
    if( fp.fail() ) {
        cout << "cannot open file" << endl;
    }
    for(i=0;i<n;i++){
        r=i*dr;
        fp << setw(15)<< r << setw(15) << R(i)*R(i) <<setw(15)<<Rho(i)<<setw(15)<<Phi(i)<< endl;
    }
    fp.close();
}
/* main function. intializes R(r) to a trial wavefunction and uses the
 shooting method to solve for the R(r) until the total energy convergese to a
 specified tolerance.
 */
int main(int argc, const char * argv[]) {
    int n=1e4, i, count = 0;
    double w = 0.5, tol = 1e-5, r, dr = RMAX/n, enew, eold, E, z=2.0;
    arrayt<double> Rho(n), RhoNew(n), RhoOld(n), Phi(n), R1(n);
    ofstream fp; fp.open("energy2.dat");
    if( fp.fail() ) {
        cout << "cannot open file" << endl;
    }
    cout << "---------Begin Testing--------"<<endl;
    testCalculateE();
    testCalculateRMax();
    testShoot();
    cout << "---------End Testing--------"<<endl;
    //initializes the trial wavefunction R1 and R2
    for(i=0;i<n;i++){
        r = i*dr;
        R1(i) = 2*sqrt(zstar/ABhor)*(zstar/ABhor)*r*r*exp(-zstar*r/ABhor);
    }
    normalize(R1);
    calculateRho(RhoOld, R1);
    calculatePhi(RhoOld, Phi);
    print(R1, RhoOld, Phi, count);
    E = calculateE(Phi, Rho, R1, z);
    enew = shoot(Phi, R1, z);
    do{
        calculateRho(RhoNew, R1);
        for(i = 0;i < n;i++){
            Rho(i) = RhoNew(i)*w + (1-w)*RhoOld(i);
            RhoOld(i) = Rho(i);
        }
        calculatePhi(Rho, Phi);
        fp <<setw(15)<< count << setw(15) << calculateE(Phi, Rho, R1, 2) << endl;
        eold = enew;
        enew = shoot(Phi,R1,z);
        count++;
        print(R1, Rho, Phi, count);
    }while(abs(eold - enew)>tol);
    fp.close();
    return EXIT_SUCCESS;
}
