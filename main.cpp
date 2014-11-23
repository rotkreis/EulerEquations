//
//  main.cpp
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include <iostream>
#include "Matrix.h"
#include "EulerEquations.h"
// Primitive
double rho(double x){
    if (x < 0.5) {
        return 1;
    }
    else {
        return  0.125;
    }
}
double p(double x){
    if (x < 0.5) {
        return 1;
    }
    else {
        return 0.1;
    }
}
double u(double x){
    return 0;
}



int main(int argc, const char * argv[]) {
    EulerSolver sol(rho, p, u);
    sol.SetCellNumber(1000);
    sol.SetGamma(1.5);
    sol.SetRange(0, 1);
    sol.SetTime(0, 0.25);
//    sol.test();
    sol.LFSolve();
    
    
    
    
    return 0;
}
