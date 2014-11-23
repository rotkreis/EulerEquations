//
//  EulerEquations.h
//  EulerEquations
//
//  Created by Li Xinrui on 11/23/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __EulerEquations__EulerEquations__
#define __EulerEquations__EulerEquations__

#include <stdio.h>
#include <cassert>
#include "Matrix.h"
typedef double (*pFunc)(double x);

class Cell : public mVector{
public:
    Cell(int n):
    mVector(n){ }
};


class Profiles{ // Mathematical Subscripts, index = 1, the first element
public:
    std::vector<mVector> data;
public:
    Profiles(int nCells){
        data = std::vector<mVector>(nCells);
        for (int i = 0; i != nCells; i++) {
            data[i] = mVector(3);
        }
    }
public:
    double& u1(int index){
        return data[index - 1][0];
    }
    double& u2(int index){
        return data[index - 1][1];
    }
    double& u3(int index){
        return data[index - 1][1];
    }
};



class EulerSolver{
public:
    int nCells;
    double xStep;
    double timeStep;
    double gamma; // ratio of specific heat...
    pFunc ivDensity;
    pFunc ivPressure;
    pFunc ivVelocity;
    double startTime;
    double finalTIme;
    double xMin;
    double xMax;
    double CFI = 0.5;
    
    // Constructor
public:
    EulerSolver(pFunc rho, pFunc p, pFunc u):
    ivDensity(rho), ivPressure(p),ivVelocity(u){}
    
    // Initiate
public:
    void SetRange(double x1, double x2){
        xMin = x1;
        xMax = x2;
    }
    void SetTime(double t1, double t2){
        startTime = t1;
        finalTIme = t2;
    }
    void SetGamma(double y){
        gamma = y;
    }
    void SetCellNumber(int n){
        nCells = n;
    }
    // Preparations
public:
    void ComputeSpatialStep(){
        xStep = (xMax - xMin) / nCells;
    }
    double IVAverage(int index, double x1, double x2, int n){
        assert(x2 > x1);
        double h = (x2 - x1) / (n - 1);
        double sum = 0;
        switch (index) {
            case 1:
                for (int i = 0; i <= n - 1; i++) {
                    sum += ivDensity(x1 + i * h);
                }
                break;
            case 2:
                for (int i = 0; i <= n - 1; i++) {
                    sum += ivDensity(x1 + i * h) * ivVelocity(x1 + i * h);
                }
                break;
            case 3:
                for (int i = 0; i <= n - 1; i++) {
                    sum += 0.5 * ivDensity(x1 + i * h) * pow(ivVelocity(x1 + i * h), 2)
                    + ivPressure(x1 + i * h) / (gamma - 1);
                }
            default:
                break;
        }
        return sum / n;
    }
    void InitiateValues(Profiles& values){
        for (int i = 1; i <= nCells; i++) {
            values.u1(i) = IVAverage(1, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
            values.u2(i) = IVAverage(2, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
            values.u3(i) = IVAverage(3, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
        }
    }
    void test(){
        ComputeSpatialStep();
        Profiles values(nCells);
        InitiateValues(values);
        for (int i = 1; i <= nCells; i++) {
            std::cout << values.u1(i) << endl;
        }
    }
    // Computations
public:
    
    
    
    
    
    
    
};

#endif /* defined(__EulerEquations__EulerEquations__) */
