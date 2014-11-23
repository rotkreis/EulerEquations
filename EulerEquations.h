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
    Cell():
    mVector() {}
    Cell(int n):
    mVector(n) {}
public: // Inheritance
    using mVector::operator=;
    using mVector::operator*=;
};

class Profiles{ // Mathematical Subscripts, index = 1, the first element
public:
    std::vector<Cell> data;
public:
    Profiles(int nCells){
        data = std::vector<Cell>(nCells);
        for (int i = 0; i != nCells; i++) {
            data[i] = Cell(3);
        }
    }
public: // Inheritance
    Cell& operator[](int index){
        return data[index - 1];
    }
public:
    double& u1(int index){
        return data[index - 1][0];
    }
    double& u2(int index){
        return data[index - 1][1];
    }
    double& u3(int index){
        return data[index - 1][2];
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
    double CFL = 0.5;
    
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
    void SetCFL(double x){
        CFL = x;
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
                break;
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
    mVector Flux(Cell& u){
        mVector temp(3);
        temp[0] = u[1];
        temp[1] = 0.5 * (3 - gamma) * pow(u[1],2) / u[0] + (gamma - 1) * u[2];
        temp[2] = gamma * u[1] * u[2] / u[0] - 0.5 * (gamma - 1) * pow(u[1],3) / pow(u[0],2);
        return temp;
    }
    void test(){
        ComputeSpatialStep();
        Profiles values(nCells);
        InitiateValues(values);
        for (int i = 1; i <= nCells; i++) {
            std::cout << values[i] << endl;
        }
    }
    double GetDensity(Profiles& profile, int index){
        return profile.u1(index);
    }
    double GetVelocity(Profiles& profile, int index) {
        return profile.u2(index) / profile.u1(index);
    }
    double GetPressure(Profiles& profile, int index){
        return (gamma - 1) * (profile.u3(index) - 0.5 * pow(profile.u2(index),2) / profile.u1(index));
    }
    // Computations
public:
    mVector RightFlux(Profiles& u, int index, double dt){
        mVector temp(3);
        if (index == nCells) { // Boundary
            temp = Flux(u[index]);
            return temp;
        }
        temp = 0.5 * (Flux(u[index]) + Flux(u[index + 1])) - (u[index + 1] - u[index]) / (2 * dt / xStep);
        return temp;
    }
    mVector LeftFlux(Profiles& u, int index, double dt){
        mVector temp(3);
        if (index == 1) { // Boundary
            temp = Flux(u[index]);
            return temp;
        }
        temp = 0.5 * (Flux(u[index - 1]) + Flux(u[index])) - (u[index] - u[index - 1]) / (2 * dt / xStep);
        return temp;
    }
    
    void LFComputeForward(Profiles& uPre, Profiles& uPost, double dt){
        for (int i = 1; i <= nCells; i++) {
            uPost[i] = uPre[i] - dt / xStep * (RightFlux(uPre, i, dt) - LeftFlux(uPre, i, dt));
        }
    }
    double LFComputeTimeStep(Profiles& profile, double atTime){
        double tMin = 0.1;
        for (int i = 1 ; i != nCells; i++) {
            double rho = GetDensity(profile, i);
            double u = GetVelocity(profile, i);
            double p = GetPressure(profile, i);
            double a = sqrt(gamma * p / rho);
            double t1 = CFL * xStep / std::abs(u);
            double t2 = CFL * xStep / std::abs(u + a);
            double t3 = CFL * xStep / std::abs(u - a);
            if (t1 < tMin) {
                tMin = t1;
            }
            if (t2 < tMin) {
                tMin = t2;
            }
            if (t3 < tMin) {
                tMin = t3;
            }
        }
        if (atTime + tMin > finalTIme) {
            tMin = finalTIme - atTime;
        }
        return tMin;
    }
    void LFSolve(){
        ComputeSpatialStep();
        Profiles uPre(nCells);
        Profiles uPost(nCells);
        InitiateValues(uPost);
        double tNow = 0;
        while (startTime + tNow < finalTIme) {
            uPre = uPost;
            double dt = LFComputeTimeStep(uPre, tNow);
            LFComputeForward(uPre, uPost, dt);
            tNow += dt;
        }
        for (int i = 1; i <= nCells; i++) {
            std::cout << uPost.u1(i) <<", ";
        }
        std::cout << "finished!";
    }
    
    
    
    
    
};

#endif /* defined(__EulerEquations__EulerEquations__) */
