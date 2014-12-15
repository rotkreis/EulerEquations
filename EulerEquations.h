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
    mVector ComputeEigenvalues(Cell& u){
        mVector eigenvalues(3);
        double rho = GetDensity(profile, index);
        double u = GetVelocity(profile, index);
        double p = GetPressure(profile, index);
        double a = sqrt(gamma * p / rho);
        eigenvalues[0] = u;
        eigenvalues[1] = u - a;
        eigenvalues[2] = u + a;
        return eigenvalues;
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
    double min(mVector& vec){
        double min = vec[0];
        for (int i = 1; i != vec.dim(); i++) {
            if(vec[i] < min){
                min = vec[i];
            }
        }
        return min;
    }
    double min(double x1, double x2){
        if (x1 < x2) {
            return x1;
        }
        else {
            return x2;
        }
    }
    double max(double x1, double x2){
        if (x1 > x2) {
            return x1;
        }
        else{
            return x2;
        }
    }

    
    double max(mVector& vec){
        double max = vec[0];
        for (int i = 1; i != vec.dim(); i++) {
            if (vec[i] > max) {
                max = vec[i];
            }
        }
        return max;
    }
// HLL
//    mVector HLLRightFlux(Profiles& profile, int index, double dt){
//        mVector temp(3);
//        mVector FL(3);
//        FL = Flux(profile[index]);
//        if (index == nCells) {
//            return FL;
//        }
//        else {
//            mVector eigenvaluesL(3);
//            mVector eigenvaluesR(3);
//            eigenvaluesL = ComputeEigenvalues(profile, index);
//            eigenvaluesR = ComputeEigenvalues(profile, index + 1);
//            mVector FR(3);
//            double SL, SR;
//            double min1 = min(eigenvaluesL);
//            double min2 = min(eigenvaluesR);
//            SL = min(min1, min2);
//            double max1 = max(eigenvaluesL);
//            double max2 = max(eigenvaluesR);
//            SR = max(max1,max2);
//            FR = Flux(profile[index + 1]);
//            if (SL > 0) {
//                temp = FL;
//            }
//            if (SL < 0 && SR > 0) {
//                temp = (SR * FL - SL * FR + SL * SR * (profile[index + 1] - profile[index]))
//                     / (SR - SL);
//            }
//            if (SR < 0) {
//                temp = FR;
//            }
//        }
//        return temp;
//    }
//    mVector HLLLeftFlux(Profiles& profile, int index, double dt){
//        mVector temp(3);
//        mVector FR(3);
//        FR = Flux(profile[index]);
//        if (index == 1) {
//            return FR;
//        }
//        else {
//            mVector FL(3);
//            FL = Flux(profile[index - 1]);
//            mVector eigenvaluesL(3);
//            mVector eigenvaluesR(3);
//            eigenvaluesL = ComputeEigenvalues(profile, index - 1);
//            eigenvaluesR = ComputeEigenvalues(profile, index);
//            double SL, SR;
//            double min1 = min(eigenvaluesL);
//            double min2 = min(eigenvaluesR);
//            SL = min(min1, min2);
//            double max1 = max(eigenvaluesL);
//            double max2 = max(eigenvaluesR);
//            SR = max(max1, max2);
//            if (SL > 0) {
//                temp = FL;
//            }
//            if (SL < 0 && SR > 0) {
//                temp = (SR * FL - SL * FR + SL * SR * (profile[index] - profile[index - 1]))
//                / (SR - SL);
//            }
//            if (SR < 0) {
//                temp = FR;
//            }
//        }
//        return temp;
//    }
//    void HLLComputeForward(Profiles& uPre, Profiles& uPost, double dt){
//        for (int i = 1; i <= nCells; i++) {
//            uPost[i] = uPre[i] - dt / xStep * (HLLRightFlux(uPre, i, dt) - HLLLeftFlux(uPre, i, dt));
//        }
//    }
//    void HLLSolve(Profiles& res) {
//        ComputeSpatialStep();
//        Profiles uPre(nCells);
//        Profiles uPost(nCells);
//        InitiateValues(uPost);
//        double tNow = 0;
//        while (startTime + tNow < finalTIme) {
//            uPre = uPost;
//            double dt = ComputeTimeStep(uPre, tNow);
//            HLLComputeForward(uPre, uPost, dt);
//            tNow += dt;
//        }
//        GetOutput(uPost, res);
//        std::cout << "HLL finished!" << std::endl;
//    }
// LW
//    mVector LWRightFlux(Profiles& profile, int index, double dt){
//        mVector temp(3);
//        if (index == nCells) {
//            temp = Flux(profile[index]);
//            return temp;
//        }
//        return temp;
//    }
//    mVector LWLeftFlux(Profiles& profile, int index, double dt){
//        mVector temp(3);
//        if (index == 1) {
//            temp = Flux(profile[index]);
//            return temp;
//        }
//        return temp;
//    }
//    void LWComputeForward(Profiles& uPre, Profiles& uPost, double dt){
//        for (int i = 1; i <= nCells; i++) {
//            uPost[i] = uPre[i] - dt / xStep * (LWRightFlux(uPre, i, dt) - LWLeftFlux(uPre, i, dt));
//        }
//    }
//    void LWSolve(Profiles& res){
//        ComputeSpatialStep();
//        Profiles uPre(nCells);
//        Profiles uPost(nCells);
//        InitiateValues(uPost);
//        double tNow = 0;
//        while (startTime + tNow < finalTIme) {
//            uPre = uPost;
//            double dt = ComputeTimeStep(uPre, tNow);
//            LWComputeForward(uPre, uPost, dt);
//            tNow += dt;
//        }
//        GetOutput(uPost, res);
//        std::cout << "LW Finished!" << std::endl;
//    }
//LF
    mVector LFFlux(Cell& left, Cell& right, double dt){
        return 0.5 * (Flux(left) + Flux(right))-(right - left) / (2 * dt / xStep);
    }
    
    void LFComputeForward(Profiles& uPre, Profiles& uPost, double dt){
        for (int i = 2; i <= nCells - 1; i++) {
            uPost[i] = uPre[i] - dt / xStep * (LFFlux(uPre[i], uPre[i+1], dt) - LFFlux(uPre[i - 1], uPre[i], dt));
        }
        uPost[1] = uPre[1] - dt / xStep * (LFFlux(uPre[1], uPre[2], dt) - LFFlux(uPre[1], uPre[1], dt));
        uPost[nCells] = uPre[1] - dt/xStep * (LFFlux(uPre[nCells - 1], uPre[nCells], dt) - LFFlux(uPre[nCells], uPre[nCells], dt));
    }
    double ComputeTimeStep(Profiles& profile, double atTime){
        double tMin = 0.1;
        for (int i = 1 ; i <= nCells; i++) {
            mVector eigenvalues(3);
            eigenvalues = ComputeEigenvalues(profile, i);
            double t1 = CFL * xStep / std::abs(eigenvalues[0]);
            double t2 = CFL * xStep / std::abs(eigenvalues[1]);
            double t3 = CFL * xStep / std::abs(eigenvalues[2]);
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
    void LFSolve(Profiles& res){
        ComputeSpatialStep();
        Profiles uPre(nCells);
        Profiles uPost(nCells);
        InitiateValues(uPost);
        double tNow = 0;
        while (startTime + tNow < finalTIme) {
            uPre = uPost;
            double dt = ComputeTimeStep(uPre, tNow);
            LFComputeForward(uPre, uPost, dt);
            tNow += dt;
        }
        GetOutput(uPost, res);
        std::cout << "LF Finished!" <<endl;
    }
    //FORCE
//    mVector FORCERightFlux(Profiles& u, int index, double dt){
//        mVector temp(3);
//        temp = .5 * (LFRightFlux(u, index, dt) + LWRightFlux(u, index, dt));
//        return temp;
//    }
//    mVector FORCELeftFlux(Profiles& u, int index, double dt){
//        mVector temp(3);
//        temp = .5 * (LFLeftFlux(u, index, dt) + LWLeftFlux(u, index, dt));
//        return temp;
//    }
//    void FORCEComputeForward(Profiles& uPre, Profiles& uPost, double dt){
//        for (int i = 1; i <= nCells; i++) {
//            uPost[i] = uPre[i] - dt / xStep * (FORCERightFlux(uPre, i, dt) - FORCELeftFlux(uPre, i, dt));
//        }
//    }
//    void FORCESolve(Profiles& res){
//        ComputeSpatialStep();
//        Profiles uPre(nCells);
//        Profiles uPost(nCells);
//        InitiateValues(uPost);
//        double tNow = 0;
//        while (startTime + tNow < finalTIme) {
//            uPre = uPost;
//            double dt = ComputeTimeStep(uPre, tNow);
//            FORCEComputeForward(uPre, uPost, dt);
//            tNow += dt;
//        }
//        GetOutput(uPost, res);
//        std::cout << "FORCE Finished!" <<endl;
//    }
    void GetOutput(Profiles& uPost, Profiles& res){
        for (int i = 1; i <= nCells; i++) {
            res.u1(i) = GetDensity(uPost, i);
            res.u2(i) = GetVelocity(uPost, i);
            res.u3(i) = GetPressure(uPost, i);
            }
        }
    };



    

#endif /* defined(__EulerEquations__EulerEquations__) */
