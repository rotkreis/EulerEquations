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
#include "Vector.h"
typedef double (*pFunc)(double x);
template<class T>
class Cell: public mVector<T>{
public:
    Cell():
    mVector<T>() {}
    Cell(int n):
    mVector<T>(n) {}
public: // Inheritance
    using mVector<T>::operator=;
    using mVector<T>::operator*=;
};

template<class T>
class Profiles{ // Mathematical Subscripts, index = 1, the first element
public:
    std::vector<Cell<T> > data;
public:
    Profiles(int nCells){
        data = std::vector<Cell<T> >(nCells);
        for (int i = 0; i != nCells; i++) {
            data[i] = Cell<T> (3);
        }
    }
public: // Inheritance
    Cell<T> operator[](int index){
        return data[index - 1];
    }
public:
    T& u1(int index){
        return data[index - 1][0];
    }
    T& u2(int index){
        return data[index - 1][1];
    }
    T& u3(int index){
        return data[index - 1][2];
    }
};


template<class T>
class EulerSolver{
public:
    int nCells;
    T xStep;
    T timeStep;
    T gamma; // ratio of specific heat...
    pFunc ivDensity;
    pFunc ivPressure;
    pFunc ivVelocity;
    T startTime;
    T finalTIme;
    T xMin;
    T xMax;
    T CFL = 0.5;
    
    // Constructor
public:
    EulerSolver(pFunc rho, pFunc p, pFunc u):
    ivDensity(rho), ivPressure(p),ivVelocity(u){}
   
    // Initiate
public:
    void SetRange(T x1, T x2){
        xMin = x1;
        xMax = x2;
    }
    void SetTime(T t1, T t2){
        startTime = t1;
        finalTIme = t2;
    }
    void SetGamma(T y){
        gamma = y;
    }
    void SetCellNumber(int n){
        nCells = n;
    }
    void SetCFL(T x){
        CFL = x;
    }
    // Preparations
public:
    void ComputeSpatialStep(){
        xStep = (xMax - xMin) / nCells;
    }
    T IVAverage(int index, T x1, T x2, int n){
        assert(x2 > x1);
        T h = (x2 - x1) / (n - 1);
        T sum = 0;
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
    void InitiateValues(Profiles<T>& values){
        for (int i = 1; i <= nCells; i++) {
            values.u1(i) = IVAverage(1, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
            values.u2(i) = IVAverage(2, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
            values.u3(i) = IVAverage(3, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
        }
    }
    mVector<T> Flux(Cell<T> u){
        mVector<T> temp(3);
        temp[0] = u[1];
        temp[1] = 0.5 * (3 - gamma) * pow(u[1],2) / u[0] + (gamma - 1) * u[2];
        temp[2] = gamma * u[1] * u[2] / u[0] - 0.5 * (gamma - 1) * pow(u[1],3) / pow(u[0],2);
        return temp;
    }
    mVector<T> ComputeEigenvalues(Profiles<T>& profile, int index){
        mVector<T> eigenvalues(3);
        T rho = GetDensity(profile, index);
        T u = GetVelocity(profile, index);
        T p = GetPressure(profile, index);
        T a = sqrt(gamma * p / rho);
        eigenvalues[0] = u;
        eigenvalues[1] = u - a;
        eigenvalues[2] = u + a;
        return eigenvalues;
    }
    Matrix<T> ComputeJacobian(Profiles<T>& profile, int index){
        Matrix<T> jacobian(3,3);
        T u1 = profile.u1(index);
        T u2 = profile.u2(index);
        T u3 = profile.u3(index);
        jacobian(0,0) = 0;
        jacobian(0,1) = 1;
        jacobian(0,2) = 0;
        jacobian(1,0) = -0.5 * (gamma - 3) * pow(u2 / u1, 2);
        jacobian(1,1) = (3 - gamma) * (u2 / u1);
        jacobian(1,2) = gamma - 1;
        jacobian(2,0) = - gamma * u2 * u3 / pow(u1, 2) + (gamma - 1) * pow(u2 / u1, 3) * gamma * u3 / u1;
        jacobian(2,1) = - 3.0 / 2 * (gamma - 1) * pow(u2 / u1, 2);
        jacobian(2,2) = gamma * (u2 / u1);
        return jacobian;
    }
    
    T GetDensity(Profiles<T>& profile, int index){
        return profile.u1(index);
    }
    T GetVelocity(Profiles<T>& profile, int index) {
        return profile.u2(index) / profile.u1(index);
    }
    T GetPressure(Profiles<T>& profile, int index){
        return (gamma - 1) * (profile.u3(index) - 0.5 * pow(profile.u2(index),2) / profile.u1(index));
    }
    // Computations
public:
    T min(mVector<T>& vec){
        T min = vec[0];
        for (int i = 1; i != vec.dim(); i++) {
            if(vec[i] < min){
                min = vec[i];
            }
        }
        return min;
    }
    T min(T x1, T x2){
        if (x1 < x2) {
            return x1;
        }
        else {
            return x2;
        }
    }
    T max(T x1, T x2){
        if (x1 > x2) {
            return x1;
        }
        else{
            return x2;
        }
    }

    
    T max(mVector<T>& vec){
        T max = vec[0];
        for (int i = 1; i != vec.dim(); i++) {
            if (vec[i] > max) {
                max = vec[i];
            }
        }
        return max;
    }
// HLL
    mVector<T> HLLRightFlux(Profiles<T>& profile, int index, T dt){
        mVector<T> temp(3);
        mVector<T> FL(3);
        FL = Flux(profile[index]);
        if (index == nCells) {
            return FL;
        }
        else {
            mVector<T> eigenvaluesL(3);
            mVector<T> eigenvaluesR(3);
            eigenvaluesL = ComputeEigenvalues(profile, index);
            eigenvaluesR = ComputeEigenvalues(profile, index + 1);
            mVector<T> FR(3);
            T SL, SR;
            T min1 = min(eigenvaluesL);
            T min2 = min(eigenvaluesR);
            SL = min(min1, min2);
            T max1 = max(eigenvaluesL);
            T max2 = max(eigenvaluesR);
            SR = max(max1,max2);
            FR = Flux(profile[index + 1]);
            if (SL > 0) {
                temp = FL;
            }
            if (SL < 0 && SR > 0) {
                temp = (SR * FL - SL * FR + SL * SR * (profile[index + 1] - profile[index]))
                     / (SR - SL);
            }
            if (SR < 0) {
                temp = FR;
            }
        }
        return temp;
    }
    mVector<T> HLLLeftFlux(Profiles<T>& profile, int index, T dt){
        mVector<T> temp(3);

        mVector<T> FR(3);
        FR = Flux(profile[index]);
        if (index == 1) {
            return FR;
        }
        else {
            mVector<T> FL(3);
            FL = Flux(profile[index - 1]);
            mVector<T> eigenvaluesL(3);
            mVector<T> eigenvaluesR(3);
            eigenvaluesL = ComputeEigenvalues(profile, index - 1);
            eigenvaluesR = ComputeEigenvalues(profile, index);
            T SL, SR;
            T min1 = min(eigenvaluesL);
            T min2 = min(eigenvaluesR);
            SL = min(min1, min2);
            T max1 = max(eigenvaluesL);
            T max2 = max(eigenvaluesR);
            SR = max(max1, max2);
            if (SL > 0) {
                temp = FL;
            }
            if (SL < 0 && SR > 0) {
                temp = (SR * FL - SL * FR + SL * SR * (profile[index] - profile[index - 1]))
                / (SR - SL);
            }
            if (SR < 0) {
                temp = FR;
            }
        }
        return temp;
    }
    T HLLComputeTimeStep(Profiles<T>& profile, T atTime){
        return LFComputeTimeStep(profile, atTime);
    }
    void HLLComputeForward(Profiles<T>& uPre, Profiles<T>& uPost, T dt){
        for (int i = 1; i <= nCells; i++) {
            uPost[i] = uPre[i] - dt / xStep * (HLLRightFlux(uPre, i, dt) - HLLLeftFlux(uPre, i, dt));
        }
    }
    void HLLSolve(Profiles<T>& res) {
        ComputeSpatialStep();
        Profiles<T> uPre(nCells);
        Profiles<T> uPost(nCells);
        InitiateValues(uPost);
        T tNow = 0;
        while (startTime + tNow < finalTIme) {
            uPre = uPost;
            T dt = HLLComputeTimeStep(uPre, tNow);
            HLLComputeForward(uPre, uPost, dt);
            tNow += dt;
        }
        GetOutput(uPost, res);
        std::cout << "HLL finished!" << std::endl;
    }
// LW
    mVector<T> LWRightFlux(Profiles<T>& profile, int index, T dt){
        mVector<T> temp(3);
        if (index == nCells) {
            temp = Flux(profile[index]);
            return temp;
        }
        temp = 0.5 * (Flux(profile[index]) + Flux(profile[index + 1])) - 0.5 * dt / xStep * ComputeJacobian(profile, index) * ComputeJacobian(profile, index + 1) *(profile[index + 1] - profile[index]);
        return temp;
    }
    mVector<T> LWLeftFlux(Profiles<T>& profile, int index, T dt){
        mVector<T> temp(3);
        if (index == 1) {
            temp = Flux(profile[index]);
            return temp;
        }
        temp = 0.5 * (Flux(profile[index - 1]) + Flux(profile[index])) - 0.5 * dt / xStep * ComputeJacobian(profile, index) * ComputeJacobian(profile, index - 1) *(profile[index] - profile[index - 1]);
        return temp;
    }
    T LWComputeTimeStep(Profiles<T>& profile, T atTime){
        return LFComputeTimeStep(profile, atTime);
    }
    void LWComputeForward(Profiles<T>& uPre, Profiles<T>& uPost, T dt){
        for (int i = 1; i <= nCells; i++) {
            uPost[i] = uPre[i] - dt / xStep * (LWRightFlux(uPre, i, dt) - LWLeftFlux(uPre, i, dt));
        }
    }
    void LWSolve(Profiles<T>& res){
        ComputeSpatialStep();
        Profiles<T> uPre(nCells);
        Profiles<T> uPost(nCells);
        InitiateValues(uPost);
        T tNow = 0;
        while (startTime + tNow < finalTIme) {
            uPre = uPost;
            T dt = LWComputeTimeStep(uPre, tNow);
            LWComputeForward(uPre, uPost, dt);
            tNow += dt;
        }
        GetOutput(uPost, res);
        std::cout << "LW Finished!" << std::endl;
    }
//LF
    mVector<T> LFRightFlux(Profiles<T>& u, int index, T dt){
        mVector<T> temp(3);
        if (index == nCells) { // Boundary
            temp = Flux(u[index]);
            return temp;
        }
        temp = 0.5 * (Flux(u[index]) + Flux(u[index + 1])) - (u[index + 1] - u[index]) / (2 * dt / xStep);
        return temp;
    }
    mVector<T> LFLeftFlux(Profiles<T>& u, int index, T dt){
        mVector<T> temp(3);
        if (index == 1) { // Boundary
            temp = Flux(u[index]);
            return temp;
        }
        temp = 0.5 * (Flux(u[index - 1]) + Flux(u[index])) - (u[index] - u[index - 1]) / (2 * dt / xStep);
        return temp;
    }
    
    void LFComputeForward(Profiles<T>& uPre, Profiles<T>& uPost, T dt){
        for (int i = 1; i <= nCells; i++) {
            uPost[i] = uPre[i] - dt / xStep * (LFRightFlux(uPre, i, dt) - LFLeftFlux(uPre, i, dt));
        }
    }
    T LFComputeTimeStep(Profiles<T>& profile, T atTime){
        T tMin = 0.1;
        for (int i = 1 ; i <= nCells; i++) {
            mVector<T> eigenvalues(3);
            eigenvalues = ComputeEigenvalues(profile, i);
            T t1 = CFL * xStep / std::abs(eigenvalues[0]);
            T t2 = CFL * xStep / std::abs(eigenvalues[1]);
            T t3 = CFL * xStep / std::abs(eigenvalues[2]);
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
    void LFSolve(Profiles<T>& res){
        ComputeSpatialStep();
        Profiles<T> uPre(nCells);
        Profiles<T> uPost(nCells);
        InitiateValues(uPost);
        T tNow = 0;
        while (startTime + tNow < finalTIme) {
            uPre = uPost;
            T dt = LFComputeTimeStep(uPre, tNow);
            LFComputeForward(uPre, uPost, dt);
            tNow += dt;
        }
        GetOutput(uPost, res);
        std::cout << "LF Finished!" <<endl;
    }
    //FORCE
    mVector<T> FORCERightFlux(Profiles<T>& u, int index, T dt){
        mVector<T> temp(3);
        temp = .5 * (LFRightFlux(u, index, dt) + LWRightFlux(u, index, dt));
        return temp;
    }
    mVector<T> FORCELeftFlux(Profiles<T>& u, int index, T dt){
        mVector<T> temp(3);
        temp = .5 * (LFLeftFlux(u, index, dt) + LWLeftFlux(u, index, dt));
        return temp;
    }
    void FORCEComputeForward(Profiles<T>& uPre, Profiles<T>& uPost, T dt){
        for (int i = 1; i <= nCells; i++) {
            uPost[i] = uPre[i] - dt / xStep * (FORCERightFlux(uPre, i, dt) - FORCELeftFlux(uPre, i, dt));
        }
    }
    T FORCEComputeTimeStep(Profiles<T>& profile, T atTime){
        return LFComputeTimeStep(profile, atTime);
    }
    void FORCESolve(Profiles<T>& res){
        ComputeSpatialStep();
        Profiles<T> uPre(nCells);
        Profiles<T> uPost(nCells);
        InitiateValues(uPost);
        T tNow = 0;
        while (startTime + tNow < finalTIme) {
            uPre = uPost;
            T dt = FORCEComputeTimeStep(uPre, tNow);
            FORCEComputeForward(uPre, uPost, dt);
            tNow += dt;
        }
        GetOutput(uPost, res);
        std::cout << "FORCE Finished!" <<endl;
    }

    
    void GetOutput(Profiles<T>& uPost, Profiles<T>& res){
        for (int i = 1; i <= nCells; i++) {
            res.u1(i) = GetDensity(uPost, i);
            res.u2(i) = GetVelocity(uPost, i);
            res.u3(i) = GetPressure(uPost, i);
        }
    }
    
};



    

#endif /* defined(__EulerEquations__EulerEquations__) */
