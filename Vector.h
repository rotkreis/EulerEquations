//
//  Vector.h
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __EulerEquations__Vector__
#define __EulerEquations__Vector__

#include <stdio.h>
#include "Matrix.h"


template <class T>
class mVector : public Matrix<T> {
public:
    using Matrix<T>::operator=;

public:
    int dim(){
        return this -> mNumRows;
    }
    int dim() const{
        return this -> mNumRows;
    }
public:
    T& operator[](int index){
        return Matrix<T>::data[0][index];
    }
    const T& operator[](int index) const{
        return Matrix<T>::data[0][index];
    }
public:
    mVector() {}
    mVector(int dim){
        this->mNumRows = dim;
        this->mNumCols = 1;
        this->data = new T*[1];
        this->data[0] = new T[dim];
        for (int i = 0; i != dim; i++) {
            this->data[0][i] = 0;
        }
    }
    mVector(const mVector& rhs){
        this->mNumRows = rhs.dim();
        this->mNumCols = 1;
        this->data = new T*[1];
        this->data[0] = new T[rhs.dim()];
        for (int i = 0; i != rhs.dim(); i++) {
            this->data[0][i] = rhs[i];
        }
    }
    T NormInf(){
        T max = 0;
        for (int i = 0; i != this->mNumRows; i++) {
            if (std::abs((double) this->data[0][i]) > max ) {
                max = std::abs((double) this->data[0][i]);
            }
        }
        return max;
    }
    T NormInf(int& position){
        T max = 0;
        for (int i = 0; i != this->mNumRows; i++) {
            if (std::abs((double) this->data[0][i]) > max ) {
                max = std::abs((double) this->data[0][i]);
                position = i;
            }
        }
        return max;
    }
    T Norm_2(){
        T sum = 0;
        for (int i = 0; i != Matrix<T>::mNumRows; i++) {
            sum += this -> pow(this->data[0][i],2);
        }
        return sum;
    }
    T Norm_1(){
        T max = 0;
        for (int i = 0; i!= this->mNumRows; i++) {
            if (std::abs(this->data[0][i])> max) {
                max = std::abs(this -> data[0][i]);
            }
        }
        return max;
    }
    
};
template <class T>
T InnerProduct(mVector<T>& v1, mVector<T>& v2){}


#endif /* defined(__EulerEquations__Vector__) */
