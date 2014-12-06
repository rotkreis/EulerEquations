//
//  Matrix.h
//  MatrixClass
//
//  Created by Li Xinrui on 10/1/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __MatrixClass__Matrix__
#define __MatrixClass__Matrix__
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <cassert>

using std::vector;
using std::cout;
using std::endl;

template<class T>
class Matrix{
//  Members
protected:
    int mNumRows = 0;
    int mNumCols = 0;
    T **data;
public:
// Constructors
    Matrix():
    mNumRows(0), mNumCols(0), data(nullptr) {
    }
    Matrix(int row, int col):
    mNumRows(row), mNumCols(col), data(new T*[mNumCols]) {
        int i,j;
        for(i = 0; i != mNumCols; ++i){
            data[i] = new T[mNumRows];
        }
        for (i = 0; i != mNumCols; ++i) {
            for (j = 0; j != mNumRows; ++j) {
                data[i][j] = 0;
            }
        }
    }

    Matrix(int dim):
    mNumRows(dim), mNumCols(dim), data(new T*[dim]) {
        int i,j;
        for(i = 0; i != dim; ++i){
            data[i] = new T[dim];
        }
        for (i = 0; i != mNumCols; ++i) {
            for (j = 0; j != mNumRows; ++j) {
                data[i][j] = 0;
            }
        }
    }
    Matrix(const Matrix& m):
    mNumRows(m.mNumRows),mNumCols(m.mNumCols),data(new T*[mNumCols]){
        for(int i = 0; i != mNumCols; ++i){
            data[i] = new T[mNumRows];
        }
        for(int i = 0; i != mNumCols; ++i)
            for(int j = 0; j!= mNumRows; ++j){
                data[i][j] = m[i][j];
            }
    }

    
// Deconstructor
    ~Matrix(){
        for (int i = 0; i != this -> mNumCols; i++) {
            delete [] data[i];
        }
        delete [] data;
    }
    
public:
    // Operators
    T* operator[](int index){
        return data[index];
    }
    const T* operator[](int index) const{
        return data[index];
    }
    T& operator()(int rowIndex, int colIndex){
        return data[colIndex][rowIndex];
    }
    const T& operator() (int rowIndex, int colIndex) const{
        return data[colIndex][rowIndex];
    }
    // Members
    int GetNumRows(){
        return mNumRows;
    }
    int GetNumRows() const{
        return mNumRows;
    }
    int GetNumCols() const{
        return mNumCols;
    }
    int GetNumCols(){
        return mNumCols;
    }
    void ResizeRow(int row){
        mNumRows = row;
    }
    void ResizeCol(int col) {
        mNumCols = col;
    }

    // Operations
    Matrix<T>& operator+=(const Matrix<T>& rhs){
        if(mNumCols != rhs.mNumCols || mNumRows != rhs.mNumRows){
            std::cout << "Unequal Dimensions" << std::endl;
        }
        for(size_t i = 0; i != rhs.mNumCols; ++i){
            for(size_t j = 0; j != rhs.mNumRows; ++j){
                data[i][j] += rhs.data[i][j];
            }
        }
        return *this;
    }
    Matrix<T>& operator-=(const Matrix<T>& rhs){
        return *this += rhs * -1.0;
    }
    Matrix<T>& operator*=(const T& rhs){
        for(int i = 0; i != mNumCols; i++)
            for(int j = 0; j != mNumRows; j++) {
                data[i][j] *= rhs;
            }
        return *this;
    }
    Matrix<T>& operator/=(const T& rhs){
        return *this *= 1/rhs;
    }

    Matrix<T>& operator=(const T* rhs){
        int k = 0;
        for(int i = 0; i != mNumCols; i++)
            for(int j = 0; j != mNumRows; j++) {
                data[i][j] = rhs[k];
                k++;
            }
        return *this;
    }
    Matrix<T>& operator=(const vector<T>& rhs){
        int k = 0;
        for(int i = 0; i != mNumCols; i++)
            for(int j = 0; j != mNumRows; j++) {
                data[i][j] = rhs[k];
                k++;
            }
        return *this;
    }
    Matrix<T>& operator=(const Matrix<T>& rhs){
        int i,j;
        if(data == nullptr){
            mNumRows = rhs.mNumRows;
            mNumCols = rhs.mNumCols;
            data = new T*[rhs.mNumCols];
            for(int i = 0; i != rhs.mNumCols; ++i){
                data[i] = new T[rhs.mNumRows];
            }
        }
        if(mNumCols != rhs.mNumCols || mNumRows != rhs.mNumRows){
            cout << "unequal dimensions cannot assign matrix" << endl;
        }
        else{
            for(i = 0; i != rhs.mNumCols; i++)
                for(j = 0; j != rhs.mNumRows; j++){
                    data[i][j] = rhs[i][j];
                }
        }
        return *this;
    }


public:
    Matrix<T>& rowSwap(int i, int j){ // Mathematical Subscripts
        T temp;
        for(int k = 0; k != mNumCols; k++){
            temp = data[k][i - 1];
            data[k][i - 1] = data[k][j - 1];
            data[k][j - 1] = temp;
        }
        return *this;
    }
    Matrix<T>& colSwap(int i, int j){
        if(i != j){
            T *p = data[i-1];
            data[i-1] = data[j-1];
            data[j-1] = p;
        }
        return *this;
    }
    Matrix<T>& transpose(){
        int i, j;
        T temp;
        // This changes the Matrix itself
        if(mNumRows != mNumCols){
            Matrix copy(*this);
            for (i = 0; i != mNumCols; i++) {
                delete [] data[i];
            }
            delete [] data;
            mNumCols = copy.mNumRows;
            mNumRows = copy.mNumCols;
            data = new T*[mNumCols];
            for(i = 0; i != mNumCols; i++){
                data[i] = new T[mNumRows];
            }
            for(i = 1; i <= mNumCols; i++)
                for(j = 1; j <= mNumRows; j++){
                    data[i - 1][j - 1] = copy[j - 1][i - 1];
                }
            
            return *this;
        }
        for(i = 1; i <= mNumCols; i++)
            for(j = 1; j < i; j++){
                if(data[i - 1][j - 1] != data[j - 1][i - 1]){
                    temp = data[i - 1][j - 1];
                    data[i - 1][j - 1] = data[j - 1][i - 1];
                    data[j - 1][i - 1] = temp;
                }
            }
        return *this;
    }
    double Norm2Vec(){
        int i,j;
        double sum = 0;
        for (i = 0; i != mNumCols; i++) {
            for (j = 0; j != mNumRows; j++) {
                sum += data[i][j] * data[i][j];
            }
        }
        return sqrt(sum);
    }
    int transpose(Matrix<T>& res){
        int i, j;
        for(i = 1; i <= mNumCols; i++)
            for(j = 1; j <= mNumRows; j++){
                res[j - 1][i - 1] = data[i - 1][j - 1];
            }
        return 0;
    }
    bool QSquare(){
        assert(mNumCols == mNumRows);
        if(mNumRows == mNumCols){
            return 1;
        }
        else return 0;
    }
    void SetIdentity(){
        if( !QSquare() ){
            return;
        }
        for (int i = 0; i != this -> mNumCols; i++) {
            data[i][i] = 1;
        }
    }

public:
    friend std::ostream& operator<<(std::ostream& out, const Matrix<T>& m){
        int i, j;
        if(m.data == nullptr) cout << "Empty Matrix";
        for (i = 0; i != m.mNumRows; i++){
            for(j = 0; j!= m.mNumCols; j++){
                out << m[j][i];
            }
            out << std::endl;
        }
        return out;
    }
};
// Operations
template<class T>
Matrix<T> operator+(const Matrix<T>& m1, const Matrix<T>& m2){
    Matrix<T> sum(m1.GetNumRows(), m1.GetNumCols());
    sum = m1;
    sum += m2;
    return sum;
}
template<class T>
Matrix<T> operator-(const Matrix<T>& m1, const Matrix<T>& m2){
    Matrix<T> sum(m1.GetNumRows(), m1.GetNumCols());
    sum = m1;
    sum -= m2;
    return sum;
}
template<class T>
Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2){
    if(m1.GetNumCols() != m2.GetNumRows()) std::cout << "Illegal Multiplication" << std::endl;
    Matrix<T> res(m1.GetNumRows(), m2.GetNumCols());
    for(int i = 0; i != m1.GetNumRows(); ++i)
        for(int j = 0; j != m2.GetNumCols(); ++j)
            for(int k = 0; k != m1.GetNumCols(); ++k){
                res(i,j) += (m1(i,k) * m2(k,j));
            }
    return res;
}
template<class T>
Matrix<T> operator*(const Matrix<T>& m1, const T& num){
    Matrix<T> res(m1.GetNumRows(),m1.GetNumCols());
    res = m1;
    res *= num;
    return res;
}
template<class T>
Matrix<T> operator*(const T& num, const Matrix<T>& m1 ){
    Matrix<T> res(m1.GetNumRows(),m1.GetNumCols());
    res = m1;
    res *= num;
    return res;
}



template<class T>
Matrix<T> operator/(const Matrix<T>& m1, const T& num){
    Matrix<T> res(m1.GetNumRows(),m1.GetNumCols());
    res = m1;
    res /= num;
    return res;
}



// Algorithms

template<class T>
int lowerSolve(const Matrix<T>& l,const Matrix<T>& rhs, Matrix<T>& res) // Array or Vector? I love vector.. Without Check on Dim
{
    int n = l.GetNumRows();
    res(0,0) = rhs(0,0) / l(0,0);
    for(int i = 1; i != n; ++i){
        T temp = 0;
        for(int j = 0; j <= i - 1; ++j){
            temp += l(i,j) * res(j,0);
        }
        res(i,0) = (rhs(i,0) - temp) / l(i,i);
    }
    return 0;
}

template<class T>
int upperSolve(const Matrix<T>& u, const Matrix<T>& rhs, Matrix<T>& res){
    assert(u.GetNumRows() == u.GetNumCols());
    int n = u.GetNumRows();
    res(n - 1,0) = rhs(n - 1,0) / u(n - 1,n - 1);
    for(int i = n - 2; i >= 0; --i){
        T temp = 0;
        for(int j = i + 1; j != n; ++j){
            temp += u(i,j) * res(j,0);
        }
        res(i,0) = (rhs(i,0) - temp) / u(i,i);
    }
    return 0;
}


template<class T>
int LUDecomp(const Matrix<T>& m, Matrix<T>& l, Matrix<T>& u){ // l, u blank Matrices
    int n = m.GetNumRows();
    if(m.GetNumRows() != m.GetNumCols()){
        std::cout << "Not Square" << std::endl;
        return 0;
    }
    int i, j, k; // Copy
    for(int i = 0; i != n; ++i)
        for(int j = 0; j != n; ++j){
            u(i,j) = m(i,j);
        }
    
    for(i = 0; i != n; ++i) l(i,i) = 1;
    
    for(k = 0; k != n - 1; ++k){ // L
        for(i = k + 1; i != n; ++i){
            l(i,k) = u(i,k) / u(k,k);
        }
        for(i = k + 1; i != n; ++i){ // U
            u(i,k) = 0;
            for(j = k + 1; j != n; ++j){
                u(i,j) = u(i,j) - l(i,k) * u(k,j);
            }
        }
    }
    return 0;
}

template<class T>
int LUPivotDecomp(const Matrix<T>& m, Matrix<T>& p, Matrix<T>& l, Matrix<T>& u){ // Mathematical Subscripts
    
    int n = m.GetNumRows();
    int i, j, k;
    int *per = new int[n - 1]; // at k, exchange GetNumRows() k with GetNumRows() per[k - 1]
    T max = 0;
    T temp = 0;
    //Test
    if(m.GetNumRows() != m.GetNumCols()){
        std::cout << "Not Square" << std::endl;
        return 0;
    }
    //Copy M into U
    for(i = 0; i != n; ++i)
        for(j = 0; j != n; ++j){
            u(i,j) = m(i,j);
        }
    //Initiate L
    for(i = 0; i != n; ++i) l(i,i) = 1;
    
    //Computation
    for(k = 1; k <= n; k++){
        // Find Permutations
        max = std::abs(u(k-1,k-1));
        per[k-1] = k;
        for(i = k; i <= n; i++){
            if(std::abs(u(i-1,k-1)) > max){
                max = std::abs(u(i-1,k-1));
                per[k-1] = i;
            }
        }
        u.rowSwap(k,per[k - 1]);
        // Compute L
        for(i = k + 1; i <= n; i++){
            l(i-1,k-1) = u(i-1,k-1) / u(k-1,k-1);
        }
        // Update U
        for(i = k + 1; i <= n; i++){
            u(i - 1,k - 1) = 0;
            for(j = k + 1; j <= n; ++j){
                u(i - 1,j - 1) = u(i - 1,j - 1) - l(i - 1,k - 1) * u(k - 1,j - 1);
            }
        }
    }
    // Compute P
    for(i = 1; i <= n; i++) p(i-1,i-1) = 1;
    for(k = 1; k <= n-1; k++) p.rowSwap(k, per[k-1]);
    
    // Compute L
    for(k = 1; k != n - 1; k++){
        if(per[k] != k + 1 ){
            for(j = 1; j != k + 1; j++){  // Exchange GetNumRows() k+1 and GetNumRows() per[k] in the Lower Left Matrix<T>
                temp = l(k,j-1);
                l(k,j-1) = l( (per[k] - 1) ,j-1);
                l( (per[k] - 1) ,j-1) = temp;
            }
        }
    }
    
    return 0;
}

template<class T>
int CholeskyDecomp(const Matrix<T>& m, Matrix<T>& l){
    int i, p, k;
    int n = m.GetNumRows();
    T temp;
    
    l(0,0) = sqrt(m(0,0));
    for(i = 1; i <= n; i++){
        l(i - 1,0) = m(i - 1,0) / l(0,0);
    }
    for(k = 2; k <= n; k++){
        temp = 0;
        for(p = 1; p <= k - 1; p++){
            temp += l(k - 1,p - 1) * l(k - 1,p - 1);
        }
        if(m(k - 1,k - 1) - temp < 0){ // Negative Diagonal
            cout << "Cholesky: Diagonal cannot Sqrt at " << k << endl;
            return 1;
        }
        else{
            l(k - 1,k - 1) = sqrt(m(k - 1,k - 1) - temp);
        }
        for(i = k + 1; i <= n; i++){
            temp = 0;
            for(p = 1; p <= k - 1; p++){
                temp += l(i - 1,p - 1) * l(k - 1,p - 1);
            }
            l(i - 1,k - 1) = (m(i - 1,k - 1) - temp) / l(k-1,k-1);
        }
    }
    
    return 0;
}

template<class T>
int LDLDecomp(const Matrix<T>& m, Matrix<T>& l, Matrix<T>& d){
    int i, j, k;
    int n = m.GetNumRows();
    T temp;
    for(i = 1; i <= n; i++){
        l(i - 1,i - 1) = 1;
    }
    for(j = 1; j <= n; j++){
        T* v = new T[j - 1];
        temp = 0;
        for(k = 1; k <= j - 1; k++){
            v[k - 1] = d(k - 1,k - 1) * l(j - 1,k - 1);
            temp += l(j - 1,k - 1) * v[k - 1];
        }
        d(j - 1,j - 1) = m(j - 1,j - 1) - temp;
        for(i = j + 1; i <= n; i++){
            temp = 0;
            for(k = 1; k <= j - 1; k++){
                temp += l(i -1,k - 1) * v[k - 1];
                
            }
            l(i - 1,j - 1) = (m(i-1,j-1) - temp) / d(j-1,j-1);
        }
        delete v;
    }
    return 0;
}


template<class T>
int LUSolve(const Matrix<T>& m, const Matrix<T>& rhs, Matrix<T>& res){ // NEED RECYCLE
    if(rhs.GetNumCols() == 1){
        int n = m.GetNumRows();
        Matrix<T> l(n,n), u(n,n);
        LUDecomp(m, l, u);
        Matrix<T> temp(n,1);
        lowerSolve(l, rhs, temp);
        upperSolve(u, temp, res);
    }
    else{
        int n = m.GetNumRows();
        Matrix<T> l(n,n), u(n,n);
        LUDecomp(m, l, u);
        for(int i = 1; i <= rhs.GetNumCols(); i++){
            Matrix<T> temp1(n,1);
            Matrix<T> temp2(n,1);
            Matrix<T> right(n,1);
            for(int j = 1; j <= rhs.GetNumRows(); j++){
                right(j - 1,0) = rhs(j - 1,i - 1);
            }
            lowerSolve(l, right, temp1);
            upperSolve(u, temp1, temp2);
            for(int j = 1; j <= rhs.GetNumRows(); j++){
                res(j - 1,i - 1) = temp2(j - 1,0);
            }
        }
    }
    return 0;
}

template<class T>
int LUPivotSolve(const Matrix<T>& m, const Matrix<T>& rhs, Matrix<T>& res){ // NEED
    if(rhs.GetNumCols() == 1){
        int n = m.GetNumRows();
        Matrix<T> p(n,n), l(n,n),u(n,n);
        LUPivotDecomp(m, p, l, u);
        Matrix<T> newb(n,1);
        newb = p * rhs;
        Matrix<T> temp(n,1);
        lowerSolve(l, newb, temp);
        upperSolve(u, temp, res);
    }
    else{
        int n = m.GetNumRows();
        Matrix<T> p(n,n), l(n,n),u(n,n);
        LUPivotDecomp(m, p, l, u);
        for(int i = 1; i <= rhs.GetNumCols(); i++){
            Matrix<T> right(n,1);
            for(int j = 1; j <= rhs.GetNumRows(); j++){
                right(j-1,0) = rhs(j - 1,i - 1);
            }
            right = p * right;
            Matrix<T> temp1(n,1);
            Matrix<T> temp2(n,1);
            lowerSolve(l, right, temp1);
            upperSolve(u, temp1, temp2);
            for(int j = 1; j <= rhs.GetNumRows(); j++){
                res(j - 1,i - 1) = temp2(j - 1,0);
            }
        }
    }
    return 0;
}

template<class T>
int CholeskySolve(const Matrix<T>&m, const Matrix<T>& rhs, Matrix<T>& res){
    int n = m.GetNumRows();
    Matrix<T> y(n,1);
    Matrix<T> l(n,n);
    Matrix<T> trans(n,n);
    if(CholeskyDecomp(m, l)){
        cout << "Cholesky Solve: Cannot Decomp" << endl;
        return 1;
    }
    lowerSolve(l, rhs, y);
    l.transpose(trans);
    upperSolve(trans, y, res);
    return 0;
}
template<class T>
int LDLSolve(const Matrix<T>& m, const Matrix<T>& rhs, Matrix<T>& res){
    int n = m.GetNumRows();
    Matrix<T> y(n,1);
    Matrix<T> l(n,n), d(n,n);
    Matrix<T> trans(n,n);
    LDLDecomp(m, l, d);
    l.transpose(trans);
    lowerSolve(l, rhs, y);
    upperSolve(d * trans, y, res);
    return 0;
}
















#endif /* defined(__Matrix<T>Class__Matrix<T>__) */
