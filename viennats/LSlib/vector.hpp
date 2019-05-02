#ifndef VECTOR_HPP_
#define VECTOR_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <cmath>
#include <algorithm>
#include <functional>
#include "math.hpp"
#include "message.h"

namespace lvlset {

    template <class T,int D> class vec {
        T x[D];
    public:

        typedef T value_type;

        vec(const vec<T,D>&);
        template <class T1, int D1> explicit vec(const vec<T1,D1>&);
        explicit vec(const T[]);

        vec();

        explicit vec(T x0,T x1,T x2) {
            x[0]=x0;
            x[1]=x1;
            if(D==3) x[2]=x2;
        }
        explicit vec(T x0,T x1) {
            x[0]=x0;
            x[1]=x1;
        }

        explicit vec(T);
        template<class X> explicit vec(X);

        template <class V> vec& operator-=(const V&);
        template <class V> vec& operator+=(const V&);
        vec& operator*=(T);
        vec& operator/=(T);


        template <class V> vec<T,D> operator-(const V&) const;
        vec<T,D> operator-() const;

        template <class V> vec<T,D> operator+(const V&) const;

        template <class V> bool operator==(const V&) const;
        template <class V> bool operator!=(const V&) const;
        template <class V> bool operator<(const V&) const;
        template <class V> bool operator<=(const V&) const;
        template <class V> bool operator>(const V&) const;
        template <class V> bool operator>=(const V&) const;

        vec<T,D> operator* (T) const;

        vec<T,D> operator/(T) const;

        T& operator[] (int);
        const T& operator[] (int) const;

        template<class V> vec<T,D>& operator=(const V&);

        T size() const;
        T element_max() const;
        void sort();
        void reverse_sort();

    };

    template <class T, int D> inline vec<T,D>::vec() {}

    template <class T, int D> template <class T1, int D1> inline vec<T,D>::vec(const vec<T1,D1>& v) {
        const int k=std::min(D,D1);
        for (int i=0;i<k;i++) x[i]=v[i];
        for (int i=k;i<D;i++) x[i]=T(0);
    }

    template <class T,int D> inline vec<T,D>::vec(const T v[]) {
        for (int i=0;i<D;i++) x[i]=v[i];
    }

    template <class T, int D> template <class V> inline  vec<T,D>& vec<T,D>::operator=(const V& v) {
        for (int i=0;i<D;i++) x[i]=v[i];
        return *this;
    }

    template <class T, int D> inline  vec<T,D>::vec(const vec<T,D>& v) {
        for (int i=0;i<D;i++) x[i]=v[i];
    }

    template <class T, int D> inline T& vec<T,D>::operator[] (int i) {
        return x[i];
    }

    template <class T, int D> inline const T& vec<T,D>::operator[] (int i) const {
        return x[i];
    }

    template <class T, int D> template <class V> inline vec<T,D>& vec<T,D>::operator-=(const V& v) {
        for (int i=0;i<D;i++) x[i]-=v[i];
        return *this;
    }

    template <class T, int D> template <class V> inline vec<T,D>& vec<T,D>::operator+=(const V& v) {
        for (int i=0;i<D;i++) x[i]+=v[i];
        return *this;
    }

    template <class T, int D> inline vec<T,D>& vec<T,D>::operator*=(T d) {
        for (int i=0;i<D;i++) x[i]*=d;
        return *this;
    }

    template <class T, int D> inline vec<T,D>& vec<T,D>::operator/=(T d) {
        for (int i=0;i<D;i++) x[i]/=d;
        return *this;
    }

    template <class T, int D> template<class V> inline vec<T,D> vec<T,D>::operator-(const V& v) const{
        vec<T,D> w;
        for (int i=0;i<D;i++) w.x[i]=x[i]-v[i];
        return w;
    }

    template <class T, int D> template<class V> inline vec<T,D> vec<T,D>::operator+(const V& v) const{
        vec<T,D> w;
        for (int i=0;i<D;i++) w.x[i]=x[i]+v[i];
        return w;
    }

    template <class T, int D> inline vec<T,D> vec<T,D>::operator-() const {
        vec<T,D> v;
        for (int i=0;i<D;i++) v.x[i]=-x[i];
        return v;
    }

    template <class T, int D> inline vec<T,D> vec<T,D>::operator*(T d) const{
        vec<T,D> v;
        for (int i=0;i<D;i++) v.x[i]=x[i]*d;
        return v;
    }


    template <class T, int D> inline vec<T,D> vec<T,D>::operator/(T d) const{
        vec<T,D> v;
        for (int i=0;i<D;i++) v.x[i]=x[i]/d;
        return v;
    }

    template <class T, int D> inline vec<T,D>::vec(T d) {
        for (int i=0;i<D;i++) x[i]=d;
    }

    template <class T, int D> template <class X> inline vec<T,D>::vec(X d) {
        for (int i=0;i<D;i++) x[i]=d[i];
    }

    template <class T, int D> template<class V> inline bool vec<T,D>::operator==(const V& v) const{
        for (int i=D-1;i>=0;i--) if (x[i]!=v[i]) return false;
        return true;
    }

    template <class T, int D> template<class V> inline bool vec<T,D>::operator!=(const V& v) const{
        for (int i=D-1;i>=0;i--) if (x[i]!=v[i]) return true;
        return false;
    }


    template <class T,int D> template<class V> inline bool vec<T,D>::operator<(const V& v) const {
        for (int i=D-1;i>=0;i--) {
            if (x[i]<v[i])
                return true;
            else if (x[i]>v[i])
                return false;
        }
        return false;
    }

    template <class T,int D> template<class V> inline bool vec<T,D>::operator>(const V& v) const {
        for (int i=D-1;i>=0;i--) {
            if (x[i]>v[i])
                return true;
            else if (x[i]<v[i])
                return false;
        }
        return false;
    }

    template <class T,int D> template<class V> inline bool vec<T,D>::operator<=(const V& v) const {
        return !(*this>v);
    }

    template <class T,int D> template<class V> inline bool vec<T,D>::operator>=(const V& v) const {
        return !(*this<v);
    }

    template <class T, int D> inline T vec<T,D>::size() const{
        return (sizeof(x)/sizeof(*x));
    }

   template <class T, int D> inline T vec<T,D>::element_max() const{
        return *std::max_element(std::begin(x), std::end(x));
    }

    //###################################################################

    template <class T, int D> inline vec<T,D> Min(const vec<T,D>& v1,const vec<T,D>& v2) {
        vec<T,D> v;
        for (int i=0;i<D;i++) v[i]=std::min(v1[i],v2[i]);
        return v;
    }

    template <class T, int D> inline vec<T,D> Max(const vec<T,D>& v1,const vec<T,D>& v2) {
        vec<T,D> v;
        for (int i=0;i<D;i++) v[i]=std::max(v1[i],v2[i]);
        return v;
    }

    template <class T, int D> inline T Element_Max(const vec<T,D>& v1,const vec<T,D>& v2) {
        vec<T,D> v;
        for (int i=0;i<D;i++) v[i]=std::max(v1[i],v2[i]);
        return v;
    }


    template <class T, int D> inline T dot(const vec<T,D>& v1, const vec<T,D>& v2) {
        T sum(0);
        for (int i=0;i<D;i++) sum+=v1[i]*v2[i];
        return sum;
    }

    // template <class V, class T, int D> inline bool operator==(const V& v0,const vec<T,D>& v1) {
    //     return (v1.operator==(v0));
    // }

    // template <class V, class T, int D> inline bool operator!=(const V& v0,const vec<T,D>& v1) {
    //     return (v1.operator!=(v0));
    // }

    template <class T, int D> inline vec<T,D> operator*(T d, const vec<T,D>& v) {
        return v*d;
    }

    template <class T> inline vec<T,2> RotateLeft(const vec<T,2>& v) {
        return vec<T,2>(-v[1],v[0]);
    }

    template <class T> inline vec<T,2> RotateRight(const vec<T,2>& v) {
        return vec<T,2>(v[1],-v[0]);
    }

    template <class T> inline vec<T,3> cross(const vec<T,3>& v1,const vec<T,3>& v2) {
        return vec<T,3>(    v1[1]*v2[2]-v1[2]*v2[1],
                            v1[2]*v2[0]-v1[0]*v2[2],
                            v1[0]*v2[1]-v1[1]*v2[0]);
    }

    template <class T> inline T cross(const vec<T,2>& v1,const vec<T,2>& v2) {
        return v1[0]*v2[1]-v1[1]*v2[0];
    }

    template <int D, class T> inline vec<T,D> Normalize(const vec<T,D>& v) {
        T n=Norm(v);
        if (n<=0.) return vec<T,D>(T(0));
        return v/n;
    }

    template <int D, class T> inline T Norm(const vec<T,D>& v) {
        T max=math::abs(v[0]);
        for (int i=1;i<D;i++) max=std::max(max,math::abs(v[i]));
        if (max==0.) return T(0.);
        T d=0.;
        for (int i=0;i<D;i++) {
            T t(v[i]/max);
            d+=t*t;
        }
        return max*std::sqrt(d);
    }

    template <int D, class T> inline T Norm2(const vec<T, D>& v) { //squared l2 norm TODO name is misleading
        return dot(v,v);
    }

    template <int D, class T> inline T NormL2(const vec<T, D>& v) {//l2 norm
        return std::sqrt(Norm2(v));
    }

    template <int D, class T> inline int ManhattanNorm(const vec<T,D>& v) {
        T k(0);
        for (int i=0;i<D;i++) k+=math::abs(v[i]);
        return k;
    }

    // template <class V, class T,int D> inline bool operator<(const V& v0,const vec<T,D>& v1) {
    //     return (v1.operator>(v0));
    // }

    // template <class V, class T,int D> inline bool operator>(const V& v0,const vec<T,D>& v1) {
    //     return (v1.operator<(v0));
    // }

    // template <class V, class T,int D> inline bool operator<=(const V& v0,const vec<T,D>& v1) {
    //     return (v1.operator>=(v0));
    // }

    // template <class V, class T,int D> inline bool operator>=(const V& v0,const vec<T,D>& v1) {
    //     return (v1.operator<=(v0));
    // }

    // template <class V, class T,int D> inline vec<T,D> operator+(const V& v0,const vec<T,D>& v1) {
    //     return v1.operator+(v0);
    // }

    // template <class V, class T,int D> inline vec<T,D> operator-(const V& v0,const vec<T,D>& v1) {
    //     return (-v1).operator+(v0);
    // }

    template <class T> inline T Volume(const vec<T,2>* p) {
        return  ((p[1]-p[0])%(p[2]-p[0]));
    }

    template <class T> inline T Volume(const vec<T,3>* p) {
        return ((p[1]-p[0])*((p[2]-p[0])%(p[3]-p[0])));
    }

    /*template <class T, int D> inline T Volume(const vec<T,D>* p) {
        return Det(mat::mat<T,D>(p));
    }*/

    /*template <class T> inline vec<T,2> NormalVector(const vec<T,2>* v) {
        return RotateLeft(Normalize(v[1]-v[0]));
    }

    template <class T> inline vec<T,3> NormalVector(const vec<T,3>* c) {        //TODO

        int max_e(0);
        T max_edge(0);

        for (int e=0;e<3;e++) {
            for (int coord=0;coord<3;coord++) {
                T edge=c[(e+1)%3][coord]-c[e][coord];
                if (max_edge<edge) {
                    max_edge=edge;
                    max_e=e;
                }
            }
        }

        if (max_edge==T(0)) return vec<T,3>(T(0));

        return Normalize(cross(
                        (c[max_e]      -c[(max_e+2)%3])/max_edge,
                        (c[(max_e+1)%3]-c[(max_e+2)%3])
                    ));
    }*/


    template <class T, int D> inline void vec<T,D>::sort() {
        std::sort(x,x+D);
    }

    template <class T, int D> inline void vec<T,D>::reverse_sort() {
        std::sort(x,x+D,std::greater<T>());
    }


    template<int D> inline int  Parity(const vec<int,D>& v) {
        return Norm(v)%2;
    }

    template<int D, class T> inline vec<T,D> get_corner(unsigned int i) {
        vec<T,D> tmp(T(0));
        for (unsigned int k=0;k<D;k++) {
            if (((1<<k) & i) != 0) tmp[k]++;
        }
        return tmp;
    }

    /*template <class T, int D> inline int vec<T,D>::GetDirection() const {
        int dir=-1;
        int num0=0;
        for (int i=0;i<D;i++) {
            if (x[i]!=T(0)) {
                num0++;
                dir=i;
            }
        }
        if (num0==1) return dir; else return -1;
    }*/

    template <class T, int D> inline int MinIndex(const vec<T,D>& v) {
        int idx=0;
        for (int i=1;i<D;i++) {
            if (v[i]<v[idx]) idx=i;
        }
        return idx;
    }

    template <class T, int D> inline int MaxIndex(const vec<T,D>& v) {
        int idx=0;
        for (int i=1;i<D;i++) {
            if (v[i]>v[idx]) idx=i;
        }
        return idx;
    }


    template <class T> inline bool Orientation(const vec<T,3>* v) {
        return (dot(cross(v[1]-v[0],v[2]-v[0]),v[3]-v[0])>=-0.);
    }

    template <class T> inline bool Orientation(const vec<T,2>* v) {
        return (dot(RotateLeft(v[1]-v[0]),v[2]-v[0])>=-0.);
    }

    template <class S, class T, int D> S& operator<<(S& s, const vec<T,D> & v) {
        s << "[" << (v[0]);
        for (int i=1;i<D;++i) s << "," << (v[i]);
        s << "]";
        return s;
    }

    template <class T, int D> int compare(const vec<T,D>& v1, const vec<T,D>& v2) {
        for (int i=D-1;i>=0;--i) {
            if (v1[i]>v2[i]) return 1;
            if (v1[i]<v2[i]) return -1;
        }
        return 0;
    }

    template <class T, int D> bool any(const vec<T,D>& v1, const vec<T,D>& v2){
        for(int i=0; i<D-1; ++i){
            if(v1[i]==v2[i]) return true;
        }
        return false;
    }

    /*template <class T> inline const vec<T,3>& Normal(const vec<T,3>* v) {
        return dot(cross(v[1]-v[0],v[2]-v[0]),v[3]-v[0])
    }*/


    /*template <class T, int D> inline vec<T,D> vec<T,D>::ReplaceNth(int i,const T& val) const {
        vec<T,D> v(x);
        v[i]=val;
        return v;
    }*/
}

#endif /*VECTOR_HPP_*/
