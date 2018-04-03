#ifndef DEF_STATISTICS
#define DEF_STATISTICS

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include "Math.h"
#include <cmath>
#include <cstdlib>
#include <omp.h>
//#include "sprng/sprng.h"
#include <chrono>
#include <random>

#include "message.h"
#include <vector>
#include <fstream>
#include <iostream>

namespace my {
  ///Contains Random number generation algorythms and other statistical tools.
  namespace stat {

        static const double epsilon=1e-10;
  unsigned int ClockSEED = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(ClockSEED);
  std::uniform_real_distribution<double> distribution(0.0,1.0);

    using namespace math;

    int* rng;
    #pragma omp threadprivate (rng)

    inline double RandomNumber() {
        return distribution(generator);
    }

       inline void PickRandomPointOnUnitCircle(double& a, double& b) {        //better on AMD
            double x,y,x2,y2,x2py2;
            do {
                x=RandomNumber()-0.5;
                x2=x*x;
                y=RandomNumber()-0.5;
                y2=y*y;
                x2py2=x2+y2;
            } while ((x2py2>=0.25) || (x2py2<=epsilon));
            a=(x2-y2)/x2py2;
            b=2*((x*y)/x2py2);
        }

        inline void PickRandomPointOnUnitCircle2(double& a, double& b) {
            double phi=RandomNumber()*Pi2;
            a=std::cos(phi);
            b=std::sin(phi);
        }

        inline void PickRandomPointOnUnitSphere(double& x, double& y, double& z) {     //better
            double x2,y2,x2py2;
            do {
                x=2*RandomNumber()-1.;
                x2=x*x;
                y=2*RandomNumber()-1.;
                y2=y*y;
                x2py2=x2+y2;
            } while (x2py2>=1.);
            double tmp=2*std::sqrt(1.-x2py2);
            x*=tmp;
            y*=tmp;
            z=1.-2*x2py2;
        }

        inline double PowerCosineSineDistributionReturnCosTheta(const double N) {
            return std::pow(RandomNumber(),1./(N+1.));
        }


        inline double ConeCosineSineDistributionReturnTheta(const double cone_angle) {

            double u, sqrt_1m_u;
            double angle;
            do {
                u=std::sqrt(RandomNumber());
                sqrt_1m_u=std::sqrt(1.-u);
                angle=cone_angle*sqrt_1m_u;
            } while (RandomNumber()*angle*u>std::cos(Pi1_2*sqrt_1m_u)*std::sin(angle));
            return angle;
        }

        inline double ConeCosineSineDistributionReturnCosTheta(const double cone_angle) {

            double u, sqrt_1m_u;
            double cosine;
            double left, right;
            do {
                u=std::sqrt(RandomNumber());
                sqrt_1m_u=std::sqrt(1.-u);
                cosine=std::cos(cone_angle*sqrt_1m_u);
                left=RandomNumber()*cone_angle*sqrt_1m_u*u;
                left*=left;
                right=std::cos(Pi1_2*sqrt_1m_u);
                right*=right;
                right*=(1.-cosine*cosine);
            } while (left>right);
            return cosine;
        }


        inline double ConeCosineSineDistributionReturnTheta2(const double cone_angle) {

            double cosine;
            double _1_m_cos_cone_angle=1.-std::cos(cone_angle);
            double angle;
            double a=Pi1_2/cone_angle;
            do {
                cosine=1-RandomNumber()*_1_m_cos_cone_angle;
                angle=std::acos(cosine);
            } while (RandomNumber()>std::cos(a*angle));

            return angle;
        }

        inline double ConeCosineSineDistributionReturnTheta3(const double cone_angle) {
            double angle;
            double sqrt;
            do {
                sqrt=std::sqrt(RandomNumber());
                angle=sqrt*cone_angle;
            } while (RandomNumber()*angle>std::cos(Pi1_2*sqrt)*std::sin(angle));

            return angle;
        }



        template<class VecType, class VecType2>
        inline void Rotate(const VecType& AverageDirection, VecType2& RandomDirection, const double sinphi, const double cosphi, double costheta, const double r2=1.) {

          costheta=std::min(costheta,1.);

            double a0;
            double a1;

            if (std::fabs(AverageDirection[0])<=std::fabs(AverageDirection[1])) {
                a0=AverageDirection[0];
                a1=AverageDirection[1];
            } else {
                a0=AverageDirection[1];
                a1=AverageDirection[0];
            }

            const double a0_a0_m1=1.-a0*a0;
            const double tmp=std::sqrt(std::max(1.-costheta*costheta,0.)/(r2*a0_a0_m1));
            const double tmp_sinphi=tmp*sinphi;
            const double tmp_cosphi=tmp*cosphi;
            const double costheta_p_a0_tmp_sinphi=costheta+a0*tmp_sinphi;

            RandomDirection[0]=a0*costheta-a0_a0_m1*tmp_sinphi;
            RandomDirection[1]=a1                 *costheta_p_a0_tmp_sinphi+AverageDirection[2]*tmp_cosphi;
            RandomDirection[2]=AverageDirection[2]*costheta_p_a0_tmp_sinphi-a1                 *tmp_cosphi;

            if (a0!=AverageDirection[0]) std::swap(RandomDirection[0],RandomDirection[1]);
        }


        template<class VecType, class VecType2>
        inline void RandomAzimuthalRotation3(const VecType& AverageDirection, VecType2& RandomDirection, const double costheta) {

            double cosphi, sinphi;
            PickRandomPointOnUnitCircle(cosphi, sinphi);

            Rotate(AverageDirection, RandomDirection, sinphi, cosphi, costheta);

        }

        template<class VecType, class VecType2>
        inline void RandomAzimuthalRotation(const VecType& AverageDirection, VecType2& RandomDirection, const double costheta) {

            double cosphi, sinphi;
            double r2;

            do {
                cosphi=RandomNumber()-0.5;
                sinphi=RandomNumber()-0.5;
                r2=cosphi*cosphi+sinphi*sinphi;
            } while (r2>=0.25 || r2<=epsilon) ;

            Rotate(AverageDirection, RandomDirection, sinphi, cosphi, costheta, r2);
        }

        template<class VecType, class VecType2>
        inline void RandomAzimuthalRotation2(const VecType& AverageDirection, VecType2& RandomDirection, const double costheta) {

            double tmp[3];
            double dot;
            do {
                PickRandomPointOnUnitSphere(tmp[0],tmp[1], tmp[2]);
                dot=AverageDirection[0]*tmp[0]+AverageDirection[1]*tmp[1]+AverageDirection[2]*tmp[2];
            } while (dot>=1.);

            double r=std::sqrt((1-costheta*costheta)/(1-dot*dot));

            RandomDirection[0]=AverageDirection[0]*costheta+(tmp[0]-AverageDirection[0]*dot)*r;
            RandomDirection[1]=AverageDirection[1]*costheta+(tmp[1]-AverageDirection[1]*dot)*r;
            RandomDirection[2]=AverageDirection[2]*costheta+(tmp[2]-AverageDirection[2]*dot)*r;
        }

        template<class VecType, class VecType2>
        inline void CosineNDistributedRandomDirection(double N, const VecType& AverageDirection, VecType2& RandomDirection) {

            double costheta=PowerCosineSineDistributionReturnCosTheta(N);

            RandomAzimuthalRotation(AverageDirection, RandomDirection, costheta);

        }

        template<class VecType, class VecType2>
        inline void CosineNDistributedRandomDirection(double N, const VecType& AverageDirection, VecType2& RandomDirection, double cos_cutoff_angle) {

            double costheta;

            do {
                costheta=PowerCosineSineDistributionReturnCosTheta(N);
            } while (costheta<cos_cutoff_angle);

            RandomAzimuthalRotation(AverageDirection, RandomDirection, costheta);
        }

        template<class VecType, class VecType2>
        inline void Cosine1DistributedRandomDirection(const VecType& AverageDirection, VecType2& RandomDirection, double twice_cos_cutoff_angle=0.) {
            double tmp;
            do {
                PickRandomPointOnUnitSphere(RandomDirection[0], RandomDirection[1], RandomDirection[2]);
                RandomDirection[0]+=AverageDirection[0];
                RandomDirection[1]+=AverageDirection[1];
                RandomDirection[2]+=AverageDirection[2];
                tmp=std::sqrt(RandomDirection[0]*RandomDirection[0]+RandomDirection[1]*RandomDirection[1]+RandomDirection[2]*RandomDirection[2]);
            } while (tmp<=twice_cos_cutoff_angle);
            RandomDirection[0]/=tmp;
            RandomDirection[1]/=tmp;
            RandomDirection[2]/=tmp;
        }

        template<class VecType, class VecType2>
        inline void CosAngleDistributedRandomDirection(double Angle,const VecType& AverageDirection, VecType2& RandomDirection) {

            double costheta=std::cos(ConeCosineSineDistributionReturnTheta(Angle));
            RandomAzimuthalRotation(AverageDirection, RandomDirection, costheta);

        }

    template<class VecType1, class VecType2, class VecType3, class ValType>
        inline void NormalDistributedStartPosition(const VecType1& center, const VecType2& dir, ValType FWHM, VecType3& position) {

            static const double fac=-std::log(16.);
            double v0,v1,rsq;

            do {
                v0=2.0*RandomNumber()-1.0;
                v1=2.0*RandomNumber()-1.0;
                rsq=v0*v0+v1*v1;
            } while (rsq>=1.0 || rsq<1e-20);

            Rotate(dir, position,v0,v1, 0.,rsq);

            rsq=std::sqrt(std::log(rsq)/(rsq*fac))*FWHM;

            for (int k=0;k<3;k++) {
              position[k]*=rsq;
              position[k]+=center[k];
            }
    }

    template<class VecType1, class VecType2, class VecType3, class ValType, class ParameterType>
        inline void NormalDistributedStartPosition2(const VecType1& center, const VecType2& dir, ValType FWHM, VecType3& position, const ParameterType& Parameter) {

            static const double fac=-std::log(16.);
            double v0,v1,rsq;

            do {
                v0=2.0*RandomNumber()-1.0;
                v1=2.0*RandomNumber()-1.0;
                rsq=v0*v0+v1*v1;
            } while (rsq>=1.0 || rsq<1e-20);

            Rotate(dir, position,v0,v1, 0.,rsq);

            rsq=std::sqrt(std::log(rsq)/(fac))*FWHM;

            for (int k=0;k<3;k++) {
              position[k]*=rsq;
              position[k]+=center[k];
            }

    }

    template<class VecType1, class VecType2, class VecType3, class ValType, class ParameterType>
        inline void LorentzDistributedStartPosition(const VecType1& center, const VecType2& dir, ValType FWHM, VecType3& position, const ParameterType& Parameter) {

      int x=4;
//      int y;

      double radius;
      radius=FWHM/2*std::sqrt(exp(my::math::Pi*RandomNumber())-1);

      double theta;
            theta = 2*my::math::Pi*RandomNumber();

      for (int i=0;i<3;i++) {
        if (Parameter.open_boundary==i) {
          position[i]=center[i];
        } else {
          if (x==4) {
            x=i;
            position[i]=radius*std::cos(theta)+center[i];
            //position[i]=radius*std::cos(theta)+center[i];
          } else {
//            y=i;
            position[i]=radius*std::sin(theta)+center[i];
            //position[i]=radius*std::sin(theta)+center[i];
          }
        }
      }

    }

    template<class VecType1, class VecType2, class VecType3, class ValType, class ParameterType>
        inline void SurfaceChargeDensityDistributedStartPosition(const VecType1& center, const VecType2& dir, ValType FWHM, VecType3& position, const ParameterType& Parameter, double distance) {

      int x=4, y=4;

      double radius, randnum;
      randnum=RandomNumber();

      radius=2*FWHM*std::sqrt((1/(1-randnum))*(1/(1-randnum))-1)/3;
      radius*=5e-9/distance;

      double theta;
            theta = 2*my::math::Pi*RandomNumber();

      for (int i=0;i<3;i++) {
        if (Parameter.open_boundary==i) {
          position[i]=center[i];
        } else {
          if (x==4) {
            x=i;
            position[i]=radius*std::cos(theta)+center[i];
          } else {
            y=i;
            position[i]=radius*std::sin(theta)+center[i];
          }
        }
      }
    }

    template<class VecType1, class VecType2, class VecType3, class ValType, class ParameterType>
        inline double SurfaceChargeDensityDistributedStartPosition(const VecType1& center, const VecType2& dir, ValType FWHM, VecType3& position, const ParameterType& Parameter) {

      int x=4;
//      int y;

      double radius, randnum;

        randnum=RandomNumber();

        radius=2*FWHM*std::sqrt((1/(1-randnum))*(1/(1-randnum))-1)/3;


      double theta;
            theta = 2*my::math::Pi*RandomNumber();

      for (int i=0;i<3;i++) {
        if (Parameter.open_boundary==i) {
          position[i]=center[i];
        } else {
          if (x==4) {
            x=i;
            position[i]=radius*std::cos(theta)+center[i];
            //position[i]=radius*std::cos(theta)+center[i];
          } else {
//            y=i;
            position[i]=radius*std::sin(theta)+center[i];
            //position[i]=radius*std::sin(theta)+center[i];
          }
        }
      }

      return radius;

    }

    template<class VecType1, class VecType2, class VecType3, class ValType, class ParameterType, class PartitionType>
        inline void SurfaceChargeDensityDistribution(const VecType1& center, const VecType2& dir, ValType voltage, VecType3& position, const ParameterType& Parameter, const PartitionType& Partition) {

          bool keep;
        double d;
        for (int i=0;i<3;i++) d=(Parameter.open_boundary==i)?(center[i]-(Partition.Max(i)-1)*Parameter.grid_delta):d;

            do {
              for (int i=0;i<3;i++) position[i]=(Parameter.open_boundary==i)?center[i]:RandomNumber()*Partition.Extension(i)*Parameter.grid_delta;

            double PositionSumSquares=0;
            for (int i=0;i<3;i++) PositionSumSquares+=(Parameter.open_boundary!=i)?(position[i]-center[i])*(position[i]-center[i]):0;

            double scd;
        scd = (d*1.6e-19)/(2*(my::math::Pi)*(PositionSumSquares+d*d)*sqrt(PositionSumSquares+d*d));

        double scd_max;
        scd_max = (1.6e-19)/(2*(my::math::Pi)*(d*d));

        double randomnu;
        randomnu=(scd_max)*RandomNumber();

        keep=scd<randomnu;

            } while (keep);
    }

    template<class VecType1, class VecType2, class ValType, class ParameterType, class PartitionType>
        inline void SurfaceChargeDensityDistribution(const VecType1& dir, ValType voltage, VecType2& position, const ParameterType& Parameter, const PartitionType& Partition, const std::vector<double>& positions, const std::vector<double>& charges) {

      //const std::vector<double, std::allocator<double> >&, const std::vector<double, std::allocator<double> >&

          bool keep;
        double d=3.7;//e-9;

        double position_max[3];
        position_max[0]=21.5;//25;
        position_max[1]=43.7;
        position_max[2]=20.;


            do {
            double scd=0;
            double scd_max=0;
            for (int i=0;i<3;i++) position[i]=(Parameter.open_boundary==i)?((Partition.Max(i)-1)*Parameter.grid_delta):RandomNumber()*Partition.Extension(i)*Parameter.grid_delta;

//            position_max[0]=23;//25;
//            position_max[1]=43.7;
//            position_max[2]=20;
              for (unsigned int i=0;i<charges.size();++i){

                double PositionSumSquares=0;
                for (int j=0;j<3;j++) PositionSumSquares+=(Parameter.open_boundary!=j)?(position[j]-positions[3*i+j])*(position[j]-positions[3*i+j]):0;


                double PositionSumSquaresMax=0;
                for (int j=0;j<3;j++) PositionSumSquaresMax+=(Parameter.open_boundary!=j)?(position_max[j]-positions[3*i+j])*(position_max[j]-positions[3*i+j]):0;

                scd   += (d*charges[i]*1.6e-19)/(2*(my::math::Pi)*(PositionSumSquares+d*d)*sqrt(PositionSumSquares+d*d));
                scd_max += (d*charges[i]*1.6e-19)/(2*(my::math::Pi)*(PositionSumSquaresMax+d*d)*sqrt(PositionSumSquaresMax+d*d));
              }

        double randomnu;
        randomnu=(scd_max)*RandomNumber()*1.2;

        keep=scd<randomnu;

            } while (keep);
    }


    template<class VecType1, class VecType2, class VecType3, class VecType4, class ValType, class ParameterType>
        inline double NanowireSurfaceCharge(const VecType1& StartPosition, const VecType2 EndPosition, const VecType3& dir, ValType FWHM, ValType Length, ValType Height, ValType Angle, VecType4& position, const ParameterType& Parameter) {

            double volume_sides, volume_line;
            volume_sides = 1;//FWHM/3;//2*FWHM/3;
          volume_line = Length*3/(my::math::Pi*FWHM);

      double line_or_sides;
      line_or_sides=(volume_sides+volume_line)*RandomNumber();

            double X,Y;
            int x=4, y=4;

      for (int i=0;i<3;i++) {
        if (Parameter.open_boundary!=i) {
          if (x==4) {
            x=i;
          } else {
            y=i;
          }
        }
      }
      //double d=3.7e-9;
      //double d=106.286545976e-9;

      double radius=0;
            if (line_or_sides <= volume_line) { // generate a particle in the line

          double randx, randy;

          randx=RandomNumber()-0.5;

//                double radius;
          radius=2*FWHM*randx/(3*sqrt(1-4*randx*randx));
          //radius=2*2*FWHM*randx/(3*sqrt(1-4*randx*randx));
                X = radius;

                randy=RandomNumber();
                Y = Length*randy;
            } else {  // generate a particle on the sides

              double randnum=RandomNumber();

//          double radius;
          radius=FWHM*std::sqrt((1/(1-randnum))*(1/(1-randnum))-1)/3;
          //radius=2*FWHM*std::sqrt((1/(1-randnum))*(1/(1-randnum))-1)/3;

          double theta;
                theta = 2*my::math::Pi*RandomNumber();

              X=radius*std::cos(theta);
              Y=radius*std::sin(theta);
              Y += Y>0?Length:0;

//              radius=2*radius*radius;
            }

            for (int i=0;i<3;i++) {
              if (i==x) {
                position[i] = cos(-Angle)*(X) - sin(-Angle)*(Y) + StartPosition[i];
              } else if (i==y) {
                position[i] = sin(-Angle)*(X) + cos(-Angle)*(Y) + StartPosition[i];
              } else {
                position[i] = StartPosition[i];
              }
            }
            return radius;
    }

    template<class VecType1, class ParameterType, class PartitionType>
        inline void Junctionless(VecType1& position, const ParameterType& Parameter, const PartitionType& Partition) {

      double where = 77.22*RandomNumber();
      double X,Z;

      if (where < 25){ //Box 1
        X = 5000*RandomNumber()+10000;
        Z = 5000*RandomNumber()+5000;
      } else if (where < 50) { //Box 2
        X = 5000*RandomNumber();
        Z = 5000*RandomNumber();
      } else if (where < 75) { //Box 3
        X = 5000*RandomNumber();
        Z = 5000*RandomNumber()+10000;
      } else { //Gate line
        X = 7400*RandomNumber()+2600;
        Z = 300*RandomNumber()+7350;
      }

      position[0]=X;
      position[1]=100;
      position[2]=Z;
    }

    template<class VecType1, class VecType2, class VecType3, class VecType4, class ValType, class ParameterType, class PartitionType>
        inline void NanowireSurfaceChargeDistribution(const VecType1& StartPosition, const VecType2 EndPosition, const VecType3& dir, ValType FWHM, ValType Length, ValType Angle, VecType4& position, const ParameterType& Parameter, const PartitionType& Partition) {

          bool keep;
        double d=0;
        int x=4;
        int y=4;

        for (int i=0;i<3;i++) d=(Parameter.open_boundary==i)?(StartPosition[i]-(Partition.Max(i)-1)*Parameter.grid_delta):d;

        double scd_max;
      scd_max = (1.6e-19)/(2*(my::math::Pi)*(d*d));

            do {
            for (int i=0;i<3;i++) {
              if (Parameter.open_boundary==i){
                position[i]=StartPosition[i];
              } else {
                position[i]=RandomNumber()*Partition.Extension(i)*Parameter.grid_delta;
                if (x==4) {  //if x has not yet been assigned
                  x=i;
                } else {  //x has been assigned, now assign y
                  y=i;
                }
              }
            }

            double alpha; //apha is the angle the line makes with the horizontal direction
            if ((EndPosition[x]-StartPosition[x])==0) {
              if ((EndPosition[y]-StartPosition[y]) > 0) {
                alpha=my::math::Pi/2;
              } else {//if ((EndPosition[2]-StartPosition[2]) < 0) {
                alpha=-my::math::Pi/2;
              }
            } else {
              alpha=std::atan((EndPosition[y]-StartPosition[y])/(EndPosition[x]-StartPosition[x]));
              if ((EndPosition[x]-StartPosition[x])<0) alpha+=my::math::Pi;
            }

            double beta; //beta is the angle between start position and position();
            if ((position[x]-StartPosition[x])==0) {
              if ((position[y]-StartPosition[y]) > 0) {
                beta=my::math::Pi/2;
              } else if ((position[y]-StartPosition[y]) < 0) {
                beta=-my::math::Pi/2;
              } else { //(position[0]-StartPosition[0]) and (position[2]-StartPosition[2])
                beta=alpha;
              }
            } else {
              beta=std::atan((position[y]-StartPosition[y])/(position[x]-StartPosition[x]));
              if ((position[x]-StartPosition[x])<0) beta+=my::math::Pi;
            }

            double hypotenuse_sq;
            hypotenuse_sq = (position[y]-StartPosition[y])*(position[y]-StartPosition[y])+(position[x]-StartPosition[x])*(position[x]-StartPosition[x]);

            double PositionSumSquares;

            if (beta<(alpha-my::math::Pi/2) || beta>(alpha+my::math::Pi/2)) {
              PositionSumSquares=hypotenuse_sq;
            } else {
              double theta;
              theta=beta-alpha;

              double DirectionalLength;
              DirectionalLength=std::sqrt(hypotenuse_sq)*std::cos(theta);

              if (DirectionalLength > Length) {
                PositionSumSquares = (position[y]-EndPosition[y])*(position[y]-EndPosition[y])+(position[x]-EndPosition[x])*(position[x]-EndPosition[x]);
              } else {
                PositionSumSquares = hypotenuse_sq*std::sin(theta)*std::sin(theta);
              }
            }

            double scd;
        scd = (d*1.6e-19)/(2*(my::math::Pi)*(PositionSumSquares+d*d)*sqrt(PositionSumSquares+d*d));

        double randomnu;
        randomnu=(scd_max)*RandomNumber();

        keep=scd<randomnu;

            } while (keep);
    }

    template<class VecType1, class VecType2, class VecType3, class VecType4, class ValType, class ParameterType, class PartitionType>
        inline void NanowireLorentzDistribution(const VecType1& StartPosition, const VecType2 EndPosition, const VecType3& dir, ValType FWHM, ValType Length, ValType Angle, VecType4& position, const ParameterType& Parameter, const PartitionType& Partition) {

          bool keep;
        int x=4;
        int y=4;

          double sigma;
          sigma=FWHM/2;

        double cauchy_max;
        cauchy_max=1/(my::math::Pi*sigma);

            do {
            for (int i=0;i<3;i++) {
              if (Parameter.open_boundary==i){
                position[i]=StartPosition[i];
              } else {
                position[i]=RandomNumber()*Partition.Extension(i)*Parameter.grid_delta;
                if (x==4) {  //if x has not yet been assigned
                  x=i;
                } else {  //x has been assigned, now assign y
                  y=i;
                }
              }
            }

            double alpha; //apha is the angle the line makes with the horizontal direction
            if ((EndPosition[x]-StartPosition[x])==0) {
              if ((EndPosition[y]-StartPosition[y]) > 0) {
                alpha=my::math::Pi/2;
              } else {//if ((EndPosition[2]-StartPosition[2]) < 0) {
                alpha=-my::math::Pi/2;
              }
            } else {
              alpha=std::atan((EndPosition[y]-StartPosition[y])/(EndPosition[x]-StartPosition[x]));
              if ((EndPosition[x]-StartPosition[x])<0) alpha+=my::math::Pi;
            }

            double beta; //beta is the angle between start position and position();
            if ((position[x]-StartPosition[x])==0) {
              if ((position[y]-StartPosition[y]) > 0) {
                beta=my::math::Pi/2;
              } else if ((position[y]-StartPosition[y]) < 0) {
                beta=-my::math::Pi/2;
              } else { //(position[0]-StartPosition[0])==0 and (position[2]-StartPosition[2])==0
                beta=alpha;
              }
            } else {
              beta=std::atan((position[y]-StartPosition[y])/(position[x]-StartPosition[x]));
              if ((position[x]-StartPosition[x])<0) beta+=my::math::Pi;
            }

            double hypotenuse_sq;
            hypotenuse_sq = (position[y]-StartPosition[y])*(position[y]-StartPosition[y])+(position[x]-StartPosition[x])*(position[x]-StartPosition[x]);

            double rc;

            if (beta<(alpha-my::math::Pi/2) || beta>(alpha+my::math::Pi/2)) {
              rc=hypotenuse_sq;
            } else {
              double theta;
              theta=beta-alpha;

              double DirectionalLength;
              DirectionalLength=std::sqrt(hypotenuse_sq)*std::cos(theta);

              if (DirectionalLength > Length) {
                rc = (position[y]-EndPosition[y])*(position[y]-EndPosition[y])+(position[x]-EndPosition[x])*(position[x]-EndPosition[x]);
              } else {
                rc = hypotenuse_sq*std::sin(theta)*std::sin(theta);
              }
            }

              double cauchy;
            cauchy=1/(my::math::Pi*sigma*(1+(rc/(sigma*sigma))));

        double randomnu;
        randomnu=(cauchy_max)*RandomNumber();

        keep=cauchy<randomnu;

            } while (keep);
    }

    //using alternative Monte Carlo
    template<class VecType1, class VecType2, class VecType3, class VecType4, class ValType, class ParameterType>
        inline void NanowireLorentzDistribution2(const VecType1& StartPosition, const VecType2 EndPosition, const VecType3& dir, ValType FWHM, ValType Length, ValType Height, ValType Angle, VecType4& position, const ParameterType& Parameter) {

            double volume_sides, volume_line;
            volume_sides = 2*Height*(FWHM/2)*(FWHM/2)*(my::math::Pi)*(my::math::Pi);
          volume_line = Height*(my::math::Pi)*(FWHM/2)*Length;

      double line_or_sides;
      line_or_sides=(volume_sides+volume_line)*RandomNumber();

            double X,Y;
            int x=4, y=4;

      for (int i=0;i<3;i++) {
        if (Parameter.open_boundary!=i) {
          if (x==4) {
            x=i;
          } else {
            y=i;
          }
        }
      }

            if (line_or_sides <= volume_line) { // generate a particle in the line

          double randx, randy;

          randx=RandomNumber()-0.5;

                double radius;
                radius=std::tan(my::math::Pi*randx);

                X = FWHM/2*radius;

                randy=RandomNumber();
                Y = Length*randy;
            } else {  // generate a particle on the sides

//          double rannum;
//           rannum=RandomNumber();

          double radius;
          //radius=std::sqrt(exp(my::math::Pi*RandomNumber())-1);
          radius=FWHM/2*std::sqrt(std::abs(exp(2*my::math::Pi*RandomNumber()))-1);

                double theta;
                theta = 2*my::math::Pi*RandomNumber();

              X=radius*std::cos(theta);
              Y=radius*std::sin(theta);
              Y += Y>0?Length:0;
            }

            for (int i=0;i<3;i++) {
              if (i==x) {
                position[i] = cos(-Angle)*(X) - sin(-Angle)*(Y) + StartPosition[i];
              } else if (i==y) {
                position[i] = sin(-Angle)*(X) + cos(-Angle)*(Y) + StartPosition[i];
              } else {
                position[i] = StartPosition[i];
              }
            }
    }

    template<class DropletType, class VecType1, class VecType2, class ParameterType, class PartitionType>
    inline void ESDDistribution(const DropletType& d, const VecType1& StartPosition, VecType2& Position, double& r, double& q, long double* Velocity, const ParameterType& Parameter, const PartitionType& Partition){

      double d_test_v; double d_test_r;

      //----------------Find the radius distribution----------------------------------
      double volume_fraction=0.42*RandomNumber()+0.58;
      double r_min_inv=1/2.5e-6; double r_max_inv=1/55e-6;
      double r_max_inv_third=exp(log(r_max_inv)/3);
      double r_min_inv_third=exp(log(r_min_inv)/3);
      double radius=1/pow(volume_fraction*(r_max_inv_third-r_min_inv_third)+r_min_inv_third,3.);

      double eta_a = 2.2e-5;  //Ns/m2
      double rho_d = 789;    // kg/m3

      //------------------------------------------------------------------------------

      //----------------Calculate the charge given r----------------------------------
      double gamma_d = 0.022;
      double permittivity = 8.854187817e-12;
      double qd=0.58*8*my::math::Pi*sqrt(gamma_d*permittivity*radius*radius*radius);
      //------------------------------------------------------------------------------

      //----------------Find the initial droplet position (cylidrical)----------------
      long double theta_i=0.5*my::math::Pi*RandomNumber()*0.5; // 45 degree spray cone
      long double phi=2*my::math::Pi*RandomNumber();
      double location_radius=0.5;
      double z_star1  = 1-location_radius*cos(theta_i);  // height
      double r_star1  = location_radius*sin(theta_i);  // radius

      double mass  = 4*my::math::Pi*rho_d*radius*radius*radius/3;


      //---FIND THE INITIAL AND THERMALDROPLET ELECTRICAL FORCES AND ACCELLERATIONS---

      double H      = StartPosition[Parameter.open_boundary]; //270 mm or 0.27 m
      double Phi_0    = 10e3;          //V
      double R      = 1e-3;          //outer radius of the nozzle (guess)
      double K_V      = 1-exp(-0.021*H/R);  //non-dimensional related to H/R ratio
      double Phi_star = K_V/(log(4*H/R));
      double E_e    = Phi_0*Phi_star/H;

      //-----------Find expected initial electric force when E-field acts alone--------
      double plusz1     =  1+z_star1;//1+z_star;
      double minusz1     =  1-z_star1;//1-z_star;
      double rootplusz1  =  std::sqrt(r_star1*r_star1+plusz1*plusz1);
      double rootminusz1 =  std::sqrt(r_star1*r_star1+minusz1*minusz1);
      double E_v1  = E_e*(1/rootminusz1+1/rootplusz1);
      double E_r1  = E_e*(plusz1/rootplusz1-minusz1/rootminusz1)/r_star1;
      double theta=atan(E_r1/E_v1);

      //-----------Find expected final electric force where E-field acts alone---------
      double t_heat     =   10e-3;
      double r_star2     =  (E_r1*(H-t_heat))/E_v1;
      double z_star2     =  0;//t_heat/H;
      double plusz2     =  1+z_star2;//1+z_star;
      double minusz2     =  1-z_star2;//1-z_star;
      double rootplusz2  =  std::sqrt(r_star2*r_star2+plusz2*plusz2);
      double rootminusz2 =  std::sqrt(r_star2*r_star2+minusz2*minusz2);
      double Eth_v2     =   E_e*(1/rootminusz2+1/rootplusz2);
      double Eth_r2     =   E_e*(plusz2/rootplusz2-minusz2/rootminusz2)/r_star2;

      //-----------Use the first and second to come up with linear dependence----------
      double E_v     = (E_v1-Eth_v2)/(H-t_heat);
      double E_r     = (theta<1e-20)?0:(E_r1-Eth_r2)/((H-t_heat)*tan(theta));

      double Fe_v    = qd*E_v;
      double Fe_r    = qd*E_r;

      double ae_v    = Fe_v/mass;    // initial dependent component
      double ae_r   = Fe_r/mass;    // initial dependent component

      //------------Calculate the required constant component of electric force--------
      double Fe_v1 = qd*E_v1;
      double Fe_r1 = qd*E_r1;
      double ae_v1 = Fe_v1/mass;    // Initial constant component
      double ae_r1 = (theta<1e-20)?0:Fe_r1/mass;    // Initial constant component

      //-------------------------------------------------------------------------------

      //------------------------------------------------------------------------------

      //--------------Find all forces acting on the droplet---------------------------
      // Find the droplet mass - knowns: radius, rho_d

      // Gravity force component acceleration
      double g     = 9.81;  //m/s

      // Stokes force component acceleration - knowns: rho_d, eta_a
      double s_f   = (4.5*eta_a)/(rho_d*radius*radius);

      // Electric force component acceleration - knowns: q_d, r_star, z_star;

      double d0_v  = 0;
      double d0_r  = 0;
      double v_0  = 0;
      double v0_v  = v_0*cos(theta);
      double v0_r  = v_0*sin(theta);


      // This is for both dimensions
      //------------------------------------------------------------------------------

      //--------------Find new position of droplet after t----------------------------
      //First, separate the forces due to velocity/displacement dependences

      double a_v = g+ae_v1;//+2*a_e;  // independent acceleration (mainly gravity) - vertical only
      double b   = s_f;        // velocity dependent acceleration  - Stokes force - vert and rad
      double c_v = ae_v;      // displacement dependent acceleration - vertical E-force
      double a_r = ae_r1;  //+2*a_e;  // independent acceleration - a_e component
      double c_r = ae_r;  // displacement dependent acceleration - radial E-force

      //--------------Start vertical onlgoings ---------------------------------------
      long double t_drop;

      double v0th_v;//=0;
      double v0th_r;//=0;
      if (b*b-4*c_v<0) {
        t_drop=0;
        d_test_v  = H-t_heat;
        d_test_r  = r_star1*H;
        v0th_v=1;
        v0th_r=0;
      } else {
        double r1_v  = (-b+std::sqrt(b*b-4*c_v))/(2);
        double r2_v  = (-b-std::sqrt(b*b-4*c_v))/(2);

        double B1_v = d0_v*(r2_v)/(r2_v-r1_v);
        double A1_v  = d0_v-B1_v;
        double B2_v = (v0_v+d0_v)/(r2_v-r1_v);
        double A2_v = -B2_v;;
        double C_v  = a_v/(r1_v*r2_v);
        double B3_v  = -a_v/(r2_v*(r2_v-r1_v));
        double A3_v  = -(B3_v+C_v);

        double iteration_stop=10000;
        double t_low  = 0;
        double t_high = 1;
        double t_check = (t_low+t_high)/2;
        while (true){
          d_test_v = (A1_v+A2_v+A3_v)*exp(r1_v*t_check)+(B1_v+B2_v+B3_v)*exp(r2_v*t_check)+C_v;
          if (t_high-t_low<1e-100) break;
          if (d_test_v > (H-t_heat)) {
            t_high=t_check;
            t_check=(t_high+t_low)/2;
          } else if (d_test_v < (H-t_heat)) {
            t_low=t_check;
            t_check=(t_high+t_low)/2;
          }
          if (iteration_stop==0) break;
          iteration_stop--;
        }

        t_drop=t_check;
        v0th_v    = (A1_v+A2_v+A3_v)*(exp(r1_v*t_drop)-1)/r1_v+(B1_v+B2_v+B3_v)*(exp(r2_v*t_drop)-1)/r2_v+C_v*t_drop+v0_v;

        //---------Now time required to reach heat zone is known------------------------

        //---------Can now calculate the radial displacement----------------------------

        if (b*b-4*c_r<0) {
          t_drop=0;
          d_test_r  = 0;
          v0th_r  = 0;
        } else {
          double r1_r  = (-b+std::sqrt(b*b-4*c_r))/(2);
          double r2_r  = (-b-std::sqrt(b*b-4*c_r))/(2);

          double B1_r = d0_r*(r2_r)/(r2_r-r1_r);
          double A1_r  = d0_r-B1_r;
          double B2_r = (v0_r+d0_r)/(r2_r-r1_r);
          double A2_r = -B2_r;;
          double C_r  = a_r/(r1_r*r2_r);
          double B3_r  = -a_r/(r2_r*(r2_r-r1_r));
          double A3_r  = -(B3_r+C_r);

          d_test_r = (A1_r+A2_r+A3_r)*exp(r1_r*t_drop)+(B1_r+B2_r+B3_r)*exp(r2_r*t_drop)+C_r;
          v0th_r   = (A1_r+A2_r+A3_r)*(exp(r1_r*t_drop)-1)/r1_r+(B1_r+B2_r+B3_r)*(exp(r2_r*t_drop)-1)/r2_r+C_r*t_drop+v0_r;
        }

      }
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------      THERMAL ZONE CALCULATIONS      -----------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
//      //---------Now we know z and r on the cusp of the thermal zone------------------
//      //---------Reset parameters to include thermal components and re-run -----------

      //--------CALCULATE THERMAL EFFECTS IN DROPLET SIZE REDUCTION-------------------
      double dth_test_v;
      double dth_test_r;
      double q0       = 373e-12;    // (373 um^2)    for water: 88e-12 (m^2)
      double q1       = 89.1;      // (8.91e-5 /um) for water: 4.3e3  (/m)
      double del_T    = 100000;
      double dK       = q0*del_T*(1+2*q1*radius);
      double r_new    = radius-radius*t_drop*exp(log(dK)/3);
      r=r_new;
      double qd_new   = 0.58*8*my::math::Pi*sqrt(gamma_d*permittivity*r_new*r_new*r_new);
      double mass_new  = 4*my::math::Pi*rho_d*r_new*r_new*r_new/3;
//
//      //---------CALCULATE THE ENERGIES ASSUMING LINEAR REDUCTION TO FINAL------------
//      //---------Electric force F_e is constant in this region and small--------------
//
      double rth_star  =   d_test_r/H;    //d_r/H;
      double zth_star  =   1-d_test_v/H;  //1-d_v/H;

      double d0th_v    = 0;//d_test_v;
      double d0th_r    = 0;//d_test_r;
      double tth_drop    = t_drop;

      double thplusz     =  1+zth_star;//1+z_star;
      double thminusz     =  1-zth_star;//1-z_star;
      double throotplusz  =  std::sqrt(rth_star*rth_star+thplusz*thplusz);
      double throotminusz =  std::sqrt(rth_star*rth_star+thminusz*thminusz);

      double Eth_v  = E_e*(1/throotminusz+1/throotplusz);
      double Eth_r  = E_e*(thplusz/throotplusz-thminusz/throotminusz)/rth_star;
      double theta_th  = atan(Eth_r/Eth_v);

      double Feth_v  = qd_new*Eth_v;
      double Feth_r  = qd_new*Eth_r;

      double aeth_v  = Feth_v/mass_new;
      double aeth_r   = Feth_r/mass_new;

      //--------------Find all other forces acting on the droplet---------------------
      // Find the droplet mass - knowns: radius, rho_d

      // Gravity force component acceleration

      // Stokes force component acceleration - knowns: rho_d, eta_a
      double sth_f = (4.5*eta_a)/(rho_d*r_new*r_new);

      // Thermophoretic force
      double kappa_a  = 0.025;
      double kappa_d  = 0.19;
      double grad_T  = 100000;
      double T    = 523;
      double rho_a  = 1.29;

      double F_th   = 3*my::math::Pi*eta_a*eta_a*r_new*3*kappa_a*grad_T/(rho_a*T*(2*kappa_a+kappa_d));
      double a_th   = F_th/mass_new;

      // Electric force component acceleration - knowns: q_d, r_star, z_star;

      double ath_v = g+aeth_v-a_th;//+2*a_e;  // independent acceleration (gravity, initial e-force, thermal force)
      double bth    = sth_f;          // velocity dependent acceleration    - Stokes force
      double cth_v = aeth_v/t_heat;      // displacement dependent acceleration   - vertical E-force
      double ath_r = (theta_th<1e-20)?0:aeth_r;                // independent acceleration - Initial electric force
      double cth_r = (theta_th<1e-20)?0:aeth_r/(t_heat*tan(theta_th));  // displacement dependent acceleration - Radial E-force

      double v_final_v=0;
      double v_final_r=0;

      if (bth*bth-4*cth_v<0) {
        tth_drop=0;
        dth_test_v  = t_heat;
        dth_test_r  = d_test_r;
        v_final_v  = 1;
        v_final_r  = 0;
      } else {
        double r1th_v  = (-bth+std::sqrt(bth*bth-4*cth_v))/(2);
        double r2th_v  = (-bth-std::sqrt(bth*bth-4*cth_v))/(2);

        double B1th_v   = d0th_v*(r2th_v)/(r2th_v-r1th_v);
        double A1th_v  = d0th_v-B1th_v;
        double B2th_v   = (v0th_v+d0th_v)/(r2th_v-r1th_v);
        double A2th_v   = -B2th_v;;
        double Cth_v  = ath_v/(r1th_v*r2th_v);
        double B3th_v  = -ath_v/(r2th_v*(r2th_v-r1th_v));
        double A3th_v  = -(B3th_v+Cth_v);

        double iteration_stopth  =10000;
        double tth_low      = 0;
        double tth_high      = t_drop;
        double tth_check    = (tth_low+tth_high)/2;
        while (true){
          dth_test_v = (A1th_v+A2th_v+A3th_v)*exp(r1th_v*tth_check)+(B1th_v+B2th_v+B3th_v)*exp(r2th_v*tth_check)+Cth_v;
          if (tth_high-tth_low<1e-10) break;
          if (dth_test_v > t_heat) {
            tth_high=tth_check;
            tth_check=(tth_high+tth_low)/2;
          } else if (dth_test_v < t_heat) {
            tth_low=tth_check;
            tth_check=(tth_high+tth_low)/2;
          }
          if (iteration_stopth==0) break;
          iteration_stopth--;
        }
        tth_drop  = tth_check;
        v_final_v = (A1th_v+A2th_v+A3th_v)*(exp(r1th_v*tth_drop)-1)/r1th_v+(B1th_v+B2th_v+B3th_v)*(exp(r2th_v*tth_drop)-1)/r2th_v+Cth_v*tth_drop+v0th_v;


        //---------Now time required to reach surface is known--------------------------

        //---------Can now calculate the radial displacement----------------------------

        if (bth*bth-4*cth_r<0) {
          tth_drop=0;
          dth_test_r  = 0;
          v_final_r  = 0;
        } else {
          double r1th_r = (-bth+std::sqrt(bth*bth-4*cth_r))/(2);
          double r2th_r = (-bth-std::sqrt(bth*bth-4*cth_r))/(2);

          double B1th_r   = d0th_r*(r2th_r)/(r2th_r-r1th_r);
          double A1th_r   = d0th_r-B1th_r;
          double B2th_r   = (v0th_r+d0th_r)/(r2th_r-r1th_r);
          double A2th_r   = -B2th_r;;
          double Cth_r    = ath_r/(r1th_r*r2th_r);
          double B3th_r   = -ath_r/(r2th_r*(r2th_r-r1th_r));
          double A3th_r   = -(B3th_r+Cth_r);

          dth_test_r = (A1th_r+A2th_r+A3th_r)*exp(r1th_r*tth_drop)+(B1th_r+B2th_r+B3th_r)*exp(r2th_r*tth_drop)+Cth_r;
          v_final_r  = (A1th_r+A2th_r+A3th_r)*(exp(r1th_r*tth_drop)-1)/r1th_r+(B1th_r+B2th_r+B3th_r)*(exp(r2th_r*tth_drop)-1)/r2th_r+Cth_r*tth_drop+v0th_r;
        }
      }


      Velocity[0] = v_final_r*cos(phi);
      Velocity[1] = -v_final_v;
      Velocity[2] = v_final_r*sin(phi);

      Position[0] = std::sqrt(dth_test_r+d_test_r)*cos(phi)+StartPosition[0];
      Position[1] = 0;
      Position[2] = std::sqrt(dth_test_r+d_test_r)*sin(phi)+StartPosition[2];

      r  = r_new;
      q  = qd_new;

    }

    template<class DropletType, class VecType1, class VecType2, class ParameterType, class PartitionType>
    inline void EvenlyDistributed(const DropletType& d, const VecType1& StartPosition, VecType2& Position, double& r, double& q, long double* Velocity, const ParameterType& Parameter, const PartitionType& Partition){


      Velocity[0] = 0;
      Velocity[1] = -1;
      Velocity[2] = 0;

      double volume_fraction=0.22*RandomNumber()+0.58;
      double r_min_inv=1/2.5e-6; double r_max_inv=1/55e-6;
      double r_max_inv_third=exp(log(r_max_inv)/3);
      double r_min_inv_third=exp(log(r_min_inv)/3);
      double radius=1/pow(volume_fraction*(r_max_inv_third-r_min_inv_third)+r_min_inv_third,3.);

      r  = radius;

      Position[0]=RandomNumber()*4*(Partition.Max(0)*Parameter.grid_delta)-2*(Partition.Max(0)*Parameter.grid_delta);
      Position[1]=0;
      Position[2]=RandomNumber()*4*(Partition.Max(2)*Parameter.grid_delta)-2*(Partition.Max(2)*Parameter.grid_delta);

      double gamma_d = 0.022;
      double permittivity = 8.854187817e-12;
      double qd=0.58*8*my::math::Pi*sqrt(gamma_d*permittivity*radius*radius*radius);
      q  = qd;

    }

    template<class VecType1, class VecType2, class ParameterType, class PartitionType>
    inline void EvenlyDistributed(const VecType1& StartPosition, VecType2& Position, const ParameterType& Parameter, const PartitionType& Partition){

        Position[0]=RandomNumber()*(Partition.Max(0)*Parameter.grid_delta-Partition.Min(0)*Parameter.grid_delta)+Partition.Min(0)*Parameter.grid_delta;
        Position[1]=0;
        Position[2]=RandomNumber()*(Partition.Max(2)*Parameter.grid_delta-Partition.Min(2)*Parameter.grid_delta)+Partition.Min(2)*Parameter.grid_delta;
    }


    template<class DropletType, class VecType1>
        inline void DiskDistribution(const DropletType d, VecType1& position) {

      double v0, v1, rsq;
            do {
                v0=2.0*RandomNumber()-1.0;
                v1=2.0*RandomNumber()-1.0;
                rsq=v0*v0+v1*v1;
            } while (rsq>=1.0 || rsq<1e-20);

      position[0] = 2*d.Radius*v0+d.Position[0];
      position[1] = 0;//d.Position[1];
      position[2] = 2*d.Radius*v1+d.Position[2];
    }

    template<class DataType> bool AnyElement(typename std::vector<DataType> vec, DataType check){
      typename std::vector<DataType>::iterator first = vec.begin(), last = vec.end();
      while(first!=last){
        if(*first==check) return true;
        ++first;
      }
      return false;
    }

  }
}

#endif //DEF_STATISTICS
