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
#include "sprng/sprng.h"
#include "message.h"
//#include "calc.h"
//#include "Partition/PartitionNeighborLinksArrays.h"
#include "LambertW.hpp"

#include <fstream>
#include <iostream>

namespace my {
        namespace stat {

        static const double epsilon=1e-10;

                using namespace math;

		int* rng;
		#pragma omp threadprivate (rng)

		inline double RandomNumber() {
		    return sprng(rng);
		}

		inline void InitRandomGenerator(unsigned int my_rank, unsigned int size, int seed, int rng_type, int rng_par) {

            #pragma omp parallel
            {
                int nstreams=size*MAX_NUM_THREADS;
                int streamnum=MAX_NUM_THREADS*my_rank;

                #ifdef _OPENMP
                streamnum+=omp_get_thread_num();
                #endif

                rng= init_sprng(rng_type,streamnum,nstreams,seed,rng_par);

            }
        }

        inline void FreeRandomGenerator() {
            #pragma omp parallel
            {
                free_sprng(rng);
            }
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

        /*inline void PickRandomPointOnUnitSphere2(double& a, double& b, double& c) {     //worse
            c=2*RandomNumber()-1;
            double tmp=std::sqrt(1-c*c);
            PickRandomPointOnUnitCircle2(a,b);
            a*=tmp;
            b*=tmp;
        }*/

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

                //std::cout << "IN ROTATE:" << endl;
                //std::cout << "costheta: " << costheta << endl;
                costheta=std::min(costheta,1.);

                //std::cout << "costheta: " << costheta << endl;
                //std::cout << "sinphi: " << sinphi << endl;
                //std::cout << "cosphi: " << sinphi << endl;
                //std::cout << "costheta: " << sinphi << endl;
                //std::cout << "r2: " << r2 << endl;

            double a0;
            double a1;

            if (std::fabs(AverageDirection[0])<=std::fabs(AverageDirection[1])) {
                a0=AverageDirection[0];
                a1=AverageDirection[1];
            } else {
                a0=AverageDirection[1];
                a1=AverageDirection[0];
            }
                //std::cout << "a0: " << a0 << endl;
                //std::cout << "a1: " << a1 << endl;

            const double a0_a0_m1=1.-a0*a0;
            //std::cout << "a0_a0_m1: " << a0_a0_m1 << endl;
            const double tmp=std::sqrt(std::max(1.-costheta*costheta,0.)/(r2*a0_a0_m1));
            //std::cout << "tmp: " << tmp << endl;
            const double tmp_sinphi=tmp*sinphi;
            //std::cout << "tmp_sinphi: " << tmp_sinphi << endl;
            const double tmp_cosphi=tmp*cosphi;
            //std::cout << "tmp_cosphi: " << tmp_cosphi << endl;
            const double costheta_p_a0_tmp_sinphi=costheta+a0*tmp_sinphi;
            //std::cout << "costheta_p_a0_tmp_sinphi: " << costheta_p_a0_tmp_sinphi << endl;

            RandomDirection[0]=a0*costheta-a0_a0_m1*tmp_sinphi;
                //std::cout << "RandomDirection[0]: " << RandomDirection[0] << endl;
            RandomDirection[1]=a1                 *costheta_p_a0_tmp_sinphi+AverageDirection[2]*tmp_cosphi;
                //std::cout << "RandomDirection[1]: " << RandomDirection[1] << endl;
            RandomDirection[2]=AverageDirection[2]*costheta_p_a0_tmp_sinphi-a1                 *tmp_cosphi;
            //std::cout << "RandomDirection[2]: " << RandomDirection[2] << endl;

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

//        template<class VecType, class VecType2>
//        inline void MyRandomAzimuthalRotation(const VecType& AverageDirection, VecType2& RandomDirection, const double costheta) {

//            double cosphi, sinphi;
//            double r2;

////            do {
//                cosphi=0.4;//(RandomNumber()-0.5)*.1;
//                sinphi=0.2;//(RandomNumber()-0.5)*.1;
//                r2=cosphi*cosphi+sinphi*sinphi;
////            } while (r2>=0.25 || r2<=epsilon) ;

//            Rotate(AverageDirection, RandomDirection, sinphi, cosphi, costheta, r2);
//        }

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

        /*template<class VecType, class VecType2>
        inline void RandomAzimuthalRotation4(const VecType& AverageDirection, VecType2& RandomDirection, const double costheta) {

            double tmp[3];
            double q;
            double dot;
            double dot2;
            do {
                do {
                    tmp[0]=2*RandomNumber()-1.;
                    tmp[1]=2*RandomNumber()-1.;
                    tmp[2]=2*RandomNumber()-1.;
                    q=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
                } while (q>=1.);
                dot=AverageDirection[0]*tmp[0]+AverageDirection[1]*tmp[1]+AverageDirection[2]*tmp[2];
                dot2=dot*dot;
            } while (dot2>=q);

            double r=std::sqrt((1-costheta*costheta)/(q-dot2));

            RandomDirection[0]=AverageDirection[0]*costheta+(tmp[0]-AverageDirection[0]*dot)*r;
            RandomDirection[1]=AverageDirection[1]*costheta+(tmp[1]-AverageDirection[1]*dot)*r;
            RandomDirection[2]=AverageDirection[2]*costheta+(tmp[2]-AverageDirection[2]*dot)*r;
        }*/

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

        /*template<class VecType, class VecType2>
        inline void Cosine3DistributedRandomDirection(const VecType& AverageDirection, VecType2& RandomDirection) {

            double tmp;
            do {
                do {
                    RandomDirection[0]=2*RandomNumber()-1;
                    RandomDirection[1]=2*RandomNumber()-1;
                    RandomDirection[2]=2*RandomNumber()-1;
                } while (RandomDirection[0]*RandomDirection[0]+RandomDirection[1]*RandomDirection[1]+RandomDirection[2]*RandomDirection[2]>=1.);
                RandomDirection[0]+=AverageDirection[0];
                RandomDirection[1]+=AverageDirection[1];
                RandomDirection[2]+=AverageDirection[2];
                tmp=std::sqrt(RandomDirection[0]*RandomDirection[0]+RandomDirection[1]*RandomDirection[1]+RandomDirection[2]*RandomDirection[2]);
            } while (tmp==0.);
            RandomDirection[0]/=tmp;
            RandomDirection[1]/=tmp;
            RandomDirection[2]/=tmp;
        }*/

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

/*			int x=4, y=4;

			double numb;
			numb = RandomNumber();

			double sig;
			sig = FWHM/(2*std::log(2));

			double radius;
			radius=sig*std::sqrt(std::log(1/((1-(std::sqrt(2*my::math::Pi)*numb))*(1-(std::sqrt(2*my::math::Pi)*numb)))));

                        double theta;
            theta = 2*my::math::Pi*RandomNumber();

			for (int i=0;i<3;i++) {
				if (Parameter.open_boundary_direction==i) {
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
*/

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
//			int y;

			double radius;
			radius=FWHM/2*std::sqrt(exp(my::math::Pi*RandomNumber())-1);
			//radius=FWHM/2*std::tan(std::sqrt(RandomNumber())*my::math::Pi*0.5);
			//radius=FWHM/2*std::sqrt(exp(my::math::Pi*RandomNumber())-1);
//			std::cout << "radius: " << radius << std::endl;

                        double theta;
            theta = 2*my::math::Pi*RandomNumber();

			for (int i=0;i<3;i++) {
				if (Parameter.open_boundary_direction==i) {
					position[i]=center[i];
				} else {
					if (x==4) {
						x=i;
						position[i]=radius*std::cos(theta)+center[i];
						//position[i]=radius*std::cos(theta)+center[i];
					} else {
//						y=i;
						position[i]=radius*std::sin(theta)+center[i];
						//position[i]=radius*std::sin(theta)+center[i];
					}
				}
			}

/*    		bool keep;

                        do {
                for (int i=0;i<3;i++) position[i]=(Parameter.open_boundary_direction==i)?center[i]:RandomNumber()*Partition.Extension(i)*Parameter.GridDelta;

                double rc=0;
                for (int i=0;i<3;i++) rc+=(Parameter.open_boundary_direction==i)?0:(position[i]-center[i])*(position[i]-center[i]);

                double sigma;
                sigma=FWHM/2;

                double cauchy;
                        cauchy=1/(my::math::Pi*sigma*(1+(rc/(sigma*sigma))));

                        double cauchy_max;
                        cauchy_max=1/(my::math::Pi*sigma);

                        double randomnu;
                        randomnu=(cauchy_max)*RandomNumber();

                        keep=cauchy<randomnu;

                        } while (keep);

*/		}

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
				if (Parameter.open_boundary_direction==i) {
					position[i]=center[i];
				} else {
					if (x==4) {
						x=i;
						position[i]=radius*std::cos(theta)+center[i];
						//position[i]=radius*std::cos(theta)+center[i];
					} else {
						y=i;
						position[i]=radius*std::sin(theta)+center[i];
						//position[i]=radius*std::sin(theta)+center[i];
					}
				}
			}
		}

                template<class VecType1, class VecType2, class VecType3, class ValType, class ParameterType>
        inline double SurfaceChargeDensityDistributedStartPosition(const VecType1& center, const VecType2& dir, ValType FWHM, VecType3& position, const ParameterType& Parameter) {

                        int x=4;
//			int y;

                        double radius, randnum;
//			while(true) {

				randnum=RandomNumber();

				radius=2*FWHM*std::sqrt((1/(1-randnum))*(1/(1-randnum))-1)/3;
//				radius=FWHM*std::sqrt((1/(1-randnum))*(1/(1-randnum))-1)/5;
				//radius*=FWHM*1e9;
//				if (radius < FWHM) break;
//			}


                        double theta;
            theta = 2*my::math::Pi*RandomNumber();

			for (int i=0;i<3;i++) {
				if (Parameter.open_boundary_direction==i) {
					position[i]=center[i];
				} else {
					if (x==4) {
						x=i;
						position[i]=radius*std::cos(theta)+center[i];
						//position[i]=radius*std::cos(theta)+center[i];
					} else {
//						y=i;
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
                for (int i=0;i<3;i++) d=(Parameter.open_boundary_direction==i)?(center[i]-(Partition.Max(i)-1)*Parameter.GridDelta):d;

            do {
                for (int i=0;i<3;i++) position[i]=(Parameter.open_boundary_direction==i)?center[i]:RandomNumber()*Partition.Extension(i)*Parameter.GridDelta;

			double PositionSumSquares=0;
			for (int i=0;i<3;i++) PositionSumSquares+=(Parameter.open_boundary_direction!=i)?(position[i]-center[i])*(position[i]-center[i]):0;

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

//    		double PositionSumSquaresMax;
//    		PositionSumSquaresMax=position_max[0]*position_max[0];

//    		double scd_max;
//    		scd_max = 100*(d*charges[0]*1.6e-19)/(2*(my::math::Pi)*(PositionSumSquaresMax+d*d)*sqrt(PositionSumSquaresMax+d*d));

//    		double scd;
//    		double scd_max;
                //for (int i=0;i<3;i++) d=(Parameter.open_boundary_direction==i)?(center[i]-(Partition.Max(i)-1)*Parameter.GridDelta):d;

            do {
                        double scd=0;
                        double scd_max=0;
                        for (int i=0;i<3;i++) position[i]=(Parameter.open_boundary_direction==i)?((Partition.Max(i)-1)*Parameter.GridDelta):RandomNumber()*Partition.Extension(i)*Parameter.GridDelta;

//        		position_max[0]=23;//25;
//        		position_max[1]=43.7;
//        		position_max[2]=20;
                for (unsigned int i=0;i<charges.size();++i){
//            		std::cout << "i: " << i << std::endl;
//            		std::cout << "charges: " << charges[i] << std::endl;

                        double PositionSumSquares=0;
                        for (int j=0;j<3;j++) PositionSumSquares+=(Parameter.open_boundary_direction!=j)?(position[j]-positions[3*i+j])*(position[j]-positions[3*i+j]):0;


                        double PositionSumSquaresMax=0;
                        for (int j=0;j<3;j++) PositionSumSquaresMax+=(Parameter.open_boundary_direction!=j)?(position_max[j]-positions[3*i+j])*(position_max[j]-positions[3*i+j]):0;

//            		std::cout << "PositionSumSquares: " << PositionSumSquares << std::endl;
//            		for (int j=0;j<3;j++) {
//            			d=(Parameter.open_boundary_direction==j)?(positions[3*i+j]-(Partition.Max(j)-1)*Parameter.GridDelta):d;
//            		}
//            		std::cout << "d: " << d << std::endl;

                        scd 	+= (d*charges[i]*1.6e-19)/(2*(my::math::Pi)*(PositionSumSquares+d*d)*sqrt(PositionSumSquares+d*d));
                        scd_max += (d*charges[i]*1.6e-19)/(2*(my::math::Pi)*(PositionSumSquaresMax+d*d)*sqrt(PositionSumSquaresMax+d*d));
//            		scd_max += (charges[i]*1.6e-19)/(2*(my::math::Pi)*(d*d));
                }
//            	std::cout << "scd: " << scd << std::endl;
//            	std::cout << "scd_max: " << scd_max << std::endl;

				double randomnu;
				randomnu=(scd_max)*RandomNumber()*1.2;

				keep=scd<randomnu;

	    } while (keep);
		}


                template<class VecType1, class VecType2, class VecType3, class VecType4, class ValType, class ParameterType>
        inline double NanowireSurfaceCharge(const VecType1& StartPosition, const VecType2 EndPosition, const VecType3& dir, ValType FWHM, ValType Length, ValType Height, ValType Angle, VecType4& position, const ParameterType& Parameter) {

//			double d=20e-9;
//			std::cout << "FWHM: " << FWHM << std::endl;
            double volume_sides, volume_line;
            volume_sides = 1;//FWHM/3;//2*FWHM/3;
//        	volume_line = Length*3/(my::math::Pi*FWHM*2);
                volume_line = Length*3/(my::math::Pi*FWHM);

//            volume_sides = Height*(FWHM/2)*(FWHM/2)*(my::math::Pi)*2*(my::math::Pi);
//        	volume_line = Height*(my::math::Pi)*(FWHM/2)*Length;

//        	std::cout << "volume_sides: " << volume_sides << ", volume_line: " << volume_line << std::endl;
			double line_or_sides;
			line_or_sides=(volume_sides+volume_line)*RandomNumber();

            double X,Y;
            int x=4, y=4;

			for (int i=0;i<3;i++) {
				if (Parameter.open_boundary_direction!=i) {
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
            } else {	// generate a particle on the sides

                double randnum=RandomNumber();

//    			double radius;
                        radius=FWHM*std::sqrt((1/(1-randnum))*(1/(1-randnum))-1)/3;
                        //radius=2*FWHM*std::sqrt((1/(1-randnum))*(1/(1-randnum))-1)/3;

                        double theta;
                theta = 2*my::math::Pi*RandomNumber();

                X=radius*std::cos(theta);
                Y=radius*std::sin(theta);
                Y += Y>0?Length:0;

//            	radius=2*radius*radius;
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

                for (int i=0;i<3;i++) d=(Parameter.open_boundary_direction==i)?(StartPosition[i]-(Partition.Max(i)-1)*Parameter.GridDelta):d;

		double scd_max;
			scd_max = (1.6e-19)/(2*(my::math::Pi)*(d*d));

	    do {
			for (int i=0;i<3;i++) {
				if (Parameter.open_boundary_direction==i){
					position[i]=StartPosition[i];
				} else {
					position[i]=RandomNumber()*Partition.Extension(i)*Parameter.GridDelta;
					if (x==4) {	//if x has not yet been assigned
						x=i;
					} else {	//x has been assigned, now assign y
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
/*
			static const double fac=-std::log(16.);
			double line_or_sides;
			line_or_sides=(vol_line+vol_oth)*RandomNumber();

	    double X,Y;

            if (line_or_sides < vol_line) { // generate a particle in the line
                double v0, v1, rsq;

                do {
                    v0=2.0*RandomNumber()-1.0;
                    v1=2.0*RandomNumber()-1.0;
                    rsq=v0*v0+v1*v1;
                } while (rsq>=1.0 || rsq<1e-20);

                //Distance = std::sqrt((EndPosition[0]-StartPosition[0])*(EndPosition[0]-StartPosition[0])+(EndPosition[2]-StartPosition[2])*(EndPosition[2]-StartPosition[2]));
                rsq=std::sqrt(std::log(rsq)/(fac))*FWHM;

                X = v0*std::sqrt(std::log(rsq)/(fac))*FWHM;

                Y = 1.0*RandomNumber();
                Y = Y*Length;

                //position[0] = X; //X + StartPosition[0];
                //position[1] = StartPosition[1];
                //position[2] = Y; //Y + StartPosition[2];

                //Rotate:

                //Rotate(dir, position,v0,v1, 0.,rsq);

                //rsq=std::sqrt(std::log(rsq)/(fac))*FWHM;

                //for (int k=0;k<3;k++) {
                //	position[k]*=rsq;
                //	position[k]+=center[k];
                //}
            } else {	// generate a particle on the sides
                static const double fac=-std::log(16.);
                double v0,v1,rsq;

                do {
                    v0=2.0*RandomNumber()-1.0;
                    v1=2.0*RandomNumber()-1.0;
                    rsq=v0*v0+v1*v1;
                } while (rsq>=1.0 || rsq<1e-20);

                //Rotate(dir, position,v0,v1, 0.,rsq);

                rsq=std::sqrt(std::log(rsq)/(fac))*FWHM;

                X = v0*std::sqrt(std::log(rsq)/(fac))*FWHM;
                Y = v1*std::sqrt(std::log(rsq)/(fac))*FWHM;
//                for (int k=0;k<3;k++) {
//                	position[k]*=rsq;
//                	position[k]+=center[k];
//                }
                if (v0>0) position[2]+=Length;
            }

                position[0] = cos(-Angle)*(X) - sin(-Angle)*(Y) + StartPosition[0];
            position[2] = sin(-Angle)*(X) + cos(-Angle)*(Y) + StartPosition[2];
            position[1] = StartPosition[1];
*/
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
				if (Parameter.open_boundary_direction==i){
					position[i]=StartPosition[i];
				} else {
					position[i]=RandomNumber()*Partition.Extension(i)*Parameter.GridDelta;
					if (x==4) {	//if x has not yet been assigned
						x=i;
					} else {	//x has been assigned, now assign y
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
				if (Parameter.open_boundary_direction!=i) {
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
            } else {	// generate a particle on the sides

//    			double rannum;
//   				rannum=RandomNumber();

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

                template<class VecType1, class VecType2, class VecType3, class ValType, class ParameterType, class PartitionType>
        inline void SprayPyrolysisAmbient(const VecType1& StartPosition, const ValType& min_radius, const ValType& max_radius, const int& flow_rate, const int& SprayDirection, const double& direction_phi, const double& direction_theta, const int& AdditionalDirection, const VecType2 direction, VecType3& position, double& FluxFactor, const ParameterType& Parameter, const PartitionType& Partition) {

			double theta_min = atan(min((Partition.Min(AdditionalDirection)*Parameter.GridDelta-StartPosition[AdditionalDirection])/(Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection]), (Partition.Min(AdditionalDirection)*Parameter.GridDelta-StartPosition[AdditionalDirection])/(Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])));
			double theta_max = atan(max((Partition.Max(AdditionalDirection)*Parameter.GridDelta-StartPosition[AdditionalDirection])/(Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection]), (Partition.Max(AdditionalDirection)*Parameter.GridDelta-StartPosition[AdditionalDirection])/(Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])));

			//if flow rate is 30ml/h:			//if flow rate is 120ml/h:
		//0-5um  - 0.57						//0-5um  - 0.55
		//5-10um - 0.32						//5-10um - 0.32
		//10-15  - 0.05						//10-15  - 0.08
		//15-20  - 0.025					//15-20  - 0.03
		//20-25  - 0.02						//20-25  - 0.011
		//25-30  - 0.01						//25-30  - 0.006
		//30-35  - 0.003					//30-35  - 0.002
		//35-40  - 0.001					//35-40  - 0.001
		//40-45  - 0.001
			double StartPosition_temp[3];
		long double v0[3];
		double factor=1;//(max_radius-min_radius)/45e-6;
		double s;
		long double t;

		double a_gr = -9.81;
		double air_viscosity = 2.2e-5;
		double density = 780;
		double air_density=1.29;
		double air_thermal_cond=0.025;
		double droplet_thermal_cond=0.19;
		double gradient=25000; //was 100000
		double temperature=250+273;//was 400+273 //523;

                long double phi;
                long double velocity;
                double rad;

                while(true){

                        for (int i=0;i<3;++i) StartPosition_temp[i]=StartPosition[i];
                        std::cout << "Pos_temp_a: (" << StartPosition_temp[0] << ", " << StartPosition_temp[1] << ", " << StartPosition_temp[2] << ")\n";

			while (true) {
				break;
				rad = (5e-6)*RandomNumber();
				double flow = RandomNumber();

					if (flow_rate==30) {
						flow*=11;
						if (flow<5) {
							rad+=10e-6;
						} else if ((flow>5) && (flow<7.5)) {
							rad+=15e-6;
						} else if ((flow>7.5) && (flow<9.5)) {
							rad+=20e-6;
						} else if ((flow>9.5) && (flow<10.5)) {
							rad+=25e-6;
						} else if ((flow>10.5) && (flow<10.8)) {
							rad+=30e-6;
						} else if ((flow>10.8) && (flow<10.9)) {
							rad+=35e-6;
						} else {
							rad+=40e-6;
						}
					} else {
						flow*=13;
						if (flow<8) {
							rad+=10e-6;
						} else if ((flow>8) && (flow<11)) {
							rad+=15e-6;
						} else if ((flow>11) && (flow<12.1)) {
							rad+=20e-6;
						} else if ((flow>12.1) && (flow<12.7)) {
							rad+=25e-6;
						} else if ((flow>12.7) && (flow<12.9)) {
							rad+=30e-6;
						} else {
							rad+=35e-6;
						}
				}
				if ((rad>min_radius) && (rad<max_radius)) break;
			}
			phi = my::math::Pi*(10*RandomNumber()-5)/180+direction_phi; //narrow nozzle

			rad=(15*RandomNumber()+30)*1e-6;

			s=4.5*air_viscosity/(density*rad*rad);

			long double min_velocity = 0.99999*(Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s/cos(phi);
			long double max_velocity = 1.00001*(Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s/cos(phi);

			factor*=(max_velocity-min_velocity)/(15.);

			velocity = (max_velocity-min_velocity)*RandomNumber()+min_velocity;

			long double theta = (theta_max-theta_min)*RandomNumber()+theta_min+direction_theta;
			factor*=(theta_max-theta_min)/(70.*my::math::Pi/180.);

			v0[0] = velocity*cos(phi)*cos(theta);
			v0[1] = velocity*cos(phi)*sin(theta);
			v0[2] = velocity*sin(phi);
			//std::cout << "v0_a: (" << v0[0] << ", " << v0[2] << ", " << v0[1] << ")\n";

			double d = StartPosition_temp[Parameter.open_boundary_direction]-((Partition.Max(Parameter.open_boundary_direction))*Parameter.GridDelta+6e-3);
			d = -d;

			double C1=(v0[2]-a_gr/s)/s;
			double a = C1;// -C1;
			double b = a_gr/s;
			double c = C1;
			double w = (v0[2]*s/a_gr-1)*exp((a-d)*s*s/a_gr); //a*s/b*exp(s*(a-d)/b);
			if (w>=std::numeric_limits<double>::max()) {
				t=0;
				} else {
					//if (w<-1/exp(1)){
					t=-(a-d)/b; // t = -(A-d)/B
					if (w>=-1/exp(1)) {
						long double lamb0 = LambertW<0>(w)/s;
						t += lamb0;
					} else {
						t = -t;
					}
			}
			//std::cout << "t: " << t << std::endl;
			//std::cout << "s: " << s << std::endl;

			double x, y;
			x=1/s*v0[0]*(1-exp(-s*t));
			y=1/s*v0[1]*(1-exp(-s*t));

			StartPosition_temp[SprayDirection]+=x;
			StartPosition_temp[AdditionalDirection]+=y;
			StartPosition_temp[Parameter.open_boundary_direction]+=d;
			//std::cout << "Pos_temp_b: (" << StartPosition_temp[0] << ", " << StartPosition_temp[1] << ", " << StartPosition_temp[2] << ")\n";

			v0[0]*=exp(-s*t);
			v0[1]*=exp(-s*t);
			v0[2]=(v0[2]-a_gr/s)*exp(-s*t) + a_gr/s;

			//std::cout << "v0_b: (" << v0[0] << ", " << v0[2] << ", " << v0[1] << ")\n";

			// Enter thermal region
			rad-=10e-6;

			double F_th = (3*my::math::Pi*air_viscosity*air_viscosity*rad/air_density)*(3*air_thermal_cond/(2*air_thermal_cond+droplet_thermal_cond))*gradient/temperature;
			double a_th = F_th*3/(4*my::math::Pi*density*rad*rad*rad);

			d = StartPosition_temp[Parameter.open_boundary_direction]-(Partition.Max(Parameter.open_boundary_direction))*Parameter.GridDelta;
			d = -d;

                        double C2=(v0[2]-(a_gr+a_th)/s)/s;
                        a = C2;
                        b = (a_gr+a_th)/s;
                        c = C2;

                        w = (v0[2]*s/(a_gr+a_th)-1)*exp((a-d)*s*s/(a_gr+a_th)); //a*s/b*exp(s*(a-d)/b);

			if (w>=std::numeric_limits<double>::max()) {
				t=0;
				} else {
					//if (w<-1/exp(1)){
					t=-(a-d)/b; // t = -(A-d)/B
					if (w>=-1/exp(1)) {
						long double lamb0 = LambertW<0>(w)/s;
						t += lamb0;
					} else {
						t = -t;
					}
			}

			x=1/s*v0[0]*(1-exp(-s*t));
			y=1/s*v0[1]*(1-exp(-s*t));

			StartPosition_temp[SprayDirection]+=x;
			StartPosition_temp[AdditionalDirection]+=y;
			StartPosition_temp[Parameter.open_boundary_direction]+=d;
			//std::cout << "Pos_temp_c: (" << StartPosition_temp[0] << ", " << StartPosition_temp[1] << ", " << StartPosition_temp[2] << ")\n";

		for (int i=0;i<3;++i) position[i]=StartPosition_temp[i];

                if ((position[0] < (Partition.Min(0))*Parameter.GridDelta) || (position[0] > (Partition.Max(0))*Parameter.GridDelta) || (position[2] < (Partition.Min(2))*Parameter.GridDelta) || (position[2] > (Partition.Max(2))*Parameter.GridDelta)) {
                } else {
                        break;
                }
                }

                direction[SprayDirection] = v0[0]*exp(-s*t);
                direction[AdditionalDirection] = v0[1]*exp(-s*t);
                direction[Parameter.open_boundary_direction] = (v0[2]-a_gr/s)*exp(-s*t) + a_gr/s;
                //std::cout << "v0_c: (" << direction[0] << ", " << direction[1] << ", " << direction[2] << ")\n";

                double v_end=0;
                for (int i=0;i<3;++i) v_end+=direction[i]*direction[i];
                for (int i=0;i<3;++i) direction[i]/=sqrt(v_end);
                FluxFactor=factor;
                }

		template<class DropletType, class VecType1, class VecType2, class ParameterType, class PartitionType>
		inline void ESDDistribution(const DropletType& d, const VecType1& StartPosition, VecType2& Position, double& r, double& q, long double* Velocity, const ParameterType& Parameter, const PartitionType& Partition){

//			ofstream out;
//			out.open("macroscopic.txt");
//			out << "d_test_v" << ", " << "d_test_r" << "\n";
			double d_test_v; double d_test_r;
//			for (int qq=0;qq<10000;qq++) {
			// Do this for T=310C = 310+573 K
//			for (int i=0;i<3;i++) Position[i]=StartPosition[i];

			//----------------Find the radius distribution----------------------------------
			double volume_fraction=0.42*RandomNumber()+0.58;
			double r_min_inv=1/2.5e-6; double r_max_inv=1/55e-6;
			double r_max_inv_third=exp(log(r_max_inv)/3);
			double r_min_inv_third=exp(log(r_min_inv)/3);
//			r=1/pow(volume_fraction*(r_max_inv_third-r_min_inv_third)+r_min_inv_third,3.);
			double radius=1/pow(volume_fraction*(r_max_inv_third-r_min_inv_third)+r_min_inv_third,3.);
//			std::cout << "radius: " << radius << "\n";

//			double Rc     = 8.3144621;		// mN/k mol
//			double T	 = 310+273;			// K
//			double N_A	 = 6.02214129e23;	// /mol
//			double rho_a = 1.29;	//kg/m3
			double eta_a = 2.2e-5;	//Ns/m2
//			double gamma = 67;		// Pa/K
//			double M_L	 = 0.0461;	// kg/mol
			double rho_d = 789;		// kg/m3
//			double del_T = 5;		// temperature difference
//			double s0	 = -1.117e-3;	// ms(-0.5)
//			double r0	 = 64.65;		//s(-0.5)
//			double q0	 = 373e-12;		//100e-12 for water
//			double q1	 = 89.1;		//400 for water
//			double t_dr	 = 1;
//			double K	 = q0*del_T*(1+2*q1*r);
//			double t_life = 4*radius*radius/(q0*del_T);

//			double c1	 = 4*gamma*M_L*D_vf*del_T/(rho_d*Rc*T);
//			double c1	 = 939;		// 110 for water;
//			double c2    = 0.276*(exp(log(rho_a/(eta_a*D_vf*D_vf))/6));// rho_a/(eta_a*D_vf*D_vf));
//			double c2	 = 8.91e-5;	// 6.452e-5 for water;
//			double D_vf  = Rc*T/(N_A*6*my::math::Pi*0.00116*radius);
//			std::cout << "t_life: " << t_life << "\n";
			//------------------------------------------------------------------------------

			//----------------Calculate the charge given r----------------------------------
			double gamma_d = 0.022;
			double permittivity = 8.854187817e-12;
			double qd=0.58*8*my::math::Pi*sqrt(gamma_d*permittivity*radius*radius*radius);
//			std::cout << "charge: " << q_d << "\n";
			//------------------------------------------------------------------------------

			//----------------Find the initial droplet position (cylidrical)----------------
			long double theta_i=0.5*my::math::Pi*RandomNumber()*0.5; // 45 degree spray cone
			long double phi=2*my::math::Pi*RandomNumber();
			double location_radius=0.5;
			double z_star1	= 1-location_radius*cos(theta_i);	// height
			double r_star1	= location_radius*sin(theta_i);	// radius

			double mass  = 4*my::math::Pi*rho_d*radius*radius*radius/3;


			//---FIND THE INITIAL AND THERMALDROPLET ELECTRICAL FORCES AND ACCELLERATIONS---

			double H	    = StartPosition[Parameter.open_boundary_direction]; //270 mm or 0.27 m
			double Phi_0    = 10e3;					//V
			double R	    = 1e-3;					//outer radius of the nozzle (guess)
			double K_V	    = 1-exp(-0.021*H/R);	//non-dimensional related to H/R ratio
			double Phi_star = K_V/(log(4*H/R));
			double E_e		= Phi_0*Phi_star/H;

			//-----------Find expected initial electric force when E-field acts alone--------
			double plusz1	   =	1+z_star1;//1+z_star;
			double minusz1	   =	1-z_star1;//1-z_star;
			double rootplusz1  =	std::sqrt(r_star1*r_star1+plusz1*plusz1);
			double rootminusz1 =	std::sqrt(r_star1*r_star1+minusz1*minusz1);
			double E_v1	= E_e*(1/rootminusz1+1/rootplusz1);
			double E_r1	= E_e*(plusz1/rootplusz1-minusz1/rootminusz1)/r_star1;
			double theta=atan(E_r1/E_v1);

			//-----------Find expected final electric force where E-field acts alone---------
			double t_heat	   = 	10e-3;
			double r_star2	   =	(E_r1*(H-t_heat))/E_v1;
			double z_star2	   =	0;//t_heat/H;
//			std:cout << "r_star2 = " << r_star2 << ", z_star2 = " << z_star2 << "\n";
			double plusz2	   =	1+z_star2;//1+z_star;
			double minusz2	   =	1-z_star2;//1-z_star;
			double rootplusz2  =	std::sqrt(r_star2*r_star2+plusz2*plusz2);
			double rootminusz2 =	std::sqrt(r_star2*r_star2+minusz2*minusz2);
			double Eth_v2	   = 	E_e*(1/rootminusz2+1/rootplusz2);
			double Eth_r2	   = 	E_e*(plusz2/rootplusz2-minusz2/rootminusz2)/r_star2;

			//-----------Use the first and second to come up with linear dependence----------
			double E_v 		= (E_v1-Eth_v2)/(H-t_heat);
			double E_r 		= (theta<1e-20)?0:(E_r1-Eth_r2)/((H-t_heat)*tan(theta));

			double Fe_v		= qd*E_v;
			double Fe_r		= qd*E_r;

			double ae_v		= Fe_v/mass;		// initial dependent component
			double ae_r 	= Fe_r/mass;		// initial dependent component

			//------------Calculate the required constant component of electric force--------
			double Fe_v1 = qd*E_v1;
			double Fe_r1 = qd*E_r1;
			double ae_v1 = Fe_v1/mass;		// Initial constant component
			double ae_r1 = (theta<1e-20)?0:Fe_r1/mass;		// Initial constant component

//			std::cout << "initial guess at ae_th_v: " << qd*Eth_v2/mass << "\n";
//			std::cout << "initial guess at ae_th_r: " << qd*Eth_r2/mass << "\n";
			//-------------------------------------------------------------------------------


//			std::cout << "ae_r = " << ae_r << ", ae_v = " << ae_v << "\n";


			//------------------------------------------------------------------------------

			//--------------Find all forces acting on the droplet---------------------------
			// Find the droplet mass - knowns: radius, rho_d
//			double mass  = 4*my::math::Pi*rho_d*radius*radius*radius/3;

			// Gravity force component acceleration
			double g     = 9.81;	//m/s
//			std::cout << "mass = " << mass << ", Fg = " << g*mass << ", ";

			// Stokes force component acceleration - knowns: rho_d, eta_a
			double s_f	 = (4.5*eta_a)/(rho_d*radius*radius);
//			std::cout << "Fs = " << s_f*mass << "\n";

			// Electric force component acceleration - knowns: q_d, r_star, z_star;
//			double H	    = StartPosition[Parameter.open_boundary_direction]; //270 mm or 0.27 m
//			double Phi_0    = 10e3;					//V
//			double R	    = 1e-3;					//outer radius of the nozzle (guess)
//			double K_V	    = 1-exp(-0.021*H/R);	//non-dimensional related to H/R ratio

			double d0_v	= 0;
			double d0_r	= 0;
			double v_0	= 0;
			double v0_v	= v_0*cos(theta);
			double v0_r	= v_0*sin(theta);


			// This is for both dimensions
//			double Phi_star = K_V/(log(4*H/R));
//			double E_e		= Phi_0*Phi_star/H;
//			double F_e		= q_d*E_e;
//			double a_e		= F_e/mass;
//
//			// Separate the vertical from radial
//			double a_e_v	= 2*H*a_e;
//			double a_e_r	= 0.5*H*a_e;
//			std::cout << "Fe_r = " << a_e_r*F_e << ", Fe_v = " << a_e_v*F_e <<"\n";
			//------------------------------------------------------------------------------

			//--------------Find new position of droplet after t----------------------------
			//First, separate the forces due to velocity/displacement dependences

			double a_v = g+ae_v1;//+2*a_e;	// independent acceleration (mainly gravity) - vertical only
			double b   = s_f;				// velocity dependent acceleration	- Stokes force - vert and rad
			double c_v = ae_v;			// displacement dependent acceleration - vertical E-force
			double a_r = ae_r1;	//+2*a_e;	// independent acceleration - a_e component
			double c_r = ae_r;	// displacement dependent acceleration - radial E-force

//			std::cout << "a_v = " << a_v << ", b = " << b << ", c_v = " << c_v << "\n";

			//--------------Start vertical onlgoings ---------------------------------------
			long double t_drop;

			double v0th_v;//=0;
			double v0th_r;//=0;
			if (b*b-4*c_v<0) {
				t_drop=0;
				d_test_v	= H-t_heat;
				d_test_r	= r_star1*H;
				v0th_v=1;
				v0th_r=0;
			} else {
				double r1_v	= (-b+std::sqrt(b*b-4*c_v))/(2);
				double r2_v	= (-b-std::sqrt(b*b-4*c_v))/(2);

//				std::cout << "r1_v = " << r1_v << ", r2_v = " << r2_v << "\n";

				double B1_v = d0_v*(r2_v)/(r2_v-r1_v);
				double A1_v	= d0_v-B1_v;
				double B2_v = (v0_v+d0_v)/(r2_v-r1_v);
				double A2_v = -B2_v;;
				double C_v	= a_v/(r1_v*r2_v);
				double B3_v	= -a_v/(r2_v*(r2_v-r1_v));
//				double B3_v	= a_v*(r1_v+r2_v-1)/(r1_v*r2_v*(1-r1_v));
				double A3_v	= -(B3_v+C_v);
//				double A3_v	= -a_v/(r1_v*r2_v)-B3_v;

//				d_test_v = (A1_v+A2_v+A3_v)*exp(r1_v*t_check)+(B1_v+B2_v+B3_v)*exp(r2_v*t_check)+C_v;
//				std::cout << "t_check: " << t_check << ", d_test_v: " << d_test_v << "\n";

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
//				std::cout << "Iterations: " << iteration_stop << ", t_check: " << t_check << ", d_test_v: " << d_test_v << "\n";

				t_drop=t_check;
				v0th_v		= (A1_v+A2_v+A3_v)*(exp(r1_v*t_drop)-1)/r1_v+(B1_v+B2_v+B3_v)*(exp(r2_v*t_drop)-1)/r2_v+C_v*t_drop+v0_v;
//				std::cout << "Fist: " << (A1_v+A2_v+A3_v)*(exp(r1_v*t_drop)-1)/r1_v << "\n";
//				std::cout << "second: " << (B1_v+B2_v+B3_v)*(exp(r2_v*t_drop)-1)/r2_v << "\n";
//				std::cout << "third: " << C_v*t_drop << "\n";

				//---------Now time required to reach heat zone is known------------------------

				//---------Can now calculate the radial displacement----------------------------
//				double a_r = ae_r1;	//+2*a_e;	// independent acceleration - a_e component
//				double c_r = ae_r;	// displacement dependent acceleration - radial E-force
	//			std::cout << "a_r = " << a_r << ", b = " << b << ", c_r = " << c_r << "\n";

				if (b*b-4*c_r<0) {
					t_drop=0;
//					dth_test_v	= t_heat;
					d_test_r	= 0;
//					v_final_v	= 1;
					v0th_r	= 0;
				} else {
					double r1_r	= (-b+std::sqrt(b*b-4*c_r))/(2);
					double r2_r	= (-b-std::sqrt(b*b-4*c_r))/(2);

		//			std::cout << "r1_r = " << r1_r << ", r2_r = " << r2_r << "\n";

					double B1_r = d0_r*(r2_r)/(r2_r-r1_r);
					double A1_r	= d0_r-B1_r;
					double B2_r = (v0_r+d0_r)/(r2_r-r1_r);
					double A2_r = -B2_r;;
					double C_r	= a_r/(r1_r*r2_r);
					double B3_r	= -a_r/(r2_r*(r2_r-r1_r));
		//				double B3_v	= a_v*(r1_v+r2_v-1)/(r1_v*r2_v*(1-r1_v));
					double A3_r	= -(B3_r+C_r);
		//				double A3_v	= -a_v/(r1_v*r2_v)-B3_v;

					d_test_r = (A1_r+A2_r+A3_r)*exp(r1_r*t_drop)+(B1_r+B2_r+B3_r)*exp(r2_r*t_drop)+C_r;
					v0th_r	 = (A1_r+A2_r+A3_r)*(exp(r1_r*t_drop)-1)/r1_r+(B1_r+B2_r+B3_r)*(exp(r2_r*t_drop)-1)/r2_r+C_r*t_drop+v0_r;
				}

			}
//			std::cout << "d_test_r = " << d_test_r << ", d_test_v = " << d_test_v << "\n";
//			std::cout << "v0th_r = " << v0th_r << ", v0th_v = " << v0th_v << "\n";
//			double d_r	= C1_r*exp(r1_r*t_drop) + C2_r*exp(r2_r*t_drop) + C3_r;
//			std::cout << "t = " << t_drop << ", d = " << d_r << "\n\n";
			//------------------------------------------------------------------------------
//			}
//			std::cout << "t_drop: " << t_drop << ", d_test_r: " << d_test_r << ", d_test_v: " << d_test_v << "\n";
			//------------------------------------------------------------------------------
			//------------------------------------------------------------------------------
			//------------------      THERMAL ZONE CALCULATIONS      -----------------------
			//------------------------------------------------------------------------------
			//------------------------------------------------------------------------------
//			//---------Now we know z and r on the cusp of the thermal zone------------------
//			//---------Reset parameters to include thermal components and re-run -----------

			//--------CALCULATE THERMAL EFFECTS IN DROPLET SIZE REDUCTION-------------------
			double dth_test_v;
			double dth_test_r;
			double q0 	  	= 373e-12;		// (373 um^2) 	 for water: 88e-12 (m^2)
			double q1 	  	= 89.1;			// (8.91e-5 /um) for water: 4.3e3  (/m)
			double del_T  	= 100000;
			double dK 	  	= q0*del_T*(1+2*q1*radius);
//			std::cout << "dK = " << dK << "\n";
//			double K	  	= dK*t_drop;
			double r_new  	= radius-radius*t_drop*exp(log(dK)/3);
//			std::cout << "dr = " << t_drop*exp(log(dK)/3) << "\n";
			r=r_new;
			double qd_new 	= 0.58*8*my::math::Pi*sqrt(gamma_d*permittivity*r_new*r_new*r_new);
			double mass_new	= 4*my::math::Pi*rho_d*r_new*r_new*r_new/3;
//			std::cout << "dK = " << dK << ", K = " << K << "\n";
//			std::cout << "r_new= " << r_new << ", q_new= " << qd_new << ", mass= " << mass_new << "\n";
//
//			//---------CALCULATE THE ENERGIES ASSUMING LINEAR REDUCTION TO FINAL------------
//			//---------Electric force F_e is constant in this region and small--------------
//
			double rth_star	= 	d_test_r/H;		//d_r/H;
			double zth_star	= 	1-d_test_v/H;	//1-d_v/H;
//			std::cout << "rth_star = " << rth_star << ", zth_star = " << zth_star << "\n";

//			double v0th_v		= (A1_v+A2_v+A3_v)*(exp(r1_v*t_drop)-1)/r1_v+(B1_v+B2_v+B3_v)*(exp(r2_v*t_drop)-1)/r2_v+C_v*t_drop+v0_v;
//			double v0th_r		= (A1_r+A2_r+A3_r)*(exp(r1_r*t_drop)-1)/r1_r+(B1_r+B2_r+B3_r)*(exp(r2_r*t_drop)-1)/r2_r+C_r*t_drop+v0_r;
			double d0th_v		= 0;//d_test_v;
			double d0th_r		= 0;//d_test_r;
			double tth_drop		= t_drop;

			double thplusz	 	=	1+zth_star;//1+z_star;
			double thminusz	 	=	1-zth_star;//1-z_star;
			double throotplusz  =	std::sqrt(rth_star*rth_star+thplusz*thplusz);
			double throotminusz =	std::sqrt(rth_star*rth_star+thminusz*thminusz);

//			double H	    = StartPosition[Parameter.open_boundary_direction]; //270 mm or 0.27 m
//			double Phi_0    = 10e3;					//V
//			double R	    = 1e-3;					//outer radius of the nozzle (guess)
//			double K_V	    = 1-exp(-0.021*H/R);	//non-dimensional related to H/R ratio
//			double Phi_star = K_V/(log(4*H/R));
//			double E_e		= Phi_0*Phi_star/H;

			double Eth_v	= E_e*(1/throotminusz+1/throotplusz);
			double Eth_r	= E_e*(thplusz/throotplusz-thminusz/throotminusz)/rth_star;
			double theta_th	= atan(Eth_r/Eth_v);
//			std::cout << "Eth_r = " << Eth_r << ", Eth_v = " << Eth_v << ", theta_th = " << theta_th << "\n";

			double Feth_v	= qd_new*Eth_v;
			double Feth_r	= qd_new*Eth_r;
//			std::cout << "Feth_r = " << Feth_r << ", Feth_v = " << Feth_v << "\n";

			double aeth_v	= Feth_v/mass_new;
			double aeth_r 	= Feth_r/mass_new;
//			std::cout << "aeth_r = " << aeth_r << ", aeth_v = " << aeth_v << "\n";

			//--------------Find all other forces acting on the droplet---------------------
			// Find the droplet mass - knowns: radius, rho_d
//			double mass  = 4*my::math::Pi*rho_d*radius*radius*radius/3;

			// Gravity force component acceleration
//			double g     = 9.81;	//m/s
//			std::cout << "mass = " << mass << ", Fg = " << g*mass << ", ";

			// Stokes force component acceleration - knowns: rho_d, eta_a
			double sth_f = (4.5*eta_a)/(rho_d*r_new*r_new);
//			std::cout << "Fs = " << s_f*mass << "\n";

			// Thermophoretic force
//			double eta_a;
			double kappa_a	= 0.025;
			double kappa_d	= 0.19;
			double grad_T	= 100000;
			double T		= 523;
			double rho_a	= 1.29;

			double F_th	 = 3*my::math::Pi*eta_a*eta_a*r_new*3*kappa_a*grad_T/(rho_a*T*(2*kappa_a+kappa_d));
			double a_th	 = F_th/mass_new;

			// Electric force component acceleration - knowns: q_d, r_star, z_star;
//			double H	    = StartPosition[Parameter.open_boundary_direction]; //270 mm or 0.27 m
//			double Phi_0    = 10e3;					//V
//			double R	    = 1e-3;					//outer radius of the nozzle (guess)
//			double K_V	    = 1-exp(-0.021*H/R);	//non-dimensional related to H/R ratio

			double ath_v = g+aeth_v-a_th;//+2*a_e;	// independent acceleration (gravity, initial e-force, thermal force)
			double bth 	 = sth_f;					// velocity dependent acceleration		- Stokes force
			double cth_v = aeth_v/t_heat;			// displacement dependent acceleration 	- vertical E-force
			double ath_r = (theta_th<1e-20)?0:aeth_r;						  	// independent acceleration - Initial electric force
//			double b is known from vertical calculation
			double cth_r = (theta_th<1e-20)?0:aeth_r/(t_heat*tan(theta_th));	// displacement dependent acceleration - Radial E-force
//			std::cout << "ath_v = " << ath_v << ", bth = " << bth << ", cth_v = " << cth_v << "\n";

			double v_final_v=0;
			double v_final_r=0;

			if (bth*bth-4*cth_v<0) {
				tth_drop=0;
				dth_test_v	= t_heat;
				dth_test_r	= d_test_r;
				v_final_v	= 1;
				v_final_r	= 0;
			} else {
				double r1th_v	= (-bth+std::sqrt(bth*bth-4*cth_v))/(2);
				double r2th_v	= (-bth-std::sqrt(bth*bth-4*cth_v))/(2);

//				std::cout << "r1_v = " << r1_v << ", r2_v = " << r2_v << "\n";

				double B1th_v 	= d0th_v*(r2th_v)/(r2th_v-r1th_v);
				double A1th_v	= d0th_v-B1th_v;
				double B2th_v 	= (v0th_v+d0th_v)/(r2th_v-r1th_v);
				double A2th_v 	= -B2th_v;;
				double Cth_v	= ath_v/(r1th_v*r2th_v);
				double B3th_v	= -ath_v/(r2th_v*(r2th_v-r1th_v));
//				double B3_v	= a_v*(r1_v+r2_v-1)/(r1_v*r2_v*(1-r1_v));
				double A3th_v	= -(B3th_v+Cth_v);
//				double A3_v	= -a_v/(r1_v*r2_v)-B3_v;

//				d_test_v = (A1_v+A2_v+A3_v)*exp(r1_v*t_check)+(B1_v+B2_v+B3_v)*exp(r2_v*t_check)+C_v;
//				std::cout << "t_check: " << t_check << ", d_test_v: " << d_test_v << "\n";

				double iteration_stopth	=10000;
				double tth_low			= 0;
				double tth_high			= t_drop;
				double tth_check		= (tth_low+tth_high)/2;
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
//				std::cout << "Iteration_th: " << iteration_stopth << ", tth_check: " << tth_check << ", dth_test_v: " << dth_test_v << "\n";

				tth_drop  = tth_check;
				v_final_v = (A1th_v+A2th_v+A3th_v)*(exp(r1th_v*tth_drop)-1)/r1th_v+(B1th_v+B2th_v+B3th_v)*(exp(r2th_v*tth_drop)-1)/r2th_v+Cth_v*tth_drop+v0th_v;


				//---------Now time required to reach surface is known--------------------------

				//---------Can now calculate the radial displacement----------------------------
//				double ath_r = aeth_r;						  	// independent acceleration - Initial electric force
//	//			double b is known from vertical calculation
//				double cth_r = aeth_r/(t_heat*tan(theta_th));	// displacement dependent acceleration - Radial E-force
//				std::cout << "ath_r = " << ath_r << ", bth = " << bth << ", cth_r = " << cth_r << "\n";

				if (bth*bth-4*cth_r<0) {
					tth_drop=0;
//					dth_test_v	= t_heat;
					dth_test_r	= 0;
//					v_final_v	= 1;
					v_final_r	= 0;
				} else {
					double r1th_r = (-bth+std::sqrt(bth*bth-4*cth_r))/(2);
					double r2th_r = (-bth-std::sqrt(bth*bth-4*cth_r))/(2);

		//			std::cout << "r1_r = " << r1_r << ", r2_r = " << r2_r << "\n";

					double B1th_r 	= d0th_r*(r2th_r)/(r2th_r-r1th_r);
					double A1th_r 	= d0th_r-B1th_r;
					double B2th_r 	= (v0th_r+d0th_r)/(r2th_r-r1th_r);
					double A2th_r 	= -B2th_r;;
					double Cth_r  	= ath_r/(r1th_r*r2th_r);
					double B3th_r 	= -ath_r/(r2th_r*(r2th_r-r1th_r));
	//				double B3_v	  	= a_v*(r1_v+r2_v-1)/(r1_v*r2_v*(1-r1_v));
					double A3th_r 	= -(B3th_r+Cth_r);
	//				double A3_v	  	= -a_v/(r1_v*r2_v)-B3_v;

					dth_test_r = (A1th_r+A2th_r+A3th_r)*exp(r1th_r*tth_drop)+(B1th_r+B2th_r+B3th_r)*exp(r2th_r*tth_drop)+Cth_r;
					v_final_r  = (A1th_r+A2th_r+A3th_r)*(exp(r1th_r*tth_drop)-1)/r1th_r+(B1th_r+B2th_r+B3th_r)*(exp(r2th_r*tth_drop)-1)/r2th_r+Cth_r*tth_drop+v0th_r;
	//				std::cout << "tth_check: " << tth_check << ", dth_test_r: " << dth_test_r << "\n";
				}
			}

//			std::cout << "v_final_r = " << v_final_r << ", v_final_v = " << v_final_v << "\n";

//			Velocity[0] = ((bth*bth-4*cth_v)<0)?0:v_final_r*cos(phi);
//			Velocity[1] = ((bth*bth-4*cth_v)<0)?-1:-v_final_v;
//			Velocity[2] = ((bth*bth-4*cth_v)<0)?0:v_final_r*sin(phi);

			Velocity[0] = v_final_r*cos(phi);
			Velocity[1] = -v_final_v;
			Velocity[2] = v_final_r*sin(phi);

			Position[0] = std::sqrt(dth_test_r+d_test_r)*cos(phi)+StartPosition[0];
			Position[1] = 0;
			Position[2] = std::sqrt(dth_test_r+d_test_r)*sin(phi)+StartPosition[2];

			r	= r_new;
			q	= qd_new;

//			double norm=0;
//            for (int i=0;i<3;i++) norm+=Velocity[i]*Velocity[i];
//            if (!(norm>0)) {
//            	std::cout << "dth_test_r: " << dth_test_r << ", v_final_r: " << v_final_r << ", " << Velocity[2] << "\n";
//            	std::cout << "Velocity: " << Velocity[0] << ", " << Velocity[1] << ", " << Velocity[2] << "\n";
//            	std::cout << "Position: " << Position[0] << ", " << Position[1] << ", " << Position[2] << "\n";
//            	std::cout << "radius = " << radius << ", charge = " << qd_new << "\n";
//            	std::cout << "tth_drop = " << tth_drop << ", v_final_v = " << v_final_v << "\n";
//            	std::cout << "t_drop = " << t_drop << ", d_test_v = " << d_test_v << "\n";
//            	std::cout << "r1_r = " << r1_r << ", r2_r = " << r2_r << "\n";
//            	std::cout << "ath_v = " << ath_v << ", bth = " << bth << ", cth_v = " << cth_v << "\n";
//            	std::cout << "ath_r = " << ath_r << ", cth_r = " << cth_r << "\n";
//
//            	std::cout << "d_test_r = " << d_test_r << ", v0th_r = " << v0th_r << ", (b*b-4*c_v<0) = " << (b*b-4*c_v<0) << "\n";
//
//
//            	std::cout << "a_v = " << a_v << ", b = " << b << ", c_v = " << c_v << "\n";
//            	std::cout << "a_r = " << a_r << ", c_r = " << c_r << "\n";
//            	std::cout << "theta_th = " << theta_th << "\n";
//
//            }

//			if (bth*bth-4*cth_v>0) {
//				std::cout << "here!\n";
//				double vtot=0;
//				for (int i=0;i<3;i++) vtot+=Velocity[i]*Velocity[i];
//				for (int i=0;i<3;i++) Velocity[i]/=std::sqrt(vtot);
//			}

//			double r_dist = std::sqrt(dth_test_r);
//			Position[0] = std::sqrt(dth_test_r+d_test_r)*cos(phi)+StartPosition[0];
//			Position[1] = 0;
//			Position[2] = std::sqrt(dth_test_r+d_test_r)*sin(phi)+StartPosition[2];
//
//			r	= r_new;
//			q	= qd_new;

//			std::cout << std::sqrt(dth_test_r + d_test_r) << ", " << phi << "\n";
//			std::cout << Position[0] << ", " << Position[2] << "\n";

//			std::cout << "Position: " << Position[0] << ", " << Position[1] << ", " << Position[2] << "\n";
//			std::cout << "Velocity: " << Velocity[0] << ", " << Velocity[1] << ", " << Velocity[2] << "\n";
//			std::cout << "radius:   " << r << "\n";
//			std::cout << "charge:   " << q << "\n";

		// This is for both dimensions
//		double Phi_star = K_V/(log(4*H/R));
//		double E_e		= Phi_0*Phi_star/H;
//		double F_e		= q_d*E_e;
//		double a_e		= F_e/mass;
//
//		// Separate the vertical from radial
//		double a_e_v	= 2*H*a_e;
//		double a_e_r	= 0.5*H*a_e;
//		std::cout << "Fe_r = " << a_e_r*F_e << ", Fe_v = " << a_e_v*F_e <<"\n";
			//std::cout << "\n";
		}

		template<class DropletType, class VecType1, class VecType2, class ParameterType, class PartitionType>
		inline void PSDDistribution(const DropletType& d, const VecType1& StartPosition, VecType2& Position, double& r, double& q, long double* Velocity, const ParameterType& Parameter, const PartitionType& Partition){


				ofstream out;
				out.open("timestuff.txt");

			out << "radius, v0, distance \n";

				for (int i=0;i<100000;i++) {
			//--------- SOME CONSTANTS -----------------------------------------------------
			//			double Rc     = 8.3144621;		// mN/k mol
			//			double T	 = 310+273;			// K
			//			double N_A	 = 6.02214129e23;	// /mol
			//			double rho_a = 1.29;	//kg/m3
						double eta_a = 2.2e-5;	//Ns/m2
			//			double gamma = 67;		// Pa/K
			//			double M_L	 = 0.0461;	// kg/mol
						double rho_d = 998;		// kg/m3
			//			double del_T = 5;		// temperature difference
			//			double s0	 = -1.117e-3;	// ms(-0.5)
			//			double r0	 = 64.65;		//s(-0.5)
			//			double q0	 = 373e-12;		//100e-12 for water
			//			double q1	 = 89.1;		//400 for water
			//			double t_dr	 = 1;
			//			double K	 = q0*del_T*(1+2*q1*r);
			//			double t_life = 4*radius*radius/(q0*del_T);

			//			double c1	 = 4*gamma*M_L*D_vf*del_T/(rho_d*Rc*T);
			//			double c1	 = 939;		// 110 for water;
			//			double c2    = 0.276*(exp(log(rho_a/(eta_a*D_vf*D_vf))/6));// rho_a/(eta_a*D_vf*D_vf));
			//			double c2	 = 8.91e-5;	// 6.452e-5 for water;
			//			double D_vf  = Rc*T/(N_A*6*my::math::Pi*0.00116*radius);
			//			std::cout << "t_life: " << t_life << "\n";
			//------------------------------------------------------------------------------

			//----------------Find the radius distribution----------------------------------
			double radius=4.5e-6*RandomNumber()+1.5e-6;
//			std::cout << "radius: " << radius << "\n";
			double mass  = 4*my::math::Pi*rho_d*radius*radius*radius/3;
			//------------------------------------------------------------------------------

			//----------------Calculate the charge given r----------------------------------
			double gamma_d = 0.022;
			double permittivity = 8.854187817e-12;
			double qd=0.58*8*my::math::Pi*sqrt(gamma_d*permittivity*radius*radius*radius);
//			std::cout << "charge: " << qd << "\n";
			//------------------------------------------------------------------------------

			//----------------Find the initial droplet direction(Cartesian)-----------------
			long double lateral_min = Partition.Min(2)*Parameter.GridDelta;
			long double lateral_max = Partition.Max(2)*Parameter.GridDelta;
			long double lateral = RandomNumber()*(lateral_max-lateral_min)-(lateral_max+lateral_min)/2;
			long double vertical = RandomNumber()*(lateral_max-lateral_min)/10;
			long double length = 280e-3;

					//normalized starting direction
			lateral  = lateral /(std::sqrt(lateral*lateral+vertical*vertical+length*length));
			vertical = vertical/(std::sqrt(lateral*lateral+vertical*vertical+length*length));
			length   = length/(std::sqrt(lateral*lateral+vertical*vertical+length*length));

//			std::cout << "direction: (" << length << ", " << vertical << ", " << lateral << "\n";
			//------------------------------------------------------------------------------

			//----------------Find the initial velocity ------------------------------------
			double v0 = 10000*RandomNumber();
			double v0_lateral  = v0*lateral;
			double v0_vertical = 0;//v0*vertical;
			double v0_length   = v0;//*length;
//			std::cout << "velocity: (" << v0_length << ", " << v0_vertical << ", " << v0_lateral << "\n";
			//------------------------------------------------------------------------------

			//----------------Find time required to reach thermal zone----------------------
			double g     = 9.81;
			double eta_d = 0.01;
			double F_s	 = 6*my::math::Pi*eta_a*radius;
			double a_s	 = F_s/mass;
//			std::cout << "g: " << g << ", a_s: " << a_s << "\n";

                        double height = 10e-3;
                double w = (g-a_s*v0_vertical)*exp(1-((a_s*(v0_vertical+a_s*height))/g))/g;
                double time_v;

                if (w>=std::numeric_limits<double>::max()) {
                        time_v=0;
                        } else {
                                time_v = (a_s*v0_vertical+a_s*a_s*height-g+g*LambertW<0>(w))/(g*a_s);
                }

//        	std::cout << "height = " << height << ", time = " << time_v << "\n";
                        //------------------------------------------------------------------------------

                        //----------------Find all directions after initial drop -----------------------
                double lat_length = v0_lateral*(1-exp(-a_s*time_v))/a_s;
                double sp_length  = v0_length*(1-exp(-a_s*time_v))/a_s;

                if ((time_v!=0)&&(sp_length>0.28)&&(sp_length<(0.28+Partition.Max(0)*Parameter.GridDelta))) {
                        double v_length   = v0_length*(exp(-a_s*time_v));
                        double v_vertical = (v0_vertical-g/a_s)*(exp(-a_s*time_v))+g/a_s;
                        out << radius << ", " << v_length << ", " << v_vertical << "\n";
                } else {
                        --i;
                }
//        	std::cout << "time = " << time_v << ", height = " << height << ", length = " << sp_length << ", lateral = " << lat_length << "\n";
                        //------------------------------------------------------------------------------

				}
			out.close();
			std::cout << "ready!";

			while (true) {}

//			double d_test_v; double d_test_r;
//			long double theta_i=0.5*my::math::Pi*RandomNumber()*0.5; // 45 degree spray cone
//			long double phi=2*my::math::Pi*RandomNumber();
//			double location_radius=0.5;
//			double z_star1	= 1-location_radius*cos(theta_i);	// height
//			double r_star1	= location_radius*sin(theta_i);	// radius
//
////			double mass  = 4*my::math::Pi*rho_d*radius*radius*radius/3;
//
//
//			//---FIND THE INITIAL AND THERMALDROPLET ELECTRICAL FORCES AND ACCELLERATIONS---
//
//			double H	    = StartPosition[Parameter.open_boundary_direction]; //270 mm or 0.27 m
//			double Phi_0    = 10e3;					//V
//			double R	    = 1e-3;					//outer radius of the nozzle (guess)
//			double K_V	    = 1-exp(-0.021*H/R);	//non-dimensional related to H/R ratio
//			double Phi_star = K_V/(log(4*H/R));
//			double E_e		= Phi_0*Phi_star/H;
//
//			//-----------Find expected initial electric force when E-field acts alone--------
//			double plusz1	   =	1+z_star1;//1+z_star;
//			double minusz1	   =	1-z_star1;//1-z_star;
//			double rootplusz1  =	std::sqrt(r_star1*r_star1+plusz1*plusz1);
//			double rootminusz1 =	std::sqrt(r_star1*r_star1+minusz1*minusz1);
//			double E_v1	= E_e*(1/rootminusz1+1/rootplusz1);
//			double E_r1	= E_e*(plusz1/rootplusz1-minusz1/rootminusz1)/r_star1;
//			double theta=atan(E_r1/E_v1);
//
//			//-----------Find expected final electric force where E-field acts alone---------
//			double t_heat	   = 	10e-3;
//			double r_star2	   =	(E_r1*(H-t_heat))/E_v1;
//			double z_star2	   =	0;//t_heat/H;
////			std:cout << "r_star2 = " << r_star2 << ", z_star2 = " << z_star2 << "\n";
//			double plusz2	   =	1+z_star2;//1+z_star;
//			double minusz2	   =	1-z_star2;//1-z_star;
//			double rootplusz2  =	std::sqrt(r_star2*r_star2+plusz2*plusz2);
//			double rootminusz2 =	std::sqrt(r_star2*r_star2+minusz2*minusz2);
//			double Eth_v2	   = 	E_e*(1/rootminusz2+1/rootplusz2);
//			double Eth_r2	   = 	E_e*(plusz2/rootplusz2-minusz2/rootminusz2)/r_star2;
//
//			//-----------Use the first and second to come up with linear dependence----------
//			double E_v 		= (E_v1-Eth_v2)/(H-t_heat);
//			double E_r 		= (theta<1e-20)?0:(E_r1-Eth_r2)/((H-t_heat)*tan(theta));
//
//			double Fe_v		= qd*E_v;
//			double Fe_r		= qd*E_r;
//
//			double ae_v		= Fe_v/mass;		// initial dependent component
//			double ae_r 	= Fe_r/mass;		// initial dependent component
//
//			//------------Calculate the required constant component of electric force--------
//			double Fe_v1 = qd*E_v1;
//			double Fe_r1 = qd*E_r1;
//			double ae_v1 = Fe_v1/mass;		// Initial constant component
//			double ae_r1 = (theta<1e-20)?0:Fe_r1/mass;		// Initial constant component
//
////			std::cout << "initial guess at ae_th_v: " << qd*Eth_v2/mass << "\n";
////			std::cout << "initial guess at ae_th_r: " << qd*Eth_r2/mass << "\n";
//			//-------------------------------------------------------------------------------
//
//
////			std::cout << "ae_r = " << ae_r << ", ae_v = " << ae_v << "\n";
//
//
//			//------------------------------------------------------------------------------
//
//			//--------------Find all forces acting on the droplet---------------------------
//			// Find the droplet mass - knowns: radius, rho_d
////			double mass  = 4*my::math::Pi*rho_d*radius*radius*radius/3;
//
//			// Gravity force component acceleration
////			double g     = 9.81;	//m/s
////			std::cout << "mass = " << mass << ", Fg = " << g*mass << ", ";
//
//			// Stokes force component acceleration - knowns: rho_d, eta_a
//			double s_f	 = (4.5*eta_a)/(rho_d*radius*radius);
////			std::cout << "Fs = " << s_f*mass << "\n";
//
//			// Electric force component acceleration - knowns: q_d, r_star, z_star;
////			double H	    = StartPosition[Parameter.open_boundary_direction]; //270 mm or 0.27 m
////			double Phi_0    = 10e3;					//V
////			double R	    = 1e-3;					//outer radius of the nozzle (guess)
////			double K_V	    = 1-exp(-0.021*H/R);	//non-dimensional related to H/R ratio
//
//			double d0_v	= 0;
//			double d0_r	= 0;
//			double v_0	= 0;
//			double v0_v	= v_0*cos(theta);
//			double v0_r	= v_0*sin(theta);
//
//
//			// This is for both dimensions
////			double Phi_star = K_V/(log(4*H/R));
////			double E_e		= Phi_0*Phi_star/H;
////			double F_e		= q_d*E_e;
////			double a_e		= F_e/mass;
////
////			// Separate the vertical from radial
////			double a_e_v	= 2*H*a_e;
////			double a_e_r	= 0.5*H*a_e;
////			std::cout << "Fe_r = " << a_e_r*F_e << ", Fe_v = " << a_e_v*F_e <<"\n";
//			//------------------------------------------------------------------------------
//
//			//--------------Find new position of droplet after t----------------------------
//			//First, separate the forces due to velocity/displacement dependences
//
//			double a_v = g+ae_v1;//+2*a_e;	// independent acceleration (mainly gravity) - vertical only
//			double b   = s_f;				// velocity dependent acceleration	- Stokes force - vert and rad
//			double c_v = ae_v;			// displacement dependent acceleration - vertical E-force
//			double a_r = ae_r1;	//+2*a_e;	// independent acceleration - a_e component
//			double c_r = ae_r;	// displacement dependent acceleration - radial E-force
//
////			std::cout << "a_v = " << a_v << ", b = " << b << ", c_v = " << c_v << "\n";
//
//			//--------------Start vertical onlgoings ---------------------------------------
//			long double t_drop;
//
//			double v0th_v;//=0;
//			double v0th_r;//=0;
//			if (b*b-4*c_v<0) {
//				t_drop=0;
//				d_test_v	= H-t_heat;
//				d_test_r	= r_star1*H;
//				v0th_v=1;
//				v0th_r=0;
//			} else {
//				double r1_v	= (-b+std::sqrt(b*b-4*c_v))/(2);
//				double r2_v	= (-b-std::sqrt(b*b-4*c_v))/(2);
//
////				std::cout << "r1_v = " << r1_v << ", r2_v = " << r2_v << "\n";
//
//				double B1_v = d0_v*(r2_v)/(r2_v-r1_v);
//				double A1_v	= d0_v-B1_v;
//				double B2_v = (v0_v+d0_v)/(r2_v-r1_v);
//				double A2_v = -B2_v;;
//				double C_v	= a_v/(r1_v*r2_v);
//				double B3_v	= -a_v/(r2_v*(r2_v-r1_v));
////				double B3_v	= a_v*(r1_v+r2_v-1)/(r1_v*r2_v*(1-r1_v));
//				double A3_v	= -(B3_v+C_v);
////				double A3_v	= -a_v/(r1_v*r2_v)-B3_v;
//
////				d_test_v = (A1_v+A2_v+A3_v)*exp(r1_v*t_check)+(B1_v+B2_v+B3_v)*exp(r2_v*t_check)+C_v;
////				std::cout << "t_check: " << t_check << ", d_test_v: " << d_test_v << "\n";
//
//				double iteration_stop=10000;
//				double t_low  = 0;
//				double t_high = 1;
//				double t_check = (t_low+t_high)/2;
//				while (true){
//					d_test_v = (A1_v+A2_v+A3_v)*exp(r1_v*t_check)+(B1_v+B2_v+B3_v)*exp(r2_v*t_check)+C_v;
//					if (t_high-t_low<1e-100) break;
//					if (d_test_v > (H-t_heat)) {
//						t_high=t_check;
//						t_check=(t_high+t_low)/2;
//					} else if (d_test_v < (H-t_heat)) {
//						t_low=t_check;
//						t_check=(t_high+t_low)/2;
//					}
//					if (iteration_stop==0) break;
//					iteration_stop--;
//				}
////				std::cout << "Iterations: " << iteration_stop << ", t_check: " << t_check << ", d_test_v: " << d_test_v << "\n";
//
//				t_drop=t_check;
//				v0th_v		= (A1_v+A2_v+A3_v)*(exp(r1_v*t_drop)-1)/r1_v+(B1_v+B2_v+B3_v)*(exp(r2_v*t_drop)-1)/r2_v+C_v*t_drop+v0_v;
////				std::cout << "Fist: " << (A1_v+A2_v+A3_v)*(exp(r1_v*t_drop)-1)/r1_v << "\n";
////				std::cout << "second: " << (B1_v+B2_v+B3_v)*(exp(r2_v*t_drop)-1)/r2_v << "\n";
////				std::cout << "third: " << C_v*t_drop << "\n";
//
//				//---------Now time required to reach heat zone is known------------------------
//
//				//---------Can now calculate the radial displacement----------------------------
////				double a_r = ae_r1;	//+2*a_e;	// independent acceleration - a_e component
////				double c_r = ae_r;	// displacement dependent acceleration - radial E-force
//	//			std::cout << "a_r = " << a_r << ", b = " << b << ", c_r = " << c_r << "\n";
//
//				if (b*b-4*c_r<0) {
//					t_drop=0;
////					dth_test_v	= t_heat;
//					d_test_r	= 0;
////					v_final_v	= 1;
//					v0th_r	= 0;
//				} else {
//					double r1_r	= (-b+std::sqrt(b*b-4*c_r))/(2);
//					double r2_r	= (-b-std::sqrt(b*b-4*c_r))/(2);
//
//		//			std::cout << "r1_r = " << r1_r << ", r2_r = " << r2_r << "\n";
//
//					double B1_r = d0_r*(r2_r)/(r2_r-r1_r);
//					double A1_r	= d0_r-B1_r;
//					double B2_r = (v0_r+d0_r)/(r2_r-r1_r);
//					double A2_r = -B2_r;;
//					double C_r	= a_r/(r1_r*r2_r);
//					double B3_r	= -a_r/(r2_r*(r2_r-r1_r));
//		//				double B3_v	= a_v*(r1_v+r2_v-1)/(r1_v*r2_v*(1-r1_v));
//					double A3_r	= -(B3_r+C_r);
//		//				double A3_v	= -a_v/(r1_v*r2_v)-B3_v;
//
//					d_test_r = (A1_r+A2_r+A3_r)*exp(r1_r*t_drop)+(B1_r+B2_r+B3_r)*exp(r2_r*t_drop)+C_r;
//					v0th_r	 = (A1_r+A2_r+A3_r)*(exp(r1_r*t_drop)-1)/r1_r+(B1_r+B2_r+B3_r)*(exp(r2_r*t_drop)-1)/r2_r+C_r*t_drop+v0_r;
//				}
//
//			}
////			std::cout << "d_test_r = " << d_test_r << ", d_test_v = " << d_test_v << "\n";
////			std::cout << "v0th_r = " << v0th_r << ", v0th_v = " << v0th_v << "\n";
////			double d_r	= C1_r*exp(r1_r*t_drop) + C2_r*exp(r2_r*t_drop) + C3_r;
////			std::cout << "t = " << t_drop << ", d = " << d_r << "\n\n";
//			//------------------------------------------------------------------------------
////			}
////			std::cout << "t_drop: " << t_drop << ", d_test_r: " << d_test_r << ", d_test_v: " << d_test_v << "\n";
//			//------------------------------------------------------------------------------
//			//------------------------------------------------------------------------------
//			//------------------      THERMAL ZONE CALCULATIONS      -----------------------
//			//------------------------------------------------------------------------------
//			//------------------------------------------------------------------------------
////			//---------Now we know z and r on the cusp of the thermal zone------------------
////			//---------Reset parameters to include thermal components and re-run -----------
//
//			//--------CALCULATE THERMAL EFFECTS IN DROPLET SIZE REDUCTION-------------------
//			double dth_test_v;
//			double dth_test_r;
//			double q0 	  	= 373e-12;		// (373 um^2) 	 for water: 88e-12 (m^2)
//			double q1 	  	= 89.1;			// (8.91e-5 /um) for water: 4.3e3  (/m)
//			double del_T  	= 100000;
//			double dK 	  	= q0*del_T*(1+2*q1*radius);
////			std::cout << "dK = " << dK << "\n";
////			double K	  	= dK*t_drop;
//			double r_new  	= radius-radius*t_drop*exp(log(dK)/3);
////			std::cout << "dr = " << t_drop*exp(log(dK)/3) << "\n";
//			r=r_new;
//			double qd_new 	= 0.58*8*my::math::Pi*sqrt(gamma_d*permittivity*r_new*r_new*r_new);
//			double mass_new	= 4*my::math::Pi*rho_d*r_new*r_new*r_new/3;
////			std::cout << "dK = " << dK << ", K = " << K << "\n";
////			std::cout << "r_new= " << r_new << ", q_new= " << qd_new << ", mass= " << mass_new << "\n";
////
////			//---------CALCULATE THE ENERGIES ASSUMING LINEAR REDUCTION TO FINAL------------
////			//---------Electric force F_e is constant in this region and small--------------
////
//			double rth_star	= 	d_test_r/H;		//d_r/H;
//			double zth_star	= 	1-d_test_v/H;	//1-d_v/H;
////			std::cout << "rth_star = " << rth_star << ", zth_star = " << zth_star << "\n";
//
////			double v0th_v		= (A1_v+A2_v+A3_v)*(exp(r1_v*t_drop)-1)/r1_v+(B1_v+B2_v+B3_v)*(exp(r2_v*t_drop)-1)/r2_v+C_v*t_drop+v0_v;
////			double v0th_r		= (A1_r+A2_r+A3_r)*(exp(r1_r*t_drop)-1)/r1_r+(B1_r+B2_r+B3_r)*(exp(r2_r*t_drop)-1)/r2_r+C_r*t_drop+v0_r;
//			double d0th_v		= 0;//d_test_v;
//			double d0th_r		= 0;//d_test_r;
//			double tth_drop		= t_drop;
//
//			double thplusz	 	=	1+zth_star;//1+z_star;
//			double thminusz	 	=	1-zth_star;//1-z_star;
//			double throotplusz  =	std::sqrt(rth_star*rth_star+thplusz*thplusz);
//			double throotminusz =	std::sqrt(rth_star*rth_star+thminusz*thminusz);
//
////			double H	    = StartPosition[Parameter.open_boundary_direction]; //270 mm or 0.27 m
////			double Phi_0    = 10e3;					//V
////			double R	    = 1e-3;					//outer radius of the nozzle (guess)
////			double K_V	    = 1-exp(-0.021*H/R);	//non-dimensional related to H/R ratio
////			double Phi_star = K_V/(log(4*H/R));
////			double E_e		= Phi_0*Phi_star/H;
//
//			double Eth_v	= E_e*(1/throotminusz+1/throotplusz);
//			double Eth_r	= E_e*(thplusz/throotplusz-thminusz/throotminusz)/rth_star;
//			double theta_th	= atan(Eth_r/Eth_v);
////			std::cout << "Eth_r = " << Eth_r << ", Eth_v = " << Eth_v << ", theta_th = " << theta_th << "\n";
//
//			double Feth_v	= qd_new*Eth_v;
//			double Feth_r	= qd_new*Eth_r;
////			std::cout << "Feth_r = " << Feth_r << ", Feth_v = " << Feth_v << "\n";
//
//			double aeth_v	= Feth_v/mass_new;
//			double aeth_r 	= Feth_r/mass_new;
////			std::cout << "aeth_r = " << aeth_r << ", aeth_v = " << aeth_v << "\n";
//
//			//--------------Find all other forces acting on the droplet---------------------
//			// Find the droplet mass - knowns: radius, rho_d
////			double mass  = 4*my::math::Pi*rho_d*radius*radius*radius/3;
//
//			// Gravity force component acceleration
////			double g     = 9.81;	//m/s
////			std::cout << "mass = " << mass << ", Fg = " << g*mass << ", ";
//
//			// Stokes force component acceleration - knowns: rho_d, eta_a
//			double sth_f = (4.5*eta_a)/(rho_d*r_new*r_new);
////			std::cout << "Fs = " << s_f*mass << "\n";
//
//			// Thermophoretic force
////			double eta_a;
//			double kappa_a	= 0.025;
//			double kappa_d	= 0.19;
//			double grad_T	= 100000;
//			double T		= 523;
//			double rho_a	= 1.29;
//
//			double F_th	 = 3*my::math::Pi*eta_a*eta_a*r_new*3*kappa_a*grad_T/(rho_a*T*(2*kappa_a+kappa_d));
//			double a_th	 = F_th/mass_new;
//
//			// Electric force component acceleration - knowns: q_d, r_star, z_star;
////			double H	    = StartPosition[Parameter.open_boundary_direction]; //270 mm or 0.27 m
////			double Phi_0    = 10e3;					//V
////			double R	    = 1e-3;					//outer radius of the nozzle (guess)
////			double K_V	    = 1-exp(-0.021*H/R);	//non-dimensional related to H/R ratio
//
//			double ath_v = g+aeth_v-a_th;//+2*a_e;	// independent acceleration (gravity, initial e-force, thermal force)
//			double bth 	 = sth_f;					// velocity dependent acceleration		- Stokes force
//			double cth_v = aeth_v/t_heat;			// displacement dependent acceleration 	- vertical E-force
//			double ath_r = (theta_th<1e-20)?0:aeth_r;						  	// independent acceleration - Initial electric force
////			double b is known from vertical calculation
//			double cth_r = (theta_th<1e-20)?0:aeth_r/(t_heat*tan(theta_th));	// displacement dependent acceleration - Radial E-force
////			std::cout << "ath_v = " << ath_v << ", bth = " << bth << ", cth_v = " << cth_v << "\n";
//
//			double v_final_v=0;
//			double v_final_r=0;
//
//			if (bth*bth-4*cth_v<0) {
//				tth_drop=0;
//				dth_test_v	= t_heat;
//				dth_test_r	= d_test_r;
//				v_final_v	= 1;
//				v_final_r	= 0;
//			} else {
//				double r1th_v	= (-bth+std::sqrt(bth*bth-4*cth_v))/(2);
//				double r2th_v	= (-bth-std::sqrt(bth*bth-4*cth_v))/(2);
//
////				std::cout << "r1_v = " << r1_v << ", r2_v = " << r2_v << "\n";
//
//				double B1th_v 	= d0th_v*(r2th_v)/(r2th_v-r1th_v);
//				double A1th_v	= d0th_v-B1th_v;
//				double B2th_v 	= (v0th_v+d0th_v)/(r2th_v-r1th_v);
//				double A2th_v 	= -B2th_v;;
//				double Cth_v	= ath_v/(r1th_v*r2th_v);
//				double B3th_v	= -ath_v/(r2th_v*(r2th_v-r1th_v));
////				double B3_v	= a_v*(r1_v+r2_v-1)/(r1_v*r2_v*(1-r1_v));
//				double A3th_v	= -(B3th_v+Cth_v);
////				double A3_v	= -a_v/(r1_v*r2_v)-B3_v;
//
////				d_test_v = (A1_v+A2_v+A3_v)*exp(r1_v*t_check)+(B1_v+B2_v+B3_v)*exp(r2_v*t_check)+C_v;
////				std::cout << "t_check: " << t_check << ", d_test_v: " << d_test_v << "\n";
//
//				double iteration_stopth	=10000;
//				double tth_low			= 0;
//				double tth_high			= t_drop;
//				double tth_check		= (tth_low+tth_high)/2;
//				while (true){
//					dth_test_v = (A1th_v+A2th_v+A3th_v)*exp(r1th_v*tth_check)+(B1th_v+B2th_v+B3th_v)*exp(r2th_v*tth_check)+Cth_v;
//					if (tth_high-tth_low<1e-10) break;
//					if (dth_test_v > t_heat) {
//						tth_high=tth_check;
//						tth_check=(tth_high+tth_low)/2;
//					} else if (dth_test_v < t_heat) {
//						tth_low=tth_check;
//						tth_check=(tth_high+tth_low)/2;
//					}
//					if (iteration_stopth==0) break;
//					iteration_stopth--;
//				}
////				std::cout << "Iteration_th: " << iteration_stopth << ", tth_check: " << tth_check << ", dth_test_v: " << dth_test_v << "\n";
//
//				tth_drop  = tth_check;
//				v_final_v = (A1th_v+A2th_v+A3th_v)*(exp(r1th_v*tth_drop)-1)/r1th_v+(B1th_v+B2th_v+B3th_v)*(exp(r2th_v*tth_drop)-1)/r2th_v+Cth_v*tth_drop+v0th_v;
//
//
//				//---------Now time required to reach surface is known--------------------------
//
//				//---------Can now calculate the radial displacement----------------------------
////				double ath_r = aeth_r;						  	// independent acceleration - Initial electric force
////	//			double b is known from vertical calculation
////				double cth_r = aeth_r/(t_heat*tan(theta_th));	// displacement dependent acceleration - Radial E-force
////				std::cout << "ath_r = " << ath_r << ", bth = " << bth << ", cth_r = " << cth_r << "\n";
//
//				if (bth*bth-4*cth_r<0) {
//					tth_drop=0;
////					dth_test_v	= t_heat;
//					dth_test_r	= 0;
////					v_final_v	= 1;
//					v_final_r	= 0;
//				} else {
//					double r1th_r = (-bth+std::sqrt(bth*bth-4*cth_r))/(2);
//					double r2th_r = (-bth-std::sqrt(bth*bth-4*cth_r))/(2);
//
//		//			std::cout << "r1_r = " << r1_r << ", r2_r = " << r2_r << "\n";
//
//					double B1th_r 	= d0th_r*(r2th_r)/(r2th_r-r1th_r);
//					double A1th_r 	= d0th_r-B1th_r;
//					double B2th_r 	= (v0th_r+d0th_r)/(r2th_r-r1th_r);
//					double A2th_r 	= -B2th_r;;
//					double Cth_r  	= ath_r/(r1th_r*r2th_r);
//					double B3th_r 	= -ath_r/(r2th_r*(r2th_r-r1th_r));
//	//				double B3_v	  	= a_v*(r1_v+r2_v-1)/(r1_v*r2_v*(1-r1_v));
//					double A3th_r 	= -(B3th_r+Cth_r);
//	//				double A3_v	  	= -a_v/(r1_v*r2_v)-B3_v;
//
//					dth_test_r = (A1th_r+A2th_r+A3th_r)*exp(r1th_r*tth_drop)+(B1th_r+B2th_r+B3th_r)*exp(r2th_r*tth_drop)+Cth_r;
//					v_final_r  = (A1th_r+A2th_r+A3th_r)*(exp(r1th_r*tth_drop)-1)/r1th_r+(B1th_r+B2th_r+B3th_r)*(exp(r2th_r*tth_drop)-1)/r2th_r+Cth_r*tth_drop+v0th_r;
//				}
//			}
//
//			Velocity[0] = v_final_r*cos(phi);
//			Velocity[1] = -v_final_v;
//			Velocity[2] = v_final_r*sin(phi);
//
//			Position[0] = std::sqrt(dth_test_r+d_test_r)*cos(phi)+StartPosition[0];
//			Position[1] = 0;
//			Position[2] = std::sqrt(dth_test_r+d_test_r)*sin(phi)+StartPosition[2];
//
//			r	= r_new;
//			q	= qd_new;
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
//			r=1/pow(volume_fraction*(r_max_inv_third-r_min_inv_third)+r_min_inv_third,3.);
			double radius=1/pow(volume_fraction*(r_max_inv_third-r_min_inv_third)+r_min_inv_third,3.);

			r	= radius;

			Position[0]=RandomNumber()*4*(Partition.Max(0)*Parameter.GridDelta)-2*(Partition.Max(0)*Parameter.GridDelta);
			Position[1]=0;
			Position[2]=RandomNumber()*4*(Partition.Max(2)*Parameter.GridDelta)-2*(Partition.Max(2)*Parameter.GridDelta);

//			Position[0]=RandomNumber()*2*(Partition.Max(0)*Parameter.GridDelta-Partition.Min(0)*Parameter.GridDelta+4*radius)+Partition.Min(0)*Parameter.GridDelta-2*radius;
//			Position[1]=0;
//			Position[2]=RandomNumber()*(Partition.Max(2)*Parameter.GridDelta-Partition.Min(2)*Parameter.GridDelta+4*radius)+Partition.Min(2)*Parameter.GridDelta-2*radius;

			double gamma_d = 0.022;
			double permittivity = 8.854187817e-12;
			double qd=0.58*8*my::math::Pi*sqrt(gamma_d*permittivity*radius*radius*radius);
			q	= qd;

		}

		template<class VecType1, class VecType2, class ParameterType, class PartitionType>
		inline void EvenlyDistributed(const VecType1& StartPosition, VecType2& Position, const ParameterType& Parameter, const PartitionType& Partition){

//			while (true) {
				Position[0]=RandomNumber()*(Partition.Max(0)*Parameter.GridDelta-Partition.Min(0)*Parameter.GridDelta)+Partition.Min(0)*Parameter.GridDelta;
				Position[1]=0;
				Position[2]=RandomNumber()*(Partition.Max(2)*Parameter.GridDelta-Partition.Min(2)*Parameter.GridDelta)+Partition.Min(2)*Parameter.GridDelta;

//	            double radius=0;
//	            for (int i=0;i<3;++i) radius+=(Parameter.open_boundary_direction==i)?0:(std::abs(Position[i]-50))*(std::abs(Position[i]-50));
//	            radius=std::sqrt(radius);
//	            if (radius <= 50) break;
//			}
                }


                template<class DropletType, class VecType1>
        inline void DiskDistribution(const DropletType d, VecType1& position) {

//			std::cout << "d.Position: " << d.Position[0] << ", " << d.Position[1] << ", " << d.Position[2] << "\n";
//			std::cout << "d.radius:   " << d.Radius << "\n";
//			std::cout << "d.charge:   " << d.Charge << "\n";

//			std::cout << d.Position[0] << ", " << d.Position[2] << "\n";

//			double randomNum	= RandomNumber();
//			double radius		= std::sqrt(randomNum)*(2*d.Radius);
////			std::cout << "radius: " << radius << "\n";
//			double theta		= 2*my::math::Pi*RandomNumber();
////			std::cout << radius << ", " << theta << "\n";

//			position[0] = radius*cos(theta)+d.Position[0];
//			position[1] = 0;//d.Position[1];
//			position[2] = radius*sin(theta)+d.Position[2];

                        double v0, v1, rsq;
            do {
                v0=2.0*RandomNumber()-1.0;
                v1=2.0*RandomNumber()-1.0;
                rsq=v0*v0+v1*v1;
            } while (rsq>=1.0 || rsq<1e-20);

			position[0] = 2*d.Radius*v0+d.Position[0];
			position[1] = 0;//d.Position[1];
			position[2] = 2*d.Radius*v1+d.Position[2];


//			ofstream out;
//			out.open("distribution.txt");
//
//			out << "x, y\n";
//
//			for (int i=0;i<1000;i++) {
//				randomNum=RandomNumber();
//				radius=std::sqrt(randomNum);
//				theta=2*my::math::Pi*RandomNumber();
//				out << radius*cos(theta) << ", " << radius*sin(theta) << "\n";
//			}
//
//			out.close();




                }

//        inline void getRandomN(double& rando) {
//
//        	rando=RandomNumber();
//
//		}


/*
                template<class VecType1, class VecType2, class VecType3, class ValType, class ParameterType, class PartitionType>
        inline void SprayPyrolysisAmbient(const VecType1& StartPosition, const ValType& min_radius, const ValType& max_radius, const int& flow_rate, const int& SprayDirection, const double& direction_phi, const double& direction_theta, const int& AdditionalDirection, const VecType2 direction, VecType3& position, double& FluxFactor, const ParameterType& Parameter, const PartitionType& Partition) {

			double theta_min = atan(min((Partition.Min(AdditionalDirection)*Parameter.GridDelta-StartPosition[AdditionalDirection])/(Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection]), (Partition.Min(AdditionalDirection)*Parameter.GridDelta-StartPosition[AdditionalDirection])/(Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])));
			double theta_max = atan(max((Partition.Max(AdditionalDirection)*Parameter.GridDelta-StartPosition[AdditionalDirection])/(Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection]), (Partition.Max(AdditionalDirection)*Parameter.GridDelta-StartPosition[AdditionalDirection])/(Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])));
//        	std::cout << "StartPosition: (" << StartPosition[0] << ", " << StartPosition[1] << ", " << StartPosition[2] << ")\n";

//			std::cout << "theta_min: " << theta_min << std::endl;
//			std::cout << "theta_max: " << theta_max << std::endl;
//			double phi_min = atan((Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])/((Partition.Max(Parameter.open_boundary_direction)*Parameter.GridDelta-StartPosition[Parameter.open_boundary_direction])));
//			double phi_max = atan((Partition.Max()))
		//if flow rate is 30ml/h:			//if flow rate is 120ml/h:
		//0-5um  - 0.57						//0-5um  - 0.55
		//5-10um - 0.32						//5-10um - 0.32
		//10-15  - 0.05						//10-15  - 0.08
		//15-20  - 0.025					//15-20  - 0.03
		//20-25  - 0.02						//20-25  - 0.011
		//25-30  - 0.01						//25-30  - 0.006
		//30-35  - 0.003					//30-35  - 0.002
		//35-40  - 0.001					//35-40  - 0.001
		//40-45  - 0.001
			double StartPosition_temp[3];
		long double v0[3];
		double factor=1;//(max_radius-min_radius)/45e-6;
		double s;
//	        long double t0;
//	        long double t1;
		long double t;
//	        double w;

		double a_gr = -9.81;
		double air_viscosity = 2.2e-5;
		double density = 780;
		double air_density=1.29;
		double air_thermal_cond=0.025;
		double droplet_thermal_cond=0.19;
		double gradient=25000; //was 100000
		double temperature=250+273;//was 400+273 //523;

                long double phi;
                long double velocity;
                double rad;

		//std::cout << "mid_radius: " << min_radius << std::endl;
		//while (true) {
		while(true){

			for (int i=0;i<3;++i) StartPosition_temp[i]=StartPosition[i];

//        		rad = (max_radius-min_radius)*RandomNumber()+min_radius;//0.1e-6;//(5e-6)*RandomNumber();
			//rad = 45e-6;//min_radius+1e-10;
			//double rad;
			while (true) {
				break;
				rad = (5e-6)*RandomNumber();
				double flow = RandomNumber();

					if (flow_rate==30) {
						flow*=11;
						if (flow<5) {
							rad+=10e-6;
						} else if ((flow>5) && (flow<7.5)) {
							rad+=15e-6;
						} else if ((flow>7.5) && (flow<9.5)) {
							rad+=20e-6;
						} else if ((flow>9.5) && (flow<10.5)) {
							rad+=25e-6;
						} else if ((flow>10.5) && (flow<10.8)) {
							rad+=30e-6;
						} else if ((flow>10.8) && (flow<10.9)) {
							rad+=35e-6;
						} else {
							rad+=40e-6;
						}
					} else {
						flow*=13;
						if (flow<8) {
							rad+=10e-6;
						} else if ((flow>8) && (flow<11)) {
							rad+=15e-6;
						} else if ((flow>11) && (flow<12.1)) {
							rad+=20e-6;
						} else if ((flow>12.1) && (flow<12.7)) {
							rad+=25e-6;
						} else if ((flow>12.7) && (flow<12.9)) {
							rad+=30e-6;
						} else {
							rad+=35e-6;
						}
				}
				if ((rad>min_radius) && (rad<max_radius)) break;
			}
//	        	rad = (8*RandomNumber()+2)*1e-6;
//	        	std::cout << "rad, " << rad;
			//double rad = (45e-6-min_radius)*RandomNumber()+min_radius;
			//double

//	        	phi = my::math::Pi*(70*RandomNumber()-35)/180+direction_phi;
			phi = my::math::Pi*(10*RandomNumber()-5)/180+direction_phi; //narrow nozzle
//	        	std::cout << "phi: " << phi*180/my::math::Pi << std::endl;

			rad=(15*RandomNumber()+30)*1e-6;
//	        	std::cout << "rad: " << rad << std::endl;;

			s=4.5*air_viscosity/(density*rad*rad);
//	        	std::cout << "s: " << s << std::endl;

			long double min_velocity = 0.99999*(Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s/cos(phi);
			long double max_velocity = 1.00001*(Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s/cos(phi);
//	        	std::cout << "(min_vel, max_vel) = (" << min_velocity << ", " << max_velocity << ")\n";

			factor*=(max_velocity-min_velocity)/(15.);

			velocity = (max_velocity-min_velocity)*RandomNumber()+min_velocity;
//	        	velocity = 15.*RandomNumber()+35.;
//	        	std::cout << "velocity: " << velocity << std::endl;;
//	        	velocity=5;
//	        	std::cout << ", phi,  " << phi << std::endl;
//	        	double phi = my::math::Pi*(70*RandomNumber()-35)/180;
			long double theta = (theta_max-theta_min)*RandomNumber()+theta_min+direction_theta;
			factor*=(theta_max-theta_min)/(70.*my::math::Pi/180.);

//	        	velocity=20.;
//	        	phi=atan(150./260.);
//	        	theta=0.;
			v0[0] = velocity*cos(phi)*cos(theta);
			v0[1] = velocity*cos(phi)*sin(theta);
			v0[2] = velocity*sin(phi);
//	        	v0[2] = -20.;
//	        	std::cout << "v0: (" << v0[0] << ", " << v0[1] << ", " << v0[2] << ")\n";


			// Minimum velocity with phi known
//	        	long double min_velocity = (Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s*(1/(exp((a_gr+StartPosition[Parameter.open_boundary_direction]*s*s)/(a_gr*a_gr))-1)+1);
//	        	long double max_velocity = (Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s*(1/(exp((a_gr+StartPosition[Parameter.open_boundary_direction]*s*s)/(a_gr*a_gr))-1)+1);
			//std::cout << "phi: " << phi*180/my::math::Pi << std::endl;
			//std::cout << "min_velocity: " << min_velocity << std::endl;
			//std::cout << "max_velocity: " << max_velocity << std::endl;
			//std::cout << "dif_velocity: " << max_velocity-min_velocity << std::endl;

//	        	// Maximum velocity when phi = 35
//
//	        	long double min_velocity3 = (Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s/cos(0);
//	        	long double max_velocity3 = (Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s/cos(0);
//	        	std::cout << "min_velocity3: " << min_velocity3 << std::endl;
//	        	std::cout << "max_velocity3: " << max_velocity3 << std::endl;
//	        	std::cout << "dif_velocity3: " << max_velocity3-min_velocity3 << std::endl;
//
//	        	long double min_velocity2 = (Partition.Min(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s/cos(phi*my::math::Pi/180);
//	        	long double max_velocity2 = (Partition.Max(SprayDirection)*Parameter.GridDelta-StartPosition[SprayDirection])*s/cos(phi*my::math::Pi/180);
////	        	std::cout << "s: " << s << std::endl;
////	        	std::cout << "a_gr: " << a_gr << std::endl;
////	        	std::cout << "W: " << -exp((0.5*s*s/a_gr)+(0.15*s*tan(35*my::math::Pi/180)/a_gr)-1) << std::endl;
////	        	std::cout << "LambertW: " << LambertW<0>(-exp((0.5*s*s/a_gr)+(0.15*s*tan(35*my::math::Pi/180)/a_gr)-1)) << std::endl;
//	        	std::cout << "min_velocity2: " << min_velocity2 << std::endl;
//	        	std::cout << "max_velocity2: " << max_velocity2 << std::endl;
//	        	std::cout << "dif_velocity2: " << max_velocity2-min_velocity2 << std::endl;

			double d = StartPosition_temp[Parameter.open_boundary_direction]-((Partition.Max(Parameter.open_boundary_direction))*Parameter.GridDelta+6e-3);
			d = -d;
//	        	std::cout << "d: " << d << std::endl;


			double C1=(v0[2]-a_gr/s)/s;
			double a = C1;// -C1;
//	        	std::cout << "g: " << a_gr << std::endl;
//	        	std::cout << "a: " << a << std::endl;
			double b = a_gr/s;
//	        	std::cout << "b: " << b << std::endl;
			double c = C1;
			double w = (v0[2]*s/a_gr-1)*exp((a-d)*s*s/a_gr); //a*s/b*exp(s*(a-d)/b);
//	        	std::cout << "w: " << w << std::endl;
			if (w>=std::numeric_limits<double>::max()) {
				t=0;
				} else {
					//if (w<-1/exp(1)){
					t=-(a-d)/b; // t = -(A-d)/B
					if (w>=-1/exp(1)) {
						long double lamb0 = LambertW<0>(w)/s;
//						std::cout << "LambertW<0>(w)/s: " << lamb0 << std::endl;

						//t=(C1-d)/b;
						t += lamb0;
//						double distance0;
//						distance0 = (v0[2]-a_gr/s)/s*(1-exp(-s*t))+a_gr*t/s;
//			        	std::cout << "t0: " << t << std::endl;
//			        	std::cout << "distance0: " << distance0 << std::endl;

//						long double lamb1 = LambertW<-1>(w)/s;
//						std::cout << "LambertW<-1>(w)/s: " << lamb1 << std::endl;
//
//						t += lamb1;
//						double distance1;
//						distance1 = (v0[2]-a_gr/s)/s*(1-exp(-s*t))+a_gr*t/s;
//			        	std::cout << "t1: " << t << std::endl;
//			        	std::cout << "distance1: " << distance1 << std::endl;
					} else {
						t = -t;
					}
			}
			//std::cout << "t1: " << d*s/a_gr+1/s-v0[2]/a_gr << std::endl;
// 	        	std::cout << "t: " << t << std::endl;


			double x, y;
			x=1/s*v0[0]*(1-exp(-s*t));
			y=1/s*v0[1]*(1-exp(-s*t));

//	        	std::cout << "(x,y) = (" << x << ", " << y << ")\n";

			StartPosition_temp[SprayDirection]+=x;
			StartPosition_temp[AdditionalDirection]+=y;
//	        	std::cout << "StartPosition_temp: " << StartPosition_temp[Parameter.open_boundary_direction] << std::endl;
			StartPosition_temp[Parameter.open_boundary_direction]+=d;
//	        	std::cout << "StartPosition_temp: " << StartPosition_temp[Parameter.open_boundary_direction] << std::endl;
//	        	std::cout << "StartPosition_temp: (" << StartPosition_temp[0] << ", " << StartPosition_temp[1] << ", " << StartPosition_temp[2] << ")\n";

			//Updating the velocity as particle enters heat zone:
			v0[0]*=exp(-s*t);
			v0[1]*=exp(-s*t);
			v0[2]=(v0[2]-a_gr/s)*exp(-s*t) + a_gr/s;
//	        	std::cout << "v2: (" << v0[0] << ", " << v0[1] << ", " << v0[2] << ")\n";

//	        	std::cout << "Vy: " << v0[2] << std::endl;

			//Now accessing the heat zone:

			rad-=10e-6;

//	        	s=4.5*air_viscosity/(density*rad*rad);

			double F_th = (3*my::math::Pi*air_viscosity*air_viscosity*rad/air_density)*(3*air_thermal_cond/(2*air_thermal_cond+droplet_thermal_cond))*gradient/temperature;
			double a_th = F_th*3/(4*my::math::Pi*density*rad*rad*rad);
//	        	std::cout << "a_th: " << a_th << ", a_gr: " << a_gr << std::endl;
//
//	        	std::cout << "Partition.Max: " << Partition.Max(Parameter.open_boundary_direction) << std::endl;
//	        	std::cout << "Partition.Max*Delta: " << (Partition.Max(Parameter.open_boundary_direction))*Parameter.GridDelta << std::endl;
			d = StartPosition_temp[Parameter.open_boundary_direction]-(Partition.Max(Parameter.open_boundary_direction))*Parameter.GridDelta;
			d = -d;
//	        	std::cout << "d2: " << d << std::endl;

                        double C2=(v0[2]-(a_gr+a_th)/s)/s;
                        a = C2;
//	        	std::cout << "a: " << a << std::endl;
                        b = (a_gr+a_th)/s;

//        		std::cout << "b: " << b << std::endl;
			c = C2;

			//double
			w = (v0[2]*s/(a_gr+a_th)-1)*exp((a-d)*s*s/(a_gr+a_th)); //a*s/b*exp(s*(a-d)/b);
//	        	std::cout << "w: " << w << std::endl;

//        		w = -a*s*exp(s*(c-d)/b)/b;
			if (w>=std::numeric_limits<double>::max()) {
				t=0;
				} else {
					//if (w<-1/exp(1)){
					t=-(a-d)/b; // t = -(A-d)/B
					if (w>=-1/exp(1)) {
						long double lamb0 = LambertW<0>(w)/s;
//						std::cout << "LambertW<0>(w)/s: " << lamb0 << std::endl;

						//t=(C1-d)/b;
						t += lamb0;
//						double distance0;
//						distance0 = (v0[2]-a_gr/s)/s*(1-exp(-s*t))+a_gr*t/s;
//			        	std::cout << "t0: " << t << std::endl;
//			        	std::cout << "distance0: " << distance0 << std::endl;

//						long double lamb1 = LambertW<-1>(w)/s;
//						std::cout << "LambertW<-1>(w)/s: " << lamb1 << std::endl;
////
//						long double t1 = lamb1 -(a-d)/b;
//						double distance1;
//						distance1 = (v0[2]-a_gr/s)/s*(1-exp(-s*t1))+a_gr*t1/s;
//			        	std::cout << "t1: " << t1 << std::endl;
//			        	std::cout << "distance1: " << distance1 << std::endl;
					} else {
						t = -t;
					}
			}
			//std::cout << "t1: " << d*s/a_gr+1/s-v0[2]/a_gr << std::endl;
// 	        	std::cout << "t2: " << t << std::endl;
			//std::cout << "t2: " << d*s/(a_gr-a_th)+1/s-v0[2]/(a_gr-a_th) << std::endl;


				//t is the time to reach the feature scale region
			x=1/s*v0[0]*(1-exp(-s*t));
			y=1/s*v0[1]*(1-exp(-s*t));
//	        	std::cout << "(x2,y2) = (" << x << ", " << y << ")\n";

			StartPosition_temp[SprayDirection]+=x;
			StartPosition_temp[AdditionalDirection]+=y;
			StartPosition_temp[Parameter.open_boundary_direction]+=d;
//	        	std::cout << "StartPosition_temp: (" << StartPosition_temp[0] << ", " << StartPosition_temp[1] << ", " << StartPosition_temp[2] << ")\n";


//	    		if ((StartPosition_temp[SprayDirection]>Partition.Min(SprayDirection)*Parameter.GridDelta) && (StartPosition_temp[SprayDirection]<Partition.Max(SprayDirection)*Parameter.GridDelta)) {
//	        		if ((StartPosition_temp[AdditionalDirection]>Partition.Min(AdditionalDirection)*Parameter.GridDelta) && (StartPosition_temp[AdditionalDirection]<Partition.Max(AdditionalDirection)*Parameter.GridDelta)) {
//					   break;
//	        		}
//	        	}
//		        std::cout << "Position (" << StartPosition_temp[0] << ", " << StartPosition_temp[1] << ", " << StartPosition_temp[2] << ")" << std::endl;

		//	if (((StartPosition_temp[SprayDirection]-Partition.Min(SprayDirection)*Parameter.GridDelta) > 1e-10) && ((StartPosition_temp[SprayDirection]-Partition.Max(SprayDirection)*Parameter.GridDelta) < -1e-10)) {
		//		break;
		//	}
		//}
//	        std::cout << "Right Pos (" << StartPosition_temp[0] << ", " << StartPosition_temp[1] << ", " << StartPosition_temp[2] << ")" << std::endl;

//	        std::cout << ", " << rad;
//	        std::cout << ", " << velocity;
//	        std::cout << ", " << StartPosition_temp[SprayDirection] << std::endl;
//	        std::cout << ", " << phi << std::endl;// *180./my::math::Pi << std::endl;


                for (int i=0;i<3;++i) position[i]=StartPosition_temp[i];
//	        std::cout << "position (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;

                if ((position[0] < (Partition.Min(0))*Parameter.GridDelta) || (position[0] > (Partition.Max(0))*Parameter.GridDelta) || (position[2] < (Partition.Min(2))*Parameter.GridDelta) || (position[2] > (Partition.Max(2))*Parameter.GridDelta)) {
//        		std::cout << "position (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
//        		std::cout << "1\n";
                } else {
//        		std::cout << "pos (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
//        		std::cout << "2\n";
                        break;
                }
                }
                if ((position[0] < (Partition.Min(0))*Parameter.GridDelta) || (position[0] > (Partition.Max(0))*Parameter.GridDelta) || (position[2] < (Partition.Min(2))*Parameter.GridDelta) || (position[2] > (Partition.Max(2))*Parameter.GridDelta)) {
                        //        		std::cout << "position (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
                        //        		std::cout << "1\n";
                } else {
                        //        		std::cout << "pos (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
                        //        		std::cout << "2\n";
                }

                direction[SprayDirection] = v0[0]*exp(-s*t);
                direction[AdditionalDirection] = v0[1]*exp(-s*t);
                direction[Parameter.open_boundary_direction] = (v0[2]-a_gr/s)*exp(-s*t) + a_gr/s;

                double v_end=0;
                for (int i=0;i<3;++i) v_end+=direction[i]*direction[i];
                for (int i=0;i<3;++i) direction[i]/=sqrt(v_end);
//	        std::cout << "direction (" << direction[0] << ", " << direction[1] << ", " << direction[2] << ")" << std::endl;
                FluxFactor=factor;

                }
*/
        }
}

#endif //DEF_STATISTICS
