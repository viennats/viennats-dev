#ifndef MODELFIB_H_
#define MODELFIB_H_

#include <stack>
#include "../Statistics.h"
#include <cassert>
#include <string>
#include <fstream>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include "../parser_actors.h"

namespace model {
///Focused ion beam model

	class FIB {

       static const double e0;

        double StartDirection[3];
        double StartPosition[3];
        double FWHM;
        std::string SputterYieldFile;
        double current;
        double target_density;

        int num_reemitted_particles;
        int num_simulated_particles;

        double end_probability;
        double sticking_probability;

        bool calculate_redeposition;
        bool use_yamamura;

        double yield_max;
        double yield_0;
        double angle_max;

        double c0;
        double c1;
        double c2;

        std::vector<double> SputterYieldX, SputterYieldY;

        static double YamamuraFunction(double c1, double c2, double cos) {
        	return std::pow(cos,-c1)*exp(c2*(1.-1./cos));
        }

    public:

        class ParticleType {
        public:
            double Direction[3];
            double Flux;
            double probability;
            bool direct;
        };

        unsigned int NumberOfParticleClusters[1];

        static const bool OutputFluxes=false;
        static const int RatesStorageSize=1;
        static const int CoverageStorageSize=0;
        static const bool ReemissionIsMaterialDependent=false;
        static const bool CalculateNormalVectors=true;
        static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;
        static const unsigned int NumberOfParticleTypes=1;
        static const bool SpatiallyEqualDistributedFlux=false;

        FIB(const std::string & Parameters) {

            using namespace boost::spirit::classic;
            using namespace parser_actors;

            calculate_redeposition=false;
            sticking_probability=0;
            num_reemitted_particles=1;
            end_probability=0;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("FWHM")  >> '='  >> real_p[assign_a(FWHM)]  >> ';') |
                            (str_p("yield_0")  >> '='  >> real_p[assign_a(yield_0)]  >> ';') |
                            (str_p("yield_max")  >> '='  >> real_p[assign_a(yield_max)]  >> ';') |
                            (str_p("angle_max")  >> '='  >> real_p[assign_a(angle_max)]  >> ';') |
                            (str_p("current")  >> '='  >> real_p[assign_a(current)]  >> ';') |
                            (str_p("num_reemitted_particles")  >> '='  >> int_p[assign_a(num_reemitted_particles)]  >> ';') |
                            (str_p("target_density")  >> '='  >> real_p[assign_a(target_density)]  >> ';') |
							(str_p("stop_criterion")  >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("sticking_probability")  >> '='  >> real_p[assign_a(sticking_probability)]  >> ';') |
                            (str_p("position")  >> '='  >> '{' >> real_p[assign_a(StartPosition[0])]  >> "," >> real_p[assign_a(StartPosition[1])] >> "," >> real_p[assign_a(StartPosition[2])] >> '}' >> ';') |
                            (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("sputter_yield_file") >> '='  >>  '\"' >> *((~ch_p('\"'))[push_back_a(SputterYieldFile)]) >> '\"' >> ';') |
                            (str_p("num_simulated_particles")  >> '='  >>int_p[assign_a(num_simulated_particles)]  >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

            double norm=std::sqrt(StartDirection[0]*StartDirection[0]+StartDirection[1]*StartDirection[1]+StartDirection[2]*StartDirection[2]);
            assert(norm>0.);
            StartDirection[0]/=norm;
            StartDirection[1]/=norm;
            StartDirection[2]/=norm;

            std::cout << StartDirection[0] <<  " "<< StartDirection[1] <<  " "<< StartDirection[2] <<  std::endl;

            calculate_redeposition=(sticking_probability>0.);
            use_yamamura=SputterYieldFile.empty();

            if (use_yamamura) {
            	double cos_max=std::cos(angle_max/180.*my::math::Pi);
            	c0=yield_0;
            	if (cos_max<=0) {
            		c1=0.;
            	} else {
            		c1=std::log(yield_max/yield_0)/(cos_max-1-std::log(cos_max));
            	}
            	c2=c1*cos_max;

            	std::cout << "c0=" << c0 << " c1=" << c1 << " c2=" << c2 << std::endl;

			} else {
				std::ifstream f(SputterYieldFile.c_str());
				assert(f);
				while (!f.eof()) {
					double x,y;
					f >> x >> y;
					SputterYieldX.push_back(std::cos(x/180.*my::math::Pi));
					SputterYieldY.push_back(y);
				}
				f.close();
			}

            NumberOfParticleClusters[0]=num_simulated_particles;
        }

        template <class NormVecType>
        void CalculateVelocity(double &Velocity, const NormVecType NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
            Velocity=std::min(0.,Rates[0]/(e0*target_density));
            //Velocity=-Rates[0]/(e0*target_density);
            //std::cout << "velocity: " << Velocity << std::endl;
        }

		template<class VecType>
		void CalculateVectorVelocity(
				double *Velocity,
				const VecType& NormalVector,
				const double *Coverages,
				const double *Rates,
				int Material,
				bool connected,
				bool visible) const {}


        static void UpdateCoverage(double *Coverages, const double *Rates) {
        }
        /*
        template <class ParameterType, class PartitionType>
        void DefineAngles(const ParameterType& Parameter, const PartitionType& Partition, double &theta_min, double &theta_max) const {
        }

        template <class PT, class ParameterType, class PartitionType>
        void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position, const ParameterType& Parameter, const PartitionType& Partition, double theta_min, double theta_max) const {
        }
*/
		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {
/*
            p.Direction[0]=0;
            p.Direction[1]=-1;
            p.Direction[2]=0;

            Position[0]=400e-9;
            Position[1]=1000e-9;
            Position[2]=400e-9;

            p.direct=true;
            p.probability=1.;
            p.Flux=current;
*/
           for (int i=0;i<3;++i) p.Direction[i]=StartDirection[i];


           //PickRandomPointOnUnitCircle(cosphi, sinphi);
             //      	r2=1;

           my::stat::NormalDistributedStartPosition(StartPosition,StartDirection, FWHM, Position);

          // double t = (StartPosition[1]-Position[1])/p.Direction[1];
           //Position[0]+=t*p.Direction[0];
           //Position[2]+=t*p.Direction[2];



           //Position[0]=my::stat::RandomNumber()*100;
           //Position[2]=my::stat::RandomNumber()*10;

           //my::stat::CosineNDistributedRandomDirection(1000.,StartDirection,p.Direction);

           p.direct=true;
           p.probability=1.;
           p.Flux=current;

        }



        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {

            if (p.direct) {
            	double dot=NormalVector[0]*p.Direction[0]+NormalVector[1]*p.Direction[1]+NormalVector[2]*p.Direction[2];
            	double cos=std::min(std::max(-dot,0.),1.);

            	if (use_yamamura) {
//            		Rates[0]-=p.Flux*c0*YamamuraFunction(c1,c2,cos);
            		Rates[0]+=p.Flux*c0*YamamuraFunction(c1,c2,cos);
            	} else {
//            		Rates[0]-=p.Flux*my::math::Interpolate(&SputterYieldX[0],&SputterYieldY[0],SputterYieldX.size(),cos);
            		Rates[0]+=p.Flux*my::math::Interpolate(&SputterYieldX[0],&SputterYieldY[0],SputterYieldX.size(),cos);
            	}
            } else {
                Rates[0]+=sticking_probability*p.probability*p.Flux;
            }
        }


		template <class PT, class VecType> void ParticleReflexion(
							const PT& p,
							std::stack<PT>& particle_stack,
							const VecType& NormalVector,
							const double* Coverages,
							int Material//,
//                            int D,
//                            double dot // dot product between the incoming particle direction and the normal vector
        					) const {

            if (calculate_redeposition) {
                if (p.direct) {
                    double dot=NormalVector[0]*p.Direction[0]+NormalVector[1]*p.Direction[1]+NormalVector[2]*p.Direction[2];
                    double cos=std::min(std::max(-dot,0.),1.);
                    double new_probability;
                    if (use_yamamura) {
                    	new_probability=c0*YamamuraFunction(c1,c2,cos);
                    } else {
                    	new_probability=my::math::Interpolate(&SputterYieldX[0],&SputterYieldY[0],SputterYieldX.size(),cos);
                    }
                    if (new_probability>end_probability) {
						for (int k=0;k<num_reemitted_particles;++k) {
							particle_stack.push(p);
							PT& p_new=particle_stack.top();
							my::stat::CosineNDistributedRandomDirection(1,NormalVector,p_new.Direction);
							p_new.Flux=p.Flux/num_reemitted_particles;
							p_new.probability=new_probability;
							p_new.direct=false;
						}
                    }
                } else {
                    double new_probability=p.probability*(1.-sticking_probability);
                    if (new_probability>end_probability) {
                        particle_stack.push(p);
                        PT& p_new=particle_stack.top();
                        my::stat::CosineNDistributedRandomDirection(1,NormalVector,p_new.Direction);
                        p_new.Flux=p.Flux;
                        p_new.probability=new_probability;
                        p_new.direct=false;
                    }
                }
            }
        }




    };

   double const FIB::e0=-1.602176487e-19;

    //const double FIB::StartDirection[3]={0.,-1.,0.};    //AXIS
}






#endif /*MODELFIB_H_*/
