#ifndef MODELTIN_ALD_H_
#define MODELTIN_ALD_H_

#include <stack>
#include "../Statistics.h"
#include <cassert>
#include <string>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

namespace model {

///TiN ALD model

	class TiN_ALD {

		static const double end_probability;
		static const double reaction_order;

		double StartDirection[3];

		double sticking_probability;
		double Flux_TDMAT;
		double Flux_NH3;
		double TDMAT_time;
		double NH3_time;
		double Purge_time;
		double Cycle_time;
		double Flux;

	public:

		class ParticleType {
		public:
			double Direction[3];
			double Probability;
			double Energy;
			double Flux;
			int Type;
		};

		unsigned int NumberOfParticleClusters[3];

		static const int CoverageStorageSize=2;
		static const int RatesStorageSize=6;

		static const unsigned int NumberOfParticleTypes=2;

        static const bool OutputFluxes=false;
		static const bool SpatiallyEqualDistributedFlux=true;

		static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;

        static const bool CalculateNormalVectors=true;

        static const bool ReemissionIsMaterialDependent=true;

        TiN_ALD(const std::string & Parameters) {

		    using namespace boost::spirit::classic;

		    double Accuracy;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                    		(str_p("direction")				>> '='  >> '{' >> real_p[assign_a(StartDirection[0])] >> ","
                    												>> real_p[assign_a(StartDirection[1])] >> ","
                    												>> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("TDMAT_time")			>> '='  >> real_p[assign_a(TDMAT_time)]  >> ';') |
                            (str_p("NH3_time")				>> '='  >> real_p[assign_a(NH3_time)]  >> ';') |
                            (str_p("Purge_time")			>> '='  >> real_p[assign_a(Purge_time)]  >> ';') |
                            (str_p("Flux_TDMAT")  			>> '='  >> real_p[assign_a(Flux_TDMAT)]  >> ';') |
                            (str_p("Flux_NH3")    			>> '='  >> real_p[assign_a(Flux_NH3)]  >> ';') |
                            (str_p("sticking_probability")  >> '='  >> real_p[assign_a(sticking_probability)]  >> ';') |
                            (str_p("statistical_accuracy")	>> '='  >>real_p[assign_a(Accuracy)]  >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

		    unsigned int num_particles=static_cast<unsigned int>(Accuracy);

		    NumberOfParticleClusters[0]=(Flux_TDMAT>0.)?num_particles:0;	//Flux
		    NumberOfParticleClusters[1]=(Flux_NH3>0.)?num_particles:0;		//Flux
		    Cycle_time = TDMAT_time+NH3_time+2.*Purge_time;
		}

		template <class VecType>
		void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
			if (Material==0) {	//Si
				Velocity=0.*Coverages[1];
			} else if (Material==1) {
				Velocity=0.*Coverages[1];										//in cm/s     //TODO
			} else Velocity=0.;
		}

		template<class VecType>
		void CalculateVectorVelocity(
				double *Velocity,
				const VecType& NormalVector,
				const double *Coverages,
				const double *Rates,
				int Material,
				bool connected,
				bool visible
				) const {}

		static void UpdateCoverage(double *Coverages, const double *Rates) {

			Coverages[0]=1;
			Coverages[1]=1;
		}

		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {

			double RelativeTimeTemp=ProcessTime / Cycle_time;
			int RTI=(int) RelativeTimeTemp-0.5;
			double RelativeTime=RelativeTimeTemp-(double) RTI; 
			if (RelativeTime <= TDMAT_time) {
				p.Type = 0;
			// Generate TDMAT particles
				my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
				p.Probability=1.;
				p.Flux=Flux;
			} else if ((RelativeTime > (TDMAT_time+Purge_time)) && (RelativeTime <= (TDMAT_time+Purge_time+NH3_time))) {
			// Generate NH3 particles
				p.Type = 1;
				my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
				p.Probability=1.;
				p.Flux=Flux;
			} else {
			// Purge time
				p.Type = 2;
				p.Probability=0.;
				p.Flux=0.;
			}
		}

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {

			Rates[p.Type]+=p.Probability;
/*			switch(p.Type) {
				case 0: {	//TDMAT
					Coverages[0]-=
					break;
				}
				case 1: {	//NH3

					break;
				}
				default: {
					assert(0);
				}
			}
*/
		}


		template <class PT, class VecType> void ParticleReflexion(
			const PT& p,
			std::stack<PT>& particle_stack,
			const VecType& NormalVector,
			const double* Coverages,
			int Material,
			int D,
			double dot // dot product between the incoming particle direction and the normal vector
			) const {

			double new_probability;
			if (Coverages[p.Type]>0.) {
				new_probability=p.Probability*(1.-sticking_probability*std::pow(Coverages[p.Type],reaction_order-1));
			} else {
				if (reaction_order<1.) {
					new_probability=0.;
				} else if (reaction_order==1.) {
					new_probability=p.Probability*(1.-sticking_probability);
				} else {
					new_probability=p.Probability;
				}
			}

			if (new_probability>=end_probability) {
				particle_stack.push(p);
				PT& p_new=particle_stack.top();
				p_new.Probability=new_probability;
				my::stat::CosineNDistributedRandomDirection(1.,NormalVector,p_new.Direction);
			}

/*			switch(p.Type) {
				case 0: {	//TDMAT


					break;
				}
				case 1: {	//NH3

					break;
				}
				default: {
					assert(0);
				}
			}
*/

		}

	};

	double const TiN_ALD::end_probability=0.001;
	double const TiN_ALD::reaction_order=1.0;
}

#endif /*MODELTIN_ALD_H_*/
