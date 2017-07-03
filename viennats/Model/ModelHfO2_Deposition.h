#ifndef MODELHFO2DEPOSITION_H_
#define MODELHFO2DEPOSITION_H_

#include <stack>
#include "../Statistics.h"
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

// This model was presented in:
// T.S. Yang et al, 
//"Chemical Vapor Deposition of HfO2 Thin Films Using the Novel Single Precursor Hafnium 3-Methyl-3-pentoxide, Hf(mp)4", 
//Chem. Mater. 2005, 17, 6713-6718

namespace model {

	class HfO2_Deposition {

		double deposition_rate;
                double temperature;
                double pressure;
		double sticking_probability;
		double reaction_order;

		double end_probability;
		double StartDirection[3];

	public:
            
        static constexpr double rate0=6514.006606266484;
        static constexpr double E=68.14;
        static constexpr double kB=0.0083144621;
        static constexpr double rho_HfO2 = 9680;
        static constexpr double M_HfO2 = 0.21049;
        static constexpr double M_gas = 0.58315728;

        static const bool OutputFluxes=false;
		static const bool SpatiallyEqualDistributedFlux=true;
		static const bool ReemissionIsMaterialDependent=true;
		static const bool CalculateConnectivities=false;
		static const bool CalculateVisibilities=false;
		static const bool CalculateNormalVectors=false;
        static const int CoverageStorageSize=1;
        static const int RatesStorageSize=1;
        static const unsigned int NumberOfParticleTypes=1;
        unsigned int NumberOfParticleClusters[NumberOfParticleTypes];

		class ParticleType {
		public:
			double Probability;
			double Direction[3];
			double Flux;
		};


		HfO2_Deposition(const std::string & Parameters) {

		    double Accuracy;

            using namespace boost::spirit::classic;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                    	    (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("temperature")  >> '='  >> real_p[assign_a(temperature)]  >> ';') |
                            (str_p("pressure")  >> '='  >> real_p[assign_a(pressure)]  >> ';') |
                            (str_p("stop_criterion")  >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';')

                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");
            temperature += 273.15;
            deposition_rate = rate0*std::exp(-E/(kB*temperature));
            std::cout << "deposition rate = " << deposition_rate << "\n";
            sticking_probability = 1e-10*deposition_rate*rho_HfO2/(M_HfO2*pressure)*std::sqrt(2000.*kB*temperature*M_gas);
            std::cout << "sticking probability = " << sticking_probability << "\n";
            reaction_order = 1.0;

             NumberOfParticleClusters[0]=static_cast<unsigned int>(Accuracy);
        }

		template<class VecType>
        void CalculateVelocity(
            double &Velocity,
            const VecType& NormalVector,
            const double *Coverages,
            const double *Rates,
            int Material, bool Connected, bool Visible) const {

		    Velocity=deposition_rate*std::pow(Rates[0],reaction_order);
                    std::cout << "Velocity = " << Velocity << "\n";
                    std::cout << "reaction_order = " << reaction_order << "\n";
                    std::cout << "Rates[0] = " << Rates[0] << "\n";
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
			Coverages[0]=Rates[0];
		};

		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position)  const {
			//my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
			my::stat::CosineNDistributedRandomDirection(1.,StartDirection,p.Direction);

			p.Probability=1.;
//			p.Flux=1.;
		}


        template <class PT, class NormVecType> void ParticleCollision(
                            const PT& p,
                            const NormVecType& NormalVector,
                            double* Rates,
                            const double* Coverages,
                            double RelTime) const {
//                                    std::cout << "ParticleCollision inside model \n";
                        std::cout << "p.Probability = " << p.Probability << "\n";
//			Rates[0]+=p.Flux*p.Probability;
                        Rates[0]+=p.Probability;
                        
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

//                    std::cout << "ParticleReflexion inside model \n";

			double new_probability;

			if (Coverages[0]>0.) {
				new_probability=p.Probability*(1.-sticking_probability*std::pow(Coverages[0],reaction_order-1));
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
				//my::stat::Cosine1DistributedRandomDirection(NormalVector,p_new.Direction);
				my::stat::CosineNDistributedRandomDirection(1.,NormalVector,p_new.Direction);
			}

		}
	};


}




#endif /*MODELHFO2DEPOSITION_H_*/
