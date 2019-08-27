#ifndef MODELNONLINEARDEPOSITION_H_
#define MODELNONLINEARDEPOSITION_H_

#include <stack>
#include "../Statistics.h"
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

// This model was presented in:
// T.S. Cale and G.B. Raupp,
// "A unified line-of-sight model of deposition in rectangular trenches",
// J. Vac. Sci. technol. B 8 (6), 1990

namespace model {

	class NonlinearDeposition {

		std::vector<double> deposition_rates;
		double distribution_order=1;
		double sticking_probability;
		double reaction_order;

		double end_probability;
		double StartDirection[3];

	public:

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


		NonlinearDeposition(const std::string & Parameters) {

		    double Accuracy;

            using namespace boost::spirit::classic;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                    	    (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("deposition_rates")  >> '='  >> '{' >> (real_p[push_back_a(deposition_rates)] % ',') >> '}'  >> ';') |
                            (str_p("sticking_probability")  >> '='  >> real_p[assign_a(sticking_probability)]  >> ';') |
                            (str_p("reaction_order")  >> '='  >> real_p[assign_a(reaction_order)]  >> ';') |
                            (str_p("stop_criterion")  >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';') |
														(str_p("distribution_order")  >> '='  >>real_p[assign_a(distribution_order)]  >> ';')

                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");





             NumberOfParticleClusters[0]=static_cast<unsigned int>(Accuracy);
        }

		template<class VecType>
        void CalculateVelocity(
            double &Velocity,
            const VecType& NormalVector,
            const double *Coverages,
            const double *Rates,
            int Material, bool Connected, bool Visible) const {
				if(Material<deposition_rates.size())
		    	Velocity=deposition_rates[Material]*std::pow(Rates[0],reaction_order);
				else Velocity = 0;
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
			my::stat::CosineNDistributedRandomDirection(distribution_order,StartDirection,p.Direction);

			p.Probability=1.;
			p.Flux=1.;
		}


        template <class PT, class NormVecType> void ParticleCollision(
                            const PT& p,
                            const NormVecType& NormalVector,
                            double* Rates,
                            const double* Coverages,
                            double RelTime) const {
//                                    std::cout << "ParticleCollision inside model \n";

			Rates[0]+=p.Flux*p.Probability;
//                        Rates[0]+=p.Probability;
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




#endif /*MODELNONLINEARDEPOSITION_H_*/
