#ifndef MODELSIMPLEDEPOSITION_H_
#define MODELSIMPLEDEPOSITION_H_

#include <stack>
#include <vector>
#include "../Statistics.h"
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

#include <iostream>

namespace model {

///Simple deposition model
	class SimpleDeposition {

		double TotalFlux;
		std::vector<double> StickingProbabilities;
		std::vector<double> Yields;
		double end_probability;
		double StartAngleDistribution;
		double ReemittedAngleDistribution;

		double StartDirection[3];

	public:

        static const bool OutputFluxes=false;
		static const bool SpatiallyEqualDistributedFlux=true;
		static const bool ReemissionIsMaterialDependent=true;
		static const bool CalculateConnectivities=false;
		static const bool CalculateVisibilities=false;
		static const bool CalculateNormalVectors=false;
        static const int CoverageStorageSize=0;
        static const int RatesStorageSize=1;
        static const unsigned int NumberOfParticleTypes=1;
        unsigned int NumberOfParticleClusters[NumberOfParticleTypes];

		class ParticleType {
		public:
			double Probability;
			double Direction[3];
			double Flux;
		};


		SimpleDeposition(const std::string & Parameters):StartAngleDistribution(1.),ReemittedAngleDistribution(1.) {

		    double Accuracy;

            using namespace boost::spirit::classic;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("flux")  >> '='  >> real_p[assign_a(TotalFlux)]  >> ';') |
                            (str_p("sticking_probabilities")  >> '='  >>  '{' >> (real_p[push_back_a(StickingProbabilities)] % ',') >> '}'  >> ';') |
                            (str_p("start_angle_distribution")  >> '='  >> real_p[assign_a(StartAngleDistribution)]  >> ';') |
                            (str_p("reemitted_angle_distribution")  >> '='  >> real_p[assign_a(ReemittedAngleDistribution)]  >> ';') |
                            (str_p("stop_criterion")  >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("yields") >> '='  >>  '{' >> (real_p[push_back_a(Yields)] % ',') >> '}'  >> ';') |
                            (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';')

                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

			//Check if sticking is between 0 and 1
			for (auto sticking : StickingProbabilities) {
				if (sticking < 0.0) {
					msg::print_warning("Negative sticking probability changed to zero");
					sticking=0.0;
				}
				else if (sticking > 1.0) {
					msg::print_warning("Sticking probability too high, changed to one");
					sticking=1.0;
				}
			}
/*

             end_probability=0.01;
             IonFlux=3.125e15;          // atoms/(cm²s)
             CFxFlux=2e18;              // atoms/(cm²s)
             yDPolyI=10e-24;            // cm³/atom
             yDPolyN=0.5e-24;           // cm³/atom
             StickingProbabilityCFxonPolymer=0.1;
             StickingProbabilityCFxonSi=0.1;
             StickingProbabilityCFxonMask=0.1;
             IonAngleDistribution=410.;
             CFxAngleDistribution=1.;*/

             NumberOfParticleClusters[0]=(TotalFlux>0.)?static_cast<unsigned int>(Accuracy):0;
        }

		template<class VecType>
        void CalculateVelocity(
            double &Velocity,
            const VecType& NormalVector,
            const double *Coverages,
            const double *Rates,
            int Material, bool Connected, bool Visible) const {
			
			double Yield=(Material < static_cast<int>(Yields.size()))?Yields[Material]:0; 
		    Velocity=Rates[0]*Yield;
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


		static void UpdateCoverage(double *Coverages, const double *Rates) {}

		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {
			my::stat::CosineNDistributedRandomDirection(StartAngleDistribution,StartDirection,p.Direction);
			p.Probability=1.;
			p.Flux=TotalFlux;
		}

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {
        	Rates[0]+=p.Flux*p.Probability;
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

			double StickingProbability = (Material < static_cast<int>(StickingProbabilities.size()))?StickingProbabilities[Material]:1.0; 
			double new_probability=p.Probability*(1.-StickingProbability);
			if (new_probability>=end_probability) {
				particle_stack.push(p);
				PT& p_new=particle_stack.top();
				p_new.Probability=new_probability;
				my::stat::CosineNDistributedRandomDirection(ReemittedAngleDistribution,NormalVector,p_new.Direction);
			}

		}
	};


}




#endif /*MODELSIMPLEDEPOSITION_H_*/
