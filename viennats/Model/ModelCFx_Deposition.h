#ifndef MODELCFX_DEPOSITION_H_
#define MODELCFX_DEPOSITION_H_

#include <stack>
#include "../Statistics.h"
#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include <string>
#include <cassert>

namespace model {
///Model for the deposition of

	class CFx_Deposition {

		double end_probability;
		double CFxFlux;           // atoms/(cm²s)
		double IonFlux;           // atoms/(cm²s)
		double yDPolyI;
        double yDPolyN;
        double StickingProbabilityCFx;
        double IonAngleDistribution;
        double CFxAngleDistribution;

		double StartDirection[3];

	public:

		static const bool OutputFluxes=false;
        static const bool SpatiallyEqualDistributedFlux=true;
        static const int CoverageStorageSize=0;
        static const bool ReemissionIsMaterialDependent=false;
        static const bool CalculateNormalVectors=false;
        static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;
        static const int RatesStorageSize=2;
        static const unsigned int NumberOfParticleTypes=2;
        unsigned int NumberOfParticleClusters[NumberOfParticleTypes];

		class ParticleType {
		public:
			double Direction[3];
			double Flux;
			double Probability;
			int Type;
		};

		CFx_Deposition(const std::string & Parameters) {

		    using namespace boost::spirit::classic;

            double Accuracy;
            yDPolyI=10e-24;            // cm³/atom
            yDPolyN=0.5e-24;           // cm³/atom
            StickingProbabilityCFx=0.01;
            IonAngleDistribution=410.;
            CFxAngleDistribution=1.;
            end_probability=0.0001;

            bool b = parse(                 //TODO improve parser
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("end_probability")  >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("CFxFlux")  >> '='  >> real_p[assign_a(CFxFlux)]  >> ';') |
                            (str_p("IonFlux")  >> '='  >> real_p[assign_a(IonFlux)]  >> ';') |
                            (str_p("statistical_accuracy")    >> '='  >>real_p[assign_a(Accuracy)]  >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

             NumberOfParticleClusters[0]=(IonFlux>0.)?static_cast<unsigned int>(Accuracy):0;
             NumberOfParticleClusters[1]=(CFxFlux>0.)?static_cast<unsigned int>(Accuracy):0;
        }


		template <class VecType>
				void CalculateVelocity(
						double &Velocity,
						const VecType& NormalVector,
						const double *Coverages,
						const double *Rates,
						int Material,
						bool Connected,
						bool Visible) const {
            Velocity=yDPolyI*Rates[0]+yDPolyN*Rates[1];
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

		template <class PT> void ParticleGeneration(
				PT& p,
				int ParticleType,
				double ProcessTime,
				double* Position) const
		{
        	switch (ParticleType) {
				case 0: {       //Ion
					my::stat::CosineNDistributedRandomDirection(IonAngleDistribution,StartDirection,p.Direction);
					p.Probability=1.;
					p.Flux=IonFlux;
					p.Type=0;
					break;
				}
				case 1: {           //CFx
					my::stat::CosineNDistributedRandomDirection(CFxAngleDistribution,StartDirection,p.Direction);
					p.Probability=1.;
					p.Flux=CFxFlux;
					p.Type=1;
					break;
				}
			}
		}

		template <class PT, class NormVecType> void ParticleCollision(
				const PT& p,
				const NormVecType& NormalVector,
				double* Rates,
				const double* Coverages,
				double RelTime) const {
			Rates[p.Type]+=p.Flux*p.Probability;
                        //std::cout << "ParticleCollision inside model \n";
		}

		template <class PT, class VecType> void ParticleReflexion(
				const PT& p,
				std::stack<PT>& particle_stack,
				const VecType& NormalVector,
				const double* Coverages,
				int Material//,
//				int D,
//				double dot
                ) const
		{
                    //std::cout << "ParticleReflexion inside model \n";
			switch(p.Type) {
				case 0: {       //ION
					break;
				}
				case 1: {
                    const double StickingProbability=StickingProbabilityCFx;

					double new_probability=p.Probability*(1.-StickingProbability);
					if (new_probability>=end_probability) {
						particle_stack.push(p);
						PT& p_new=particle_stack.top();
						p_new.Probability=new_probability;
						my::stat::CosineNDistributedRandomDirection(1.,NormalVector,p_new.Direction);
					}
					break;
				}
			}

		}
	};
}
#endif /*MODELCFX_DEPOSITION_H_*/
