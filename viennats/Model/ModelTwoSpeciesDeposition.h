#ifndef MODELTWOSPECIESDEPOSITION_H_
#define MODELTWOSPECIESDEPOSITION_H_

#include <stack>
#include "../Statistics.h"
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

// This model was presented in:
// M.M. IslamRaja, C. Chang, J.P. McVittie, M.A. Cappelli, and K.C. Saraswat
// "Two precursor model for low pressure chemical vapor deposition of silicon dioxide from tetraethylorthosilicate",
// J. Vac. Sci. technol. B 11 (3), 1993

namespace model {

	class TwoSpeciesDeposition {

            	static const double rho_SiO2; 	//in g/cm^3
                static const double N_A;   //in /mol
                static const double M_SiO2;     //molar mass in g/mol

//		double rate_main;
//		double rate_intermediate;
		double flux_main;
		double flux_intermediate;
		double sticking_main;
		double sticking_intermediate;
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
        static const int CoverageStorageSize=0;
        static const int RatesStorageSize=2;
        static const unsigned int NumberOfParticleTypes=2;
        unsigned int NumberOfParticleClusters[NumberOfParticleTypes];

		class ParticleType {
		public:
			double Probability;
			double Direction[3];
			double Flux;
                        int Type;
		};


		TwoSpeciesDeposition(const std::string & Parameters) {

		    double Accuracy;
                    reaction_order = 1.;

            using namespace boost::spirit::classic;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                    	    (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
//                            (str_p("rate_main")  >> '='  >> real_p[assign_a(rate_main)]  >> ';') |
//                            (str_p("rate_intermediate")  >> '='  >> real_p[assign_a(rate_intermediate)]  >> ';') |
                            (str_p("flux_main")  >> '='  >> real_p[assign_a(flux_main)]  >> ';') |
                            (str_p("flux_intermediate")  >> '='  >> real_p[assign_a(flux_intermediate)]  >> ';') |
                            (str_p("sticking_main")  >> '='  >> real_p[assign_a(sticking_main)]  >> ';') |
                            (str_p("sticking_intermediate")  >> '='  >> real_p[assign_a(sticking_intermediate)]  >> ';') |
                            (str_p("reaction_order")  >> '='  >> real_p[assign_a(reaction_order)]  >> ';') |
                            (str_p("stop_criterion")  >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';')

                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");


//            if (flux_main > 0) rate_main = flux_main/rho_SiO2;
//            if (flux_intermediate > 0) rate_intermediate = flux_intermediate/rho_SiO2;

             end_probability *= std::min(sticking_main,sticking_intermediate);
             NumberOfParticleClusters[0]=static_cast<unsigned int>(Accuracy);
             NumberOfParticleClusters[1]=static_cast<unsigned int>(Accuracy);
        }

		template<class VecType>
        void CalculateVelocity(
            double &Velocity,
            const VecType& NormalVector,
            const double *Coverages,
            const double *Rates,
            int Material, bool Connected, bool Visible) const {

//		    Velocity=deposition_rate*std::pow(Rates[0],reaction_order);
//?                    Velocity = rate_main*Rates[0] + rate_intermediate*Rates[1];
                    Velocity = M_SiO2*(Rates[0] + Rates[1])/(rho_SiO2 * N_A);
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


		static void UpdateCoverage(double *Coverages, const double *Rates) {};

		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position)  const {
			//my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
                        switch (ParticleType) {
                            case 0: {       //main
                                my::stat::CosineNDistributedRandomDirection(1.,StartDirection,p.Direction);
                                p.Probability=1.;
                                p.Flux=flux_main;//1.;
                                p.Type=0;
                                break;
                            }
                            case 1: {           //intermediate
                                my::stat::CosineNDistributedRandomDirection(1.,StartDirection,p.Direction);
                                p.Probability=1.;
                                p.Flux=flux_intermediate;//1.;
                                p.Type=1;
                                break;
                            }
			}

//			my::stat::CosineNDistributedRandomDirection(1.,StartDirection,p.Direction);
//
//			p.Probability=1.;
//                        p.Type=ParticleType;
//			p.Flux=1.;
		}


        template <class PT, class NormVecType> void ParticleCollision(
                            const PT& p,
                            const NormVecType& NormalVector,
                            double* Rates,
                            const double* Coverages,
                            double RelTime) const {
//                                    std::cout << "ParticleCollision inside model \n";

//			Rates[0]+=p.Flux*p.Probability;
//                        Rates[0]+=p.Probability;
                        Rates[p.Type]+=p.Flux*p.Probability;

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

                    double new_probability=0.;

                        switch (p.Type) {
                            case 0: {           //main
                                new_probability=p.Probability*(1.-sticking_main);
                                break;
                            }
                            case 1: {           //intermediate
                                new_probability=p.Probability*(1.-sticking_intermediate);
                                break;
                            }
			}

//			if (Coverages[0]>0.) {
//				new_probability=p.Probability*(1.-sticking_probability*std::pow(Coverages[0],reaction_order-1));
//			} else {
//				if (reaction_order<1.) {
//					new_probability=0.;
//				} else if (reaction_order==1.) {
//					new_probability=p.Probability*(1.-sticking_probability);
//				} else {
//					new_probability=p.Probability;
//				}
//			}

			if (new_probability>=end_probability) {
				particle_stack.push(p);
				PT& p_new=particle_stack.top();
				p_new.Probability=new_probability;
				//my::stat::Cosine1DistributedRandomDirection(NormalVector,p_new.Direction);
				my::stat::CosineNDistributedRandomDirection(1.,NormalVector,p_new.Direction);
			}

		}
	};

//	double const TwoSpeciesDeposition::rho_SiO2 = 2.3e22; 	//in atoms/cm^3 //2.6e22
        double const TwoSpeciesDeposition::rho_SiO2 = 2.2;             //in g/cm^3
        double const TwoSpeciesDeposition::N_A = 6.022140857e23;   //in /mol
        double const TwoSpeciesDeposition::M_SiO2 = 60.1;    //molar mass in g/mol

}




#endif /*MODELTWOSPECIESDEPOSITION_H_*/
