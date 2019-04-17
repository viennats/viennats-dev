#ifndef MODELTIN_PEALD_H_
#define MODELTIN_PEALD_H_

#include <stack>
#include "../Statistics.h"
#include <cassert>
#include <string>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

namespace model {

///TiN ALD model

    	class TiN_PEALD {

            const static double C[]; //{C_NH3_TDMAT, C_TiN_TDMAT, C_TDMAT_NH3, C_TiN_NH3}
            const static double n[]; //{n_NH3, n_TDMAT}

		double molecular_thickness;
                //double reduction_ratio;
                double step_size;
                int step;
                double m_Plasma;
                double m_TDMAT;
                double sticking_probability;
		double StartDirection[3];
		double end_probability;
		//double reaction_order;
//                double TDMAT_flux;
//                double NH3_flux;
                double Flux;


	public:

        class ParticleType {
	public:
		double Direction[3];
		double Probability;
		double Flux;
                int Type;
	};

	static const int CoverageStorageSize=12;
	static const int RatesStorageSize=1;

	static const unsigned int NumberOfParticleTypes=1;
	unsigned int NumberOfParticleClusters[NumberOfParticleTypes];

        static const bool OutputFluxes=false;
	static const bool SpatiallyEqualDistributedFlux=true;

	static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;

        static const bool CalculateNormalVectors=true;

        static const bool ReemissionIsMaterialDependent=false;


        TiN_PEALD(const std::string & Parameters, const int &current_step) {

		    double Accuracy;
		    using namespace boost::spirit::classic;

                    step=current_step;
                    step_size=1e-4;
                    molecular_thickness=0.433;
                    m_Plasma=13.58;
                    m_TDMAT=1.;
                    Flux=1.;
//                    if (step == 1){
////                        std::cout << "AGH\n";
//                        TDMAT_flux = 1.;
//                        NH3_flux = 0.;
//                    } else {
////                        std::cout << "AGH2\n";
//                        TDMAT_flux = 0.;
//                        NH3_flux = 1.;
//                    }

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("step_size")  >> '='  >> real_p[assign_a(step_size)] >> ';') |
                            (str_p("molecular_thickness")  >> '='  >> real_p[assign_a(molecular_thickness)] >> ';') |
//                            (str_p("m_Plasma")  >> '='  >> real_p[assign_a(m_Plasma)] >> ';') |
//                            (str_p("m_TDMAT")  >> '='  >> real_p[assign_a(m_TDMAT)] >> ';') |
//                            (str_p("reaction_order")  >> '='  >> real_p[assign_a(reaction_order)]  >> ';') |
                            (str_p("stop_criterion")  >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("sticking_probability")  >> '='  >> real_p[assign_a(sticking_probability)] >> ';') |
                            (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';')


                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");
//       		  NumberOfParticleClusters[0]=1;
                    NumberOfParticleClusters[0]=static_cast<unsigned int>(Accuracy);

		}

		template <class VecType>
		void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
                    Velocity = molecular_thickness*std::max(Coverages[11]/m_TDMAT,0.);
//                    Velocity = molecular_thickness*std::max(Coverages[5]/m_Plasma + Coverages[11]/m_TDMAT,0.);
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

//		static void UpdateCoverage(double *Coverages, const double *Rates) {
//		static void UpdateCoverage(double *Coverages, const double *Rates, double &delta_time, double &CurrentTime) {
                    void UpdateCoverage(double *Coverages, const double *Rates, double &delta_time) const {//, double &CurrentTime) const {

                    int Type = step-1;

                        //TDMAT on NH3 : NH3 on TDMAT
                        double Coverages0_6 = Coverages[6*Type+0];
                        //TDMAT on TiN : NH3 on TiN
                        double Coverages1_7 = Coverages[6*Type+1];

                        delta_time = step_size;
//                        std::cout << "step_size = " << step_size << "\n";
//                        std::cout << "Coverages1_7 = " << Coverages1_7 << "\n";
//                        std::cout << "C[2*Type+1] = " << C[2*Type+1] << "\n";
//                        std::cout << "Coverages[(6*Type+10)%12] = " << Coverages[(6*Type+10)%12] << "\n";
//                        std::cout << "Coverages0_6 = " << Coverages0_6 << "\n";
//                        std::cout << "n[Type] = " << n[Type] << "\n";
//                        delta_time = std::min(
//                                        std::max(
//                                            (step_size * Coverages1_7) /
//                                            (C[2*Type+1]*std::pow(std::max(Coverages[(6*Type+10)%12] +
//                                                            Coverages0_6 - Coverages1_7,0.),
//                                                         n[Type])),
//                                        step_size),
//                                     0.01);

                        Coverages[6*Type+0] = std::max(
                                                Coverages[(6*Type+8)%12]
                                                - std::pow(
                                                    std::max(
                                                        std::pow(
                                                            Coverages[(6*Type+8)%12] - Coverages0_6,
                                                            1.-n[Type])
                                                        - (1 - n[Type])*Rates[0]*C[2*Type+0]*delta_time,
                                                    0.),
                                                1./(1.-n[Type])),
                                            0.);

                        Coverages[6*Type+1] = Coverages1_7
                                              + Rates[0]*C[2*Type+1]
                                                    * delta_time
                                                        * std::max(
                                                            std::pow(
                                                                std::max(
                                                                    Coverages[(6*Type+10)%12] + Coverages0_6 - Coverages1_7,
                                                                0.)
                                                            ,n[Type]),
                                                        0.);

                        // NH3 or TDMAT after a NH3 or TDMAT step
                        Coverages[6*Type+2] = Coverages[(6*Type+9)%12]  + Coverages[6*Type+1];
                        // NH3 or TDMAT after a TDMAT or NH3 step
                        Coverages[6*Type+3] = Coverages[(6*Type+8)%12]  - Coverages[6*Type+0];
                        // TiN after a TDMAT or NH3 step
                        Coverages[6*Type+4] = Coverages[(6*Type+10)%12] + Coverages[6*Type+0] - Coverages[6*Type+1];
                        // Change in coverage for the relevant species
                        Coverages[6*Type+5] = std::max((Coverages[6*Type+0] - Coverages0_6)/delta_time, 0.);
		}

		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {

                    	p.Type=ParticleType;
//                        std::cout << "ParticleType = " << ParticleType << "\n";
//                        std::cout << "TDMAT_flux = " << TDMAT_flux << "\n";
//                        std::cout << "NH3_flux = " << NH3_flux << "\n";
//			switch(ParticleType) {
//				case 0:	//TDMAT
//					my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
//					p.Probability=1.;
//					p.Flux=TDMAT_flux;
//					break;
//				case 1:	//NH3
//					my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
//					p.Probability=1.;
//					p.Flux=NH3_flux;
//					break;
//				default:
//					assert(0);
//			}

                        my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
			p.Probability=1.;
			p.Flux=Flux;

		}

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {
//                        std::cout << "ParticleCollision = " << "\n";

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
                    	double new_probability;
//                        std::cout << "ParticleCollision = " << "\n";

//			if (Coverages[0]>0.) {
				new_probability=p.Probability*(1.-sticking_probability);
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

              double const TiN_PEALD::C[] = {9.0, 7.5, 6.0, 6.0};
              double const TiN_PEALD::n[] = {1.6, 1.05};
}

#endif /*MODELTIN_PEALD_H_*/
