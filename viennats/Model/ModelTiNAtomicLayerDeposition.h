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

            const static double C[]; //{C_NH3_TDMAT, C_TiN_TDMAT, C_TDMAT_NH3, C_TiN_NH3}
            const static double n[]; //{n_NH3, n_TDMAT}

		double molecular_thickness;
                double reduction_ratio;
                double step_size;
                int step;
                double m_NH3;
                double m_TDMAT;

	public:
            
        class ParticleType {
	public:
		double Direction[3];
		double Probability;
		double Flux;
                int Type;
	};

	unsigned int NumberOfParticleClusters[1];

	static const int CoverageStorageSize=12;
	static const int RatesStorageSize=1;

	static const unsigned int NumberOfParticleTypes=1;

        static const bool OutputFluxes=false;
	static const bool SpatiallyEqualDistributedFlux=true;

	static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;

        static const bool CalculateNormalVectors=true;

        static const bool ReemissionIsMaterialDependent=false;
        

        TiN_ALD(const std::string & Parameters, const int &current_step) {

		    using namespace boost::spirit::classic;

                    step=current_step;
                    step_size=1e-4;
                    molecular_thickness=0.433;
                    m_NH3=4.;
                    m_TDMAT=1.;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("step_size")  >> '='  >> real_p[assign_a(step_size)] >> ';') |
                            (str_p("molecular_thickness")  >> '='  >> real_p[assign_a(molecular_thickness)] >> ';') |
                            (str_p("m_NH3")  >> '='  >> real_p[assign_a(m_NH3)] >> ';') |
                            (str_p("m_TDMAT")  >> '='  >> real_p[assign_a(m_TDMAT)] >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");
       		    NumberOfParticleClusters[0]=1;
		}

		template <class VecType>
		void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
                    Velocity = molecular_thickness*std::max(Coverages[5]/m_NH3 + Coverages[11]/m_TDMAT,0.);
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

                        double Coverages0_6 = Coverages[6*Type+0];
                        double Coverages1_7 = Coverages[6*Type+1];

                        delta_time = std::min(
                                        std::max(
                                            (step_size * Coverages1_7) / 
                                            (C[2*Type+1]*std::pow(std::max(Coverages[(6*Type+10)%12] + 
                                                            Coverages0_6 - Coverages1_7,0.), 
                                                         n[Type])), 
                                        step_size),
                                     0.01);

                        Coverages[6*Type+0] = std::max(
                                                Coverages[(6*Type+8)%12] 
                                                - std::pow(
                                                    std::max(
                                                        std::pow(
                                                            Coverages[(6*Type+8)%12] - Coverages0_6,
                                                            1.-n[Type])
                                                        - (1 - n[Type])*C[2*Type+0]*delta_time, 
                                                    0.),
                                                1./(1.-n[Type])), 
                                            0.);

                        Coverages[6*Type+1] = Coverages1_7 
                                              + C[2*Type+1] 
                                                    * delta_time 
                                                        * std::max(
                                                            std::pow(
                                                                std::max(
                                                                    Coverages[(6*Type+10)%12] + Coverages0_6 - Coverages1_7,
                                                                0.)
                                                            ,n[Type]), 
                                                        0.);

                        Coverages[6*Type+2] = Coverages[(6*Type+9)%12]  + Coverages[6*Type+1];
                        Coverages[6*Type+3] = Coverages[(6*Type+8)%12]  - Coverages[6*Type+0];
                        Coverages[6*Type+4] = Coverages[(6*Type+10)%12] + Coverages[6*Type+0] - Coverages[6*Type+1];
                        Coverages[6*Type+5] = std::max((Coverages[6*Type+0] - Coverages0_6)/delta_time, 0.);
		}

		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {

//                    	p.Type=ParticleType;
//			switch(ParticleType) {
//				case 0:	//TTIP
//					my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
//					p.Probability=1.;
//					p.Flux=1.-step;
//					break;
//				case 1:	//H2O
//					my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
//					p.Probability=1.;
//					p.Flux=step;
//					break;
//				default:
//					assert(0);
//			}
		}

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {
//                    Rates[p.Type]+=p.Flux;
		}


		template <class PT, class VecType> void ParticleReflexion(
			const PT& p,
			std::stack<PT>& particle_stack,
			const VecType& NormalVector,
			const double* Coverages,
			int Material,
			int D,
			double dot // dot product between the incoming particle direction and the normal vector
			) const {}
	};

            double const TiN_ALD::C[] = {0.75, 0.8, 3.8, 1.2};
            double const TiN_ALD::n[] = {0.9, 1.2};
}

#endif /*MODELTIN_ALD_H_*/