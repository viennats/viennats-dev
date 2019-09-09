#ifndef MODELTIO2_ALD_H_
#define MODELTIO2_ALD_H_

#include <stack>
#include "../Statistics.h"
#include <cassert>
#include <string>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

namespace model {

///TiO2 ALD model

    	class TiO2_ALD {

//            double times;

            const static double temperature;
            const static double k_B;
            const static double gas;
            const static double kJm_eV;

//            const static double s_TTIP;
//            const static double s_H2O;
            const static double kdes_TTIP;
            const static double kdes_H2O;

            const static double Cov0_TTIP;
            const static double Cov0_H2O;
            const static double k2_pyrolysis;
            const static double k_hydrolysis;

            const static double error;

                double molecular_thickness;
                double end_probability;
		double reaction_order;
                double FluxTTIP;
                double FluxH2O;

		double StartDirection[3];

		//double sticking_TTIP;
		//double sticking_H2O;
                double sticking[2];

	public:

        class ParticleType {
	public:
		double Direction[3];
		double Probability;
		double Flux;
		int Type;
	};

//	unsigned int NumberOfParticleClusters[3];

	static const int CoverageStorageSize=2;
	static const int RatesStorageSize=2;

	static const unsigned int NumberOfParticleTypes=2;
	unsigned int NumberOfParticleClusters[NumberOfParticleTypes];

        static const bool OutputFluxes=false;
	static const bool SpatiallyEqualDistributedFlux=true;

	static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;

        static const bool CalculateNormalVectors=true;

        static const bool ReemissionIsMaterialDependent=true;

        TiO2_ALD(const std::string & Parameters) {

		    using namespace boost::spirit::classic;

		    double Accuracy;
                    molecular_thickness = 0.5e-8;//0.5Angstrom (cm)
            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction")            >> '='  >> '{' >> real_p[assign_a(StartDirection[0])] >> ","
                    						          >> real_p[assign_a(StartDirection[1])] >> ","
                    						          >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("reaction_order")       >> '='  >> real_p[assign_a(reaction_order)]  >> ';') |

                            //(str_p("sticking_TTIP")        >> '='  >> real_p[assign_a(sticking[0])]  >> ';') |
                            //(str_p("sticking_H2O")         >> '='  >> real_p[assign_a(sticking[1])]  >> ';') |
                            (str_p("end_probability")      >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("TTIP_Flux")            >> '='  >> real_p[assign_a(FluxTTIP)]  >> ';') |
                            (str_p("H2O_Flux")             >> '='  >> real_p[assign_a(FluxH2O)]  >> ';') |
                            (str_p("statistical_accuracy") >> '='  >>real_p[assign_a(Accuracy)]  >> ';') |
                            (str_p("molecular_thickness")  >> '='  >> real_p[assign_a(molecular_thickness)] >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if(!b) msg::print_error("Failed interpreting process parameters!");

		    unsigned int num_particles=static_cast<unsigned int>(Accuracy);

       		    NumberOfParticleClusters[0]=(FluxTTIP>0.)?num_particles:0;		//TTIP flux
		    NumberOfParticleClusters[1]=(FluxH2O>0.)?num_particles:0;		//H2O flux
//                    std::cout<< << "\n";
//                    time_step=0.001;
//                    this->setTimes(1.);
//                    timesint=1000;
//                    std::cout << "times int = " << timesint << "\n";
                    //times=10;
//                    std::cout << "times = " << times << "\n";
		}

		template <class VecType>
		void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
//                     k2_pyrolysis*kdes_TTIP*Coverages0 + k_hydrolysis*Cov0_H2O*Coverages1
//                    Velocity = 0.;//std::max(Coverages[1]*Coverages[0]*Cov0_TTIP*Cov0_H2O/k_hydrolysis + Coverages[1]*Coverages[0]/Cov0_H2O,0.);
//                    std::cout << "Velocity = " << Velocity << "\n";
                    Velocity = molecular_thickness*std::max(Coverages[0]*Coverages[3]+Coverages[1]*Coverages[2],0.);
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

//		static void UpdateCoverage(double *Coverages, const double *Rates, double &time, double &CurrentTime) {
                void UpdateCoverage(double *Coverages, const double *Rates, double &delta_time) const {//, double &CurrentTime) const {
                    //delta_time=1e-5;
//                    std::cout << "    Coverages[0] = " << Coverages[0] << ", ";
//                    std::cout << "    Coverages[1] = " << Coverages[1] << "\n";
                    const double Coverages0 = Coverages[0]*Cov0_TTIP;
                    const double Coverages1 = Coverages[1]*Cov0_H2O;

                    const double C_TTIP = (Rates[0]*sticking[0]//*(1.-Coverages[0])
                                        - kdes_TTIP*Coverages0
                                        - k2_pyrolysis*Coverages0*Coverages0
                                        - k_hydrolysis*Coverages0*Coverages1)/Cov0_TTIP;
//                                C_TTIP /= Cov0_TTIP;
                    const double t_TTIP = std::min(0.01/std::abs(C_TTIP),0.001);
                    const double C_H2O  = (Rates[1]*sticking[1]//*(1.-Coverages[1])
                                        - kdes_H2O *Coverages1
                                        - k_hydrolysis*Coverages0*Coverages1)/Cov0_H2O;
//                                C_H2O /= Cov0_H2O;
                    const double t_H2O = std::min(0.01/std::abs(C_H2O),0.001);

                    delta_time = std::min(t_TTIP,t_H2O);

                    Coverages[0] += delta_time*(C_TTIP);
                    Coverages[1] += delta_time*(C_H2O);
                    Coverages[0] = Coverages[0]>=1.?1.:Coverages[0];
                    Coverages[1] = Coverages[1]>=1.?1.:Coverages[1];
                    Coverages[2] = C_TTIP;
                    Coverages[3] = C_H2O;
//                    std::cout << "1 - Coverages[0] = " << Coverages[0] << ", ";
//                    std::cout << "1 - Coverages[1] = " << Coverages[1] << "\n";

//                    const double A0 = Rates[0]/Cov0_TTIP;    // 2e3
//                    const double B0 = kdes_TTIP + k2_pyrolysis*kdes_TTIP*Coverages0 + k_hydrolysis*Cov0_H2O*Coverages1 + A0; // 7.803e11
//
//                    const double A1 = Rates[1]/Cov0_H2O;     // 2.427e4
//                    const double B1 = kdes_H2O + k_hydrolysis*Cov0_TTIP*Coverages0 + A1;   // 8.574e13
//                    delta_time = std::min(error/B0,error/B1);
//                    //time = 0.001;
//
//                    Coverages[0] = std::max(A0/B0 + (Coverages0 - A0/B0)*std::exp(-B0*(delta_time)),0.);   // 2.563e-9
//                    Coverages[1] = std::max(A1/B1 + (Coverages1 - A1/B1)*std::exp(-B1*(delta_time)),0.);   // 2.830e-10
		}

		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {

			p.Type=ParticleType;
			switch(ParticleType) {
				case 0:	//TTIP
					my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
					p.Probability=1.;
					p.Flux=FluxTTIP;
					break;
				case 1:	//H2O
					my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
					p.Probability=1.;
					p.Flux=FluxH2O;
					break;
				default:
					assert(0);
			}
		}

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {

            Rates[p.Type]+=p.Flux*p.Probability*std::max(1.-Coverages[p.Type],0.);
//                        Rates[p.Type]+=p.Flux*p.Probability*sticking[p.Type];//*std::max(1.-Coverages[0],0.);

/*
                    switch(p.Type) {
			case 0: {   //TTIP
//                            const double StickingProbability=std::max(s_TTIP*(1.-Coverages[0]),0.);
                            Rates[0]+=p.Flux*p.Probability;//StickingProbability;
                            break;
			}
			case 1: {   //H2O
//                            const double StickingProbability=std::max(s_H2O*(1.-Coverages[1]),0.);
                            Rates[1]+=p.Flux*p.Probability;//StickingProbability;
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
							int Material//,
//                            int D,
//                            double dot // dot product between the incoming particle direction and the normal vector
			) const {

                    double new_probability;

//                    	if (Coverages[p.Type]>0.) {
//				new_probability=p.Probability*(1.-sticking[p.Type]*std::pow(Coverages[p.Type],reaction_order-1));
//			} else {
//				if (reaction_order<1.) {
//					new_probability=0.;
//				} else if (reaction_order==1.) {
					new_probability=std::max(p.Probability*(1.-sticking[p.Type]),0.);
//				} else {
//					new_probability=p.Probability;
//				}
//			}
//                            const double new_probability=p.Probability*(1-StickingProbability);
                            if (new_probability>=end_probability) {
                                particle_stack.push(p);
                                PT& p_new=particle_stack.top();
                                p_new.Probability=new_probability;
				//my::stat::Cosine1DistributedRandomDirection(NormalVector,p_new.Direction);
				my::stat::CosineNDistributedRandomDirection(1.,NormalVector,p_new.Direction);
                            }
/*
                    switch(p.Type) {
			case 0: {   //TTIP
//                            const double StickingProbability=std::max(s_TTIP*(1.-Coverages[0]),0.);
			if (Coverages[0]>0.) {
				new_probability=p.Probability*(1.-sticking_TTIP*std::pow(Coverages[0],reaction_order-1));
			} else {
				if (reaction_order<1.) {
					new_probability=0.;
				} else if (reaction_order==1.) {
					new_probability=p.Probability*(1.-sticking_TTIP);
				} else {
					new_probability=p.Probability;
				}
			}
//                            const double new_probability=p.Probability*(1-StickingProbability);
                            if (new_probability>=end_probability) {
                                particle_stack.push(p);
                                PT& p_new=particle_stack.top();
                                p_new.Probability=new_probability;
                                my::stat::Cosine1DistributedRandomDirection(NormalVector,p_new.Direction);
                            }
                            break;
			}
			case 1: {   //H2O
//                            const double StickingProbability=std::max(s_H2O*(1.-Coverages[1]),0.);
//                            const double new_probability=p.Probability*(1-StickingProbability);
			if (Coverages[1]>0.) {
				new_probability=p.Probability*(1.-sticking_H2O*std::pow(Coverages[1],reaction_order-1));
			} else {
				if (reaction_order<1.) {
					new_probability=0.;
				} else if (reaction_order==1.) {
					new_probability=p.Probability*(1.-sticking_H2O);
				} else {
					new_probability=p.Probability;
				}
			}
                            if (new_probability>=end_probability) {
                                particle_stack.push(p);
				PT& p_new=particle_stack.top();
				p_new.Probability=new_probability;
				my::stat::Cosine1DistributedRandomDirection(NormalVector,p_new.Direction);
                            }
                            break;
                        }
			default: {
                            assert(0);
                        }
                    }
*/
		}
	};

//	class TiO2_ALD_H2O {
//
//		double end_probability;
//		double reaction_order;
//
//		double StartDirection[3];
//
//		double sticking_probability;
//                double ProcTime;
//
//	public:
//
//	class ParticleType {
//	public:
//		double Direction[3];
//		double Probability;
//		double Flux;
//	};
//
//	unsigned int NumberOfParticleClusters[3];
//
//	static const int CoverageStorageSize=10;
//	static const int RatesStorageSize=2;
//
//	static const unsigned int NumberOfParticleTypes=2;
//
//        static const bool OutputFluxes=false;
//	static const bool SpatiallyEqualDistributedFlux=true;
//
//	static const bool CalculateVisibilities=false;
//        static const bool CalculateConnectivities=false;
//
//        static const bool CalculateNormalVectors=true;
//
//        static const bool ReemissionIsMaterialDependent=true;
//
//        TiO2_ALD_H2O (const std::string & Parameters, const double ProcessTime) {
//
//		    using namespace boost::spirit::classic;
//
//		    double Accuracy;
//
//            bool b = parse(
//                    Parameters.begin(),
//                    Parameters.end(),
//                    *(
//                            (str_p("direction")            >> '='  >> '{' >> real_p[assign_a(StartDirection[0])] >> ","
//                    						          >> real_p[assign_a(StartDirection[1])] >> ","
//                    						          >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
//                            (str_p("reaction_order")       >> '='  >> real_p[assign_a(reaction_order)]  >> ';') |
//                            (str_p("sticking_probability") >> '='  >> real_p[assign_a(sticking_probability)]  >> ';') |
//                            (str_p("end_probability")      >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
//                            (str_p("statistical_accuracy") >> '='  >>real_p[assign_a(Accuracy)]  >> ';')
//                    ),
//                    space_p | comment_p("//") | comment_p("/*", "*/")).full;
//
//            if (!b) msg::print_error("Failed interpreting process parameters!");
//
//		    unsigned int num_particles=static_cast<unsigned int>(Accuracy);
//
//		    NumberOfParticleClusters[0]=num_particles;	//Flux
//                    ProcTime=ProcessTime;
//		}
//
//		template <class VecType>
//		void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
//                    Velocity=0.;//1e18*std::pow(Rates[0],reaction_order);
////                    std::cout<<"Vel = " << Velocity<<"\n";
//
//		}
//
//		template<class VecType>
//		void CalculateVectorVelocity(
//				double *Velocity,
//				const VecType& NormalVector,
//				const double *Coverages,
//				const double *Rates,
//				int Material,
//				bool connected,
//				bool visible
//				) const {}
//
//		static void UpdateCoverage(double *Coverages, const double *Rates) {
//
////                    Coverages[2]=std::max(Rates[1],Coverages[2]);
////                    Coverages[1]=std::max(Coverages[1]-Rates[0],0.);
////                    Coverages[0]=1.-Coverages[1]-Coverages[2];
//
//                    Coverages[5]=Rates[0];
//                    Coverages[6]=Rates[1];
//                    Coverages[7]=std::max(Coverages[2]-Coverages[5],0.);
//                    Coverages[8]=Coverages[3]+Coverages[6];
//                    Coverages[9]=std::max(Coverages[4]+Coverages[5]-Coverages[6],0.);
//
//		}
//
//		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {
//
//                    my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
//                    p.Probability=1.;
//                    p.Flux=1.;
//		}
//
//        template <class PT, class NormVecType> void ParticleCollision(
//                                    const PT& p,
//                                    const NormVecType& NormalVector,
//                                    double* Rates,
//                                    const double* Coverages,
//                                    double RelTime) const {
////            if (Coverages[1]>0.) std::cout<<"Coverages[1]="<<Coverages[1]<<"\n";
////			Rates[0]+=p.Flux*p.Probability*std::max(Coverages[1],0.);
////                        Rates[1]+=p.Flux*p.Probability*0.2*std::max(1-Coverages[1]-Coverages[2],0.);
//
//                    double C_T_N=3.8;
//                    double C_TiN_N=1.2;
//                    double n_N=1.2;
//
//                    Rates[0]+=std::max(C_T_N  * std::pow(Coverages[2]-Rates[0],         n_N),0.);
//                    Rates[1]+=std::max(C_TiN_N* std::pow(Coverages[4]+Rates[0]-Rates[1],n_N),0.);
//
//
//		}
//
//
//		template <class PT, class VecType> void ParticleReflexion(
//			const PT& p,
//			std::stack<PT>& particle_stack,
//			const VecType& NormalVector,
//			const double* Coverages,
//			int Material,
//			int D,
//			double dot // dot product between the incoming particle direction and the normal vector
//			) const {
//
//			double new_probability;
//                        if ((Coverages[5]+Coverages[6])>0.) {
//				new_probability=p.Probability*(1.-sticking_probability*std::pow((Coverages[5]+Coverages[6]),reaction_order-1.));
//			} else {
//				if (reaction_order<1.) {
//					new_probability=0.;
//				} else if (reaction_order==1.) {
//					new_probability=p.Probability*(1.-sticking_probability);
//				} else {
//					new_probability=p.Probability;
//				}
//			}
//
//			if (new_probability>=end_probability) {
//				particle_stack.push(p);
//				PT& p_new=particle_stack.top();
//				p_new.Probability=new_probability;
//				my::stat::CosineNDistributedRandomDirection(1.,NormalVector,p_new.Direction);
//			}
//		}
//
//	};

        double const TiO2_ALD::temperature=200+273.15; //K
        double const TiO2_ALD::k_B=8.6173303e-5;     // eV/K
        double const TiO2_ALD::kJm_eV=0.010364272301331;
        double const TiO2_ALD::gas=8.3144598;

//        double const TiO2_ALD::s_TTIP=0.3;     // a.u.
//        double const TiO2_ALD::s_H2O=0.91;     // a.u.
        double const TiO2_ALD::Cov0_TTIP=1.5e14;   // /cm2
        double const TiO2_ALD::Cov0_H2O=1.5e15;    // /cm2
//        double const TiO2_ALD::kdes_TTIP=2.0e10*std::exp(-114*TiO2_ALD::kJm_eV/(TiO2_ALD::k_B*TiO2_ALD::temperature));       // /s
//        double const TiO2_ALD::kdes_H2O=5.0e9*std::exp(-81*TiO2_ALD::kJm_eV/(TiO2_ALD::k_B*TiO2_ALD::temperature));        // /s
//        double const TiO2_ALD::k2_pyrolysis=1.3e3*std::exp(-185*TiO2_ALD::kJm_eV/(TiO2_ALD::k_B*TiO2_ALD::temperature));    // cm2/s
//        double const TiO2_ALD::k_hydrolysis=4.7e-7*std::exp(-86*TiO2_ALD::kJm_eV/(TiO2_ALD::k_B*TiO2_ALD::temperature));    // cm2/s

        double const TiO2_ALD::kdes_TTIP    = 2.0e10*std::exp(-114000/(TiO2_ALD::gas*TiO2_ALD::temperature));       // /s
        double const TiO2_ALD::kdes_H2O     = 5.0e9* std::exp(- 81000/(TiO2_ALD::gas*TiO2_ALD::temperature));        // /s
        double const TiO2_ALD::k2_pyrolysis = 1.3e3* std::exp(-185000/(TiO2_ALD::gas*TiO2_ALD::temperature));    // cm2/s
        double const TiO2_ALD::k_hydrolysis = 4.7e-7*std::exp(- 86000/(TiO2_ALD::gas*TiO2_ALD::temperature));    // cm2/s

        double const TiO2_ALD::error=0.001;     // a.u.
}

#endif /*MODELTIO2_ALD_H_*/
