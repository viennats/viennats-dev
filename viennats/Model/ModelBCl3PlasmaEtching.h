#ifndef MODELBCl3PLASMAETCHING_H_
#define MODELBCl3PLASMAETCHING_H_

#include <stack>
#include "../Statistics.h"
#include <cassert>
#include <string>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

///Namespace containing all models.
namespace model {

///SF6/O2 plasma etching model

	class BCl3PlasmaEtching {

		static const double rho_HfO2;//=5e22; 		//in atoms/cm^3
		static const double rho_SiO2; 	//in atoms/cm^3

		static const double k_Sigma_HfO2;	//in cm^-2*s^-1
		static const double beta_Sigma_HfO2;	//in cm^-2*s^-1



		static const double MinEnergy;	//TODO //in eV
		static const double end_probability;

		static const double Gamma_O;
		static const double Gamma_F;

		static const double Eref_max;
		static const double Phi_inflect;
		static const double n_r;
		static const double n_l;

		static const double A_HfO2;
		double A_O=0;
		static const double A_p;
		static const double A_SiO2;

		static const double Eth_HfO2;
		static const double Eth_O;
		static const double Eth_p;
		static const double Eth_SiO2;

		double min_ion_energy;
		double delta_ion_energy;


		//static const double min_ion_energy=100;
		//static const double delta_ion_energy=40;

		static const double ion_exponent;
		static const double Phi_min;

		double StartDirection[3];

		double FluxIon;
		double FluxO=0.;
		double FluxBCl3;




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

		static const unsigned int NumberOfParticleTypes=3;

        static const bool OutputFluxes=false;
		static const bool SpatiallyEqualDistributedFlux=true;

		static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;

        static const bool CalculateNormalVectors=true;

        static const bool ReemissionIsMaterialDependent=true;
//        int stacked_mat;

		BCl3PlasmaEtching(const std::string & Parameters) {//, int stacked_material) {

		    using namespace boost::spirit::classic;
//		    stacked_mat=stacked_material;

		    double Accuracy;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                    		(str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("flux_ion")  >> '='  >> real_p[assign_a(FluxIon)]  >> ';') |
                            (str_p("flux_BCl3")  >> '='  >> real_p[assign_a(FluxBCl3)]  >> ';') |
                            (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';') |
                            (str_p("min_ion_energy")  >> '='  >>real_p[assign_a(min_ion_energy)]  >> ';') |
                            (str_p("delta_ion_energy")  >> '='  >>real_p[assign_a(delta_ion_energy)]  >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

		    unsigned int num_particles=static_cast<unsigned int>(Accuracy);

		    NumberOfParticleClusters[0]=(FluxIon>0.)?num_particles:0;	//Ion flux
		    NumberOfParticleClusters[1]=(FluxO>0.)?num_particles:0;		//Oxygen flux
		    NumberOfParticleClusters[2]=(FluxBCl3>0.)?num_particles:0;		//Fluorine flux


		}

		template <class VecType>
		void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
			//Velocity=-1.;
			//return;
			if (Material==2) {	//HfO2
				Velocity=(-(k_Sigma_HfO2*Coverages[1]/4.+Rates[0]+Rates[1]*Coverages[1])/rho_HfO2);	//in cm/s      //TODO
				if (Coverages[1]>1) std::cout << "Coverages[1]= " << Coverages[1] << ", Rates[0]= " << Rates[0] << ", Rates[1]= " << Rates[1] << "\n";
			} else if (Material==3) {
				Velocity=(-Rates[5]*Coverages[1]/rho_SiO2);										//in cm/s     //TODO
			} else Velocity=0.;
//assert(0);

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

			Coverages[1]=Rates[4]/(
				Rates[4]+
					(k_Sigma_HfO2+2*Rates[1])
						*
					(1.+Rates[3]/
							(beta_Sigma_HfO2+Rates[2]))
			);

			Coverages[0]=Rates[3]/(
				Rates[3]+
					(beta_Sigma_HfO2+Rates[2])
						*
					(1.+Rates[4]/
							(k_Sigma_HfO2+2*Rates[1]))
			);

		}

		template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {

			p.Type=ParticleType;
			switch(ParticleType) {
				case 0: //ION
					my::stat::CosineNDistributedRandomDirection(ion_exponent,StartDirection,p.Direction);
					p.Probability=1;
					p.Energy=min_ion_energy+delta_ion_energy*my::stat::RandomNumber();
					p.Flux=FluxIon;
					break;
				case 1:	//OXYGEN
					my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
					p.Energy=1.;
					p.Probability=1.;
					p.Flux=FluxO;
					break;
				case 2:	//FLUORINE
					my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
					p.Energy=1.;
					p.Probability=1.;
					p.Flux=FluxBCl3;	//5.5e18
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

			switch(p.Type) {
				case 0: {	//ION

					double f_HfO2_phi, f_O_phi, f_p_phi, f_SiO2_phi;

					const double dot=NormalVector[0]*p.Direction[0]+NormalVector[1]*p.Direction[1]+NormalVector[2]*p.Direction[2];
					const double cos=std::min(std::max(-dot,0.),1.);
					const double angle=std::acos(cos);

					if (cos>0.5) f_HfO2_phi=1.; else f_HfO2_phi=std::max(3.-6.*angle/my::math::Pi,0.);
					if (cos>0.5) f_O_phi=1.; else f_O_phi=std::max(3.-6.*angle/my::math::Pi,0.);
					//if (cos>0.5) f_p_phi=1.; else f_p_phi=3.-6.*angle/my::math::Pi;
					f_p_phi=1.;
					//if (cos>0.5) f_p_phi=1.; else f_p_phi=std::asin(cos)/std::asin(0.5);			//TODO
					//if (cos>0.5) f_SiO2_phi=1.; else f_SiO2_phi=3.-6.*angle/my::math::Pi;
					f_SiO2_phi=1.;

					const double sqrtE=std::sqrt(p.Energy);
					const double Y_p=A_p*std::max(sqrtE-std::sqrt(Eth_p),0.)*f_p_phi;
					const double Y_HfO2=A_HfO2*std::max(sqrtE-std::sqrt(Eth_HfO2),0.)*f_HfO2_phi;
					const double Y_O=A_O*std::max(sqrtE-std::sqrt(Eth_O),0.)*f_O_phi;
					const double Y_SiO2=A_SiO2*std::max(sqrtE-std::sqrt(Eth_SiO2),0.)*f_SiO2_phi;

					Rates[0]+=p.Flux*Y_p;
					Rates[1]+=p.Flux*Y_HfO2;
					Rates[2]+=p.Flux*Y_O;
					Rates[5]+=p.Flux*Y_SiO2;
					break;
				}
				case 1: {	//OXYGEN
					const double StickingProbability=std::max(Gamma_O*(1.-Coverages[0]-Coverages[1]),0.);
					Rates[3]+=p.Flux*StickingProbability;
					break;
				}
				case 2: {		//FLUORINE
					const double StickingProbability=std::max(Gamma_F*(1.-Coverages[0]-Coverages[1]),0.);
					Rates[4]+=p.Flux*StickingProbability;
					break;
				}
				default: {
					assert(0);
				}
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


			switch(p.Type) {
				case 0: {	//ION

					const double A=1./(1.+(n_l/n_r)*(my::stat::Pi1_2/Phi_inflect-1.));

					const double dot=(NormalVector[0]*p.Direction[0]+NormalVector[1]*p.Direction[1]+NormalVector[2]*p.Direction[2]);
					const double dot2=dot+dot;
					double ReflectedDirection[3];
					ReflectedDirection[0]=p.Direction[0]-NormalVector[0]*dot2;
					ReflectedDirection[1]=p.Direction[1]-NormalVector[1]*dot2;
					ReflectedDirection[2]=p.Direction[2]-NormalVector[2]*dot2;

					const double Phi=std::acos(std::max(std::min(-dot,1.),0.));

					double Eref_peak;
					if (Phi>=Phi_inflect) {
						Eref_peak=Eref_max*(1-(1-A)*std::pow((my::stat::Pi1_2-Phi)/(my::stat::Pi1_2-Phi_inflect),n_r));
					} else {
						Eref_peak=Eref_max*A*std::pow(Phi/Phi_inflect,n_l);
					}

					const double TempEnergy=Eref_peak*p.Energy;
					double NewEnergy;
					do {
						NewEnergy=	TempEnergy + (std::min((p.Energy-TempEnergy),TempEnergy)+0.05*p.Energy)*
											(1-2.*my::stat::RandomNumber())*
											std::sqrt(std::fabs(std::log(my::stat::RandomNumber())));
											//std::cos(my::math::Pi2*my::stat::RandomNumber())*std::sqrt(-2.*std::log(1.-my::stat::RandomNumber()));
					} while (NewEnergy>p.Energy || NewEnergy<0.);

					if (NewEnergy>=MinEnergy) {
						particle_stack.push(p);
						PT& p_new=particle_stack.top();
						p_new.Energy=NewEnergy;
						my::stat::CosAngleDistributedRandomDirection(my::stat::Pi1_2-std::min(Phi,Phi_min),ReflectedDirection,p_new.Direction);
					}

					break;
				}
				case 1: {	//OXYGEN

					const double StickingProbability=std::max(Gamma_O*(1.-Coverages[0]-Coverages[1]),0.);
					const double new_probability=p.Probability*(1-StickingProbability);
					if (new_probability>=end_probability) {
						particle_stack.push(p);
						PT& p_new=particle_stack.top();
						p_new.Probability=new_probability;
						my::stat::Cosine1DistributedRandomDirection(NormalVector,p_new.Direction);
					}
					break;
				}
				case 2: {		//FLUORINE

					const double StickingProbability=std::max(Gamma_F*(1.-Coverages[0]-Coverages[1]),0.);
					const double new_probability=p.Probability*(1-StickingProbability);
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

		}

	};
	double const BCl3PlasmaEtching::rho_HfO2 = 0.084365*6.022e23;//5e22;
	double const BCl3PlasmaEtching::rho_SiO2=2.6e22; 	//in atoms/cm^3

	double const BCl3PlasmaEtching::k_Sigma_HfO2=3e17;	//in cm^-2*s^-1
	double const BCl3PlasmaEtching::beta_Sigma_HfO2=4e13;	//in cm^-2*s^-1



	double const BCl3PlasmaEtching::MinEnergy=1.;	//TODO //in eV
	double const BCl3PlasmaEtching::end_probability=0.001;

	double const BCl3PlasmaEtching::Gamma_O=1.;
	double const BCl3PlasmaEtching::Gamma_F=0.00064;

	double const BCl3PlasmaEtching::Eref_max=1.;
	double const BCl3PlasmaEtching::Phi_inflect=1.55334303;
	double const BCl3PlasmaEtching::n_r=1.;
	double const BCl3PlasmaEtching::n_l=10.;

	double const BCl3PlasmaEtching::A_HfO2=7.;
	double const BCl3PlasmaEtching::A_p=0.0337;
	double const BCl3PlasmaEtching::A_SiO2=0.2;

	double const BCl3PlasmaEtching::Eth_HfO2=15.;
	double const BCl3PlasmaEtching::Eth_O=10.;
	double const BCl3PlasmaEtching::Eth_p=20.;
	double const BCl3PlasmaEtching::Eth_SiO2=8.;

	double const BCl3PlasmaEtching::ion_exponent=1000.;

	double const BCl3PlasmaEtching::Phi_min=1.3962634;

//	const double BCl3PlasmaEtching::Phi_min;

}

#endif /*MODELBCl3PLASMAETCHING_H_*/
