#ifndef DEF_PARAMETERS
#define DEF_PARAMETERS

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.
                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------
   Contact:         viennats@iue.tuwien.ac.at
   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <string>
#include <map>
#include <fstream>
#include <algorithm>

#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/classic_clear_actor.hpp>
#include <boost/spirit/include/classic_push_back_actor.hpp>
#include <boost/spirit/include/classic_increment_actor.hpp>
#include <boost/spirit/include/classic_decrement_actor.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_if.hpp>

#include "Partition/Partition.h"

#include "boundaries.h"

#include "message.h"

#include "parser_actors.h"

#include "sprng/sprng.h"

#define BOOST_SPIRIT_DEBUG


namespace par {



	class Parameters {
	public:

	    class ProcessParameterType {
        public:

	    	//AFM wire only!
	        double ProcessDistance;
	        double AFMStartPosition[3];
	        double AFMEndPosition[3];

	        double Masks;

            int AddLayer;                   //the number of level set layers which should be added to the geometry before the process is started

            double ProcessTime;             //the process time in seconds

            std::string ModelName;          //the name of the used process model

            std::string ModelParameters;    //the model parameters, which are passed to the model given by "ModelName"

            int IterationCycles;            //the number of iterations (process calculation) before advancing to the next time step,
                                            //may be necessary in case of an iterative problem, as is the case in the presence of coverages which themselves depend
                                            //on the flux

     //       int Iterations;

     //       int TotalTime;

            int StartIterationCycles;       //the number of iterations at the process simulation start, might be necessary to calculate initial values for coverages

            double MaxTimeStep;             //maximal allowed time step

            double smoothing_max_curvature;               //parameter for smoothing surface level set every time step
            double smoothing_min_curvature;               //parameter for smoothing surface level set every time step
            int smoothing_material_level;
            int smoothing_max_iterations;

            bool print_coverages;
            bool print_rates;
            bool print_velocities;
            bool print_materials;

            enum FiniteDifferenceSchemeType {ENGQUIST_OSHER_1ST_ORDER, ENGQUIST_OSHER_2ND_ORDER, LAX_FRIEDRICHS_1ST_ORDER, LAX_FRIEDRICHS_2ND_ORDER} ;
            FiniteDifferenceSchemeType FiniteDifferenceScheme;

            double LaxFriedrichsDissipationCoefficient;

            partition::DataStructureType partition_data_structure;
            partition::SplittingType partition_splitting_strategy;
            partition::SurfaceAreaHeuristicLambdaType partition_surface_area_heuristic_lambda;



            //Outputtimes
            std::vector<double> output_times;
            int output_times_periodicity;
            double output_times_period_length;
            bool initial_output;
            bool final_output;
            bool MaskLayer;
            bool GrowNewOxide;

            ProcessParameterType() {
                clear();
            }

            void clear() {

            	//AFM wire only!
    	        ProcessDistance=0;
    	        AFMStartPosition[0]=0;AFMStartPosition[1]=0;AFMStartPosition[1]=0;
    	        AFMEndPosition[0]=0;AFMEndPosition[1]=0;AFMEndPosition[1]=0;

                //default values
                AddLayer=0;
                Masks=0;
                ProcessTime=0;
                ModelName.clear();
                ModelParameters.clear();
                IterationCycles=0;
                StartIterationCycles=0;
     //           Iterations=0;
     //           TotalTime=0;
                MaxTimeStep=std::numeric_limits<double>::max();
                print_coverages=false;
                print_rates=false;
                print_velocities=false;
                print_materials=false;
                FiniteDifferenceScheme=ENGQUIST_OSHER_1ST_ORDER;

                output_times_periodicity=1;
                output_times_period_length=0;
                initial_output=false;
                final_output=false;
                MaskLayer=false;

                smoothing_max_curvature=std::numeric_limits<double>::max();
                smoothing_min_curvature=-std::numeric_limits<double>::max();
                smoothing_material_level=0;
                smoothing_max_iterations=100;

                partition_data_structure=partition::NEIGHBOR_LINKS_ARRAYS;
                partition_splitting_strategy=partition::SURFACE_AREA_HEURISTIC;
                partition_surface_area_heuristic_lambda=0.8;

            }
        };

		//parameters for random number generator
		static int RNG_Seed;
		static const int RNG_Type=SPRNG_MLFG;
		static const int RNG_Par=SPRNG_DEFAULT;

		//dimension
		int Dimensions;
//		int Orientation;

		//maximal time step for level set method (in grid units)
		double TimeStepRatio;

		//boundary conditions of simulation domain
		std::vector<bnc::boundary_conditions_type> boundary_conditions;

		//input geometry and output path
		std::string InputFile;
		std::vector<std::string> InputFiles;
		std::string OutputPath;

		bool surface_geometry;
		bool report_import_errors;
		//Outputtimes
		//std::vector<double> OutputTimes;


		//grid constant
		double GridDelta;

		//input transformations
		double InputScale;
		std::vector<int> MapMaterials;
		std::vector<int> InputTransformationDirections;
		std::vector<bool> InputTransformationSigns;
		std::vector<double> InputShift;
		bool change_input_parity;

		std::vector<double> DefaultDiskOrientation;
		std::vector<int> IgnoreMaterials;

		bool remove_bottom;
//		bool separate_materials;

		double snap_to_boundary_eps;

		int process_cycles;

		//maximum number of cores
		int OpenMP_threads;

		//domain extension in case of EXTENDED boundary conditions
        double DomainExtension;


        //disk properties
        double ReceptorRadius;
        double FurtherTrackingDistance;

        //output file format
        bool print_vtk;     //TODO
        bool print_dx;      //TODO


        bool print_coverages;
        bool print_rates;
        bool print_velocities;
        bool print_materials;

        bool PrintStatistics;

        int open_boundary_direction;
        bool is_open_boundary_negative;
        int max_extended_starting_position;

        int AddLayer;

        //TODO

		//static const bool PrintStatistics=true;

		Parameters(const std::string& FileName);

		std::string GetCompleteOutputFileName(const std::string& FileName) const {
			std::ostringstream oss;
			oss << par::Parameters::OutputPath<< FileName;
			return oss.str();
		}

		std::list<ProcessParameterType> ProcessParameters;

	};



	int Parameters::RNG_Seed;

	Parameters::ProcessParameterType tmp_process;
	bnc::boundary_conditions_type tmp_boundary_condition;

	class tmp_counter_check_class {
		int& z;
	 public:

		 tmp_counter_check_class(int& a) : z(a) {}

		 void clear() {
			 z=0;
		 }

		 bool operator()() const{
			return (z>0);
		 }

		 void operator++() {
			 ++z;
		 }

		 void operator--() {
			 --z;
		 }
	 };

	using namespace boost::spirit::classic;

	struct par_grammar : public grammar<par_grammar>
	{
	private:
		tmp_counter_check_class& tmp_counter_check;

	public:
		Parameters& par;

		par_grammar(Parameters& p, tmp_counter_check_class& t2) :  tmp_counter_check(t2), par(p)  {}

	    template <typename ScannerT> struct definition {

	    	 rule<ScannerT> 	rule_all,
								rule_comment,
								rule_input_file,
								rule_output_path,
								//rule_output_times,
								rule_surface_geometry,
								rule_report_import_errors,
								rule_CFL_condition,
								rule_grid_delta,
								rule_input_scale,
								rule_input_transformation,
								rule_input_shift,
								rule_default_disk_orientation,
								rule_ignore_materials,
								rule_change_input_parity,
								rule_random_seed,
								rule_dimensions,
								rule_OpenMP_threads,
								rule_Domain_Extension,
								rule_boundary_condition2A,
								rule_boundary_condition2B,
                                rule_boundary_conditions2,
                                rule_material_mapping,
                                rule_ReceptorRadius,
                                rule_FurtherTrackingDistance,
                                rule_print_vtk,
                                rule_print_dx,

                                rule_print_velocities,
                                rule_print_coverages,
                                rule_print_rates,
                                rule_print_materials,
                                rule_print_statistics,
                                rule_open_boundary,
                                rule_max_extended_starting_position,
                                rule_remove_bottom,
//                                rule_separate_materials,
                                rule_snap_to_boundary_eps,
                                rule_processing_cycles,
//                                rule_crystal_orientation,
                                rule_add_layer,

                                //AFM wire only!
                                rule_process_distance,
								rule_process_startPosition,
								rule_process_endPosition,

								rule_process_time,
								rule_process_top_mask,
								rule_process_grow_new_oxide,
								rule_process_smoothing_max_curvature,
								rule_process_smoothing_min_curvature,
								rule_process_smoothing_material_level,
								rule_process_smoothing_max_iterations,
								rule_process_partition_data_structure,
								rule_process_partition_splitting_strategy,
								rule_process_partition_surface_area_heuristic_lambda,
								rule_process_MaxTimeStep,
								rule_process,
								rule_process_addlayer,
								rule_process_IterationCycles,
								rule_process_StartIterationCycles,
								rule_process_model_name,
								rule_process_output_times_periodicity,
                                rule_process_output_times_period_length,
                                rule_process_initial_output,
								rule_process_final_output,
								rule_process_print_velocities,
								rule_process_print_coverages,
								rule_process_print_rates,
								rule_process_print_materials,
								rule_process_parameters,
								rule_processes,
								rule_process_finite_difference_scheme,
								rule_process_lax_friedrichs_dissipation_coefficient,
								rule_process_output_times;

	        definition(par_grammar const& self)	{

	            using namespace parser_actors;


	        	Parameters& p=self.par;

	            rule_input_file = (str_p("GeometryFiles") | str_p("geometry_files") | str_p("GeometryFile") | str_p("geometry_file")) >>'=' >> '\"' >> *((~ch_p('\"'))[push_back_a(p.InputFile)]) >> '\"' >> ';';
	            rule_output_path = (str_p("OutputPath") | str_p("output_path"))  >> '=' >> '\"' >> *((~ch_p('\"'))[push_back_a(p.OutputPath)]) >> '\"'  >> ';';
	            rule_surface_geometry = (str_p("surface_geometry")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.surface_geometry)]) >> ';');
	            rule_report_import_errors = (str_p("report_import_errors")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.report_import_errors)]) >> ';');
	            rule_CFL_condition   = (str_p("CFL-Condition") | str_p("cfl_condition"))  >> '='  >> real_p[assign_a(p.TimeStepRatio)]  >> ';';
				rule_grid_delta   = (str_p("GridDelta") | str_p("grid_delta"))  >> '='  >> real_p[assign_a(p.GridDelta)]  >> ';';
				rule_input_scale   = (str_p("InputScale") | str_p("input_scale"))  >> '='  >> real_p[assign_a(p.InputScale)]  >> ';';
				rule_input_transformation = (str_p("input_transformation") >> '=' >> '{' >> (('\"' >> (((ch_p('+') | '-') >> (ch_p('x') | 'y' | 'z'))[assign_input_transformation(p.InputTransformationDirections, p.InputTransformationSigns)]) >> '\"') % ',') >> '}' >> ';');
				rule_input_shift = (str_p("input_shift") >> '=' >> '{' >> (real_p[push_back_a(p.InputShift)] % ',') >> '}' >> ';');
				rule_default_disk_orientation = (str_p("default_disk_orientation") >> '=' >> '{' >> (real_p[push_back_a(p.DefaultDiskOrientation)] % ',') >> '}' >> ';');
				rule_ignore_materials = (str_p("ignore_materials") >> '=' >> '{' >> (int_p[push_back_a(p.IgnoreMaterials)] % ',') >> '}' >> ';');
				rule_change_input_parity=(str_p("change_input_parity")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.change_input_parity)]) >> ';');
				rule_random_seed   = (str_p("RandomSeed") | str_p("random_seed"))  >> '='  >> int_p[assign_a(p.RNG_Seed)]  >> ';';
				rule_dimensions   = (str_p("Dimensions") | str_p("num_dimensions"))  >> '='  >> int_p[assign_a(p.Dimensions)]  >> ';';
				rule_OpenMP_threads = (str_p("OpenMP_threads") | str_p("omp_threads")) >> '=' >> int_p[assign_a(p.OpenMP_threads)] >> ';';
				rule_Domain_Extension = (str_p("Domain_Extension") | str_p("domain_extension")) >> '=' >> real_p[assign_a(p.DomainExtension)] >> ';';
				rule_ReceptorRadius = (str_p("ReceptorRadius") | str_p("receptor_radius")) >> '=' >> real_p[assign_a(p.ReceptorRadius)] >> ';';
				rule_FurtherTrackingDistance = (str_p("FurtherTrackingDistance") | str_p("further_tracking_distance")) >> '=' >> real_p[assign_a(p.FurtherTrackingDistance)] >> ';';
				rule_print_vtk = (str_p("print_vtk")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.print_vtk)]) >> ';');
				rule_print_dx = (str_p("print_dx")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.print_dx)]) >> ';');
				rule_print_velocities = (str_p("print_velocities")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.print_velocities)]) >> ';');
                rule_print_coverages = (str_p("print_coverages")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.print_coverages)]) >> ';');
                rule_print_rates = (str_p("print_rates")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.print_rates)]) >> ';');
                rule_print_materials = (str_p("print_materials")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.print_materials)]) >> ';');
                rule_print_statistics = (str_p("print_statistics")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.PrintStatistics)]) >> ';');
                rule_max_extended_starting_position = (str_p("max_extended_starting_position") >> '='  >>  int_p[assign_a(p.max_extended_starting_position)] >>  ';');
                rule_open_boundary = (str_p("open_boundary") >> '=' >> '\"' >> (((ch_p('+') | '-') >> (ch_p('x') | 'y' | 'z'))[assign_dir(p.open_boundary_direction, p.is_open_boundary_negative)]) >> '\"' >> ';');
                rule_remove_bottom = str_p("remove_bottom")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.remove_bottom)]) >> ';';
//                rule_separate_materials = str_p("separate_materials")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(p.separate_materials)]) >> ';';
                rule_snap_to_boundary_eps  = str_p("snap_to_boundary_eps")  >> '='  >> real_p[assign_a(p.snap_to_boundary_eps)]  >> ';';
                rule_processing_cycles = str_p("process_cycles")  >> '='  >> int_p[assign_a(p.process_cycles)]  >> ';';

				rule_boundary_condition2A = (   (str_p("REFLECTIVE")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.min, bnc::REFLECTIVE_BOUNDARY)]) |
                                                (str_p("SYMMETRIC")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.min, bnc::REFLECTIVE_BOUNDARY)]) |              //synonym for reflective
                                                (str_p("INFINITE")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.min, bnc::INFINITE_BOUNDARY)]) |
                                                (str_p("PERIODIC")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.min, bnc::PERIODIC_BOUNDARY)]) |
                                                (str_p("EXTENDED")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.min, bnc::EXTENDED_BOUNDARY)])
                                            );

				rule_boundary_condition2B = (   (str_p("REFLECTIVE")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.max, bnc::REFLECTIVE_BOUNDARY)]) |
                                                (str_p("SYMMETRIC")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.max, bnc::REFLECTIVE_BOUNDARY)]) |              //synonym for reflective
                                                (str_p("INFINITE")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.max, bnc::INFINITE_BOUNDARY)]) |
                                                (str_p("PERIODIC")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.max, bnc::PERIODIC_BOUNDARY)]) |
                                                (str_p("EXTENDED")[assign_enum<bnc::boundary_condition_type>(tmp_boundary_condition.max, bnc::EXTENDED_BOUNDARY)])
                                            );


                rule_boundary_conditions2  = str_p("boundary_conditions")  >> '='  >>  '{' >>
                            list_p( (ch_p('{') >> rule_boundary_condition2A >> ',' >>  rule_boundary_condition2B >> ch_p('}'))[push_back_a(p.boundary_conditions, tmp_boundary_condition)], ',')
                            >> '}'  >> ';';
                rule_material_mapping  = str_p("material_mapping")  >> '='  >>  '{' >> (int_p[push_back_a(p.MapMaterials)] % ',') >> '}'  >> ';';
                rule_add_layer = str_p("add_layer")  >> '='  >>int_p[assign_a(p.AddLayer)]  >> ';';



				rule_process_time = str_p("process_time")  >> '='  >> real_p[assign_a(tmp_process.ProcessTime)]  >> ';';

				//AFM wire only
				rule_process_distance = str_p("process_distance")  >> '='  >> real_p[assign_a(tmp_process.ProcessDistance)]  >> ';';
				rule_process_startPosition = str_p("start_position")  >> '='  >> '{' >> real_p[assign_a(tmp_process.AFMStartPosition[0])]  >> "," >> real_p[assign_a(tmp_process.AFMStartPosition[1])] >> "," >> real_p[assign_a(tmp_process.AFMStartPosition[2])] >> '}' >> ';';
				rule_process_endPosition = str_p("end_position")  >> '='  >> '{' >> real_p[assign_a(tmp_process.AFMEndPosition[0])]  >> "," >> real_p[assign_a(tmp_process.AFMEndPosition[1])] >> "," >> real_p[assign_a(tmp_process.AFMEndPosition[2])] >> '}' >> ';';

				rule_process_top_mask = (str_p("mask_layer") >> '=' >> ((str_p("true") | str_p("false"))[assign_bool(tmp_process.MaskLayer)]) >> ';');
				rule_process_grow_new_oxide = (str_p("grow_new_oxide") >> '=' >> ((str_p("true") | str_p("false"))[assign_bool(tmp_process.GrowNewOxide)]) >> ';');
				rule_process_smoothing_max_curvature = str_p("smoothing_max_curvature")  >> '='  >>real_p[assign_a(tmp_process.smoothing_max_curvature)]  >> ';';
				rule_process_smoothing_min_curvature = str_p("smoothing_min_curvature")  >> '='  >>real_p[assign_a(tmp_process.smoothing_min_curvature)]  >> ';';
				rule_process_smoothing_material_level = str_p("smoothing_material_level")  >> '='  >> int_p[assign_a(tmp_process.smoothing_material_level)]  >> ';';
				rule_process_smoothing_max_iterations = str_p("smoothing_max_iterations")  >> '='  >> int_p[assign_a(tmp_process.smoothing_max_iterations)]  >> ';';
				rule_process_addlayer = str_p("add_layer")  >> '='  >>int_p[assign_a(tmp_process.AddLayer)]  >> ';';
				rule_process_output_times_periodicity=str_p("output_times_periodicity") >> '=' >> int_p[assign_a(tmp_process.output_times_periodicity)]  >> ';';
                rule_process_output_times_period_length=str_p("output_times_period_length") >> '=' >> real_p[assign_a(tmp_process.output_times_period_length)]  >> ';';
                rule_process_initial_output = (str_p("initial_output")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(tmp_process.initial_output)]) >> ';');
				rule_process_final_output = (str_p("final_output")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(tmp_process.final_output)]) >> ';');
                rule_process_print_velocities = (str_p("print_velocities")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(tmp_process.print_velocities)]) >> ';');
				rule_process_print_coverages = (str_p("print_coverages")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(tmp_process.print_coverages)]) >> ';');
				rule_process_print_rates = (str_p("print_rates")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(tmp_process.print_rates)]) >> ';');
				rule_process_print_materials = (str_p("print_materials")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(tmp_process.print_materials)]) >> ';');

				rule_process_MaxTimeStep = (str_p("MaxTimeStep") | str_p("max_time_step"))  >> '='  >> real_p[assign_a(tmp_process.MaxTimeStep)]  >> ';';
				rule_process_IterationCycles = (str_p("IterationCycles") | str_p("iteration_cycles"))  >> '='  >>int_p[assign_a(tmp_process.IterationCycles)]  >> ';';
				rule_process_StartIterationCycles = (str_p("StartIterationCycles") | str_p("start_iteration_cycles"))  >> '='  >>int_p[assign_a(tmp_process.StartIterationCycles)]  >> ';';
//				rule_process_Iterations = (str_p("Iterations") | str_p("iterations")) >> '=' >> int_p[assign_a(tmp_process.Iterations)] >> ';';
//				rule_process_TotalTime = (str_p("total_time") | str_p("Total_time")) >> '=' >> int_p[assign_a(tmp_process.TotalTime)] >> ';';
				rule_process_model_name = str_p("model_name") >> '=' >>  '\"' >> *((~ch_p('\"'))[push_back_a(tmp_process.ModelName)]) >> '\"' >> ';';
				rule_process_parameters=str_p("parameters") >>'=' 	>> ch_p('{')[clear_a(self.tmp_counter_check)]
				                                                    >> *(	(ch_p('{')[push_back_a(tmp_process.ModelParameters)][increment_a(self.tmp_counter_check)]) |
				                                                    		(if_p(self.tmp_counter_check)[ch_p('}')[push_back_a(tmp_process.ModelParameters)][decrement_a(self.tmp_counter_check)]].else_p[nothing_p]) |
				                                                    		((~ch_p('}'))[push_back_a(tmp_process.ModelParameters)])

				                                                    	)
				                                                	>> '}' >> ';';

				rule_process_output_times  = str_p("output_times")  >> '='  >>  '{' >> (real_p[push_back_a(tmp_process.output_times)] % ',') >> '}'  >> ';';


				typedef Parameters::ProcessParameterType::FiniteDifferenceSchemeType FiniteDifferenceSchemeType;

				rule_process_finite_difference_scheme=str_p("finite_difference_scheme") >> '=' >> (
				        str_p("ENGQUIST_OSHER_1ST_ORDER")[assign_enum<FiniteDifferenceSchemeType>(tmp_process.FiniteDifferenceScheme,tmp_process.ENGQUIST_OSHER_1ST_ORDER)] |
				        str_p("ENGQUIST_OSHER_2ND_ORDER")[assign_enum<FiniteDifferenceSchemeType>(tmp_process.FiniteDifferenceScheme,tmp_process.ENGQUIST_OSHER_2ND_ORDER)] |
				        str_p("LAX_FRIEDRICHS_1ST_ORDER")[assign_enum<FiniteDifferenceSchemeType>(tmp_process.FiniteDifferenceScheme,tmp_process.LAX_FRIEDRICHS_1ST_ORDER)] |
				        str_p("LAX_FRIEDRICHS_2ND_ORDER")[assign_enum<FiniteDifferenceSchemeType>(tmp_process.FiniteDifferenceScheme,tmp_process.LAX_FRIEDRICHS_2ND_ORDER)]
                                                                                                    ) >> ';';

				rule_process_lax_friedrichs_dissipation_coefficient=str_p("dissipation_coefficient") >> '=' >> real_p[assign_a(tmp_process.LaxFriedrichsDissipationCoefficient)] >> ';';

				rule_process_partition_data_structure=str_p("partition_data_structure") >> '=' >> (
                        str_p("NEIGHBOR_LINKS_ARRAYS")[assign_enum<partition::DataStructureType>(tmp_process.partition_data_structure,partition::NEIGHBOR_LINKS_ARRAYS)] |
                        str_p("FULL_GRID")[assign_enum<partition::DataStructureType>(tmp_process.partition_data_structure,partition::FULL_GRID)] |
                        str_p("UP_DOWN_LINKED_TREE")[assign_enum<partition::DataStructureType>(tmp_process.partition_data_structure,partition::UP_DOWN_LINKED_TREE)]
                                                                                                    ) >> ';';

				rule_process_partition_splitting_strategy=str_p("partition_splitting_strategy") >> '=' >> (
                        str_p("SPATIAL_MEDIAN")[assign_enum<partition::SplittingType>(tmp_process.partition_splitting_strategy,partition::SPATIAL_MEDIAN)] |
                        str_p("OBJECT_MEDIAN")[assign_enum<partition::SplittingType>(tmp_process.partition_splitting_strategy,partition::OBJECT_MEDIAN)] |
                        str_p("SURFACE_AREA_HEURISTIC")[assign_enum<partition::SplittingType>(tmp_process.partition_splitting_strategy,partition::SURFACE_AREA_HEURISTIC)]
                                                                                                    ) >> ';';

				rule_process_partition_surface_area_heuristic_lambda = str_p("partition_surface_area_heuristic_lambda")  >> '='  >>real_p[assign_a(tmp_process.partition_surface_area_heuristic_lambda)]  >> ';';


				rule_process = ch_p('{') >> *(
                                                rule_process_time |
                                                rule_process_distance |
                                                rule_process_startPosition |
                                                rule_process_endPosition |
                                                rule_process_top_mask |
                                                rule_process_grow_new_oxide |
                                                rule_process_smoothing_max_curvature |
                                                rule_process_smoothing_min_curvature |
                                                rule_process_smoothing_material_level |
                                                rule_process_smoothing_max_iterations |
                                                rule_process_partition_data_structure |
                                                rule_process_partition_splitting_strategy |
                                                rule_process_partition_surface_area_heuristic_lambda |
                                                rule_process_IterationCycles |
												rule_process_StartIterationCycles |
									//			rule_process_Iterations |
									//			rule_process_TotalTime |
												rule_process_addlayer |
												rule_process_MaxTimeStep |
												rule_process_model_name |
												rule_process_initial_output |
												rule_process_final_output |
												rule_process_output_times_periodicity       |
												rule_process_output_times_period_length     |
												rule_process_output_times |
                                                rule_process_parameters |
                                                rule_process_print_velocities |
                                                rule_process_print_coverages |
                                                rule_process_print_rates |
                                                rule_process_print_materials |
                                                rule_process_finite_difference_scheme |
                                                rule_process_lax_friedrichs_dissipation_coefficient
								) >> '}';


				rule_processes = (str_p("Processes") | str_p("processes"))  >> '='  >> '{' >> list_p((rule_process[push_back_a(p.ProcessParameters, tmp_process)][clear_a(tmp_process)]) , ',') >> '}'  >> ';';

				rule_all=*(
                                rule_comment                        |
                                rule_input_file                     |
                                rule_output_path                    |
								rule_surface_geometry				|
								rule_report_import_errors			|
                                rule_CFL_condition                  |
                                rule_grid_delta                     |
                                rule_input_scale                    |
                                rule_input_transformation           |
                                rule_input_shift					|
                                rule_default_disk_orientation		|
                                rule_ignore_materials				|
                                rule_change_input_parity            |
                                rule_random_seed                    |
                                rule_dimensions                     |
                                rule_OpenMP_threads                 |
                                rule_Domain_Extension               |
                                rule_boundary_conditions2           |
                                rule_material_mapping				|
                                rule_processes                      |
                                rule_ReceptorRadius                 |
                                rule_FurtherTrackingDistance        |
                                rule_print_vtk						|
                                rule_print_dx						|
                                rule_print_velocities               |
                                rule_print_coverages                |
                                rule_print_rates                    |
                                rule_print_materials                |
                                rule_print_statistics               |
                                rule_max_extended_starting_position |
                                rule_open_boundary                  |
                                rule_snap_to_boundary_eps           |
                                rule_remove_bottom                  |
                                rule_processing_cycles				|
//                                rule_crystal_orientation			|
                                rule_add_layer
                                //anychar_p
                            )>>end_p;

            }


	        rule<ScannerT> const&
	        start() const { return rule_all; }
	    };
	};


	namespace {

        void ReadFile( std::string FileName, std::string & str);

        class include_actor {
            std::string & str;
        public:
            include_actor(std::string & s) : str(s) {}

            template <class iter>
            void operator()(const iter&  a, const iter& b ) const {
                ReadFile(std::string(a,b), str);
            }
        };

        void ReadFile( std::string FileName, std::string & str) {

             //read parameters
            std::ifstream ifs(FileName.c_str());

            file_iterator<char> first(FileName.c_str());
            file_iterator<char> last(first.make_end());

             parse(
                    first,
                    last,

                    *(
                            str_p("#include") >> (*blank_p) >> '\"' >> (*blank_p) >> (*(~ch_p('\"')))[include_actor(str)] >> '\"' >>  (*blank_p) >> eol_p |
                            space_p |
                            anychar_p[push_back_a(str)]
                     ),
                     comment_p("//") | comment_p("/*", "*/")

            ).full;

            ifs.close();

        }
}

	inline Parameters::Parameters(const std::string& FileName)  {

	    //default values
	    OpenMP_threads=0;
        DomainExtension=0;

        FurtherTrackingDistance=3.;
        ReceptorRadius=0.8;

        InputScale=1;
        change_input_parity=false;

        print_dx=false;
        print_vtk=true;

        print_coverages=false;
        print_rates=false;
        print_velocities=false;
        print_materials=false;

        //PrintStatistics=true;
        PrintStatistics=false;

        open_boundary_direction=1;
        is_open_boundary_negative=false;

        remove_bottom=true;
        snap_to_boundary_eps=1e-6;

        process_cycles=1;


        max_extended_starting_position=1000;
        AddLayer=0;

        surface_geometry=false;
        report_import_errors=true;

	    std::string str;
	    ReadFile(FileName, str);

		int tmp_counter;
	    tmp_counter_check_class tmp_counter_check(tmp_counter);

	    //parse parameter file
		par_grammar my_parser(*this, tmp_counter_check);

		bool b = parse(
				str.begin(),
				str.end(),
				my_parser
				).full;

		if (!b) msg::print_error("Failed reading parameter file!");

//----- Section to parse through the file name(s) and determine input format ---
		std::size_t start = str.find("ile=\"");
		if (start>str.length()) {
			start = str.find("iles=\"");
			start++;
		}
		start+=4;
		std::size_t end = str.find("\";",start);
		std::string FNstr=str.substr(start+1,end-start-1);
	    bool done=false;
	    while (!done) {
	        std::size_t found = FNstr.find(",");
	        if (found!=std::string::npos) {
	          InputFiles.push_back (FNstr.substr(0,found));
	          FNstr=FNstr.substr(found+1);
	        } else {
	        	InputFiles.push_back(FNstr);
	        	done=true;
	        }
	    }

//----- Generate the output folder ---

	    struct stat st;
	    if(stat(OutputPath.c_str(),&st) == -1)
        	mkdir(OutputPath.c_str(),mode_t(0777));

		//test boundary condtions for periodicity, and check conformity
		for (unsigned int h=0;h<boundary_conditions.size();h++) {

		    if (boundary_conditions[h].min==bnc::PERIODIC_BOUNDARY) assert(boundary_conditions[h].max==bnc::PERIODIC_BOUNDARY);
		    if (boundary_conditions[h].max==bnc::PERIODIC_BOUNDARY) assert(boundary_conditions[h].min==bnc::PERIODIC_BOUNDARY);
		    //if (boundary_conditions[h].min==INFINITE_BOUNDARY) assert(h==static_cast<unsigned int>(open_boundary_direction));
		    if (h==static_cast<unsigned int>(open_boundary_direction)) assert(boundary_conditions[h].min==bnc::INFINITE_BOUNDARY);
		    //if (boundary_conditions[h].max==INFINITE_BOUNDARY) assert(h==static_cast<unsigned int>(open_boundary_direction));
		    if (h==static_cast<unsigned int>(open_boundary_direction)) assert(boundary_conditions[h].max==bnc::INFINITE_BOUNDARY);
		    //if (boundary_conditions[h].min==INFINITE_BOUNDARY) assert(boundary_conditions[h].max==INFINITE_BOUNDARY);
            //if (boundary_conditions[h].max==INFINITE_BOUNDARY) assert(boundary_conditions[h].min==INFINITE_BOUNDARY);

		}

		//add OutputTimes
		for(std::list<ProcessParameterType>::iterator pIter=ProcessParameters.begin();pIter!=ProcessParameters.end();++pIter) {
		    if (pIter->output_times_period_length>0 && pIter->output_times_periodicity>=0) {

		        unsigned int tmp=pIter->output_times.size();
		        pIter->output_times.resize(pIter->output_times_periodicity*tmp);

		        for (int h=1;h<pIter->output_times_periodicity;++h) {
		            for (unsigned int k=0;k<tmp;k++) {
		                pIter->output_times[tmp*h+k]=pIter->output_times[k]+h*pIter->output_times_period_length;
		            }
		        }
		    }


            if (pIter->initial_output) pIter->output_times.push_back(0);
//            int this_size_here = pIter->output_times.size();
//            msg::print_message("output_times.size() before: ", this_size_here);

            if (pIter->final_output) pIter->output_times.push_back(pIter->ProcessTime);

            std::sort(pIter->output_times.begin(),pIter->output_times.end());

            std::vector<double>::iterator it=std::unique(pIter->output_times.begin(),pIter->output_times.end());
            pIter->output_times.resize(it-pIter->output_times.begin());

//            int this_size_here2 = pIter->output_times.size();
//            msg::print_message("output_times.size() after: ", this_size_here2);

		}

		//process cycles
		{
		    std::list<ProcessParameterType> tmp;
		    for (int i=0;i<process_cycles;++i) {
		        copy(ProcessParameters.begin(), ProcessParameters.end(),back_inserter(tmp));
		    }
		    swap(ProcessParameters, tmp);
		}

	}

}
#endif //DEF_PARAMETERS
