#define BOOST_SPIRIT_QI_DEBUG

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/variant/recursive_variant.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <complex>
#include <list>
#include <vector>

#include "boundaries.h"
#include "Partition/Partition.h"

/*

TO ADD NEW PARAMETER TO THIS PARSER:
- Declare new variable in Parameters struct
- Add variable to BOOST_FUSION_ADAPT_STRUCT
- Declare new rule at the end of par_grammar
- Define the rule in par_grammar constructor
- Add new rule to "definition" rule in par_grammar constructor
- if it is a required parameter: add variable to validate() in Parameters struct

*/


namespace proc{
        enum FiniteDifferenceSchemeType {ENGQUIST_OSHER_1ST_ORDER, ENGQUIST_OSHER_2ND_ORDER, LAX_FRIEDRICHS_1ST_ORDER, LAX_FRIEDRICHS_2ND_ORDER};
}


namespace qi = boost::spirit::qi;
//namespace phoenix = boost::phoenix;

namespace client{
  // function to find line the error occurred
  template<typename Iter>
  unsigned find_error_line(Iter first_iter,
    Iter error_iter) {
      unsigned line_count = 1;
      while(first_iter != error_iter){
        ++first_iter;
        if(*first_iter == '\n') ++line_count;
      }
      return line_count;
    }


  template<typename Iter>
  void report_input_error(Iter first_iter, Iter last_iter,
    Iter error_iter, std::string error_msg) {
      std::string first(first_iter, error_iter);
      std::string last(error_iter, last_iter);
      auto first_pos = first.rfind('\n');
      auto last_pos = last.find('\n');
      auto error_line = ((first_pos == std::string::npos) ? first
      : std::string(first, first_pos + 1))
      + std::string(last, 0, last_pos);
      auto error_pos = (error_iter - first_iter) + 1;
      if (first_pos != std::string::npos) {
          error_pos -= (first_pos + 1);
      }
      std::cerr << error_msg << std::endl;
      std::cerr << "In line " << find_error_line(first_iter, error_iter)  << ", column " << error_pos << std::endl
      << error_line << std::endl
      << std::setw(error_pos) << '^' << "--- here"
      << std::endl;
  }

}

// lazy function for error reporting
struct ReportError {
    // the result type must be explicit for Phoenix
    template<typename, typename, typename, typename>
    struct result { typedef void type; };

    // function for spirit to report errors properly
    // contract the string to the surrounding new-line characters
    template<typename Iter>
    void operator()(Iter first_iter, Iter last_iter,
        Iter error_iter, const qi::info& what) const {

            std::string first(first_iter, error_iter);
            std::string last(error_iter, last_iter);
            auto first_pos = first.rfind('\n');
            auto last_pos = last.find('\n');
            auto error_line = ((first_pos == std::string::npos) ? first
            : std::string(first, first_pos + 1))
            + std::string(last, 0, last_pos);
            auto error_pos = (error_iter - first_iter) + 1;
            if (first_pos != std::string::npos) {
                error_pos -= (first_pos + 1);
            }
            std::cerr << "INPUT ERROR: Parsing error near " << what << std::endl;
            std::cerr << "In line " << client::find_error_line(first_iter, error_iter)  << ", column " << error_pos << std::endl
            << error_line << std::endl
            << std::setw(error_pos) << '^' << "--- here"
            << std::endl;
    }
};

    const boost::phoenix::function<ReportError> report_error = ReportError();

    //symbol table for true/false
    struct boolean_ : qi::symbols<char, bool>{
        boolean_(){
            add
            ("true", true)
            ("false", false);
        }
    }boolean_symbols;

    //symbol table for true/false
    struct boundary_condition_ : qi::symbols<char, bnc::boundary_condition_type>{
        boundary_condition_(){
            add
            ("INFINITE", bnc::INFINITE_BOUNDARY)
            ("REFLECTIVE", bnc::REFLECTIVE_BOUNDARY)
            ("PERIODIC", bnc::PERIODIC_BOUNDARY)
            ("EXTENDED", bnc::EXTENDED_BOUNDARY);
        }
    }boundary_condition_symbols;

    //symbol table for cartesian directions
    struct cartesian_dir_ : qi::symbols<char, int>{
        cartesian_dir_(){
            add
            ("-z", -3)
            ("-y", -2)
            ("-x", -1)
            ("+x", 1)
            ("x", 1)
            ("+y", 2)
            ("y", 2)
            ("+z", 3)
            ("z", 3);
        }
    }cartesian_dir_symbols;

    //symbols for finite difference scheme
    struct finite_difference_ : qi::symbols<char, proc::FiniteDifferenceSchemeType>{
        finite_difference_(){
            add
            ("ENGQUIST_OSHER_1ST_ORDER", proc::ENGQUIST_OSHER_1ST_ORDER)
            ("ENGQUIST_OSHER_2ND_ORDER", proc::ENGQUIST_OSHER_2ND_ORDER)
            ("LAX_FRIEDRICHS_1ST_ORDER", proc::LAX_FRIEDRICHS_1ST_ORDER)
            ("LAX_FRIEDRICHS_2ND_ORDER", proc::LAX_FRIEDRICHS_2ND_ORDER);
        }
    }finite_difference_symbols;

    //symbols for partition_data_structure
    struct partition_data_structure_ : qi::symbols<char, partition::DataStructureType>{
        partition_data_structure_(){
            add
            ("NEIGHBOR_LINKS_ARRAYS", partition::NEIGHBOR_LINKS_ARRAYS)
            ("FULL_GRID", partition::FULL_GRID)
            ("UP_DOWN_LINKED_TREE", partition::UP_DOWN_LINKED_TREE);
        }
    }partition_data_structure_symbols;

    //symbols for partition splitting strategy
    struct partition_splitting_strategy_ : qi::symbols<char, partition::SplittingType>{
        partition_splitting_strategy_(){
            add
            ("SPATIAL_MEDIAN", partition::SPATIAL_MEDIAN)
            ("OBJECT_MEDIAN", partition::OBJECT_MEDIAN)
            ("SURFACE_AREA_HEURISTIC", partition::SURFACE_AREA_HEURISTIC);
        }
    }partition_splitting_strategy_symbols;


    namespace client{
        //define shortcuts for namespace
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        //Define struct to keep input data
        struct Parameters{
            Parameters(std::string);

            struct ProcessParameterType{
                double ProcessDistance;
                double AFMStartPosition[3];
                double AFMEndPosition[3];
                int AddLayer;                   //the number of level set layers added to the geometry before the process is started
                double ProcessTime;             //the process time in seconds
                unsigned int ALDStep;
                std::string ModelName;          //the name of the used process model
                std::string ModelParameters;    //the model parameters, which are passed to the model given by "ModelName"
                int IterationCycles;            //the number of iterations (process calculation) before advancing to the next time step,
                int StartIterationCycles;       //the number of iterations at the process simulation start, might be necessary to calculate initial values for coverages
                double MaxTimeStep;             //maximal allowed time step
                double smoothing_max_curvature; //parameter for smoothing surface level set every time step
                double smoothing_min_curvature; //parameter for smoothing surface level set every time step
                int smoothing_material_level;
                int smoothing_max_iterations;
                bool print_coverages;
                bool print_rates;
                bool print_velocities;
                bool print_materials;
                proc::FiniteDifferenceSchemeType FiniteDifferenceScheme;
                double LaxFriedrichsDissipationCoefficient;
                partition::DataStructureType partition_data_structure;
                partition::SplittingType partition_splitting_strategy;
                partition::SurfaceAreaHeuristicLambdaType partition_surface_area_heuristic_lambda;

                //Outputtimes
                std::vector<double> output_times;
                std::vector<double> output_volume;
                int output_times_periodicity;
                double output_times_period_length;
                bool initial_output;
                bool final_output;
                bool MaskLayer;
                bool GrowNewOxide;

                ProcessParameterType(){
                    AddLayer=0;
                    ProcessTime=0;
                    ALDStep=1;
                    ModelName.clear();
                    ModelParameters.clear();
                    IterationCycles=0;
                    StartIterationCycles=0;
                    MaxTimeStep=std::numeric_limits<double>::max();
                    print_coverages=false;
                    print_rates=false;
                    print_velocities=false;
                    print_materials=false;
                    FiniteDifferenceScheme=proc::ENGQUIST_OSHER_1ST_ORDER;

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

            // ADD CUSTOM PARAMETERS HERE
            std::vector<std::string> geometry_files;
            std::string output_path;
            bool surface_geometry;
            bool report_import_errors;
            double cfl_condition;
            double input_scale;
            double grid_delta;
            std::vector<int> input_transformation;
            std::vector<double> input_shift;
            std::vector<double> default_disc_orientation;
            std::vector<int> ignore_materials;
            bool change_input_parity;
            int random_seed;
            int num_dimensions;
            int omp_threads;
            double domain_extention;
            double receptor_radius;
            double further_tracking_distance;
            int bits_per_distance;
            bool output_volume_extract_single_materials;
            bool print_vtk;
            bool print_vtp;
            bool print_dx;
            bool print_lvst;
            bool print_volume_tetra;
            bool print_volume_hull;
            bool print_velocities;
            bool print_coverages;
            bool print_rates;
            bool print_materials;
            bool print_statistics;
            int max_extended_starting_position;
            int open_boundary;
            bool remove_bottom;
            double snap_to_boundary_eps;
            int process_cycles;
            std::vector<int> material_mapping;
            int add_layer;
            std::vector<bnc::boundary_condition_type> boundary_condition;
            std::list<ProcessParameterType> process_parameters;

            //Extra parameters deduced from parsed ones
            std::vector<bool> input_transformation_signs;
            bool open_boundary_negative;
            std::vector<bnc::boundary_conditions_type> boundary_conditions;

        private:
            //error helper function
            std::string make_error(std::string name){
                std::ostringstream oss;
                oss << "INPUT ERROR: Required parameter \"" << name << "\" not set!" << std::endl;
                return oss.str();
            }


            //ADD "REQUIRED" PARAMETERS TO BE ADDED HERE
            bool validate(){
                std::string error;
                if(geometry_files.empty()) error+=make_error("geometry_file");
                //if(output_path.empty()) error+=make_error("output_path");
                if(cfl_condition==-1.) error+=make_error("cfl_condition");
                if(grid_delta==-1.) error+=make_error("grid_delta");
                if(default_disc_orientation.empty()) error+=make_error("default_disc_orientation");
                if(num_dimensions==-1) error+=make_error("num_dimensions");
                if(boundary_condition.empty()) error+=make_error("boundary_conditions");

                //check if there are errors
                if(error.compare("")){
                    std::cerr << error << std::endl;
                    return false;
                }else return true;
            }

            Parameters();   //make private so it cannot be called by mistake
        };
    }

    // ADD STRUCT MEMBERS THAT SHOULD BE PARSED HERE SO THEY CAN BE ACCESSED BY THE PARSER
    BOOST_FUSION_ADAPT_STRUCT(client::Parameters,
        (std::vector<std::string>, geometry_files)
        (std::string, output_path)
        (bool, surface_geometry)
        (bool, report_import_errors)
        (double, cfl_condition)
        (double, input_scale)
        (double, grid_delta)
        (std::vector<int>, input_transformation)
        (std::vector<double>, input_shift)
        (std::vector<double>, default_disc_orientation)
        (std::vector<int>, ignore_materials)
        (bool, change_input_parity)
        (int, random_seed)
        (int, num_dimensions)
        (int, omp_threads)
        (double, domain_extention)
        (double, receptor_radius)
        (double, further_tracking_distance)
        (int, bits_per_distance)
        (bool, output_volume_extract_single_materials)
        (bool, print_vtk)
        (bool, print_vtp)
        (bool, print_dx)
        (bool, print_lvst)
        (bool, print_volume_tetra)
        (bool, print_volume_hull)
        (bool, print_velocities)
        (bool, print_coverages)
        (bool, print_rates)
        (bool, print_materials)
        (bool, print_statistics)
        (int, max_extended_starting_position)
        (int, open_boundary)
        (bool, remove_bottom)
        (double, snap_to_boundary_eps)
        (int, process_cycles)
        (std::vector<int>, material_mapping)
        (int, add_layer)
        (std::vector<bnc::boundary_condition_type>, boundary_condition)
        (std::list<client::Parameters::ProcessParameterType>, process_parameters)
    )

    BOOST_FUSION_ADAPT_STRUCT(client::Parameters::ProcessParameterType,
        (double, ProcessDistance)
        (double, AFMStartPosition[3])
        (double, AFMEndPosition[3])
        (int, AddLayer)
        (double, ProcessTime)
        (unsigned int, ALDStep)
        (std::string, ModelName)
        (std::string, ModelParameters)
        (int, IterationCycles)
        (int, StartIterationCycles)
        (double, MaxTimeStep)
        (double, smoothing_max_curvature)
        (double, smoothing_min_curvature)
        (int, smoothing_material_level)
        (int, smoothing_max_iterations)
        (bool, print_coverages)
        (bool, print_rates)
        (bool, print_velocities)
        (bool, print_materials)
        (proc::FiniteDifferenceSchemeType, FiniteDifferenceScheme)
        (double, LaxFriedrichsDissipationCoefficient)
        (partition::DataStructureType, partition_data_structure)
        (partition::SplittingType, partition_splitting_strategy)
        (partition::SurfaceAreaHeuristicLambdaType, partition_surface_area_heuristic_lambda)
        (std::vector<double>, output_times)
        (std::vector<double>, output_volume)
        (int, output_times_periodicity)
        (double, output_times_period_length)
        (bool, initial_output)
        (bool, final_output)
        (bool, MaskLayer)
        (bool, GrowNewOxide)
    )

    namespace client{
        //Testing function to output parsed value directly
        std::string output_parsed(std::string value){
            std::cout << value << std::endl;
            return value;
        }


        template <typename Iterator>
        struct skipper_grammar : qi::grammar<Iterator>{
            skipper_grammar() : skipper_grammar::base_type(skip_on){
                using ascii::char_;
                using ascii::space;
                using qi::eol;

                skip_on = space | eol | ("/*" >> *(char_ - "*/") >> "*/") | ("//" >> *(char_ - eol) >> eol);
            }
            qi::rule<Iterator> skip_on;
        };

        template <typename Iterator, typename Skipper>
        struct par_grammar : qi::grammar<Iterator, Parameters(), Skipper>{

            par_grammar() : par_grammar::base_type(definition, "Definition"){
                using qi::lit;
                using qi::lexeme;
                using ascii::char_;
                using qi::int_;
                using qi::double_;
                using qi::on_error;
                using qi::fail;
                using qi::_val;
                using qi::_1;

                using boost::phoenix::val;
                using boost::phoenix::construct;
                using namespace qi::labels;

                //Define names for error parsing
                definition.name("Parameter File");
                quotes.name("name");
                listed.name("list");
                boolean.name("boolean argument");
                process_distance.name("process_distance");
                model_parameters.name("model_parameters");

                //define rules to be used for parsing ----------------------------------
                quotes %= '"' >> lexeme[+(char_ - '"') > '"'];
                listed %= '"' >> +(lexeme[+(char_ - '"' - ',')] % ',') >> '"';
                boolean %= boolean_symbols[_val=_1];
                intvec %= '{' >> +(lexeme[int_] % ',') > '}';
                doublevec %= '{' >> +(lexeme[double_] % ',') > '}';
                boundaryvec %= (('{' > +(lexeme[boundary_condition_symbols] % ',') > '}') % ',');
                model_param_end = lit(';') >> lit('}') >> lit(';');



                //ADD CUSTOM RULES HERE
                //rules have to be in the same order in which the variables are declared in the struct
                geometry_file %= (lit("geometry_files") | lit("GeometryFiles") | lit("geometry_file") | lit("GeometryFile")) > "=" > listed > ";";
                output_path %= (lit("output_path") | lit("OutputPath")) > "=" > quotes > ";";
                surface_geometry %= lit("surface_geometry") > "=" > boolean > ";";
                report_import_errors %= lit("report_import_errors") > "=" > boolean > ";";
                cfl_condition %= lit("cfl_condition") > "=" > lexeme[double_] > ";";
                input_scale %= lit("input_scale") > "=" > lexeme[double_] > ";";
                grid_delta %= lit("grid_delta") > "=" > lexeme[double_] > ";";
                input_transformation %= lit("input_transformation") > '=' > intvec > ';';
                input_shift %= lit("input_shift") > '=' > doublevec > ';';
                default_disc_orientation %= (lit("default_disc_orientation") | lit("default_disk_orientation")) > '=' > doublevec > ';';
                ignore_materials %= lit("ignore_materials") > '=' > intvec > ';';
                change_input_parity %= lit("change_input_parity") > '=' > boolean > ';';
                random_seed %= lit("random_seed") > '=' > lexeme[int_] > ';';
                num_dimensions %= lit("num_dimensions") > '=' > lexeme[int_] > ';';
                omp_threads %= lit("omp_threads") > '=' > lexeme[int_] > ';';
                domain_extension %= lit("domain_extension") > '=' > lexeme[int_] > ';';
                receptor_radius %= lit("receptor_radius") > '=' > lexeme[double_] > ';';
                further_tracking_distance %= lit("further_tracking_distance") > '=' > lexeme[double_] > ';';
                bits_per_distance %= lit("bits_per_distance") > '=' > lexeme[int_] > ';';
                output_volume_extract_single_materials %= lit("output_volume_extract_single_materials") > '=' > boolean > ';';
                print_vtk %= lit("print_vtk") > '=' > boolean > ';';
                print_vtp %= lit("print_vtp") > '=' > boolean > ';';
                print_dx %= lit("print_dx") > '=' > boolean > ';';
                print_lvst %= lit("print_lvst") > '=' > boolean > ';';
                print_volume_tetra %= lit("print_volume_tetra") > '=' > boolean > ';';
                print_volume_hull %= lit("print_volume_hull") > '=' > boolean > ';';
                print_velocities %= lit("print_velocities") > '=' > boolean > ';';
                print_coverages %= lit("print_coverages") > '=' > boolean > ';';
                print_rates %= lit("print_rates") > '=' > boolean > ';';
                print_materials %= lit("print_materials") > '=' > boolean > ';';
                print_statistics %= lit("print_statistics") > '=' > boolean > ';';
                max_extended_starting_position %= lit("max_extended_starting_postion") > '=' > lexeme[int_] > ';';
                open_boundary %= lit("open_boundary") > '=' > '"' > lexeme[cartesian_dir_symbols] > '"' > ';';
                remove_bottom %= lit("remove_bottom") > '=' > boolean > ';';
                snap_to_boundary_eps %= lit("snap_to_boundary_eps") > '=' > lexeme[double_] > ';';
                process_cycles %= lit("process_cycles") > '=' > lexeme[int_] > ';';
                material_mapping %= lit("material_mapping") > '=' > doublevec > ';';
                add_layer %= lit("add_layer") > '=' > lexeme[int_] > ';';
                boundary_conditions %= lit("boundary_conditions") > '=' > '{' > boundaryvec > '}' > ';';

                // PROCESS PARSING
                process_distance %= lit("process_distance") > '=' > lexeme[double_] > ';';
                afm_start_position %= lit("afm_start_position") > '=' > lexeme[double_] > ';';
                afm_end_position %= lit("afm_end_position") > '=' > lexeme[double_] > ';';
                process_time %= lit("process_time") > '=' > lexeme[double_] > ';';
                ald_step %= (lit("ald_step") | lit("ALD_step")) > '=' > lexeme[int_] > ';';
                model_name %= lit("model_name") > '=' > quotes > ';';
                model_parameters %= lit("parameters") > '=' > '{' > +(char_ - (model_param_end)) > model_param_end;    //do not use lexeme to allow skipping for end detection
                iteration_cycles %= lit("iteration_cycles") > '=' > lexeme[int_] > ';';
                start_iteration_cycles %= lit("start_iteration_cycles") > '=' > lexeme[int_] > ';';
                max_time_step %= lit("max_time_step") > '=' > lexeme[double_] > ';';
                smoothing_max_curvature %= lit("smoothing_max_curvature") > '=' > lexeme[double_] > ';';
                smoothing_min_curvature %= lit("smoothing_min_curvature") > '=' > lexeme[double_] > ';';
                smoothing_material_level %= lit("smoothing_material_level") > '=' > lexeme[int_] > ';';
                smoothing_max_iterations %= lit("smoothing_max_iterations") > '=' > lexeme[int_] > ';';
                finite_difference_scheme %= lit("finite_difference_scheme") > '=' > lexeme[finite_difference_symbols] > ';';
                dissipation_coefficient %= lit("dissipation_coefficient") > '=' > lexeme[double_] > ';';
                partition_data_structure %= lit("partition_data_structure") > '=' > lexeme[partition_data_structure_symbols] > ';';
                partition_splitting_strategy %= lit("partition_splitting_strategy") > '=' > lexeme[partition_splitting_strategy_symbols] > ';';
                partition_surface_area_heuristic_lambda %= lit("partition_surface_area_heuristic_lambda") > '=' > lexeme[double_] > ';';
                output_times %= (lit("output_times") - lit("output_times_periodicity") - lit("output_times_period_length")) > '=' > doublevec > ';';
                output_volume %= lit("output_volume") > '=' > doublevec > ';';
                output_times_periodicity %= lit("output_times_periodicity") > '=' > lexeme[int_] > ';';
                output_times_period_length %= lit("output_times_period_length") > '=' > lexeme[double_] > ';';
                initial_output %= lit("initial_output") > '=' > boolean > ';';
                final_output %= lit("final_output") > '=' > boolean > ';';
                mask_layer %= lit("mask_layer") > '=' > boolean > ';';
                grow_new_oxide %= lit("grow_new_oxide") > '=' > boolean > ';';

                // PROCESS PARAMETERS
                //ORDER OF RULES MUST MATCH VARIABLE DECLARATION IN STRUCT
                process_params %= process_distance ^ afm_start_position ^ afm_end_position ^ add_layer ^
                process_time ^ ald_step ^ model_name ^ model_parameters ^ iteration_cycles ^
                start_iteration_cycles ^ max_time_step ^ smoothing_max_curvature ^ smoothing_min_curvature ^
                smoothing_material_level ^ smoothing_max_iterations ^ print_coverages ^ print_rates ^
                print_velocities ^ print_materials ^ finite_difference_scheme ^ dissipation_coefficient ^
                partition_data_structure ^ partition_splitting_strategy ^ partition_surface_area_heuristic_lambda ^
                output_times ^ output_volume ^ output_times_periodicity ^ output_times_period_length ^
                initial_output ^ final_output ^ mask_layer ^ grow_new_oxide;

                process %= lit('{') > process_params > lit('}');
                processes %= lit("processes") > lit('=') > lit('{') > (process % ',') > lit('}') > lit(';');

                //ORDER OF RULES MUST MATCH VARIABLE DECLARATION IN STRUCT
                definition %= geometry_file ^ output_path ^ surface_geometry ^ report_import_errors ^
                cfl_condition ^ input_scale ^ grid_delta ^ input_transformation ^
                input_shift ^ default_disc_orientation ^ ignore_materials ^
                change_input_parity ^ random_seed ^ num_dimensions ^ omp_threads ^ domain_extension ^
                receptor_radius ^ further_tracking_distance ^ bits_per_distance ^ output_volume_extract_single_materials ^ print_vtk ^ print_vtp ^ print_dx ^ print_lvst ^
                print_volume_tetra ^ print_volume_hull ^ print_velocities ^ print_coverages ^
                print_rates ^ print_materials ^
                print_statistics ^ max_extended_starting_position ^ open_boundary ^
                remove_bottom ^ snap_to_boundary_eps ^ process_cycles ^ material_mapping ^
                add_layer ^ boundary_conditions ^ processes;

                // what to do when parsing fails
                on_error<fail>
                (
                    definition
                    , report_error(_1, _2, _3, _4)
                );
            }
        private:
            // Rules for the parser ------------------------------------------------------
            qi::rule<Iterator, Parameters(), Skipper> definition;
            qi::rule<Iterator, std::string(), Skipper> quotes;
            qi::rule<Iterator, std::vector<std::string>(), Skipper> listed;
            qi::rule<Iterator, std::vector<int>(), Skipper> intvec;
            qi::rule<Iterator, std::vector<double>(), Skipper> doublevec;
            qi::rule<Iterator, bool(), Skipper> boolean;
            qi::rule<Iterator, Parameters::ProcessParameterType(), Skipper> process;
            qi::rule<Iterator, Parameters::ProcessParameterType(), Skipper> process_params;
            qi::rule<Iterator, std::list<Parameters::ProcessParameterType>(), Skipper> processes;
            qi::rule<Iterator, Skipper> model_param_end;

            // ADD CUSTOM RULES HERE (use right rule template)
            qi::rule<Iterator, std::vector<std::string>(), Skipper> geometry_file;
            qi::rule<Iterator, std::vector<bnc::boundary_condition_type>(), Skipper> boundaryvec;
            qi::rule<Iterator, std::vector<bnc::boundary_condition_type>(), Skipper> boundary_conditions;
            qi::rule<Iterator, std::vector<int>(), Skipper> input_transformation,
                ignore_materials,
                material_mapping;
            qi::rule<Iterator, std::vector<double>(), Skipper> input_shift,
                default_disc_orientation,
                output_times,
                output_volume;
            qi::rule<Iterator, std::string(), Skipper> output_path,
                model_name,
                model_parameters;
            qi::rule<Iterator, bool(), Skipper> surface_geometry,
                report_import_errors,
                change_input_parity,
                output_volume_extract_single_materials,
                print_vtk,
                print_vtp,
                print_dx,
                print_lvst,
                print_volume_hull,
                print_volume_tetra,
                print_velocities,
                print_coverages,
                print_rates,
                print_materials,
                print_statistics,
                remove_bottom,
                initial_output,
                final_output,
                mask_layer,
                grow_new_oxide;
            qi::rule<Iterator, double(), Skipper> cfl_condition,
                input_scale,
                grid_delta,
                domain_extension,
                receptor_radius,
                further_tracking_distance,
                snap_to_boundary_eps,
                process_distance,
                afm_start_position,
                afm_end_position,
                process_time,
                max_time_step,
                smoothing_max_curvature,
                smoothing_min_curvature,
                dissipation_coefficient,
                output_times_period_length;
            qi::rule<Iterator, int(), Skipper> random_seed,
                num_dimensions,
                omp_threads,
                bits_per_distance,
                max_extended_starting_position,
                open_boundary,
                process_cycles,
                add_layer,
                ald_step,
                iteration_cycles,
                start_iteration_cycles,
                smoothing_material_level,
                smoothing_max_iterations,
                output_times_periodicity,
                partition_data_structure,
                partition_splitting_strategy,
                partition_surface_area_heuristic_lambda,
                finite_difference_scheme;
        };
    }


    client::Parameters::Parameters(std::string fileName){
        surface_geometry=false;
        report_import_errors=true;
        input_scale=1.;
        change_input_parity=false;
        omp_threads=0;
        domain_extention=0.;
        receptor_radius=0.8;
        further_tracking_distance=3.;
        bits_per_distance=8;
        output_volume_extract_single_materials=true;
        print_vtk=false;
        print_vtp=false;
        print_dx=false;
        print_lvst=true;
        print_volume_tetra=true;
        print_volume_hull=false;
        print_velocities=false;
        print_coverages=false;
        print_rates=false;
        print_materials=false;
        print_statistics=false;
        max_extended_starting_position=1000;
        open_boundary_negative=false;
        remove_bottom=true;
        snap_to_boundary_eps=1e-6;
        process_cycles=1;
        add_layer=0;

        //REQUIRED (set dummy variables to check if they were declared)
        cfl_condition=-1.;
        grid_delta=-1.;
        num_dimensions=-1;
        open_boundary=0;

        //Open file
        std::string input;
        std::ifstream inFile(fileName);
        if(inFile){
            std::ostringstream oss;
            oss << inFile.rdbuf();
            input = oss.str();
            inFile.close();
        }else abort();


        //Start parsing
        typedef client::skipper_grammar<std::string::iterator> skipper_type;
        client::par_grammar<std::string::iterator, skipper_type> g;

        std::string::iterator start=input.begin(), end=input.end();

        bool r = phrase_parse(start, end, g, skipper_type(), *this);

        //check if grammar is correct
        if(!r){
            report_input_error(input.begin(), end, start, "INPUT ERROR: Illegal parameter name or other unexpected name specified.");

            abort();
        }

        //check if file was parsed until the end
        if(start!=end){
            report_input_error(input.begin(), end, start, "INPUT ERROR: Illegal parameter name or duplicate parameter definition.");
            abort();
        }

        //check if semantics are correct
        if(!validate()) abort();


        //IMPORTANT: MODEL PARAMETERS ARE MISSING }; AT THE END, NEED TO ADD BEFORE ANYTHING ELSE
        for(std::list<client::Parameters::ProcessParameterType>::iterator it=process_parameters.begin(); it!=process_parameters.end(); ++it) it->ModelParameters.append(";");

        // fix input transformation
        for(size_t i=0; i<input_transformation.size(); ++i) input_transformation_signs[i] = input_transformation[i] >= 0;

        //set boundary conditions to boundary_conditions_type
        boundary_conditions.resize(3);
        for(int i=0; i<3; ++i){
            boundary_conditions[i].min = boundary_condition[2*i];
            boundary_conditions[i].max = boundary_condition[2*i+1];
        }

        //fix boundary parameters and deduce open boundary if necessary
        if(open_boundary==0){
            for(size_t i=0; i<boundary_condition.size(); ++i){
                if(boundary_condition[i]==bnc::INFINITE_BOUNDARY){
                    open_boundary=i/2;
                    open_boundary_negative=(i+1)%2;
                    break;
                }
            }
        }else{
            if(open_boundary < 0){
                open_boundary_negative=true;
                open_boundary*=-1;
            }
            --open_boundary;
        }

        //fix bits_per_distance so there are no overflow bits when writing to lvst file
        //meaning that bits_per_distance are mapped to 1, 2, 4, 8, 16, 24, 32, 48, 56, 64
        int bpd = 1;
        while( bpd < bits_per_distance || (bpd < CHAR_BIT ? CHAR_BIT % bpd : bpd % CHAR_BIT) ){
          if(bpd < CHAR_BIT) bpd++;
          else bpd += CHAR_BIT; //once bpd (bits per distance) reaches the amount of bits per char, the increase will be CHAR_BIT instead of 1
        }
        if(bits_per_distance == 0) bits_per_distance = CHAR_BIT;
        else bits_per_distance = bpd > 64 ? 64 : bpd;

        //Create directory, if it does not exist; make it relative to the parameter file's path
        std::string relative_path, relative_name;

        std::size_t last_slash = fileName.find_last_of('/');
        if(last_slash != std::string::npos){//if parameter file is in a different directory
          relative_path = fileName.substr(0,last_slash+1);//path of parameter file
          relative_name = fileName.substr(last_slash+1);//filename
        }
        else relative_name = fileName;


        if(output_path.front() != '/'){
          if(output_path.size() == 0){//default output path; example: par_trench.txt -> output_trench/
            output_path = relative_name.substr(0, relative_name.find_last_of('.'));//exclude file ending

            if(output_path.substr(0,3) == "par")//check if par is in the filename
              output_path = output_path.substr(3);//exclude par in file name
            output_path = relative_path + "output" + output_path;
          }
          else
            output_path = relative_path + output_path;
        }

        //make path of input files relative too
        for(unsigned int i=0; i<geometry_files.size(); i++){
          if(geometry_files[i].front() != '/') geometry_files[i] = relative_path + geometry_files[i];
        }

        //add _bits_per_distance/ to output path
        {
          std::ostringstream oss;
          if(output_path.back() == '/')
            oss << output_path.substr(0, output_path.size()-1);
          else
            oss << output_path;
          oss <<  "_" << bits_per_distance << "bit/";
          output_path = oss.str();
          if(!boost::filesystem::exists(output_path)) {
              msg::print_message("Output directory not found! Creating new directory: "+output_path+"\n");
                        boost::filesystem::path dir(output_path);
                        if(!boost::filesystem::create_directory(dir)) msg::print_error("Could not create directory!");
          }
        }

        //test boundary condtions for periodicity, and check conformity
        for (unsigned int h=0;h<boundary_conditions.size();h++) {
            if (boundary_conditions[h].min==bnc::PERIODIC_BOUNDARY) assert(boundary_conditions[h].max==bnc::PERIODIC_BOUNDARY);
            if (boundary_conditions[h].max==bnc::PERIODIC_BOUNDARY) assert(boundary_conditions[h].min==bnc::PERIODIC_BOUNDARY);
            if (h==static_cast<unsigned int>(open_boundary)) assert(boundary_conditions[h].min==bnc::INFINITE_BOUNDARY);
            if (h==static_cast<unsigned int>(open_boundary)) assert(boundary_conditions[h].max==bnc::INFINITE_BOUNDARY);
        }

        //add OutputTimes
        for(std::list<ProcessParameterType>::iterator pIter=process_parameters.begin();pIter!=process_parameters.end();++pIter) {
            // add output times for initial and final output
            if (pIter->initial_output) pIter->output_times.push_back(0);
            if (pIter->final_output) pIter->output_times.push_back(pIter->ProcessTime);

            // Create output times from period length and periodicity
            if (pIter->output_times_period_length>0 && pIter->output_times_periodicity>=0) {
                if(pIter->output_times.empty()) pIter->output_times.push_back(0.);  //if no output times, use 0 as start

                unsigned int tmp=pIter->output_times.size();
                pIter->output_times.resize(pIter->output_times_periodicity*tmp);

                for (int h=1;h<pIter->output_times_periodicity;++h) {
                    for (unsigned int k=0;k<tmp;k++) {
                        pIter->output_times[tmp*h+k]=pIter->output_times[k]+h*pIter->output_times_period_length;
                    }
                }
            }

            // sort volume output and add it to surface output, so that a surface is output each time a volume is output. this is required for the advection to stop at each volume output. It is also necessary for a downsampling of the volume being output
            // Sort Volume Output
            std::sort(pIter->output_volume.begin(), pIter->output_volume.end());
            std::vector<double>::iterator it = std::unique(pIter->output_volume.begin(), pIter->output_volume.end());
            pIter->output_volume.resize(it-pIter->output_volume.begin());
            // add volume output times to surface output times
            pIter->output_times.insert(pIter->output_times.end(), pIter->output_volume.begin(), pIter->output_volume.end());

            //make sure they are in correct order and is unique
            std::sort(pIter->output_times.begin(),pIter->output_times.end());
            it = std::unique(pIter->output_times.begin(), pIter->output_times.end());
            pIter->output_times.resize(it-pIter->output_times.begin());


        }
    }
