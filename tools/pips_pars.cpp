/*--------------------------------------------------------------------------*/
/*---------------------------- File pips_pars.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small tool for parsing the standard option maps that enable the support for
 * plain PIPS-IPM++ parameters into the SMS++ PIPS-IPM++ solver class.
 *
 * The tool generates two files, <prefix>_defs.h and <prefix>_maps.h in the
 * specified path.  The path should be the include directory of the MILPSolver
 * source tree.
 *
 * PIPS-IPM++ does not expose a C API to enumerate options like CPLEX, Gurobi,
 * SCIP, or HiGHS.  Its parameters are registered as assignments to
 * bool_options, int_options, double_options, and string_options in
 * PIPS-IPM/Core/Options.
 * Therefore, this generator parses those source files directly.
 *
 * By default bool PIPS-IPM++ options are mapped to SMS++ int parameters, as in
 * the other generators where boolean solver options are handled as integer
 * SMS++ parameters.
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Matematica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>
#include <getopt.h>

/*--------------------------------------------------------------------------*/

bool verbose = false;              ///< If the tool should be verbose
std::string path{};                ///< Path for output files
std::string source_path{};         ///< Path to PIPS-IPM/Core/Options
std::string exe{};                 ///< Name of the executable file
std::string docopt_desc{};         ///< Tool description
std::string solver_class = "PIPSMILPSolver";
std::string solver_header = "PIPSMILPSolver.h";
std::string prefix = "PIPS";

/*--------------------------------------------------------------------------*/

/// Gets the name of the executable from its full path
std::string get_filename( const std::string & fullpath )
{
 std::size_t found = fullpath.find_last_of( "/\\" );
 return( fullpath.substr( found + 1 ) );
}

/*--------------------------------------------------------------------------*/

/// Checks whether a regular file exists
bool file_exists( const std::string & filename )
{
 struct stat buffer;
 return( stat( filename.c_str() , &buffer ) == 0 && S_ISREG( buffer.st_mode ) );
}

/*--------------------------------------------------------------------------*/

/// Checks whether a directory exists
bool dir_exists( const std::string & dirname )
{
 struct stat buffer;
 return( stat( dirname.c_str() , &buffer ) == 0 && S_ISDIR( buffer.st_mode ) );
}

/*--------------------------------------------------------------------------*/

/// Removes one trailing slash, if present
std::string strip_trailing_slash( const std::string & p )
{
 if( p.size() > 1 && ( p.back() == '/' || p.back() == '\\' ) )
  return( p.substr( 0 , p.size() - 1 ) );
 return( p );
}

/*--------------------------------------------------------------------------*/

/// Returns the basename of a path
std::string basename( const std::string & p )
{
 std::string pp = strip_trailing_slash( p );
 const std::size_t found = pp.find_last_of( "/\\" );
 if( found == std::string::npos )
  return( pp );
 return( pp.substr( found + 1 ) );
}

/*--------------------------------------------------------------------------*/

/// Finds the PIPS-IPM/Core/Options directory from a user-supplied path
std::string normalize_source_path( const std::string & p )
{
 const std::string pp = strip_trailing_slash( p );

 if( file_exists( pp + "/PIPSIPMppOptions.C" ) &&
     file_exists( pp + "/Options.C" ) &&
     file_exists( pp + "/AbstractOptions.C" ) )
  return( pp );

 if( file_exists( pp + "/Core/Options/PIPSIPMppOptions.C" ) )
  return( pp + "/Core/Options" );

 if( file_exists( pp + "/PIPS-IPM/Core/Options/PIPSIPMppOptions.C" ) )
  return( pp + "/PIPS-IPM/Core/Options" );

 throw std::runtime_error( "Could not locate PIPS-IPM/Core/Options from "
  "path: " + p );
}

/*--------------------------------------------------------------------------*/

/// Finds the PIPS-IPM/Core/Options directory automatically
std::string find_source_path()
{
 const std::vector< std::string > candidates = {
  ".",
  "..",
  "../..",
  "PIPS-IPM/Core/Options",
  "../PIPS-IPM/Core/Options",
  "../../PIPS-IPM/Core/Options",
  "PIPS-IPMpp/PIPS-IPM/Core/Options",
  "../PIPS-IPMpp/PIPS-IPM/Core/Options",
  "../../PIPS-IPMpp/PIPS-IPM/Core/Options"
 };

 for( const auto & c : candidates ) {
  try {
   return( normalize_source_path( c ) );
  }
  catch( const std::exception & ) {}
 }

 throw std::runtime_error(
  "Could not automatically locate PIPS-IPM/Core/Options. "
  "Use -s or --source to specify it." );
}

/*--------------------------------------------------------------------------*/

/// Prints the tool description and usage
void docopt()
{
 // http://docopt.org
 std::cout << docopt_desc << std::endl;
 std::cout << "Usage:\n"
           << "  " << exe << " [-v] [-s <source>] [-c <class>] [-i <header>]"
           << " [-p <prefix>] <path>\n"
           << "  " << exe << " -h | --help\n"
           << std::endl
           << "Options:\n"
           << "  -v, --verbose          Make the tool verbose.\n"
           << "  -s, --source <source>  Path to PIPS-IPM/Core/Options, to "
            "PIPS-IPM,"
           << " or to the PIPS-IPM++ root.\n"
           << "                         If omitted, the tool searches common "
            "local"
           << " repository locations.\n"
           << "  -c, --class <class>    Solver class name used in the "
            "generated maps"
           << " [default: PIPSMILPSolver].\n"
           << "  -i, --include <header> Solver header included by the "
            "generated maps"
           << " [default: PIPSMILPSolver.h].\n"
           << "  -p, --prefix <prefix>  Prefix used for macros, arrays and "
            "output files"
           << " [default: PIPS].\n"
           << "  -h, --help             Print this help.\n";
}

/*--------------------------------------------------------------------------*/

/// Processes the command line arguments
void process_args( int argc , char ** argv )
{
 const char * const short_opts = "vhs:c:i:p:";
 const option long_opts[] = {
  { "verbose" , no_argument ,       nullptr , 'v' } ,
  { "help" ,    no_argument ,       nullptr , 'h' } ,
  { "source" ,  required_argument , nullptr , 's' } ,
  { "class" ,   required_argument , nullptr , 'c' } ,
  { "include" , required_argument , nullptr , 'i' } ,
  { "prefix" ,  required_argument , nullptr , 'p' } ,
  { nullptr ,   no_argument ,       nullptr , 0 }
 };

 // Options
 while( true ) {
  const auto opt = getopt_long( argc , argv , short_opts , long_opts , 
                                nullptr );

  if( -1 == opt )
   break;

  switch( opt ) {
   case( 'v' ):
    verbose = true;
    break;
   case( 'h' ):
    docopt();
    exit( 0 );
   case( 's' ):
    source_path = std::string( optarg );
    break;
   case( 'c' ):
    solver_class = std::string( optarg );
    break;
   case( 'i' ):
    solver_header = std::string( optarg );
    break;
   case( 'p' ):
    prefix = std::string( optarg );
    break;
   case( '?' ):
   default:
    std::cout << "Try " << exe << "' --help' for more information.\n";
    exit( 1 );
  }
 }

 // Last argument
 if( optind < argc )
  path = std::string( argv[ optind ] );
}

/*--------------------------------------------------------------------------*/

/// Custom terminate function to print the exception message
void smspp_terminate( void ) {
 std::cerr << "Uncaught exception in executing SMS++:\n";
 try {
  std::rethrow_exception( std::current_exception() );
 }
 catch( const std::exception & e ) {
  std::cerr << "\tException type: " << typeid( e ).name() << "\n";
  std::cerr << "\tException message: " << e.what() << "\n";
 } catch( ... ) {
  std::cerr << "\tUnknown exception" << std::endl;
 }
 std::abort(); // or exit(1)
}

/*--------------------------------------------------------------------------*/

/// Reads an entire file into a string
std::string read_file( const std::string & filename )
{
 std::ifstream file( filename );
 if( ! file )
  throw std::runtime_error( "Could not open " + filename + ": " +
                            std::strerror( errno ) );

 return( std::string( std::istreambuf_iterator< char >( file ) ,
                      std::istreambuf_iterator< char >() ) );
}


/*--------------------------------------------------------------------------*/

/// Trims leading and trailing whitespace
std::string trim( const std::string & s )
{
 const auto first = s.find_first_not_of( " \t\r\n" );
 if( first == std::string::npos )
  return( "" );
 const auto last = s.find_last_not_of( " \t\r\n" );
 return( s.substr( first , last - first + 1 ) );
}

/*--------------------------------------------------------------------------*/

/// Removes C/C++ comments while preserving string literals
std::string remove_comments( const std::string & text )
{
 std::string out;
 out.reserve( text.size() );
 bool in_string = false;
 bool in_char = false;
 bool in_line_comment = false;
 bool in_block_comment = false;
 bool escaped = false;

 for( std::size_t i = 0 ; i < text.size() ; ++i ) {
  const char c = text[ i ];
  const char n = ( i + 1 < text.size() ? text[ i + 1 ] : '\0' );

  if( in_line_comment ) {
   if( c == '\n' ) {
    in_line_comment = false;
    out.push_back( c );
   }
   else
    out.push_back( ' ' );
   continue;
  }

  if( in_block_comment ) {
   if( c == '*' && n == '/' ) {
    in_block_comment = false;
    out.push_back( ' ' );
    out.push_back( ' ' );
    ++i;
   }
   else
    out.push_back( c == '\n' ? '\n' : ' ' );
   continue;
  }

  if( in_string ) {
   out.push_back( c );
   if( escaped )
    escaped = false;
   else if( c == '\\' )
    escaped = true;
   else if( c == '"' )
    in_string = false;
   continue;
  }

  if( in_char ) {
   out.push_back( c );
   if( escaped )
    escaped = false;
   else if( c == '\\' )
    escaped = true;
   else if( c == '\'' )
    in_char = false;
   continue;
  }

  if( c == '/' && n == '/' ) {
   in_line_comment = true;
   out.push_back( ' ' );
   out.push_back( ' ' );
   ++i;
   continue;
  }

  if( c == '/' && n == '*' ) {
   in_block_comment = true;
   out.push_back( ' ' );
   out.push_back( ' ' );
   ++i;
   continue;
  }

  if( c == '"' )
   in_string = true;
  else if( c == '\'' )
   in_char = true;

  out.push_back( c );
 }

 return( out );
}

/*--------------------------------------------------------------------------*/

/// Splits a function argument list at top-level commas
std::vector< std::string > split_top_level_arguments( const std::string & args )
{
 std::vector< std::string > result;
 std::string current;
 int paren_depth = 0;
 int brace_depth = 0;
 int bracket_depth = 0;
 int angle_depth = 0;
 bool in_string = false;
 bool in_char = false;
 bool escaped = false;

 for( std::size_t i = 0 ; i < args.size() ; ++i ) {
  const char c = args[ i ];

  if( in_string ) {
   current.push_back( c );
   if( escaped )
    escaped = false;
   else if( c == '\\' )
    escaped = true;
   else if( c == '"' )
    in_string = false;
   continue;
  }

  if( in_char ) {
   current.push_back( c );
   if( escaped )
    escaped = false;
   else if( c == '\\' )
    escaped = true;
   else if( c == '\'' )
    in_char = false;
   continue;
  }

  if( c == '"' ) {
   in_string = true;
   current.push_back( c );
   continue;
  }

  if( c == '\'' ) {
   in_char = true;
   current.push_back( c );
   continue;
  }

  switch( c ) {
   case( '(' ): ++paren_depth; break;
   case( ')' ): --paren_depth; break;
   case( '{' ): ++brace_depth; break;
   case( '}' ): --brace_depth; break;
   case( '[' ): ++bracket_depth; break;
   case( ']' ): --bracket_depth; break;
   case( '<' ): ++angle_depth; break;
   case( '>' ): if( angle_depth > 0 ) --angle_depth; break;
   default: break;
  }

  if( c == ',' && paren_depth == 0 && brace_depth == 0 &&
      bracket_depth == 0 && angle_depth == 0 ) {
   result.push_back( trim( current ) );
   current.clear();
  }
  else
   current.push_back( c );
 }

 result.push_back( trim( current ) );
 return( result );
}

/*--------------------------------------------------------------------------*/

/// Finds the closing parenthesis matching open_pos
std::size_t find_matching_paren( const std::string & text ,
                                 std::size_t open_pos )
{
 int depth = 0;
 bool in_string = false;
 bool in_char = false;
 bool escaped = false;

 for( std::size_t i = open_pos ; i < text.size() ; ++i ) {
  const char c = text[ i ];

  if( in_string ) {
   if( escaped )
    escaped = false;
   else if( c == '\\' )
    escaped = true;
   else if( c == '"' )
    in_string = false;
   continue;
  }

  if( in_char ) {
   if( escaped )
    escaped = false;
   else if( c == '\\' )
    escaped = true;
   else if( c == '\'' )
    in_char = false;
   continue;
  }

  if( c == '"' ) {
   in_string = true;
   continue;
  }
  if( c == '\'' ) {
   in_char = true;
   continue;
  }

  if( c == '(' )
   ++depth;
  else if( c == ')' ) {
   --depth;
   if( depth == 0 )
    return( i );
  }
 }

 return( std::string::npos );
}

/*--------------------------------------------------------------------------*/

/// Extracts a string literal argument such as "PARAM" or std::string("PARAM")
std::string extract_string_literal( const std::string & arg )
{
 static const std::regex literal_regex( "\\\"([^\\\"]+)\\\"" );
 std::smatch match;
 if( std::regex_search( arg , match , literal_regex ) )
  return( match[ 1 ].str() );
 return( "" );
}

/*--------------------------------------------------------------------------*/

/// Parses PIPS-IPM++ add_parameter(...) registrations.
void parse_add_parameter_calls( const std::string & text ,
                            std::map< int , std::string > & int_parameters ,
                            std::map< int , std::string > & dbl_parameters ,
                            std::map< int , std::string > & str_parameters ,
                            std::set< std::string > & int_seen ,
                            std::set< std::string > & dbl_seen ,
                            std::set< std::string > & str_seen ,
                            int & int_counter ,
                            int & dbl_counter ,
                            int & str_counter )
{
 const std::string clean = remove_comments( text );
 const std::string needle = "add_parameter";
 std::size_t pos = 0;

 while( ( pos = clean.find( needle , pos ) ) != std::string::npos ) {
  // Avoid matching the add_parameter function declarations/definitions.
  const std::size_t prev = clean.find_last_not_of( " \t\r\n" , pos == 0 ? 
                              0 : pos - 1 );
  if( prev != std::string::npos &&
      ( std::isalnum( static_cast< unsigned char >( clean[ prev ] ) ) ||
        clean[ prev ] == '_' || clean[ prev ] == ':' ) ) {
   pos += needle.size();
   continue;
  }

  const std::size_t open = clean.find( '(' , pos + needle.size() );
  if( open == std::string::npos )
   break;

  const std::size_t close = find_matching_paren( clean , open );
  if( close == std::string::npos )
   break;

  const std::vector< std::string > args =
   split_top_level_arguments( clean.substr( open + 1 , 
                                close - open - 1 ) );

  if( args.size() >= 4 ) {
   const std::string name = extract_string_literal( args[ 1 ] );
   const std::string default_value = trim( args[ 3 ] );

   if( ! name.empty() ) {
    const bool is_string = default_value.find( "std::string" ) != std::string::npos ||
                           ( ! default_value.empty() && default_value[ 0 ] == '"' );
    const bool is_bool = default_value == "true" || default_value == "false";
    const bool is_double = default_value.find( '.' ) != std::string::npos ||
                           default_value.find( 'e' ) != std::string::npos ||
                           default_value.find( 'E' ) != std::string::npos;

    if( is_string ) {
     if( str_seen.insert( name ).second )
      str_parameters.insert( { str_counter++ , name } );
    }
    else if( is_bool ) {
     if( int_seen.insert( name ).second )
      int_parameters.insert( { int_counter++ , name } );
    }
    else if( is_double ) {
     if( dbl_seen.insert( name ).second )
      dbl_parameters.insert( { dbl_counter++ , name } );
    }
    else {
     if( int_seen.insert( name ).second )
      int_parameters.insert( { int_counter++ , name } );
    }
   }
  }

  pos = close + 1;
 }
}

/*--------------------------------------------------------------------------*/

/// Parses options of the requested PIPS-IPM++ map type from a source file
void parse_options( const std::string & text ,
                    const std::string & map_name ,
                    std::map< int , std::string > & parameters ,
                    std::set< std::string > & seen ,
                    int & counter )
{
 /*
  * The usual PIPS-IPM++ spelling is
  *   bool_options["PARAM"] = value;
  * but this parser deliberately accepts a few equivalent spellings as well,
  * so that small upstream formatting changes do not silently produce empty
  * generated headers.
  */
 const std::vector< std::regex > regexes = {
  std::regex( "(?:this->)?" + map_name +
              "\\s*\\[\\s*\\\"([^\\\"]+)\\\"\\s*\\]\\s*=",
              std::regex::ECMAScript ),
  std::regex( "(?:this->)?" + map_name +
              "\\s*\\.\\s*emplace\\s*\\(\\s*\\\"([^\\\"]+)\\\"",
              std::regex::ECMAScript ),
  std::regex( "(?:this->)?" + map_name +
              "\\s*\\.\\s*insert\\s*\\(.*?\\{\\s*\\\"([^\\\"]+)\\\"",
              std::regex::ECMAScript )
 };

 for( const auto & option_regex : regexes ) {
  auto begin = std::sregex_iterator( text.begin() , text.end() , option_regex );
  auto end = std::sregex_iterator();

  for( auto it = begin ; it != end ; ++it ) {
   const std::string name = ( *it )[ 1 ].str();
   if( seen.insert( name ).second )
    parameters.insert( { counter++ , name } );
  }
 }
}

/*--------------------------------------------------------------------------*/

/// Writes a std::array<string, ...> mapping from SMS++ ids to native options
void write_forward_map( std::ofstream & maps_file ,
                        const std::string & type ,
                        const std::map< int , std::string > & parameters )
{
 maps_file
  << "const std::array< std::string, " << prefix << "_NUM_"
  << type << "_PARS >"
  << " " << solver_class << "::SMSpp_to_" << prefix << "_"
  << ( type == "INT" ? "int" : type == "DBL" ? "dbl" : "str" )
  << "_pars{" << std::endl;

 for( const auto & i : parameters )
  maps_file << " \"" << i.second << "\"," << std::endl;

 maps_file << "};" << std::endl << std::endl;
}

/*--------------------------------------------------------------------------*/

/// Writes a reverse std::array<pair<string, int>, ...> mapping
void write_reverse_map( std::ofstream & maps_file ,
                        const std::string & type ,
                        const std::string & first_name ,
                        const std::map< int , std::string > & parameters )
{
 maps_file
  << "const std::array< std::pair< std::string, int >, " << prefix
  << "_NUM_" << type << "_PARS >" << std::endl
  << " " << solver_class << "::" << prefix << "_to_SMSpp_"
  << ( type == "INT" ? "int" : type == "DBL" ? "dbl" : "str" )
  << "_pars{" << std::endl
  << " {" << std::endl;

 for( const auto & i : parameters )
  maps_file << "  { \"" << i.second << "\", " << first_name << " + "
            << i.first << " }," << std::endl;

 maps_file
  << " }" << std::endl
  << "};" << std::endl
  << std::endl;
}

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 // override the default terminate handler to print the exception message
 std::set_terminate( smspp_terminate );

 // Manage options and help
 path = "../include";
 docopt_desc = "PIPS-IPM++ parameter map generator.\n";
 exe = get_filename( argv[ 0 ] );
 process_args( argc , argv );

 if( source_path.empty() )
  source_path = find_source_path();
 else
  source_path = normalize_source_path( source_path );

 source_path = strip_trailing_slash( source_path );

 if( verbose ) {
  std::cout << "PIPS-IPM++ options path is " << source_path << std::endl;
  std::cout << "Output path is " << path << std::endl;
  std::cout << "Solver class is " << solver_class << std::endl;
  std::cout << "Solver header is " << solver_header << std::endl;
  std::cout << "Prefix is " << prefix << std::endl;
 }

 std::string defs_filename = prefix + "_defs.h";
 std::string maps_filename = prefix + "_maps.h";

 auto defs_path = path + "/" + defs_filename;
 auto maps_path = path + "/" + maps_filename;

 std::map< int , std::string > int_parameters;
 std::map< int , std::string > dbl_parameters;
 std::map< int , std::string > str_parameters;

 std::set< std::string > int_seen;
 std::set< std::string > dbl_seen;
 std::set< std::string > str_seen;

 int int_counter = 0;
 int dbl_counter = 0;
 int str_counter = 0;

 const std::vector< std::string > option_sources = {
  "PIPSIPMppOptions.C"
 };

 for( const auto & file : option_sources ) {
  const std::string filename = source_path + "/" + file;
  if( ! file_exists( filename ) ) {
   if( verbose )
    std::cout << "Skipping missing file " << filename << std::endl;
   continue;
  }

  if( verbose )
   std::cout << "Parsing " << filename << std::endl;

  const std::string text = read_file( filename );

  // SMS++ represents boolean parameters as int parameters.
  parse_options( text , "bool_options" , int_parameters , int_seen ,
                 int_counter );
  parse_options( text , "int_options" , int_parameters , int_seen ,
                 int_counter );
  parse_options( text , "double_options" , dbl_parameters , dbl_seen ,
                 dbl_counter );
  parse_options( text , "string_options" , str_parameters , str_seen ,
                 str_counter );

  parse_add_parameter_calls( text , int_parameters , dbl_parameters ,
                             str_parameters , int_seen , dbl_seen , str_seen ,
                             int_counter , dbl_counter , str_counter );
 }

 if( int_counter + dbl_counter + str_counter == 0 )
  throw std::runtime_error(
   "No PIPS-IPM++ options were found. Check that -s points to the source "
   "tree containing Core/Options/PIPSIPMppOptions.C, Options.C and "
   "AbstractOptions.C, and run with -v to see which files are parsed." );

 // Generate defs file
 std::ofstream defs_file;
 defs_file.open( defs_path );
 if( ! defs_file )
  throw std::runtime_error( "Could not open output file " + defs_path + ": " +
                            std::strerror( errno ) );

 defs_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#ifndef __" << prefix << "_DEFS" << std::endl
           << "#define __" << prefix << "_DEFS" << std::endl
           << std::endl
           << "#define " << prefix << "_NUM_INT_PARS " << int_counter
           << std::endl
           << "#define " << prefix << "_NUM_DBL_PARS " << dbl_counter
           << std::endl
           << "#define " << prefix << "_NUM_STR_PARS " << str_counter
           << std::endl
           << std::endl
           << "#endif //__" << prefix << "_DEFS" << std::endl;

 defs_file.close();
 std::cout << "Defs file written on " << defs_path << std::endl;

 // Generate maps file
 std::ofstream maps_file;
 maps_file.open( maps_path );
 if( ! maps_file )
  throw std::runtime_error( "Could not open output file " + maps_path + ": " +
                            std::strerror( errno ) );

 maps_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#include <array>" << std::endl
           << "#include <string>" << std::endl
           << "#include <utility>" << std::endl
           << "#include \"" << solver_header << "\"" << std::endl
           << std::endl
           << "using namespace SMSpp_di_unipi_it;" << std::endl
           << std::endl;

 write_forward_map( maps_file , "INT" , int_parameters );
 write_forward_map( maps_file , "DBL" , dbl_parameters );
 write_forward_map( maps_file , "STR" , str_parameters );

 write_reverse_map( maps_file , "INT" , "intFirst" + prefix + "Par" ,
                    int_parameters );
 write_reverse_map( maps_file , "DBL" , "dblFirst" + prefix + "Par" ,
                    dbl_parameters );
 write_reverse_map( maps_file , "STR" , "strFirst" + prefix + "Par" ,
                    str_parameters );

 maps_file.close();
 std::cout << "Maps file written on " << maps_path << std::endl;

 if( verbose ) {
  std::cout << prefix << "_NUM_INT_PARS " << int_counter << std::endl;
  std::cout << prefix << "_NUM_DBL_PARS " << dbl_counter << std::endl;
  std::cout << prefix << "_NUM_STR_PARS " << str_counter << std::endl;
 }

 return( 0 );
}
