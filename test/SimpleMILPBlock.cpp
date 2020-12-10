/*--------------------------------------------------------------------------*/
/*----------------------- File SimpleMILPBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SimpleMILPBlock class.
 *
 * \version 0.21
 *
 * \date 31 - 01 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SimpleMILPBlock.h"

#include <iostream>
#include <numeric>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
using namespace std;

/*--------------------------------------------------------------------------*/
/*------------------------------ MACROS ------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

template< typename T>
static void read_T( istream & iStrm , T & t )
{
 iStrm >> eatcomments;
 int c = iStrm.peek();

 switch( c ) {
  case 'I' :
  case 'i' : t = Inf<T>();
             break;
  case '-' : iStrm.get();
             read_T( iStrm , t );
             t = - t;
             return;
  case 'M' :
  case 'm' : t = -Inf<T>();
             break;
  default :  iStrm >> t;
             return;
  }

 do { c = iStrm.get(); c = iStrm.peek();
  } while( ( c != iStrm.widen( ' ' ) ) &&
	   ( c != iStrm.widen( '\n' ) ) &&
	   ( c != iStrm.widen( '\t' ) ) );
 }

/*--------------------------------------------------------------------------*/

static inline int read_int( istream & iStrm )
{
 int d;
 read_T( iStrm , d );
 return( d );
 }

/*--------------------------------------------------------------------------*/

static inline double read_dbl( istream & iStrm )
{
 double d;
 read_T( iStrm , d );
 return( d );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SimpleMILPBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( SimpleMILPBlock );

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void SimpleMILPBlock::print( std::ostream &output ) const
{
 output << "SimpleMILPBlock with " << x.size() << " vars and " << A.size()
	<< " const" << std::endl;
 }

/*--------------------------------------------------------------------------*/

void SimpleMILPBlock::load( std::istream &input )
{
 if( x.size() )
  throw( std::logic_error( "loading a non-empty SimpleMILPBlock" ) );

 int n = read_int( input );
 if( n < 0 )
  throw( std::invalid_argument( "invalid number of variables" ) );

 int m = read_int( input );
 if( m < 0 )
  throw( std::invalid_argument( "invalid number of constraints" ) );

 // initialize (static) variables and box constraints - - - - - - - - - - - -

 x.resize( n );
 Bx.resize( n );
 for( unsigned int i = 0 ; i < n ; i++ ) {
  Bx[ i ].set_variable( & x[ i ] , eNoMod );
  Bx[ i ].set_lhs( read_dbl( input ) , eNoMod );
  Bx[ i ].set_rhs( read_dbl( input ) , eNoMod );
  x[ i ].set_Block( this );
  Bx[ i ].set_Block( this );
  }

 // initialize objective function - - - - - - - - - - - - - - - - - - - - - -

 int np = read_int( input );
 if( ( np < 0 ) || ( np > n ) )
  throw( invalid_argument( "invalid number of nonzeroes" ) );

 LinearFunction::v_coeff_pair p( np );
 for( int h = 0 ; h < np ; h++ ) {
  const int i = read_int( input );
  if( ( i < 0 ) || ( i > n ) )
   throw( invalid_argument( "invalid variable name" ) );

  p[ h ].first = &x[ i ];
  p[ h ].second = read_dbl( input );
  }

 f.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
 f.set_Block( this );

 // initialize the constraints- - - - - - - - - - - - - - - - - - - - - - - -

 A.resize( m );

 for( int j = 0 ; j < m ; j++ ) {
  A[ j ].set_lhs( read_dbl( input ) );
  A[ j ].set_rhs( read_dbl( input ) );

  np = read_int( input );
  if( ( np < 0 ) || ( np > n ) )
   throw( invalid_argument( "invalid number of nonzeroes" ) );

  p.resize( np );
  for( int h = 0 ; h < np ; h++ ) {
   const int i = read_int( input );
   if( ( i < 0 ) || ( i >= n ) )
    throw( invalid_argument( "invalid variable name" ) );

   p[ h ].first = &x[ i ];
   p[ h ].second = read_dbl( input );
   }

  A[ j ].set_function( new LinearFunction( std::move( p ) , 0 ) );
  A[ j ].set_Block( this );
  }

 // set the abstract representation - - - - - - - - - - - - - - - - - - - - -
 // explicitly clear static Constraint and Variable: this should not be
 // needed, but a single SimpleMILPBlock could be loaded multiple times,
 // in which case the abstract representation would pile up
 reset_static_constraints();
 reset_static_variables();
 reset_objective();

 // this is not required by the standard, but a SimpleMILPBlock only has the
 // "abstract representation", so it must always be created

 set_objective( & f , eNoMod );
 add_static_variable( x );
 add_static_constraint( A );
 add_static_constraint( Bx );

 // issue the NBModification - - - - - - - - - - - - - - - - - - - - - - - -

 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );

 }  // end( SimpleMILPBlock::load )

/*--------------------------------------------------------------------------*/

void SimpleMILPBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );  // invoke method of Block

 group.putAtt( "infinity" , netCDF::NcDouble() , Inf<double>() );

 // dimensions and variable declarations

 netCDF::NcDim ndim = group.addDim( "num_vars" , x.size() );
 netCDF::NcDim mdim = group.addDim( "num_constraints" , A.size() );

 // Matrix stored in rmatbeg/rmatval/rmatind format
 // with lower and upper row bounds
 netCDF::NcDim nz = group.addDim( "nzcount" ,
				  std::accumulate( A.begin() , A.end() ,
						   Index( 0 ) ,
                                  [&]( Index acc , auto &&lf ) {
				     return( acc + lf.get_num_active_var() );
				  } ) );

 netCDF::NcVar rmatval = group.addVar( "rmatval" , netCDF::NcDouble() , nz );
 netCDF::NcVar rmatind = group.addVar( "rmatind" , netCDF::NcUint64() , nz );
 netCDF::NcVar rmatbeg = group.addVar( "rmatbeg" , netCDF::NcUint64() ,
				       mdim );
 netCDF::NcVar xlb     = group.addVar( "x_lb" , netCDF::NcDouble() , ndim );
 netCDF::NcVar xub     = group.addVar( "x_ub" , netCDF::NcDouble() , ndim );
 netCDF::NcVar rlb     = group.addVar( "row_lb" , netCDF::NcDouble() , mdim );
 netCDF::NcVar rub     = group.addVar( "row_ub" , netCDF::NcDouble() , mdim );
 netCDF::NcVar obj     = group.addVar( "obj" , netCDF::NcDouble() , ndim );
 netCDF::NcVar objoff  = group.addVar( "obj_offset" , netCDF::NcDouble() );

 // inserts into netcdf require an index variable:
 assert( sizeof( Index ) <= sizeof( size_t ) );
 std::vector<size_t> idx( 1 );

 // add objective function
 auto fl = static_cast<const LinearFunction*>( f.get_function() );
 assert( fl );
 auto cc = fl->get_v_var();

 double val[ 1 ] = { fl->get_constant_term() };
 objoff.putVar( val );

 Index i = 0;
 for( Index j = 0 ; j < cc.size() ; ++i ) {
  idx[ 0 ] = i;
  if( std::distance( x.data() ,
		     static_cast< const ColVariable * >( cc[ j ].first ) )
      == i )
   obj.putVar( idx , cc[ j++ ].second );
  else
   obj.putVar( idx , double( 0 ) );
  }

 for( ; i < x.size() ; ++i ) {
  idx[ 0 ] = i;
  obj.putVar( idx , double( 0 ) );
  }
 
 // add variable bounds
 for( uint64_t i = 0 ; i < x.size() ; i++ ) {
  idx[ 0 ] = i;
  xlb.putVar( idx , Bx[ i ].get_lhs() );
  xub.putVar( idx , Bx[ i ].get_rhs() );
  }

 // loop over constraints, filling matrix and rhs/lhs
 // Index offset = 0;
 unsigned long long offset = 0;
 for( Index i = 0 ; i < A.size() ; i++ ) {
  idx[ 0 ] = i;
  rlb.putVar( idx , A[ i ].get_lhs() );
  rub.putVar( idx , A[ i ].get_rhs() );
  rmatbeg.putVar( idx , offset );

  // loop over row nonzeros
  auto f = dynamic_cast<const LinearFunction*>( A[ i ].get_function() );
  auto cc = f->get_v_var();

  for( Index j = 0 ; j < cc.size() ; j++ ) {
   /* FIXME: could do one update for all of the entries instead */
   idx[ 0 ] = offset++;
   rmatind.putVar( idx , std::distance( x.data() ,
		      static_cast< const ColVariable * >( cc[ j ].first ) ) );
   rmatval.putVar( idx , cc[ j ].second );
   }
  } // for each constraint of A
 }  // end( SimpleMILPBlock::serialize )

/*--------------------------------------------------------------------------*/

void SimpleMILPBlock::guts_of_deserialize( netCDF::NcGroup & group ) 
{
 Block::deserialize( group );  // invoke method of Block

 // collect details to fill 'this'
 std::vector<size_t> idx( 1 );

 double dinf = Inf<double>();
 netCDF::NcGroupAtt inf = group.getAtt( "infinity" );
 if( ! inf.isNull() )
  inf.getValues( &dinf );
 
 size_t n , m , nz;
 ::deserialize_dim( group , "num_vars" , n , false );
 ::deserialize_dim( group , "num_constraints" , m , false );
 ::deserialize_dim( group , "nzcount" , nz , false );

 netCDF::NcVar xlb     = group.getVar( "x_lb" );
 netCDF::NcVar xub     = group.getVar( "x_ub" );
 netCDF::NcVar obj     = group.getVar( "obj" );
 netCDF::NcVar objoff  = group.getVar( "obj_offset" );
 netCDF::NcVar rlb     = group.getVar( "row_lb" );
 netCDF::NcVar rub     = group.getVar( "row_ub" );
 netCDF::NcVar rmatval = group.getVar( "rmatval" );
 netCDF::NcVar rmatind = group.getVar( "rmatind" );
 netCDF::NcVar rmatbeg = group.getVar( "rmatbeg" );

 x.resize( n );
 Bx.resize( n );
 A.resize( m );

 // insert variables and box constraints
 for( size_t i = 0 ; i < n ; i++ ) {
  Bx[ i ].set_variable( &x[ i ] , eNoMod );
  double val;
  idx[ 0 ] = i;
  xlb.getVar( idx , &val );
  if( val <= -dinf )
   val = -Inf<double>();
  Bx[ i ].set_lhs( val , eNoMod );
  xub.getVar( idx , &val );
  if( val >= dinf )
   val = Inf<double>();
  Bx[ i ].set_rhs( val , eNoMod );
  Bx[ i ].set_Block( this );
  x[ i ].set_Block( this );
  }

 // insert objective. FIXME: currently objective is non-sparse
 LinearFunction::v_coeff_pair p;
 for( size_t i = 0 ; i < n ; i++ ) {
  idx[ 0 ] = i;
  double val;
  obj.getVar( idx , &val );
  if( val )
   p.push_back( LinearFunction::coeff_pair( &x[ i ] , val ) );
  }

 double objoff_val;
 objoff.getVar( &objoff_val );
 f.set_function( new LinearFunction( std::move( p ) , objoff_val ) , eNoMod );
 f.set_Block( this );

 // insert constraints

 // row bounds and rows
 for( size_t j = 0 ; j < m ; j++ ) {
  double val;
  idx[ 0 ] = j;
  rlb.getVar( idx , &val );
  if( val <= -dinf )
   val = -Inf<double>();
  A[ j ].set_lhs( val , eNoMod );
  rub.getVar( idx , &val );
  if( val >= dinf )
   val = Inf<double>();
  A[ j ].set_rhs( val , eNoMod );

  // constraint rows:
  uint64_t offset;
  rmatbeg.getVar( idx , &offset );
  size_t row_nzcnt;
  if( j == m - 1 )
   row_nzcnt = nz-offset;
  else {
   uint64_t next_offset;
   rmatbeg.getVar( { j + 1 } , &offset );
   row_nzcnt = next_offset - offset;
   }

  p.resize( row_nzcnt );
  for( int h = 0 ; h < row_nzcnt ; h++ ) {
   uint64_t i;
   double coeff;
   idx[ 0 ] = offset + h;
   rmatind.getVar( idx , &i );
   rmatval.getVar( idx , &coeff );
   if( ( i >= this->get_x().size() ) )
    throw( invalid_argument( "invalid variable idx" ) );

   p[ h ].first = &( this->get_x()[ i ] );
   p[ h ].second = coeff;
   }

  A[ j ].set_function( new LinearFunction( std::move( p  ) , 0 ) , eNoMod );
  A[ j ].set_Block( this );
  }

 // set the abstract representation - - - - - - - - - - - - - - - - - - - - -
 // explicitly clear static Constraint and Variable: this should not be
 // needed, but a single SimpleMILPBlock could be loaded multiple times,
 // in which case the abstract representation would pile up
 reset_static_constraints();
 reset_static_variables();
 reset_objective();

 // this is not required by the standard, but a SimpleMILPBlock only has the
 // "abstract representation", so it must always be created

 set_objective( & f );
 add_static_variable( x );
 add_static_constraint( A );
 add_static_constraint( Bx );

 // issue the NBModification - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );

 }  // end( SimpleMILPBlock::guts_of_deserialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File SimpleMILPBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
