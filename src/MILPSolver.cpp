/*--------------------------------------------------------------------------*/
/*-------------------------- File MILPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MILPSolver class.
 *
 * \version 0.10
 *
 * \date 22 - 12 - 2016
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <queue>
#include <functional>   // std::bind
#include <chrono>
#include <cstdlib>

#include "MILPSolver.h"
#include "Block.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

extern int t_pc;


int mycallback (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
 
 return( static_cast<MILPSolver *>( cbhandle )->Callback( env , cbdata , wherefrom , useraction_p ) );
 }


/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::count_const (LinearConstraint & lconst, int & count){
count++;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_var (ColVariable & lvar, int & count, int & count2){
count++;
//cout<<"[ 57 ; " << count << " ]";
count2 = count2 + lvar.get_num_active_const();
}

/*--------------------------------------------------------------------------*/

int MILPSolver::ind_const (LinearConstraint * p_const){

int i=0;

auto it1 = lower_bound (v_s_const_int.begin(), v_s_const_int.end(), p_const, [&](const_int pair, LinearConstraint * pconst) {
    return pair.first < pconst;
});

if (it1 == v_s_const_int.end()) {
         it1 = (v_s_const_int.rbegin()+1).base();
       }
       else if (it1 != v_s_const_int.begin() && it1->first > p_const) {
           --it1;
       }

  if (it1 != v_s_const_int.end())
	 i = it1->second + std::distance (it1->first, p_const);
  else
    throw( std::invalid_argument(
            "Current Constraint is not defined in the CPLEX Coeff Matrix" ) );

  return i;


}

/*--------------------------------------------------------------------------*/

LinearConstraint * MILPSolver::const_ind (int i){

LinearConstraint * p_const;
auto it1 = lower_bound (v_int_s_const.begin(), v_int_s_const.end(), i, [&](int_const pair, int i) {
    return pair.first < i;
});

if (it1 == v_int_s_const.end()) {
         it1 = (v_int_s_const.rbegin()+1).base();
       }
       else if (it1 != v_int_s_const.begin() && it1->first > i) {
           --it1;
       }

  if (it1 != v_int_s_const.end())
	 p_const = it1->second;
  else
    throw( std::invalid_argument(
            "Current index is not defined in the CPLEX Coeff Matrix" ) );

  return p_const;

}

/*--------------------------------------------------------------------------*/

int MILPSolver::ind_var (ColVariable * p_var){

int i;
auto it1 = lower_bound (v_s_var_int.begin(), v_s_var_int.end(), p_var, [&](var_int pair, ColVariable * pvar) {
    return pair.first < pvar;
});

if (it1 == v_s_var_int.end()) {
         it1 = (v_s_var_int.rbegin()+1).base();
       }
       else if (it1 != v_s_var_int.begin() && it1->first > p_var) {
           --it1;
       }
  if (it1 != v_s_var_int.end()){
	 i = it1->second +  std::distance (it1->first, p_var);
    }
  else
    throw( std::invalid_argument(
            "Current Variable is not defined in the CPLEX Coeff Matrix" ) );

  return i;

}

/*--------------------------------------------------------------------------*/

ColVariable * MILPSolver::var_ind (int i){

ColVariable * p_var;
auto it1 = lower_bound (v_int_s_var.begin(), v_int_s_var.end(), i, [&](int_var pair, int i) {
    return pair.first < i;
});

if (it1 == v_int_s_var.end()) {
         it1 = (v_int_s_var.rbegin()+1).base();
       }
       else if (it1 != v_int_s_var.begin() && it1->first > i) {
           --it1;
       }

  if (it1 != v_int_s_var.end())
	 p_var = it1->second; 
  else
    throw( std::invalid_argument(
            "Current index is not defined in the CPLEX Coeff Matrix" ) );

  return p_var;

}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_s_const (LinearConstraint & lconst, char * sense, double * rhs,  
int & first, int & i){

/* In order to define the type of the operation of the constraint we need to
   check and compare the lhs and the rhs of the LinearConstraint due to the 
   fact that is built in the following form:
         LHS <= ( some function from Variables to reals ) <= RHS
*/
if( lconst.get_lhs() == lconst.get_rhs() ){ //equality
sense[i] = ('E');
rhs[i] = (lconst.get_rhs());//setting the rhs
}
else if( lconst.get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
sense[i] = ('L');
rhs[i] = (lconst.get_rhs());//setting the rhs, which is equal to lhs
}
else if( lconst.get_rhs() == Inf<double>() ){  // inequality (greate/equal)
sense[i] = ('G');
rhs[i] = (lconst.get_lhs());//setting the rhs, which is equal to rhs
}

//check if we are in dynamic or static part of the problem

if(first==0) {
v_s_const_int.push_back( std::make_pair ( &lconst, i ) );
v_int_s_const.push_back( std::make_pair ( i, &lconst ) );
}

first++; //increasing the counter, since we are in static part
i++;

}

/*--------------------------------------------------------------------------*/





void MILPSolver::scan_d_const (LinearConstraint & lconst, char * sense, 
double * rhs, int & i){

/* In order to define the type of the operation of the constraint we need to
   check and compare the lhs and the rhs of the LinearConstraint due to the 
   fact that is built in the following form:
         LHS <= ( some function from Variables to reals ) <= RHS
*/
if( lconst.get_lhs() == lconst.get_rhs() ){ //equality
sense[i] = ('E');
rhs[i] = (lconst.get_rhs());//setting the rhs
}
else if( lconst.get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
sense[i] = ('L');
rhs[i] = (lconst.get_rhs());//setting the rhs, which is equal to lhs
}
else if( lconst.get_rhs() == Inf<double>() ){  // inequality (greate/equal)
sense[i] = ('G');
rhs[i] = (lconst.get_lhs());//setting the rhs, which is equal to rhs
}

v_d_const_int.push_back( std::make_pair ( &lconst, i ) );
v_int_d_const.push_back( std::make_pair ( i, &lconst ) );
i++;

}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_s_var (ColVariable & lvar, int    *matbeg, int    *matcnt, int    *matind, double *matval, /*double *objective,
    double *q_objective,*/ double *lb,  double *ub, char *xctype,  int & first, int & i){

if(first==0) {
v_s_var_int.push_back( std::make_pair ( &lvar, i ) );
v_int_s_var.push_back( std::make_pair ( i, &lvar ) );
//cout<<endl<<"For the BIEWS = " << lvar.get_ub();
}


/* setting the lower and upper bounds of the variable to the corresponing
  vector */
if (lvar.get_lb() == -Inf<double>() )
    lb[i] = -CPX_INFBOUND;
else
    lb[i] = (lvar.get_lb() );

if (lvar.get_ub() == Inf<double>() )
    ub[i] = CPX_INFBOUND;
else 
    ub[i] = (lvar.get_ub() );
// seting the type of the variable 

 switch( lvar.get_type() ) {
  case( ColVariable::integer ):      xctype[i] = ('I'); break;
  case( ColVariable::binary ):       xctype[i] = ('B'); break;
  case( ColVariable::continuous ):   xctype[i] = ('C');
  }


// seting the Variable's no. of non-zero elements in the CPLEX coeff matrix
matcnt[i] = (lvar.get_num_active_const());

// seting the index of index of the beginning of column in the CPLEX coefficient matrix
if (i==0)
matbeg[i] = 0;
else
matbeg[i] = matbeg[i-1] + matcnt[i-1];

/*we need to locate where lies the index in matval and matind for the corres-
  ponding variable
*/
if (i == 0) indexed[i] =0;
else if (i>0) indexed[i] = indexed[i-1] + matcnt[i-1];



/* seting the coefficients and the indexes of the corresponging rows in the 
   CPLEX coeff matrix */

for(int j = 0 ; j < lvar.get_num_active_const() ; j++){

//retreive each constraint that the variable is active
LinearConstraint * p_const = dynamic_cast<LinearConstraint *> ( lvar.get_active_const( j ) );

//locating the corresponding coefficient of the variable in the constraint
auto it = find_if (p_const->get_v_var()->begin(),p_const->get_v_var()->end(), [&](LinearConstraint::coeff_pair pair) {
    return pair.first == &lvar;
});


  if (it != p_const->get_v_var()->end())
    matval[ /*index*/ indexed[i] + j ] = (it->second);
  else
    throw( std::invalid_argument(
            "Variable is not active in the examined Constraint" ) );

/* Passing the index of the constraint to the corresponding vector matind, note that here
   we need to check if we are in the dynamic or static part of the problem, since we will
   have to use the vector of pairs in order to locate the index of the examined constraint
   in the CPLEX coeff matrix
*/

//locating the index of the constraint in CPLEX
 matind[/*index*/ indexed[i] + j] = ind_const(p_const);

} //for j

first++; //increasing the counter, since we are in static part
i++;

}

void MILPSolver::scan_l_of (LinearObjectiveFunction * l_obj_fun, double *objective, double *q_objective){

    int k;
    if (l_obj_fun->get_v_var() != nullptr){
    	for( auto &el : *l_obj_fun->get_v_var() ){
        	k = ind_var(el.first);
        	objective[k] = (el.second);
    	}
    }
}
/*--------------------------------------------------------------------------*/
void MILPSolver::scan_q_of (DQuadObjectiveFunction * q_obj_fun, double *objective, double *q_objective){

      int k;
    if (q_obj_fun->get_v_var() != nullptr){
        for( auto &el : *q_obj_fun->get_v_var() ){
               k = ind_var(el.first);
               objective[k] = (el.second);         
        }
    }// if case

    if (q_obj_fun->get_v_q_var() != nullptr){
        for( auto &el : *q_obj_fun->get_v_q_var() ){
            k = ind_var(el.first);
            q_objective[k] = (el.second);
        }
    }// if case
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_d_var (ColVariable & lvar,  int    *matbeg, int    *matcnt, int    *matind, double *matval, double *objective,
  double *q_objective, double *lb,  double *ub, char *xctype, int & i ){

v_d_var_int.push_back( std::make_pair ( &lvar, i ) );
v_int_d_var.push_back( std::make_pair ( i, &lvar ) );


/* setting the lower and upper bounds of the variable to the corresponing
  vector */

lb[i] = lvar.get_lb() ;
ub[i] = lvar.get_ub() ;

// seting the type of the variable 

 switch( lvar.get_type() ) {
  case( ColVariable::integer ):      xctype[i] = 'I'; break;
  case( ColVariable::binary ):       xctype[i] = 'B'; break;
  case( ColVariable::continuous ):   xctype[i] = 'C';
  }

// seting the index of variable to the CPLEX coefficient matrix
matbeg[i] = i;

// seting the Variable's no. of non-zero elements in the CPLEX coeff matrix
matcnt[i] = (lvar.get_num_active_const());


/*we need to locate where lies the index in matval and matind for the corres-
  ponding variable
*/
int index =0;
if (i>0){

for (int j = 0; j < i ; j++){

index = index + matcnt[j];

}

}


/* seting the coefficients and the indexes of the corresponging rows in the 
   CPLEX coeff matrix */
for(int j = 0 ; j < lvar.get_num_active_const() ; j++){

//retreive each constraint that the variable is active
Constraint * p_cnst = lvar.get_active_const( j );

LinearConstraint * p_const = dynamic_cast<LinearConstraint *>(p_cnst);

//locating the corresponding coefficient of the variable in the constraint
auto it = find_if (p_const->get_v_var()->begin(),p_const->get_v_var()->end(), [&](LinearConstraint::coeff_pair pair) {
    return pair.first == &lvar;
});


  if (it != p_const->get_v_var()->end())
    matval[index + j] = (it->second);
  else
    throw( std::invalid_argument(
            "Variable is not active in the examined Constraint" ) );
 
/* Passing the index of the constraint to the corresponding vector matind, note that here
   we need to check if we are in the dynamic or static part of the problem, since we will
   have to use the vector of pairs in order to locate the index of the examined constraint
   in the CPLEX coeff matrix
*/

//locating the index of the constraint in CPLEX
auto it1 = find_if (v_d_const_int.begin(), v_d_const_int.end(), [&](const_int pair) {
    return pair.first == p_const;
});

  if (it1 != v_d_const_int.end())
    matind[index +j] = (it->second);
  else
    throw( std::invalid_argument(
            "Current Constraint is not defined in the CPLEX Coeff Matrix" ) );

} //for j 


// seting the coefficient of the variable for the objective function
//Checking the type of the objective function
if( lvar.get_Block()->get_objective_function().type() ==  typeid(LinearObjectiveFunction *) ){

	LinearObjectiveFunction * obj_fun = boost::any_cast< LinearObjectiveFunction * >(lvar.get_Block()->get_objective_function());
	if(obj_fun->is_active(&lvar) == true){
		//in case variable is active in obj.function we retreive the corresponding coefficient
		auto it = find_if (obj_fun->get_v_var()->begin(),obj_fun->get_v_var()->end(), [&](LinearObjectiveFunction::coeff_pair pair) {
			return pair.first == &lvar;
	        	});
	        objective[i] = (it->second);
        }
	else { //in case variable is not part of obj.function we set coefficient to zero
		objective[i] = 0;
	}

}// if linear objective function
else if( lvar.get_Block()->get_objective_function().type() ==  typeid(DQuadObjectiveFunction *) ){

	DQuadObjectiveFunction * q_obj_fun = boost::any_cast<DQuadObjectiveFunction * >(lvar.get_Block()->get_objective_function());

	if(q_obj_fun->is_active(&lvar) == true){
		//in case variable is active in obj.function we retreive the corresponding coefficient
		auto it = find_if (q_obj_fun->get_v_var()->begin(),q_obj_fun->get_v_var()->end(), [&](LinearObjectiveFunction::coeff_pair pair) {
			return pair.first == &lvar;
	        	});
	        objective[i] = (it->second);
        }
	else { //in case variable is not part of obj.function we set coefficient to zero
		objective[i] = 0;
	}


	if(q_obj_fun->is_q_active(&lvar) == true){
	//in case variable is active in obj.function we retreive the corresponding coefficient
		auto it = find_if (q_obj_fun->get_v_q_var()->begin(),q_obj_fun->get_v_q_var()->end(), [&](LinearObjectiveFunction::coeff_pair pair) {
		    return pair.first == &lvar;
		});
        	q_objective[i] = (it->second);
	}
	else { //in case variable is not part of obj.function we set coefficient to zero
		q_objective[i] = 0;
	}

}// if dquad objectiv function		
else throw( std::invalid_argument( "Unknown type of Objective Function" ) );

i++;

}

/*--------------------------------------------------------------------------*/

void MILPSolver::set_Block( Block *block ){

  if( f_Block ) {                        // was attached to some other Block
  f_Block->unregister_Solver( this );  // no more so
  CPXfreeprob (env, &milp);
  CPXcloseCPLEX (&env);
  milp = 0;
  env = 0;
}

 f_Block = block;                      // this is the new block now
 f_Block->register_Solver( this );     // register to it

  // Initialisation the CPLEX environement and problem

  int status ;
  env = CPXopenCPLEX(&status);
  milp  = CPXcreateprob(env,&status,"MILPCPX");



  /**********Passing all the data of the Block to the CPLEX Problem**********/ 
  /* Following a Breadth First Search we proceed with scanning the received
     Block and all of each corresponding children if any, in order to populate 
     the CPLEX data needed to define the corresponding MILP. This is done by 
     the use of the following two CPLEX commands: 
     
     1)  int CPXcopylp(CPXCENVptr env, CPXLPptr lp, int numcols, int numrows, 
         int objsense, double const * objective, double const * rhs, char const * 
         sense, int const * matbeg, int const * matcnt, int const * matind, 
         double const * matval, double const * lb, double const * ub, double 
         const * rngval)
     2)  int CPXcopyctype(CPXCENVptr env, CPXLPptr lp, char const * xctype) 
   
    */

/* The following vectors are used in order to retreive from the Block all the
   information needed in order to construct and pass to CPLEX the corresponding
   coeff matrix and all the rest of the information needed 

 Where:
   -objective: refers to a vector with the total coefficients of the objective
    function
   -rhs: refers to a vector with the rhs of all the constraints
   -sense: refers to a vector with the type of all the constraints
   -matbeg: refers to a vector with the indices of all the variables
   -matcnt: refers to a vector with the non-zero elements to correspond to each
    variable, keeping track of the indicing from matbeg
   -matval: refers to a vector with the coefficients that correspond to each
    variable, keeping track of the indicing from matbeg
   -matind: refers to a vector that asociates the coefficients that correspond to each
    variable with the corresponding row that they refer to
   -lb,ub: refers to two array with the lower and upper bounds respectively 
    of each variable
   -xctype: refers to a vector that describes the type of each variable
   -first: is an integer used to denote if the examined constraint/variable is
    static or dynamic and in case is static if its the first element of the exa-
    mined type of variable.
   -obj_fun: pointer to the objective function of the Block
   -numcols,numrows,nzelements: total number of variables, constraints and non-zero
                                elements  of the Block respectively

*/

numrows = 0;
numcols = 0;
nzelements = 0;

/* we need to scan through the whole block in order to properly define the values of 
   total number of variables, constraints and non-zero elements
*/

 std::queue<Block *> Q; //creating the queue

	Q.push(f_Block);//passing the root, i.e. the father block


	while(!Q.empty())//iterating for the father block and all children
	{
		Block * q_Block = Q.front(); //block to be examined
		Q.pop(); //take out from the queue the examined block


        
		for (int i = 0; i < q_Block->get_nested_Blocks().size(); ++i)
		{
                       //inserting all the children of the examined block
			Q.push(q_Block->get_nested_Blocks()[i]); 
		}
		
		//We scan all the static constraints/rows
        
		for (int i =0; i< q_Block->get_static_constraints().size(); i++){
	           auto f1 = std::bind( &MILPSolver::count_const, this, placeholders::_1, std::ref(numrows) );
		   un_any_const_static(q_Block->get_static_constraints()[i],f1,un_any_type<LinearConstraint>() );

		}

		//We scan all the dynamic constraints/rows
		for (int i =0; i< q_Block->get_dynamic_constraints().size(); i++){
		cout<<endl<<"num_rows before = " << numrows;
		auto f1 = std::bind( &MILPSolver::count_const, this, placeholders::_1, std::ref(numrows) );
        	un_any_const_dynamic(q_Block->get_dynamic_constraints()[i],f1,un_any_type<LinearConstraint>() );
		cout<<endl<<"num_rows after = " << numrows;
		}

		//We scan all the static Variables/columns
        
		for (int i =0; i< q_Block->get_static_variables().size(); i++){
		auto f1 = std::bind( &MILPSolver::count_var, this, placeholders::_1, std::ref(numcols), std::ref(nzelements) );
		un_any_const_static(q_Block->get_static_variables()[i],f1,un_any_type<ColVariable>() );

		}

		//We scan all the dynamic Variables/columns
		for (int i =0; i< q_Block->get_dynamic_variables().size(); i++){
		auto f1 = std::bind( &MILPSolver::count_var, this, placeholders::_1, std::ref(numcols), std::ref(nzelements) );
        	un_any_const_dynamic(q_Block->get_dynamic_variables()[i],f1,un_any_type<ColVariable>() );
		}


	}//while loop

  int    *matbeg = new int[numcols];
  int    *matcnt = new int[numcols];
  int    *matind = new int[ nzelements ] ;
  double *matval = new double[ nzelements ];
         indexed.resize(numcols);
  double *rhs    = new double[ numrows ];
  char   *sense  = new char[ numrows ];
  double *objective    = new double[ numcols ];
  double *q_objective  = new double[ numcols ];
  for (int i = 0 ; i < numcols ; i++) { objective[i] = 0; q_objective[i] = 0; }
  double *lb    = new double[ numcols ];
  double *ub    = new double[ numcols ];
  char   *xctype  = new char[ numcols ];

/**Proceeding with  Breadth First Scan/Search of the Block and its children**/
  
	Q.push(f_Block);//passing the root, i.e. the father block

	int var =0;
        int col = 0; 
        /* counters used to locate the index in the arrays of constraints and
           variables that need to be populated with the Block data
	*/


   	while(!Q.empty())//iterating for the father block and all children
	{
		Block * q_Block = Q.front(); //block to be examined
		Q.pop(); //take out from the queue the examined block
		
		for (int i = 0; i < q_Block->get_nested_Blocks().size(); ++i)
		{
                       //inserting all the children of the examined block
			Q.push(q_Block->get_nested_Blocks()[i]); 
		}
        	/*Scanning and passing all the data of the examined block  
                  to the corresponding CPLEX data. We do this by first scan-
                  ning the static part of the problem and then the dynamic
                  part. This is because we want to have an order in the col-
                  muns and rows of the CPLEX coeff matrix where the static
                  part is being followed by the dynamic one */
		
		//We scan all the static constraints
		for (int i =0; i< q_Block->get_static_constraints().size(); i++){
			int first = 0; //variable used to locate the first element of each type of constraint 
			auto f1 = std::bind( &MILPSolver::scan_s_const, this, placeholders::_1, std::ref(sense), std::ref(rhs), std::ref(first), std::ref(col) );
          		un_any_const_static(q_Block->get_static_constraints()[i],f1,un_any_type<LinearConstraint>() );
		}

		//We scan all the dynamic constraints
		for (int i =0; i< q_Block->get_dynamic_constraints().size(); i++){
			auto f1 = std::bind( &MILPSolver::scan_d_const, this, placeholders::_1, std::ref(sense), std::ref(rhs), std::ref(col) );
	        	un_any_const_dynamic(q_Block->get_dynamic_constraints()[i],f1,un_any_type<LinearConstraint>() );
		}

	}//while loop



     std::sort( v_s_const_int.begin() , v_s_const_int.end() );
     std::sort( v_int_s_const.begin() , v_int_s_const.end() );
     std::sort( v_d_const_int.begin() , v_d_const_int.end() );
     std::sort( v_int_d_const.begin() , v_int_d_const.end() );


	Q.push(f_Block);//passing the root, i.e. the father block


	while(!Q.empty())//iterating for the father block and all children
	{

		Block * q_Block = Q.front(); //block to be examined
		Q.pop(); //take out from the queue the examined block
		
		for (int i = 0; i < q_Block->get_nested_Blocks().size(); ++i)
		{
                       //inserting all the children of the examined block
			Q.push(q_Block->get_nested_Blocks()[i]); 
		}
        	/*Scanning and passing all the data of the examined block  
                  to the corresponding CPLEX data. We do this by first scan-
                  ning the static part of the problem and then the dynamic
                  part. This is because we want to have an order in the col-
                  muns and rows of the CPLEX coeff matrix where the static
                  part is being followed by the dynamic one */

/*cout<<endl<<"set_block_714 = " << q_Block->get_static_variables().size();
cout<<endl<<"//";*/

		//We scan all the static Variables
		for (int i =0; i< q_Block->get_static_variables().size(); i++){
		int first = 0; //variable used to locate the first element of each type of variable
     /*   cout<<endl<<"In 720_For i = " << i ;
        cout<<endl<<"//";*/
        std::chrono::time_point<std::chrono::system_clock> t0, t1;
        t0 = std::chrono::system_clock::now();
		auto f1 = std::bind( &MILPSolver::scan_s_var, this, placeholders::_1, std::ref(matbeg), std::ref(matcnt), 
				     std::ref(matind), std::ref(matval), /*std::ref(objective), std::ref(q_objective), */
                                     std::ref(lb), std::ref(ub), std::ref(xctype), std::ref(first), std::ref(var) );
        	un_any_const_static(q_Block->get_static_variables()[i],f1,un_any_type<ColVariable>() );
        t1= std::chrono::system_clock::now();
        std::chrono::duration<double> t1_fin = t1-t0;
//        cout << " . . . process lasted: " << t1_fin.count();
		}

/*cout<<endl<<"set_block_726";
cout<<endl<<"//";*/

		//We scan all the dynamic Variables
		for (int i =0; i< q_Block->get_dynamic_variables().size(); i++){
		auto f1 = std::bind( &MILPSolver::scan_d_var,this, placeholders::_1, std::ref(matbeg), std::ref(matcnt), 
				     std::ref(matind), std::ref(matval), std::ref(objective), std::ref(q_objective), 
                                     std::ref(lb), std::ref(ub), std::ref(xctype), std::ref(var) );
        	un_any_const_dynamic(q_Block->get_dynamic_variables()[i],f1,un_any_type<ColVariable>() );
		}

/*cout<<endl<<"set_block_740";
cout<<endl<<"//";*/

	}//while loop

/*cout<<endl<<"set_block_748";
cout<<endl<<"//";*/
		obj_type=0;
		  
		//Checking the type of the objective function
		if( f_Block->get_objective_function().type() ==  typeid(LinearObjectiveFunction * ) ) {
		LinearObjectiveFunction * obj_fun = boost::any_cast<LinearObjectiveFunction * >(f_Block->get_objective_function());
		switch( obj_fun->get_type() ) { //note that we check the objective type from the father Block
		case( ObjectiveFunction::eMax ):   obj_type=-1; break;
		case( ObjectiveFunction::eMin ):   obj_type=1;
			}
		}
		else if( f_Block->get_objective_function().type() ==  typeid( DQuadObjectiveFunction * ) ){
		DQuadObjectiveFunction * q_obj_fun = boost::any_cast<DQuadObjectiveFunction * >(f_Block->get_objective_function());
		switch( q_obj_fun->get_type() ) { //note that we check the objective type from the father Block
		case( ObjectiveFunction::eMax ):   obj_type=-1; break;
		case( ObjectiveFunction::eMin ):   obj_type=1;
			}
		}
		else throw( std::invalid_argument( "Unknown type of Objective Function" ) );

        

/*cout<<endl<<"set_block_770";
cout<<endl<<"//";*/
 

   /* Order all the different vectors of pairs in ascending order based
   on the adress of the constraints/variables or the index of CPLEX
   rows/columns respectively. */

     std::sort( v_s_var_int.begin() , v_s_var_int.end() );
     std::sort( v_int_s_var.begin() , v_int_s_var.end() );
     std::sort( v_d_var_int.begin() , v_d_var_int.end() );
     std::sort( v_int_d_var.begin() , v_int_d_var.end() );
     



	Q.push(f_Block);//passing the root, i.e. the father block
     /*   cout<<endl<<"In Building of";
        cout<<endl<<"//";*/
        std::chrono::time_point<std::chrono::system_clock> t00, t11;
        t00 = std::chrono::system_clock::now();
    int sel = 0;
	while(!Q.empty())//iterating for the father block and all children
	{
		Block * q_Block = Q.front(); //block to be examined
		Q.pop(); //take out from the queue the examined block
      
		for (int i = 0; i < q_Block->get_nested_Blocks().size(); ++i)
		{
                       //inserting all the children of the examined block
			Q.push(q_Block->get_nested_Blocks()[i]); 
		}
		if( q_Block->get_objective_function().type() ==  typeid(LinearObjectiveFunction * ) ) {
		    LinearObjectiveFunction * l_obj_fun = boost::any_cast<LinearObjectiveFunction * >(q_Block->get_objective_function());
		    //We scan all objective functions
            scan_l_of (l_obj_fun, objective, q_objective);
        }

		if( q_Block->get_objective_function().type() ==  typeid(DQuadObjectiveFunction * ) ) {
		    DQuadObjectiveFunction * q_obj_fun = boost::any_cast<DQuadObjectiveFunction * >(q_Block->get_objective_function());
		    //We scan all objective functions
            scan_q_of (q_obj_fun, objective, q_objective);
        }

        sel++;
    
	}//while loop
   

        t11= std::chrono::system_clock::now();
        std::chrono::duration<double> t11_fin = t11-t00;
     //   cout << " . . . process lasted: " << t11_fin.count();

/* Due to the fact that CPLEX receives as arguments arrays and our data
   has been stored with the usage of vectors (due to the fact that they 
   can be dynamically resized) we need to create temporary arrays to which
   we will pass the vectors by pointers.
*/

  //creation of the problem
CPXcopylp(env,milp,numcols,numrows,obj_type,objective,rhs,sense,
	    matbeg,matcnt,matind,matval,lb,ub,NULL);



 // copy the objective function coeff matrix
  if( f_Block->get_objective_function().type() ==  typeid(DQuadObjectiveFunction *) )
       CPXcopyqpsep(env,milp,q_objective);
       
CPXcopyctype(env,milp,xctype);

  cuts=t_pc;

  if(cuts==1){
	
   // variables are referred according to original model in callback function
   CPXsetintparam( env , CPX_PARAM_MIPCBREDLP , CPX_OFF );
   CPXsetintparam( env , CPX_PARAM_PRELINEAR , CPX_OFF );

   status = CPXsetusercutcallbackfunc        ( env, mycallback ,  this );
   status = CPXsetlazyconstraintcallbackfunc ( env, mycallback ,  this );
  }

  delete []matbeg;
  delete []matcnt;
  delete []matind;
  delete []matval;

  delete []rhs;
  delete []sense;
  delete []objective;
  delete []q_objective;
  delete []lb;
  delete []ub;
  delete []xctype;



}


/*--------------------------------------------------------------------------*/

Solver::OFValue MILPSolver::get_lb( void ){

  OFValue lower_bound;


if(obj_type==1){ // for minimization problem

  if     (sol_status == kUnbounded)  lower_bound = - Inf<OFValue>() ;
  else if(sol_status == kInfeasible) lower_bound =   Inf<OFValue>() ;
  else                               CPXgetbestobjval (env, milp, &lower_bound);
}

else if (obj_type==-1){ // for maximization problem

  if     (sol_status == kUnbounded)  lower_bound =   Inf<OFValue>() ;
  else if(sol_status == kInfeasible) lower_bound = - Inf<OFValue>() ;
  else                               CPXgetbestobjval (env, milp, &lower_bound);

}

 return( lower_bound );


}

 /*--------------------------------------------------------------------------*/

Solver::OFValue MILPSolver::get_ub( void ){

    OFValue upper_bound;


if(obj_type==1){ // for minimization problem

  if     (sol_status == kUnbounded)  upper_bound = - Inf<OFValue>() ;
  else if(sol_status == kInfeasible) upper_bound =   Inf<OFValue>() ;
  else                               CPXgetobjval (env, milp, &upper_bound);
}

else if (obj_type==-1){ // for maximization problem

  if     (sol_status == kUnbounded)  upper_bound =   Inf<OFValue>() ;
  else if(sol_status == kInfeasible) upper_bound = - Inf<OFValue>() ;
  else                               CPXgetobjval (env, milp, &upper_bound);

}

 return( upper_bound );


}

/*--------------------------------------------------------------------------*/

void MILPSolver::set_par( const int par , const int value )
{
 switch( par ) {
  case( kLogVerb ): f_log_verb = value; CPXsetintparam(env, CPX_PARAM_SCRIND, value); break; //1
  case( kMaxSol ):  f_max_sol = value;  CPXsetintparam(env, CPXPARAM_MIP_Limits_Solutions, value) ; break;
  default:                              CPXsetintparam(env, par-kLastAlgPar,value ); break;
  }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::set_par( const int par , const double value )
{
 switch( par ) {
  case( kMaxTime ):  f_max_time = value;  CPXsetdblparam(env, CPXPARAM_TimeLimit,value ); break; //10
  case( kRelAcc ):   f_rel_acc = value;   /*CPXsetdblparam(env, ,value );*/ break;
  case( kAbsAcc ):   f_abs_acc = value;   /*CPXsetdblparam(env, ,value );*/ break;
  case( kUpCutOff ): f_up_cutoff = value; CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_UpperCutoff,value ); break;
  case( kLwCutOff ): f_lw_cutoff = value; CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_LowerCutoff,value ); break;
  case( kRAccSol ):  f_r_acc_sol = value; CPXsetdblparam(env, CPXPARAM_MIP_Pool_RelGap,value ); break;
  case( kAAccSol ):  f_a_acc_sol = value; CPXsetdblparam(env, CPXPARAM_MIP_Pool_AbsGap,value ); break;
  case( kFAccSol ):  f_f_acc_sol = value; /*CPXsetdblparam(env, ,value );*/ break;
  default:                                CPXsetdblparam(env, par-kLastAlgPar,value ); break;
  }
 }

void MILPSolver::set_par( const int par , const long value )
{
 switch( par ) {
  case( kMaxIter ): f_max_iter = value; CPXsetlongparam(env, CPXPARAM_MIP_Limits_Nodes, value); break; // 0 or max
  default:                              CPXsetlongparam(env, par-kLastAlgPar,value ); break;
  }
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::solve( ){

  int status;
/*cout<<endl<<"996 with env = " << env;
cout<<endl<<"//";

cout<<endl<<"999 with milp = " << milp;
cout<<endl<<"//";*/

  CPXwriteprob (env, milp, "mipex1.lp", NULL);



  CPXmipopt(env,milp);

  nodes = CPXgetnodecnt (env, milp);

  get_var_solution();

  status = CPXgetstat(env,milp);

// obtaining appropriate sol_status
     if( status == 101 || 102 || 104 )   sol_status = kOK;
else if( status == 118  )   sol_status = kUnbounded;
else if( status == 103  )   sol_status = kInfeasible;
//else if( status ==  ||  )   sol_status = kLowPrecision;
else if( status == 107 || 108  )   sol_status = kStopTime;
else if( status == 105 || 106  )   sol_status = kStopIter;
else if( status == 109 || 110 )   sol_status = kError;
  
sol_status = status;

     return( sol_status );

}

int MILPSolver::Callback (CPXCENVptr env, void *cbdata, int wherefrom, int *useraction_p) {

//pass to Block the solution of current node
  //cout<<endl<<"numcols = " << numcols;
  double *tmpx = new double [numcols];

   *useraction_p = CPX_CALLBACK_DEFAULT; 

  CPXgetcallbacknodex( env , cbdata , wherefrom, tmpx,0,numcols-1);//obtain solution from CPLEX
/*   for(int i = 0 ; i < numcols ; i++){
 	cout<< " [ " << tmpx[i] << " ] " ;
   }*/
  //Retrieve Solution for all the Variables with respect to the order
  int col=0;

 std::queue<Block *> Q; //creating the queue

	Q.push(f_Block);//passing the root, i.e. the father block


	while(!Q.empty())//iterating for the father block and all children
	{
		Block * q_Block = Q.front(); //block to be examined
		Q.pop(); //take out from the queue the examined block


		for (int i = 0; i < q_Block->get_nested_Blocks().size(); ++i)
		{
                       //inserting all the children of the examined block
			Q.push(q_Block->get_nested_Blocks()[i]); 
		}

  		// Passing values for all static Variables
		for (int i =0; i< q_Block->get_static_variables().size(); i++){
		auto f1 = std::bind( &MILPSolver::set_var_value, this, placeholders::_1, std::ref(tmpx), std::ref(col) );
        	un_any_const_static(q_Block->get_static_variables()[i],f1,un_any_type<ColVariable>() );
		}

		//We scan all the dynamic Variables
		for (int i =0; i< q_Block->get_dynamic_variables().size(); i++){
		auto f1 = std::bind( &MILPSolver::set_var_value, this, placeholders::_1, std::ref(tmpx), std::ref(col) );
        	un_any_const_dynamic(q_Block->get_dynamic_variables()[i],f1,un_any_type<ColVariable>() );
		}

		q_Block->generate_dynamic_constraints( ); // we call this method where is expected to be written the mechanism
							 // for constructing new constraints
	}//while loop

//checking if new constraints have been creating by looking at the list of modifications
if(v_mod.size() > 0){

	int    mynzcnt = 0;
	double myrhs = 0;
	char sense = ('L');
		cout<<endl<<"New Modification received to Solver its size is: " << v_mod.size();
		cout<<endl<<"//";
	for (auto it = v_mod.begin() ; it != v_mod.end() ; it++ ){

		 //B *b = dynamic_cast< B* >( doSomething.get() );
		
		shared_ptr<BlockModificationAD> dmod = dynamic_pointer_cast<BlockModificationAD > (*it);
		
		if (dmod){
		
		//if(dmod->f_type == BlockModificationAD::eAddConst){ // we need to add all attached new constraints to CPLEX
                    		cout<<endl<<"Beefore stuff ";
				cout<<endl<<"//";
			/*std::list<LinearConstraint>*/ auto v_cuts = boost::any_cast< std::list<LinearConstraint> * >( dmod->whc_list ); 
                    		cout<<endl<<"No. of Constraints to be added: " << v_cuts->size();
			for (auto ct = v_cuts->begin() ; ct != v_cuts->end() ; ct++){
                    		cout<<endl<<"suppa";
				//pass the new cuts to CPLEX
				mynzcnt = ct->get_num_active_var();
				if( ct->get_lhs() == ct->get_rhs() ){ //equality
				sense = ('E');
				myrhs = (ct->get_rhs());//setting the rhs
				}

				else if( ct->get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
				sense = ('L');
				myrhs = (ct->get_rhs());//setting the rhs, which is equal to lhs
				}

				else if( ct->get_rhs() == Inf<double>() ){  // inequality (greate/equal)
				sense = ('G');
				myrhs = (ct->get_lhs());//setting the rhs, which is equal to rhs
				}

         			int    *mycutind = new int   [ mynzcnt ];
         			double *mycutval = new double[ mynzcnt ];

				int i=0;
				for(auto it = ct->get_v_var()->begin() ; it !=ct->get_v_var()->end() ; ++it){
					mycutind[i] = ind_var(it->first);
					mycutval[i] = it->second;									
					i++;		
				}

				CPXcutcallbackadd (env, cbdata, wherefrom, mynzcnt, myrhs, sense, mycutind, mycutval, CPX_USECUT_PURGE);

				delete []mycutind;
				delete []mycutval;
				
				
			} // for (auto ct:v_cuts)
			
		}// if dmod check
		else{
			cout<<"PROBLEEEEEMOOOOOOOOOO";
		}
			
	}//for-loops of modifications

	int j = v_mod.size();	
	auto it = v_mod.begin();
	for (int i = 0; i < j && it != v_mod.end(); i++) {
		it = v_mod.erase(it);
		}
	cout<<endl<<"check clear = " << v_mod.size();
}//if existing modifications




 return( 0 );

 }


/*--------------------------------------------------------------------------*/

void MILPSolver::set_var_value( ColVariable & lvar, double *tmpx, int & i )
{
  
lvar.set_value (tmpx[i]);
i++;

}

/*--------------------------------------------------------------------------*/

void MILPSolver::get_var_solution()
{
  double *tmpx = new double [numcols];

  CPXgetmipx(env,milp,tmpx,0,numcols-1);//obtain solution from CPLEX

  //Retrieve Solution for all the Variables with respect to the order
  int col=0;

 std::queue<Block *> Q; //creating the queue

	Q.push(f_Block);//passing the root, i.e. the father block


	while(!Q.empty())//iterating for the father block and all children
	{
		Block * q_Block = Q.front(); //block to be examined
		Q.pop(); //take out from the queue the examined block


		for (int i = 0; i < q_Block->get_nested_Blocks().size(); ++i)
		{
                       //inserting all the children of the examined block
			Q.push(q_Block->get_nested_Blocks()[i]); 
		}

  		// Passing values for all static Variables
		for (int i =0; i< q_Block->get_static_variables().size(); i++){
		auto f1 = std::bind( &MILPSolver::set_var_value, this, placeholders::_1, std::ref(tmpx), std::ref(col) );
        	un_any_const_static(q_Block->get_static_variables()[i],f1,un_any_type<ColVariable>() );
		}

		//We scan all the dynamic Variables
		for (int i =0; i< q_Block->get_dynamic_variables().size(); i++){
		auto f1 = std::bind( &MILPSolver::set_var_value, this, placeholders::_1, std::ref(tmpx), std::ref(col) );
        	un_any_const_dynamic(q_Block->get_dynamic_variables()[i],f1,un_any_type<ColVariable>() );
		}

	}//while loop

double objval;
CPXgetobjval (env, milp, &objval);

if( f_Block->get_objective_function().type() ==  typeid(LinearObjectiveFunction *) ){

	LinearObjectiveFunction * obj_fun = boost::any_cast< LinearObjectiveFunction * >(f_Block->get_objective_function());
	obj_fun->set_value(objval);

}// if linear objective function
else if( f_Block->get_objective_function().type() ==  typeid(DQuadObjectiveFunction *) ){

	DQuadObjectiveFunction * obj_fun = boost::any_cast< DQuadObjectiveFunction * >(f_Block->get_objective_function());
    obj_fun->set_value(objval);
}
  delete []tmpx;

} 

/*--------------------------------------------------------------------------*/

void MILPSolver::add_modifications( sp_Mod &mod )
{

VariableModification* var_mod = dynamic_cast< VariableModification* >( mod.get() );
ObjFunModification* obj_mod = dynamic_cast< ObjFunModification* >( mod.get() );
ConstraintModification* const_mod = dynamic_cast< ConstraintModification* >( mod.get() );
BlockModificationAD* dyn_mod = dynamic_cast< BlockModificationAD* >( mod.get() );

//Checking if Modification refers to Obj_Function, Constraint or Variable
if( var_mod )
    var_modification (var_mod);
else if( obj_mod )
    of_modification (obj_mod);
else if( const_mod )
   const_modification (const_mod);
else if( dyn_mod )
   dynamic_modification (dyn_mod);
else 
    throw( std::invalid_argument( "Unknown type of Modification" ) ); 

}

/*--------------------------------------------------------------------------*/

void MILPSolver::dynamic_modification( BlockModificationAD* mod )
{
    switch( mod->f_type ) { 

        case( BlockModificationAD::eAddConst ): {
        //we check if we receive a single constraint or a vector of constraints
        
         if( mod->mod_list.type() == typeid( std::vector< LinearConstraint*> ) ) {

              std::vector< LinearConstraint*> v_const = boost::any_cast< std::vector<LinearConstraint * > >( mod->mod_list );
              for (int i = 0 ; i < v_const.size(); ++i)
                    add_dynamic_constraint(v_const[i]);

        }
        else if ( mod->mod_list.type() == typeid( LinearConstraint*) ) {
             
             LinearConstraint * p_const = boost::any_cast< LinearConstraint * >( mod->mod_list );
             add_dynamic_constraint(p_const);
        }
        else 
        throw( std::invalid_argument( "Received unexpected type of Constraint to be added" ) ); 

        break;
        }//case
     case( BlockModificationAD::eAddVar ): {

        //we check if we receive a single constraint or a vector of constraints
        
         if( mod->mod_list.type() == typeid( std::vector< ColVariable*> ) ) {

              std::vector< ColVariable*> v_vars = boost::any_cast< std::vector<ColVariable * > >( mod->mod_list );
              for (int i = 0 ; i < v_vars.size(); ++i)
                    add_dynamic_variable(v_vars[i]);

           }
        else if ( mod->mod_list.type() == typeid( ColVariable*) ) {
             
             ColVariable * p_var = boost::any_cast< ColVariable * >( mod->mod_list );
             add_dynamic_variable(p_var);
          }
          else 
        throw( std::invalid_argument( "Received unexpected type of Variable to be added" ) ); 

        break;
        }//case 
        case( BlockModificationAD::eDelConst ): {
        
        //we check if we receive a single constraint or a vector of constraints
        
         if( mod->mod_list.type() == typeid( std::list< LinearConstraint*> ) ) {

              std::list< LinearConstraint*> l_const = boost::any_cast< std::list<LinearConstraint * > >( mod->mod_list );
              for (auto it = l_const.begin() ; it != l_const.end(); ++it)
                    remove_dynamic_constraint(*it);

        }
        else if ( mod->mod_list.type() == typeid( LinearConstraint*) ) {
             
             LinearConstraint * p_const = boost::any_cast< LinearConstraint * >( mod->mod_list );
             remove_dynamic_constraint(p_const);
        }
        else 
        throw( std::invalid_argument( "Received unexpected type of Constraint to be removed" ) ); 

        break;

        }//case
        case( BlockModificationAD::eDelVar ): {

        //we check if we receive a single constraint or a vector of constraints
        
         if( mod->mod_list.type() == typeid( std::list< ColVariable*> ) ) {

              std::list< ColVariable*> l_var = boost::any_cast< std::list<ColVariable * > >( mod->mod_list );
              for (auto it = l_var.begin() ; it != l_var.end(); ++it)
                    remove_dynamic_variable(*it);

        }
        else if ( mod->mod_list.type() == typeid( ColVariable*) ) {
             
             ColVariable * p_const = boost::any_cast< ColVariable * >( mod->mod_list );
             remove_dynamic_variable(p_const);
        }
        else 
        throw( std::invalid_argument( "Received unexpected type of Constraint to be removed" ) ); 

        break;
        }//case

    }//switch
    
}

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_constraint(LinearConstraint * r_const){


        /* In order to add a new dynamic constraint to CPLEX we use the
           routine addrows that needs the following information: */
        
        int nz = r_const->get_num_active_var() ;// no of nonzero elements
        int matbeg[2] = { 0, nz};
        int matind[nz];
        double matval[nz]; 
        double rhs[1];
        char sense[1];
        int i = 0;
        for (auto it = r_const->get_v_var()->begin() ; it != r_const->get_v_var()->end() ; ++it ){

          /* locating the index in CPLEX searching first in the dynamic and then
             in the static part
          */  
            auto it1 = find_if (v_int_d_var.begin(),v_int_d_var.end(), [&](MILPSolver::int_var pair) {
                    return pair.second == it->first;
                    });
          if (it1 != v_int_d_var.end())
                matind[i] = (it1->first);
          else
                matind[i] = ind_var(it->first);
         
          matval[i] = it->second;
          i++;
        }

        //adding the rhs and the sense
        if( r_const->get_lhs() == r_const->get_rhs() ){ //equality
            sense[0] = ('E');
            rhs[0] = (r_const->get_rhs());//setting the rhs
        }
        else if( r_const->get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
            sense[0] = ('L');
            rhs[0] = (r_const->get_rhs());//setting the rhs, which is equal to lhs
        }
        else if( r_const->get_rhs() == Inf<double>() ){  // inequality (greate/equal)
            sense[0] = ('G');
            rhs[0] = (r_const->get_lhs());//setting the rhs, which is equal to rhs
        }

        //adding the new constraint in the corresponding vectors of pairs       
        v_d_const_int.push_back( std::make_pair ( r_const, v_d_const_int.back().second+1 ) );
        v_int_d_const.push_back( std::make_pair ( v_int_d_const.back().first+1, r_const ) );

        CPXaddrows (env, milp, 0, 1, nz, rhs, sense, matbeg, matind, matval, NULL, NULL);


}

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_variable(ColVariable * r_var){

   /* In order to add a new dynamic constraint to CPLEX we use the
   routine addrows that needs the following information: */

//setting the non-zero coeff
int nz = r_var->get_num_active_const();

//initializing the matrices
int matbeg[2] = { 0, nz};
int matind[nz];
double matval[nz]; 

        int i = 0;
        for (auto it = r_var->active_constraints().begin() ; it != r_var->active_constraints().end() ; ++it ){

          /* locating the index in CPLEX searching first in the dynamic and then
             in the static part
          */  
            LinearConstraint * p_const = dynamic_cast<LinearConstraint *>(*it);
            auto it1 = find_if (v_int_d_const.begin(),v_int_d_const.end(), [&](MILPSolver::int_const pair) {
                    return pair.second == p_const;
                    });
          if (it1 != v_int_d_const.end())
                matind[i] = (it1->first);
          else
                matind[i] = ind_const(p_const);
         
          //find the coefficient of the constraint that refers to the examined variable
         auto it2 = find_if (p_const->get_v_var()->begin(),p_const->get_v_var()->end(), [&](LinearConstraint::coeff_pair pair) {
                    return pair.first == r_var;
                    });
            
          matval[i] = it2->second;
          i++;
        }


//setting the bounds
double  lb[1];
double  ub[1];

lb [0] = r_var->get_lb();
ub [0] = r_var->get_ub();

 v_d_var_int.push_back( std::make_pair ( r_var, v_d_var_int.back().second+1 ) );
 v_int_d_var.push_back( std::make_pair ( v_int_d_var.back().first+1, r_var ) );

CPXaddcols( env, milp, 1, nz, NULL, matbeg, matind, matval, lb, ub, NULL );

}
/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_constraint(LinearConstraint * r_const)
{

//all we need is to find the index of the constraint in CPLEX matrix
int i=0;
auto it = find_if (v_int_d_const.begin(),v_int_d_const.end(), [&](MILPSolver::int_const pair) {
                    return pair.second == r_const;
                    });
          if (it != v_int_d_const.end())
                i = it->first;
          else
        throw( std::invalid_argument( "Constraint is not inside the CPLEX matrix" ) ); 

CPXdelrows (env, milp, i, i+1);

//remove constraint from the corresponding vector of pairs
v_int_d_const.erase(it);
auto it2 = find_if (v_d_const_int.begin(),v_d_const_int.end(), [&](MILPSolver::const_int pair) {
                    return pair.first == r_const;
                    });

 if (it2 != v_d_const_int.end())
                v_d_const_int.erase(it2);
          else
        throw( std::invalid_argument( "Constraint is not inside the CPLEX matrix" ) ); 

}

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_variable(ColVariable * r_var)
{

//all we need is to find the index of the constraint in CPLEX matrix
int i=0;
auto it = find_if (v_int_d_var.begin(),v_int_d_var.end(), [&](MILPSolver::int_var pair) {
                    return pair.second == r_var;
                    });
          if (it != v_int_d_var.end())
                i = it->first;
          else
        throw( std::invalid_argument( "Constraint is not inside the CPLEX matrix" ) ); 

CPXdelcols (env, milp, i, i+1);

//remove constraint from the corresponding vector of pairs
v_int_d_var.erase(it);
auto it2 = find_if (v_d_var_int.begin(),v_d_var_int.end(), [&](MILPSolver::var_int pair) {
                    return pair.first == r_var;
                    });

 if (it2 != v_d_var_int.end())
                v_d_var_int.erase(it2);
          else
        throw( std::invalid_argument( "Constraint is not inside the CPLEX matrix" ) ); 

}

/*--------------------------------------------------------------------------*/

void MILPSolver::var_modification( VariableModification* mod )
{


        ColVariable* colvar = dynamic_cast< ColVariable* >( mod->f_variable );

		switch( mod->f_type ) { 
		case( VariableModification::eFixVar ):   {
			//in case the Variable is set to be fixed we set the 
                        //  lower and upper bound equal to the Variable value

			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
			char * lu = new char [cnt] ; //array with the type of bound for each element
                        lu[0] = 'B' ;
			double * bd = new double [cnt]; //array with the value of bounds
			bd[0] = colvar->get_value();
			CPXchgbds (env, milp, cnt, indices, lu, bd);
            delete []indices; delete []lu; delete []bd;
            break;
        }//case
		case( VariableModification::eUnFixVar ):   {
			//in case the Variable is set to be unfixed we set the 
                        //lower and upper bound equal to their original value

			int cnt = 2; //total number of bounds to be changed
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
			indices[1] = ind_var(colvar);
			char * lu = new char [cnt]; //array with the type of bound for each element
                        lu[0] = 'L' ;
                        lu[1] = 'U' ;
			double * bd = new double [cnt]; //array with the value of bounds
			bd[0] = colvar->get_lb();
			bd[1] = colvar->get_ub();
			CPXchgbds (env, milp, cnt, indices, lu, bd);

            delete []indices; delete []lu; delete []bd;
                        break;
        }//case
		case( ColVariableModification::eChgLB ):   {
                        //changing the LB
			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
			char * lu = new char [cnt]; //array with the type of bound for each element
                        lu[0] = 'L' ;
			double * bd = new double [cnt]; //array with the value of bounds
			bd[0] = colvar->get_lb();
			CPXchgbds (env, milp, cnt, indices, lu, bd);

            delete []indices; delete []lu; delete []bd;
                        break;
        }//case
		case( ColVariableModification::eChgUB ):   {
                        //changing the UB

			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
			char * lu = new char [cnt]; //array with the type of bound for each element
                        lu[0] = 'U' ;
			double * bd = new double [cnt]; //array with the value of bounds
			bd[0] = colvar->get_ub();
			CPXchgbds (env, milp, cnt, indices, lu, bd);

            delete []indices; delete []lu; delete []bd;
                        break;
        }//case
		case( ColVariableModification::eChgType ):   {
                        //in case the type of Variable is set to be changed, we 
                        //need to upgrade equivalently the lower and upper bound 
                        //as well

			int cnt = 1; //total number of variables that change type
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
   			char * ctype = new char [cnt]; //array with the ne2 type of the variable
                        switch( colvar->get_type() ) {
			  case( ColVariable::integer ):      ctype[0] = ('I'); break;
			  case( ColVariable::binary ):       ctype[0] = ('B'); break;
			  case( ColVariable::continuous ):   ctype[0] = ('C');
			  }
			CPXchgctype (env, milp, cnt, indices, ctype);

			int cnt2 = 2; //total number of bounds to be changed
			int * indices2 = new int [cnt2]; //array with the CPLEX index of variable
			indices2[0] = ind_var(colvar);
			indices2[1] = ind_var(colvar);
			char * lu = new char [cnt]; //array with the type of bound for each element
                        lu[0] = 'L' ;
                        lu[1] = 'U' ;
			double * bd = new double [cnt2]; //array with the value of bounds
			bd[0] = colvar->get_lb();
			bd[1] = colvar->get_ub();
			CPXchgbds (env, milp, cnt2, indices2, lu, bd);

            delete []indices;  delete []indices2; delete []lu; delete []bd; delete []ctype;
                        
        }//case

			}//switch



}

/*--------------------------------------------------------------------------*/

void MILPSolver::of_modification( ObjFunModification* mod )
{


		switch( mod->f_type ) { 
			case( ObjFunModification::eSetMin ):  { 
			//setting the objective function to minimize
			CPXchgobjsen (env, milp, 1);
                        break;
        }//case
		case( ObjFunModification::eSetMax ):   {
			//setting the objective function to maximize
			CPXchgobjsen (env, milp, -1);
                        break;
        }//case
		case( LinearOFModification::eAddVar ):  { 
                        //adding linear coefficients
            LinearOFModification* l_mod = dynamic_cast< LinearOFModification* >( mod );
			int cnt = l_mod->v_variables->size();
			int * indices = new int [cnt]; //array with the CPLEX index of variables
			double * values = new double [cnt]; //array with the new values of variables

            int i=0;
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it){
			indices[i]= ind_var (it->first) ;
			values[i]= it->second;
            i++;
			}
			CPXchgobj(env, milp, cnt, indices, values);
            delete []indices;  delete []values;
                        break;
        }//case
		case( LinearOFModification::eRemoveVar ):   {
                        //deleting linear coefficients
            LinearOFModification* l_mod = dynamic_cast< LinearOFModification* >( mod );
			int cnt = l_mod->v_variables->size();
            int * indices = new int [cnt]; //array with the CPLEX index of variables
			double * values = new double [cnt]; //array with the new values of variables
            
            int i=0;
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it){
			indices[i]= ind_var (it->first) ;
			values[i]= 0;
            i++;
			}
			CPXchgobj(env, milp, cnt, indices, values);
            delete []indices;  delete []values;
                        break;
        }//case
		case( LinearOFModification::eModifyCoeff ):   {
                       //changin linear coefficients
            LinearOFModification* l_mod = dynamic_cast< LinearOFModification* >( mod );
			int cnt = l_mod->v_variables->size();
			int * indices = new int [cnt]; //array with the CPLEX index of variables
			double * values = new double [cnt]; //array with the new values of variables
            
            int i=0;
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it){
			indices[i]= ind_var (it->first) ;
			values[i]= 0;
            i++;
			}
			CPXchgobj(env, milp, cnt, indices, values);
            delete []indices;  delete []values;
			break;
        }//case

		case( DQOFModification::eqAddVar ):   {
                        //adding quad coefficients
            DQOFModification* q_mod = dynamic_cast< DQOFModification* >( mod );
            for (auto it = q_mod->v_variables->begin(); it< q_mod->v_variables->end(); ++it)
            	CPXchgqpcoef(env, milp,ind_var (it->first),ind_var (it->first), it->second ) ;
            break;
        }//case
		case( DQOFModification::eqRemoveVar ): {  
                        //deleting quad coefficients
            DQOFModification* q_mod = dynamic_cast< DQOFModification* >( mod );
			for (auto it = q_mod->v_variables->begin(); it< q_mod->v_variables->end(); ++it)
            	CPXchgqpcoef(env, milp,ind_var (it->first),ind_var (it->first), 0 ) ;
            break;
        }//case
		case( DQOFModification::eqModifyCoeff ):  { 
                       //changing quad coefficients
            DQOFModification* q_mod = dynamic_cast< DQOFModification* >( mod );
			for (auto it = q_mod->v_variables->begin(); it< q_mod->v_variables->end(); ++it)
            	CPXchgqpcoef(env, milp,ind_var (it->first),ind_var (it->first), it->second ) ;
            break;
		}//case
} //switch


}

/*--------------------------------------------------------------------------*/

void MILPSolver::const_modification( ConstraintModification* mod )
{

		switch( mod->f_type ) {

		case (ConstraintModification::eRelaxConst):{
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( mod->f_constraint );
            /*In order to relax the constraint all we do is transform it to an inequality with
              rhs equal to infinity
            */
  
            int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt];  //array with the CPLEX index of constraint
			double * values = new double [cnt]; //array with the new rhs of constraint
			char * sense = new char [cnt]; //array with the new sense of constraint
            sense[0] = ('G');
            values[0] =  -Inf<double>();
   			indices[0] = ind_const(r_const);
            CPXchgrhs (env, milp, cnt, indices, values);
			CPXchgsense (env, milp, cnt, indices, sense);
			delete []indices;  delete []values; delete []sense; 
			break;
            }//case

		case (ConstraintModification::eEnforceConst):{
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( mod->f_constraint );
            /*In order to enforce a relaxed constraint all we need to do is reverse the process
              of relaxing it, by changing the sense and the rhs back to the original form of the
              constraint
            */

             int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt];  //array with the CPLEX index of constraint
			double * values = new double [cnt]; //array with the new rhs of constraint
			char * sense = new char [cnt]; //array with the new sense of constraint

    		if( r_const->get_lhs() == r_const->get_rhs() ){ //equality
				sense[0] = ('E');
				values[0] = (r_const->get_rhs());//setting the rhs
			}
			else if(r_const->get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
				sense[0] = ('L');
				values[0] = (r_const->get_rhs());//setting the rhs, which is equal to lhs
			}
			else if( r_const->get_rhs() == Inf<double>() ){  // inequality (greate/equal)
				sense[0] = ('G');
				values[0] = (r_const->get_lhs());//setting the rhs, which is equal to rhs
			}		
			
			indices[0] = ind_const(r_const);
            CPXchgrhs (env, milp, cnt, indices, values);
			CPXchgsense (env, milp, cnt, indices, sense);
			delete []indices;  delete []values; delete []sense; 
			break;

            }//case

		case (RowConstraintModification::eChgLHS): {
            RowConstraintModification* l_mod = dynamic_cast< RowConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );
			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt];  //array with the CPLEX index of constraint
			double * values = new double [cnt]; //array with the new rhs of constraint
			char * sense = new char [cnt]; //array with the new sense of constraint
			
			if( r_const->get_lhs() == r_const->get_rhs() ){ //equality
				sense[0] = ('E');
				values[0] = (r_const->get_rhs());//setting the rhs
			}
			else if(r_const->get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
				sense[0] = ('L');
				values[0] = (r_const->get_rhs());//setting the rhs, which is equal to lhs
			}
			else if( r_const->get_rhs() == Inf<double>() ){  // inequality (greate/equal)
				sense[0] = ('G');
				values[0] = (r_const->get_lhs());//setting the rhs, which is equal to rhs
			}		
			
			indices[0] = ind_const(r_const);
            CPXchgrhs (env, milp, cnt, indices, values);
			CPXchgsense (env, milp, cnt, indices, sense);
			delete []indices;  delete []values; delete []sense; 
    		break;
        } //case
		case (RowConstraintModification::eChgRHS): {

            RowConstraintModification* l_mod = dynamic_cast< RowConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );

			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt];  //array with the CPLEX index of constraint
			double * values = new double [cnt]; //array with the new rhs of constraint
			char * sense = new char [cnt]; //array with the new sense of constraint
			
			if(  r_const->get_lhs() == r_const->get_rhs() ){ //equality
				sense[0] = ('E');
				values[0] = (r_const->get_rhs());//setting the rhs
			}
			else if( r_const->get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
				sense[0] = ('L');
				values[0] = (r_const->get_rhs());//setting the rhs, which is equal to lhs
			}
			else if( r_const->get_rhs() == Inf<double>() ){  // inequality (greate/equal)
				sense[0] = ('G');
				values[0] = (r_const->get_lhs());//setting the rhs, which is equal to rhs
			}		
			
			indices[0] = ind_const(r_const);
            CPXchgrhs (env, milp, cnt, indices, values);
			CPXchgsense (env, milp, cnt, indices, sense);			
            delete []indices;  delete []values; delete []sense; 
		    break;
        } //case
			
		case( LinearConstraintModification::eAddVar ):   {
                        //adding linear coefficients

            LinearConstraintModification* l_mod = dynamic_cast< LinearConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );
			
           for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it)
				 CPXchgcoef (env, milp, ind_const (r_const), ind_var (it->first), it->second);
		   break;
        } //case

		case( LinearConstraintModification::eRemoveVar ):   {

            LinearConstraintModification* l_mod = dynamic_cast< LinearConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );
                        //deleting linear coefficients
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it)
				 CPXchgcoef (env, milp, ind_const (r_const), ind_var (it->first), 0);
		   break;
        } //case

		case( LinearConstraintModification::eModifyCoeff ):   {

            LinearConstraintModification* l_mod = dynamic_cast< LinearConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );
                       //changin linear coefficients
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it)
				 CPXchgcoef (env, milp, ind_const (r_const), ind_var (it->first), it->second);
		   break;
        } //case

}//switch


}

/*--------------------------------------------------------------------------*/
/*----------------------- End File MILPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/ 
