/*
 * Copyright © 2005 The Walking Fish Group (WFG).
 *
 * This material is provided "as is", with no warranty expressed or implied.
 * Any use is at your own risk. Permission to use or copy this software for
 * any purpose is hereby granted without fee, provided this notice is
 * retained on all copies. Permission to modify the code and to distribute
 * modified code is granted, provided a notice that the code was modified is
 * included with the above copyright notice.
 *
 * http://www.wfg.csse.uwa.edu.au/
 */


/*
 * main.cpp
 *
 * This file contains a simple driver for testing the WFG problems and
 * transformation functions from the WFG test problem toolkit.
 *
 * Changelog:
 *   2005.06.01 (Simon Huband)
 *     - Corrected commments to indicate k and l are the number of position
 *       and distance parameters, respectively (not the other way around).
 */


//// Standard includes. /////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <cmath>
#include <fstream>
#include "knapsack.h"

using namespace std;
//// Toolkit includes. //////////////////////////////////////////////////////

#include "Toolkit/ExampleProblems.h"
#include "Toolkit/TransFunctions.h"


//// Used namespaces. ///////////////////////////////////////////////////////

using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;
using std::vector;
using std::string;


double *max_weight = 0; //max_weight
int *item_num = 0; // item numbering

double **original_weight = 0; //weight[o][i] o: object, i: item
double **original_profit = 0; //profit[o][i] o: object, i: item
double **original_profit_per_weight = 0; //profit[o][i] / weight[o][i]
double *original_q_j = 0; //q_j = max(profit[o][i] / weight[o][i])
int *original_x_check = 0; //desicion number
int *x_check = 0; //desicion number

double **weight = 0; //weight[o][i] o: object, i: item
double **profit = 0; //profit[o][i] o: object, i: item
double **profit_per_weight = 0; //profit[o][i] / weight[o][i]
double *q_j = 0; //q_j = max(profit[o][i] / weight[o][i])



//// Local functions. ///////////////////////////////////////////////////////

namespace
{

//** Using a uniform random distribution, generate a number in [0,bound]. ***

double next_double( const double bound = 1.0 )
{
  assert( bound > 0.0 );

  return bound * rand() / static_cast< double >( RAND_MAX );
}


//** Create a random Pareto optimal solution for WFG1. **********************

vector< double > WFG_1_random_soln( const int k, const int l )
{
  vector< double > result;  // the result vector


  //---- Generate a random set of position parameters.

  for( int i = 0; i < k; i++ )
  {
    // Account for polynomial bias.
    result.push_back( pow( next_double(), 50.0 ) );
  }


  //---- Set the distance parameters.

  for( int i = k; i < k+l; i++ )
  {
    result.push_back( 0.35 );
  }


  //---- Scale to the correct domains.

  for( int i = 0; i < k+l; i++ )
  {
    result[i] *= 2.0*(i+1);
  }


  //---- Done.

  return result;
}


//** Create a random Pareto optimal solution for WFG2-WFG7. *****************

vector< double > WFG_2_thru_7_random_soln( const int k, const int l )
{
  vector< double > result;  // the result vector


  //---- Generate a random set of position parameters.

  for( int i = 0; i < k; i++ )
  {
    result.push_back( next_double() );
  }


  //---- Set the distance parameters.

  for( int i = k; i < k+l; i++ )
  {
    result.push_back( 0.35 );
  }


  //---- Scale to the correct domains.

  for( int i = 0; i < k+l; i++ )
  {
    result[i] *= 2.0*(i+1);
  }


  //---- Done.

  return result;
}


//** Create a random Pareto optimal solution for WFG8. **********************

vector< double > WFG_8_random_soln( const int k, const int l )
{
  vector< double > result;  // the result vector


  //---- Generate a random set of position parameters.

  for( int i = 0; i < k; i++ )
  {
    result.push_back( next_double() );
  }


  //---- Calculate the distance parameters.

  for( int i = k; i < k+l; i++ )
  {
    const vector< double >  w( result.size(), 1.0 );
    const double u = TransFunctions::r_sum( result, w  );

    const double tmp1 = fabs( floor( 0.5 - u ) + 0.98/49.98 );
    const double tmp2 = 0.02 + 49.98*( 0.98/49.98 - ( 1.0 - 2.0*u )*tmp1 );

    result.push_back( pow( 0.35, pow( tmp2, -1.0 ) ));
  }


  //---- Scale to the correct domains.

  for( int i = 0; i < k+l; i++ )
  {
    result[i] *= 2.0*(i+1);
  }


  //---- Done.

  return result;
}


//** Create a random Pareto optimal solution for WFG9. **********************

vector< double > WFG_9_random_soln( const int k, const int l )
{
  vector< double > result( k+l );  // the result vector


  //---- Generate a random set of position parameters.

  for( int i = 0; i < k; i++ )
  {
    result[i] = next_double();
  }


  //---- Calculate the distance parameters.

  result[k+l-1] = 0.35;  // the last distance parameter is easy

  for( int i = k+l-2; i >= k; i-- )
  {
    vector< double > result_sub;
    for( int j = i+1; j < k+l; j++ )
    {
      result_sub.push_back( result[j] );
    }

    const vector< double > w( result_sub.size(), 1.0 );
    const double tmp1 = TransFunctions::r_sum( result_sub, w  );

    result[i] = pow( 0.35, pow( 0.02 + 1.96*tmp1, -1.0 ) );
  }


  //---- Scale to the correct domains.

  for( int i = 0; i < k+l; i++ )
  {
    result[i] *= 2.0*(i+1);
  }


  //---- Done.

  return result;
}


//** Create a random Pareto optimal solution for I1. *****************

vector< double > I1_random_soln( const int k, const int l )
{
  vector< double > result;  // the result vector


  //---- Generate a random set of position parameters.

  for( int i = 0; i < k; i++ )
  {
    result.push_back( next_double() );
  }


  //---- Set the distance parameters.

  for( int i = k; i < k+l; i++ )
  {
    result.push_back( 0.35 );
  }


  //---- Done.

  return result;
}


//** Create a random Pareto optimal solution for I2. **********************

vector< double > I2_random_soln( const int k, const int l )
{
  vector< double > result( k+l );  // the result vector


  //---- Generate a random set of position parameters.

  for( int i = 0; i < k; i++ )
  {
    result[i] = next_double();
  }


  //---- Calculate the distance parameters.

  result[k+l-1] = 0.35;  // the last distance parameter is easy

  for( int i = k+l-2; i >= k; i-- )
  {
    vector< double > result_sub;
    for( int j = i+1; j < k+l; j++ )
    {
      result_sub.push_back( result[j] );
    }

    const vector< double > w( result_sub.size(), 1.0 );
    const double tmp1 = TransFunctions::r_sum( result_sub, w  );

    result[i] = pow( 0.35, pow( 0.02 + 1.96*tmp1, -1.0 ) );
  }


  //---- Done.

  return result;
}


//** Create a random Pareto optimal solution for I3. **********************

vector< double > I3_random_soln( const int k, const int l )
{
  vector< double > result;  // the result vector


  //---- Generate a random set of position parameters.

  for( int i = 0; i < k; i++ )
  {
    result.push_back( next_double() );
  }


  //---- Calculate the distance parameters.

  for( int i = k; i < k+l; i++ )
  {
    const vector< double >  w( result.size(), 1.0 );
    const double u = TransFunctions::r_sum( result, w  );

    const double tmp1 = fabs( floor( 0.5 - u ) + 0.98/49.98 );
    const double tmp2 = 0.02 + 49.98*( 0.98/49.98 - ( 1.0 - 2.0*u )*tmp1 );

    result.push_back( pow( 0.35, pow( tmp2, -1.0 ) ));
  }


  //---- Done.

  return result;
}


//** Create a random Pareto optimal solution for I4. **********************

vector< double > I4_random_soln( const int k, const int l )
{
  return I1_random_soln( k, l );
}


//** Create a random Pareto optimal solution for I5. **********************

vector< double > I5_random_soln( const int k, const int l )
{
  return I3_random_soln( k, l );
}

//** Generate a random solution for a given problem. ************************

vector< double > problem_random_soln
(
  const int k,
  const int l,
  const std::string fn
)
{
  if ( fn == "wfg1" )
  {
    return WFG_1_random_soln( k, l );
  }
  else if
  (
    fn == "wfg2" ||
    fn == "wfg3" ||
    fn == "wfg4" ||
    fn == "wfg5" ||
    fn == "wfg6" ||
    fn == "wfg7"
  )
  {
    return WFG_2_thru_7_random_soln( k, l );
  }
  else if ( fn == "wfg8" )
  {
    return WFG_8_random_soln( k, l );
  }
  else if ( fn == "wfg9" )
  {
    return WFG_9_random_soln( k, l );
  }
  else if ( fn == "I1" )
  {
    return I1_random_soln( k, l );
  }
  else if ( fn == "I2" )
  {
    return I2_random_soln( k, l );
  }
  else if ( fn == "I3" )
  {
    return I3_random_soln( k, l );
  }
  else if ( fn == "I4" )
  {
    return I4_random_soln( k, l );
  }
  else if ( fn == "I5" )
  {
    return I5_random_soln( k, l );
  }
  else
  {
    assert( false );
    return vector< double >();
  }
}


//** Calculate the fitness for a problem given some parameter set. **********

vector< double > problem_calc_fitness
(
  const vector< double >& z,
  const int k,
  const int M,
  const std::string fn
)
{
  if ( fn == "wfg1" )
  {
    return Problems::WFG1( z, k, M );
  }
  else if ( fn == "wfg2" )
  {
    return Problems::WFG2( z, k, M );
  }
  else if ( fn == "wfg3" )
  {
    return Problems::WFG3( z, k, M );
  }
  else if ( fn == "wfg4" )
  {
    return Problems::WFG4( z, k, M );
  }
  else if ( fn == "wfg5" )
  {
    return Problems::WFG5( z, k, M );
  }
  else if ( fn == "wfg6" )
  {
    return Problems::WFG6( z, k, M );
  }
  else if ( fn == "wfg7" )
  {
    return Problems::WFG7( z, k, M );
  }
  else if ( fn == "wfg8" )
  {
    return Problems::WFG8( z, k, M );
  }
  else if ( fn == "wfg9" )
  {
    return Problems::WFG9( z, k, M );
  }
  else if ( fn == "I1" )
  {
    return Problems::I1( z, k, M );
  }
  else if ( fn == "I2" )
  {
    return Problems::I2( z, k, M );
  }
  else if ( fn == "I3" )
  {
    return Problems::I3( z, k, M );
  }
  else if ( fn == "I4" )
  {
    return Problems::I4( z, k, M );
  }
  else if ( fn == "I5" )
  {
    return Problems::I5( z, k, M );
  }
  else
  {
    assert( false );
    return vector< double >();
  }
}


//** Convert a double vector into a string. *********************************

string make_string( const vector< double >& v )
{
  std::ostringstream result;

  if( !v.empty() )
  {
    result << v.front();
  }

  for( int i = 1; i < static_cast< int >( v.size() ); i++ )
  {
    result << " " << v[i];
  }

  return result.str();
}

}  // unnamed namespace


void hoge(double y[], double fit[], char problem[],int M, int k, int l){
	//---- Get the function name.

	const std::string fn(problem);
	

	//---- Generate values for desired function.

	if
		(
		fn == "wfg1" ||
		fn == "wfg2" ||
		fn == "wfg3" ||
		fn == "wfg4" ||
		fn == "wfg5" ||
		fn == "wfg6" ||
		fn == "wfg7" ||
		fn == "wfg8" ||
		fn == "wfg9" ||
		fn == "I1" ||
		fn == "I2" ||
		fn == "I3" ||
		fn == "I4" ||
		fn == "I5"
		)
	{
	
		
		srand(0);  // seed the random number generator

		// Generate count random fitness values.
	
	

	//	const vector< double >& z = problem_random_soln(k, l, fn);
	
		vector<double> z(k + l);



		for (int i = 0; i < k + l; i++){
			z[i] = y[i];
			z[i] *= 2.0 * (i + 1);
		}

		

		const vector< double >& f = problem_calc_fitness(z, k, M, fn);

		for (int j = 0; j < f.size(); j++){
			fit[j] = f[j];
		}

			/* std::cout << make_string( f ) << std::endl;*/
		
	}
	else if
		(
		fn == "b_poly" ||
		fn == "b_flat" ||
		fn == "s_linear" ||
		fn == "s_decept" ||
		fn == "s_multi"
		)
	{
		const int count = 10000;  // the number of times (-1) to sample the function

		// Sample the transformation function count+1 times.
		for (int i = 0; i <= count; i++)
		{
			const double y = static_cast< double >(i) / count;
			double new_y;

			if (fn == "b_poly")
			{
				new_y = TransFunctions::b_poly(y, 20.0);
			}
			else if (fn == "b_flat")
			{
				new_y = TransFunctions::b_flat(y, 0.7, 0.4, 0.5);
			}
			else if (fn == "s_linear")
			{
				new_y = TransFunctions::s_linear(y, 0.35);
			}
			else if (fn == "s_decept")
			{
				new_y = TransFunctions::s_decept(y, 0.35, 0.005, 0.05);
			}
			else if (fn == "s_multi")
			{
				new_y = TransFunctions::s_multi(y, 5, 10, 0.35);
			}
			else
			{
				assert(false);
				exit(1);
			}

			std::cout << y << " " << new_y << std::endl;
		}
	}
	else if (fn == "b_param" || fn == "r_sum" || fn == "r_nonsep")
	{
		srand(0);

		const int count = 10000;

		// Randomly sample the transformation count times.
		for (int i = 0; i < count; i++)
		{
			vector< double > y;

			y.push_back(next_double());
			y.push_back(next_double());
			double new_y;

			if (fn == "b_param")
			{
				new_y = TransFunctions::b_param(y[0], y[1], 0.5, 2, 10);
			}
			else if (fn == "r_sum")
			{
				vector< double > w;

				w.push_back(1.0);
				w.push_back(5.0);

				new_y = TransFunctions::r_sum(y, w);
			}
			else if (fn == "r_nonsep")
			{
				new_y = TransFunctions::r_nonsep(y, 2);
			}
			else
			{
				assert(false);
				exit(1);
			}

			std::cout << y[0] << " " << y[1] << " " << new_y << std::endl;
		}
	}
	else
	{
		std::cout << "Invalid fn.\n";
		exit(1);
	}
	
}




void dtlz1(double y[], double fit[], int ob, int k){
	int x_m = k;
	double g_x_m = 0.0;

	for (int i = ob - 1; i < ob + k - 1; i++){
		g_x_m += (y[i] - 0.5) * (y[i] - 0.5) - cos(20.0 * M_PI * (y[i] - 0.5));
	}
	g_x_m += x_m;
	g_x_m *= 100.0;

	for (int i = 0; i < ob; i++){
		fit[i] = 0.5 * (1.0 + g_x_m);
	}

	for (int i = 0; i < ob - 1; i++){
		for (int k = 0; k < ob - i - 1; k++){
			fit[i] *= y[k];
		}
	}
	for (int i = 1; i < ob; i++){
		fit[i] *= (1.0 - y[ob - i - 1]);
	}
}



void dtlz2(double y[], double fit[], int ob,int k){
	int x_m = k;
	double g_x_m = 0.0;

	for (int i = ob - 1; i < ob + k - 1; i++){
		g_x_m += (y[i] - 0.5) * (y[i] - 0.5);
	}

	for (int i = 0; i < ob; i++){
		fit[i] = (1.0 + g_x_m);
	}

	for (int i = 0; i < ob - 1; i++){
		for (int k = 0; k < ob - i - 1; k++){
			fit[i] *= cos(0.5 * y[k] * M_PI);
		}
	}
	for (int i = 1; i < ob; i++){
		fit[i] *= sin(0.5 * y[ob - i - 1] * M_PI);
	}
}

void dtlz3(double y[], double fit[], int ob,int k){
	int x_m = k;
	double g_x_m = 0.0;

	for (int i = ob - 1; i < k + ob - 1; i++){
		g_x_m += (y[i] - 0.5) * (y[i] - 0.5) - cos(20.0 * M_PI * (y[i] - 0.5));
	}
	g_x_m += x_m;
	g_x_m *= 100.0;

	for (int i = 0; i < ob; i++){
		fit[i] = (1.0 + g_x_m);
	}

	for (int i = 0; i < ob; i++){
		for (int k = 0; k < ob - i - 1; k++){
			fit[i] *= cos(0.5 * y[k] * M_PI);
		}
	}
	for (int i = 1; i < ob; i++){
		fit[i] *= sin(0.5 * y[ob - i - 1] * M_PI);
	}
}

void dtlz4(double y[], double fit[], int ob, int k){
	int x_m = k;
	double g_x_m = 0.0;
	double alpha = 100.0;

	for (int i = ob - 1; i < ob + k - 1; i++){
		g_x_m += (y[i] - 0.5) * (y[i] - 0.5);
	}

	for (int i = 0; i < ob; i++){
		fit[i] = (1.0 + g_x_m);
	}

	for (int i = 0; i < ob - 1; i++){
		for (int k = 0; k < ob - i - 1; k++){
			fit[i] *= cos(0.5 * pow(y[k], alpha) * M_PI);
		}
	}
	for (int i = 1; i < ob; i++){
		fit[i] *= sin(0.5 * pow(y[ob - i - 1], alpha) * M_PI);
	}
}

void dtlz7(double y[], double fit[], int ob, int k){
	int x_m = k;
	double g_x_m = 0.0;
	double h = 0.0;
	double *x = new double[ob + k - 1];
	
	for (int i = 0; i < ob + k - 1; i++){
		x[i] = y[i];
	}

	for (int i = ob - 1; i < ob + k - 1; i++){
		g_x_m += x[i];
	}

	for (int i = 0; i < ob - 1; i++){
		fit[i] = x[i];
	}

	g_x_m = 1.0 + ((9.0 * g_x_m) / ((double)x_m));

	for (int i = 0; i < ob - 1; i++){
		h += (fit[i] / (1.0 + g_x_m)) * (1.0 + sin(3.0 * M_PI * fit[i]));
	}

	h = (double)ob - h;

	fit[ob - 1] = (1.0 + g_x_m) * h;
	
	delete[] x;

}

void knapsack_file_read(char *file, int &item, int &ob){

	ifstream fin(file);
	if (!fin){
		cout << "cant open the file" << endl;
		exit(1);
	}
	fin >> item;
	fin >> ob;


	item_num = new int[item];
	weight = new double*[ob];
	profit = new double*[ob];
	profit_per_weight = new double*[ob];
	max_weight = new double[ob];
	q_j = new double[item];
	x_check = new int[item];
	for (int i = 0; i < item; i++){
		x_check[i] = 1;
	}
	for (int o = 0; o < ob; o++){
		fin >> max_weight[o];
		weight[o] = new double[item];
		profit[o] = new double[item];
		profit_per_weight[o] = new double[item];
		for (int i = 0; i < item; i++){
			fin >> weight[o][i];
			fin >> profit[o][i];
			profit_per_weight[o][i] = profit[o][i] / weight[o][i];
			item_num[i] = i;
		}
	}

	ofstream check1("weight_file_check.txt");
	for (int i = 0; i < item; i++){
		for (int o = 0; o < ob; o++){
			check1 << weight[o][i] << " ";
		}
		check1 << endl;
	}
	ofstream check2("profit_file_check.txt");
	for (int i = 0; i < item; i++){
		for (int o = 0; o < ob; o++){
			check2 << profit[o][i] << " ";
		}
		check2 << endl;
	}
	check1.close();
	check2.close();


	for (int i = 0; i < item; i++){
		int max = 0;
		for (int o = 0; o < Q_J; o++){
			if (profit_per_weight[max][i] < profit_per_weight[o][i]){
				max = o;
			}
		}
		q_j[i] = profit_per_weight[max][i];
	}

	original_weight = new double*[ob];
	original_profit = new double*[ob];
	original_profit_per_weight = new double*[ob];
	original_q_j = new double[item];
	original_x_check = new int[item];

	for (int o = 0; o < ob; o++){
		original_weight[o] = new double[item];
		original_profit[o] = new double[item];
		original_profit_per_weight[o] = new double[item];
		for (int i = 0; i < item; i++){
			original_weight[o][i] = weight[o][i];
			original_profit[o][i] = profit[o][i];
			original_profit_per_weight[o][i] = profit_per_weight[o][i];
		}

	}
	for (int i = 0; i < item; i++){
		original_q_j[i] = q_j[i];
		original_x_check[i] = x_check[i];
	}

	fin.close();
}

void sorting_profit_per_weight(int item, int ob){

	double min = 0.0;
	double temp_double;
	int temp_int;
	int k = 0;
	double **temp_array;

	//profit_per_weight and item_num sorting
	for (int i = 0; i < item; i++){
		for (int j = item - 1; j > i; j--){
			if (q_j[j - 1] > q_j[j]){
				temp_double = q_j[j];
				q_j[j] = q_j[j - 1];
				q_j[j - 1] = temp_double;


				temp_int = item_num[j];
				item_num[j] = item_num[j - 1];
				item_num[j - 1] = temp_int;
			}
		}

	}

	//weight sorting
	temp_array = new double *[ob];
	for (int o = 0; o < ob; o++){
		temp_array[o] = new double[item];
		for (int i = 0; i < item; i++){
			temp_array[o][i] = weight[o][i];
		}
	}

	for (int o = 0; o < ob; o++){
		for (int i = 0; i < item; i++){
			weight[o][i] = temp_array[o][item_num[i]];
		}
	}

	//profit sorting
	for (int o = 0; o < ob; o++){
		for (int i = 0; i < item; i++){
			temp_array[o][i] = profit[o][i];
		}
	}
	for (int o = 0; o < ob; o++){
		for (int i = 0; i < item; i++){
			profit[o][i] = temp_array[o][item_num[i]];
		}
	}

	//profit_per_weight sorting
	for (int o = 0; o < ob; o++){
		for (int i = 0; i < item; i++){
			temp_array[o][i] = profit_per_weight[o][i];
		}
	}
	for (int o = 0; o < ob; o++){
		for (int i = 0; i < item; i++){
			profit_per_weight[o][i] = temp_array[o][item_num[i]];
		}
	}

	for (int n = 0; n < ob; n++){
		delete[] temp_array[n];
	}
	delete[] temp_array;

}


void fitness_knap(int y[], double fit[], int item, int ob){
	double *sum_profit = new double[ob];
	double *sum_weight = new double[ob];

	for (int o = 0; o < ob; o++){
		sum_profit[o] = 0.0;
		sum_weight[o] = 0.0;
	}

	/*int i = 0;
	for (int o = 0; o < ob; o++){
	for (i = 0; i < item; i++){
	sum_profit[o] += original_profit[o][i] * y[i];
	sum_weight[o] += original_weight[o][i] * y[i];
	}
	}
	i = 0;
	int flag = 0;
	for (int o = 0; o < Q_J; o++){
	if (sum_weight[o] > max_weight[o]){
	flag = 1;
	break;
	}
	}
	while (flag == 1){
	for (int o = 0; o < ob; o++){
	sum_profit[o] -= profit[o][i] * y[item_num[i]];
	sum_weight[o] -= weight[o][i] * y[item_num[i]];
	}
	y[item_num[i]] = 0;
	flag = 0;
	for (int o = 0; o < Q_J; o++){
	if (sum_weight[o] > max_weight[o]){
	flag = 1;
	break;
	}
	}
	i++;
	}*/
	int i = item - 1;
	for (; i >= 0; i--){
		int flag = 0;
		for (int o = 0; o < Q_J; o++){
			if (sum_weight[o] + original_weight[o][item_num[i]] * y[item_num[i]] > max_weight[o]){
				flag = 1;
				break;
			}
		}
		if (flag == 0){
			for (int o = 0; o < ob; o++){
				sum_weight[o] += original_weight[o][item_num[i]] * y[item_num[i]];
				sum_profit[o] += original_profit[o][item_num[i]] * y[item_num[i]];
			}
		}
		else{
			break;
		}
	}
	for (; i >= 0; i--){
		y[item_num[i]] = 0;
	}

	for (int o = 0; o < ob; o++){
		fit[o] = sum_profit[o];
	}
	delete[] sum_weight;
	delete[] sum_profit;

}

void repair_output(int item, int ob){
	ofstream fout_repair("repair.txt");
	fout_repair << fixed << setprecision(16);
	for (int i = 0; i < item; i++){
		fout_repair << q_j[i] << "\t" << item_num[i] << endl;
	}
	fout_repair.close();
}

void check_input_file_output(int y[], int item, int ob){

	ofstream fout_file("file.txt");
	double *fit = new double[ob];
	for (int i = 0; i < item; i++){
		fout_file << y[i] << " ";
	}
	fout_file << endl << endl;
	fitness_knap(y, fit, item, ob);
	fout_file << fit[0] << " " << fit[1] << endl;
	for (int i = 0; i < item; i++){
		fout_file << y[i] << " ";
	}
	fout_file << endl << endl;

	for (int i = 0; i < item; i++){
		y[i] = i % 2;
	}
	for (int i = 0; i < item; i++){
		fout_file << y[i] << " ";
	}
	fout_file << endl << endl;
	fitness_knap(y, fit, item, ob);
	fout_file << fit[0] << " " << fit[1] << endl;
	for (int i = 0; i < item; i++){
		fout_file << y[i] << " ";
	}
	fout_file << endl << endl;

	for (int i = 0; i < item; i++){
		y[i] = (i + 1) % 2;
	}
	for (int i = 0; i < item; i++){
		fout_file << y[i] << " ";
	}
	fout_file << endl << endl;
	fitness_knap(y, fit, item, ob);
	fout_file << fit[0] << " " << fit[1] << endl;
	for (int i = 0; i < item; i++){
		fout_file << y[i] << " ";
	}
	fout_file << endl << endl;

	delete[] fit;
	fout_file.close();

}

void fitness(double y[], double fit[], char problem[], int ob, int k, int l,int k_dtlz){
	
	if (strcmp("dtlz1", problem) == 0){
		dtlz1(y, fit, ob, k_dtlz);
		//std::cout << "dtlz1" << std::endl;
	}
	else if (strcmp("dtlz2", problem) == 0){
		dtlz2(y, fit, ob, k_dtlz);
		//std::cout << "dtlz2" << std::endl;
	}
	else if (strcmp("dtlz3", problem) == 0){
		dtlz3(y, fit, ob, k_dtlz);
		//std::cout << "dtlz3" << std::endl;
	}
	else if (strcmp("dtlz4", problem) == 0){
		dtlz4(y, fit, ob, k_dtlz);
		//std::cout << "dtlz4" << std::endl;
	}
	else if (strcmp("dtlz7", problem) == 0){
		dtlz7(y, fit, ob, k_dtlz);
		//std::cout << "dtlz7" << std::endl;
	}
	else{
		hoge(y, fit, problem, ob, k, l);
	}
}
//// Standard functions. ////////////////////////////////////////////////////