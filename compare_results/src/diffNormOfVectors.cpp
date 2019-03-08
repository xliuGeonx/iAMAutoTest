// $Id$ 
/**
 * @file diffNormOfVectors.cpp
 * @brief Barracuda checking results
 * @author Josue Barboza
 * @date July 11th 2016
 *
 * @copyright  All rights reserved - GeonX SA 2016
*/ 
//$Log$
// -----------------------------------------------------------------------------
// HEADERS
// -----------------------------------------------------------------------------
//#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/exceptions.hpp>


// Avoid conflicting declaration of min/max macros in windows headers
#if !defined(NOMINMAX) && (defined(_WIN32) || defined(_WIN32_)  || defined(WIN32) || defined(_WIN64))
# define NOMINMAX
# ifdef max
#  undef   max
#  undef   min
# endif
#endif

namespace bf=boost::filesystem;

bool is_file( const std::string& filename )
{
  std::ifstream file( filename.c_str() );
  if( !file.is_open() )
  {
    std::cout << " Could not open file : " << filename << std::endl;
    return false;
  }
  else
  {
    file.close();
  }
  return true;
}

// ----------------------------------------------------------------- global_test
bool l2_norm( std::istream& ref, std::istream& test, const double& tol )
{
  bool res = true;
  size_t size1, size2;
  ref >> size1;
  test >> size2;
  if( size1 != size2 )
  {
    std::cout << " The number of values in files does not match ( " 
      << size1 << " - " << size2 << " )"
         << std::endl;
    return false;
  }

  double norm_ref(0.0);
  double norm_test (0.0);
  double norm_diff(0.0);
  
  double ref_val(0.0);
  double test_val(0.0);
  std::string ref_str, test_str;

  for( size_t i = 0; i < size1; ++i )
  {
    ref >> ref_str;
    test >> test_str;

    try
    {
      ref_val = boost::lexical_cast<double>( ref_str );
    }
    catch( boost::exception const & e )
    {
      boost::rethrow_exception( boost::copy_exception(e) );
    }   

    try
    {
      test_val = boost::lexical_cast<double>( test_str );
    }
    catch( boost::exception const & e )
    {
      boost::rethrow_exception( boost::copy_exception(e) );
    }  
    norm_ref += ref_val*ref_val;
    norm_test += test_val*test_val;
    norm_diff += ( ref_val - test_val ) * ( ref_val - test_val );
  }
  norm_ref = std::sqrt( norm_ref );
  norm_test = std::sqrt( norm_test );
  norm_diff = std::sqrt( norm_diff );
  if( std::abs(norm_ref) > std::numeric_limits<double>::epsilon())
  {
    res = (std::abs( (norm_ref-norm_test) / norm_ref) < tol);
  }
  else
  {
    res = ( std::abs(norm_test) < tol );
  }

  return res;
}

// ------------------------------------------------------------------ local_test
bool l1_norm( std::istream& ref, std::istream& test, const double& tol )
{
  bool res = true;
  size_t size1, size2;
  size_t line_num(0);
  ref >> size1;
  test >> size2;

  if( size1 != size2 )
  {
    std::cout << " The number of values in files does not match ( " 
      << size1 << " - " << size2 << " )"
         << std::endl;
    return false;
  }

  while( size1-- > 0 && res == true )
  {
    ++line_num;
//    bool   is_ref_num(false), is_test_num(false);
    double ref_val(1), test_val(-1);
    std::string ref_str, test_str;
    
    ref  >> ref_str;
    test >> test_str;

    try
    {
      ref_val = boost::lexical_cast<double>( ref_str );
    }
    catch( boost::exception const & e )
    {
      boost::rethrow_exception( boost::copy_exception(e) );
    }   

    try
    {
      test_val = boost::lexical_cast<double>( test_str );
    }
    catch( boost::exception const & e )
    {
      boost::rethrow_exception( boost::copy_exception(e) );
    }  
    if( boost::math::isfinite( ref_val ) )        // finite value 
    {
      if( ref_val != 0.0 )                        // non zero
      {
        res = std::abs( (ref_val - test_val) / ref_val ) < tol;
      }
      else                                        // exact zero
      {
        res = std::abs( test_val ) < tol;
      }
    }
    else                                          // non finite value
    {
      res = ( ( boost::math::isinf( ref_val )  &&
                boost::math::isinf( test_val ) &&
                ( ref_val == test_val ) )      ||
              ( boost::math::isnan( ref_val )  &&
                boost::math::isnan( test_val ) ) );
    }

    if( !res ) 
    {
      std::cout << " the local values (respect to the tolerance) are different, at line "
      <<line_num<<" : test value is \""<< test_str 
      <<"\" while reference value is \""<< ref_str <<"\""<<std::endl;
    }
  }

  return res;
}


/* Function used to check that 'opt1' and 'opt2' are not specified
  at the same time. */
void conflicting_options(const boost::program_options::variables_map& vm, 
                        const char* opt1, const char* opt2)
{
  if (vm.count(opt1) && !vm[opt1].defaulted() 
      && vm.count(opt2) && !vm[opt2].defaulted())
      throw boost::program_options::error(std::string("Conflicting options '") 
                        + opt1 + "' and '" + opt2 + "'.");
}
// -----------------------------------------------------------------------------
// MAIN
// -----------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  try
  {
    std::vector<std::string> reference_files_name(4), files_name(4);
    bool is_comparable = true;
    bool l2(true), l1(false);
    bool verbosity(false);
    //bool filecheck(false);
    double tolerance(1.0e-8);

    boost::program_options::options_description prg_description("Allowed options:");
    prg_description.add_options()
      ("help,h", "print usage message")
      ("filecheck", "enable filecheck")
      ("verbose,v", boost::program_options::bool_switch(&verbosity), "enable verbosity")
      ("l1_norm", boost::program_options::bool_switch(&l1), "use relative L1 norm")
      ("l2_norm", boost::program_options::bool_switch(&l2), "use relative L2 norm")
      ("tolerance,t", boost::program_options::value<double>(&tolerance)->default_value(1.0e-8), 
                                          "relative tolerance")
      ("reference,r", boost::program_options::value< std::vector<std::string> >(&reference_files_name),
               "reference files, reference files sequence to compare")
      ("file,f", boost::program_options::value< std::vector<std::string> >(&files_name),
              "computated file to comapre")     
      ;

     boost::program_options::variables_map variablesMap;
     boost::program_options::store(
       boost::program_options::parse_command_line(argc, argv, prg_description), variablesMap );
     boost::program_options::notify(variablesMap);

    if(variablesMap.count("help"))
    {
      std::cout<<prg_description<<std::endl;
      return 0;
    }

    if (variablesMap.count("filecheck"))
    {
      for (std::vector<std::string>::const_iterator fiter = files_name.begin();
        fiter != files_name.end(); ++fiter)
      {
        bf::path pdir(*fiter);
        if (!bf::exists(pdir))
        {
          std::cout << "The result file " << *fiter <<" does not exists" << std::endl;
          exit(2);
        }
      }

      {
        std::cout << "OK" << std::endl;
        return 0;
      }     
    }

    conflicting_options(variablesMap, "l1_norm", "l2_norm");
    if( reference_files_name.size() != files_name.size())
    {
        std::cout << "The number of reference files and testing files should be egal" << std::endl;
        exit( 2 );
    }
    std::vector<std::string>::const_iterator iter_ref_filename = reference_files_name.begin();
    std::vector<std::string>::const_iterator iter_testing_filename = files_name.begin();
    for(;iter_ref_filename!=reference_files_name.end() ; ++iter_ref_filename, ++iter_testing_filename)
    {
      const std::string& ref_filename = *iter_ref_filename;
      const std::string& testing_filename = *iter_testing_filename;
//      bool is_comparable(false);
      if(is_file(ref_filename) && is_file(testing_filename))
      {
        std::ifstream reference( ref_filename.c_str() );
        std::ifstream tested( testing_filename.c_str() );
        bool is_the_same(false);
        if( l1 )
        {
          is_the_same = l1_norm( reference, tested, tolerance );
        }
        else
        {
          is_the_same = l2_norm( reference, tested, tolerance );
        }
  
        reference.close();
        tested.close();
        if( verbosity )
        {
          if( is_the_same )
          {
            std::cout << "Ok" << std::endl;
          }
          else
          {
            std::cout << "Not ok" << std::endl;
          }
        }
        is_comparable = ( is_comparable && is_the_same );  
      }
      else
      {
        exit(2);
      }
    }
    return (is_comparable?0:1 );
  }
  catch(std::exception& e) 
  {
    std::cerr << e.what() << "\n";
  }
}
// =============================================================================
// == END OF FILE ==============================================================

