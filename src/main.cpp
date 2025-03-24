/**
 *Alec Rosentrater
 *R2U2 MLTL Benchmark Generation Tool
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <tuple>
#include "formula.h"
#include <unistd.h>

using namespace std;

//input parameters
string formula_file = "";
string finitetracefile = "";
int trace_length;

//Globals for tracking size of finite trace
//int max_time = 0;
int max_atoms = 0;
int num_atoms = 0;

//Globals for storing which pattern is to be used
bool pattern_sat = false;
bool pattern_unsat = false;
bool pattern_median1 = false;//satisfiable then unsatisfiable
bool pattern_median2 = false;//unsatisfiable then satisfiable
bool pattern_random_trace = false;
bool pattern_random_oracle = false;
bool pattern_fail = false;

//Global for storing the 1st frame of the finite trace for ease of access deep in recursion
vector <int> zero_states;
vector <string> atom_names;

//Global for storing random seed
int rand_seed = 0;

///// << operator overloading for easier printing and debugging finite traces
template <typename S>
ostream& operator<<(ostream& os, const vector<S>& vector){
  for (auto element : vector) {
    os << element << "\n";
  }
  return os;
}



//////////////////////////////////////////////
void print_helper () {
  cout << "MLTL Benchmark Generator help message.Currently limited to 1 formula per file.\n";
  cout << "Usage: ./BenchmarkGenerator [option] 'formula.mltl' [length]" << endl;
  cout << "Options: --help, -S, -U, -M1, -M2, -R, -SR" << endl;
  cout << "\t -S: Almost-Satisfiable Pattern" << endl;
  cout << "\t -U: Almost-Unsatisfiable Pattern" << endl;
  cout << "\t -M1: Median-Satisfiable Pattern with Almost-Satisfiable for the first half and Almost-Unsatisfiable for the second half of the generated trace" << endl;
  cout << "\t -M2: Median-Satisfiable Pattern with Almost-Unsatisfiable for the first half and Almost-Satisfiable for the second half of the generated trace" << endl;
  cout << "\t -R: Random-Satisfiable Pattern from trace variable assignment" << endl;
  cout << "\t -SR: Semi-Random-Satisfiable Pattern" << endl;
  // cout << "Supported symbols: ! & | <- <-> U R G F" << endl;
  cout << "Length: desired integer length of generated trace" << endl;
}

///////////////////////////////////////////////
void Oracle_printer(vector<bool> Oracle){
  for(int i = 0; i < Oracle.size(); i++){
    cout << i << ", ";
    printf("%s\n",Oracle[i]?"true":"false");
  }
  printf("\n");
}

/////////////////////////////////////////////
vector<string> percentage_reporter(vector<bool> Oracle){
  float sat_count = 0;
  float unsat_count = 0;
  float total = 0;
   for(int i = 0; i < Oracle.size(); i++){
    if (Oracle[i]){
      sat_count++;
    }
    else if (!Oracle[i]){
      unsat_count++;
    }
    else{
      cerr << "Error in Oracle percentage evaluation";
      cerr << Oracle[i] << endl;
    }
    total++;
  }//end new oracle generation for loop
   float percent_sat =(sat_count / total) * 100;
   float percent_unsat =(unsat_count / total) * 100;
   string line1 = "#" + to_string(percent_sat) + "% satisfing(" + to_string(sat_count) + " frames)";
   string line2 =  "#" + to_string(percent_unsat) + "% dissatisfying (" + to_string(unsat_count) + " frames)";
   string line3 = "#" + to_string(total) + " total frames generated";
     vector<string> return_vec;
     return_vec.push_back(line1);
     return_vec.push_back(line2);
     return_vec.push_back(line3);
     return return_vec;
     // cout << percent_sat << "% satisfiable with " << sat_count << " frames satisfiable of " << total << " total frames" << endl;
     // cout << percent_unsat << "% unsatisfiable with " <<  unsat_count << " frames unsatisfiable of " << total << " total frames" << endl;
}
/////////////////////////////////////////////
void model_printer(std::tuple<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >, std::vector<bool, std::allocator<bool> >, mltl::Formula> result){
   mltl::Formula out_formula = std::get<mltl::Formula>(result);
   mltl::Formula *out_formula_point = &out_formula;
   vector<bool> Oracle = std::get<1>(result);
   vector<vector<int>> trace_prime =  std::get<0>(result);

   cout << percentage_reporter(Oracle) << endl;
   
   cout << out_formula_point->to_string() << endl;
   Oracle_printer(Oracle);
   cout << trace_prime << endl;
}
///////////////////////////////////////////
//generates the trace output csv file
void output_trace_printer(std::tuple<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >, std::vector<bool, std::allocator<bool> >, mltl::Formula> result,int number, string name){
  
  vector<vector<int>> trace =  std::get<0>(result); //extract the trace from the given result object
  //open or create the output file
  string filename = name + "_" + to_string(number);
  string outfilename = filename + "_trace.csv";
  ofstream outfile (outfilename);
  //print the trace to the output file
  // cout << "trace size: " << trace.size() << endl;
  // cout << "trace[0] size: " << trace[0].size() << endl;
  // cout << "trace[0][0]: " << trace[0][0] << endl;
  for(int i =0; i < trace.size(); i++){
    for(int j = 0; j < trace[0].size(); j++){
      outfile << trace[i][j];
      if (j != trace[0].size()-1) {
        outfile << ",";
      }
      //outfile << ",";
    }
    outfile << "\n";
  }//end outer for
  //close the output file
  outfile.close();
}

///////////////////////////////////////////
//generates the results text file
void output_result_printer(std::tuple<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >, std::vector<bool, std::allocator<bool> >, mltl::Formula> result, int formula_number,string pattern,string filename){
  mltl::Formula out_formula = std::get<mltl::Formula>(result);
  mltl::Formula *out_formula_point = &out_formula;
  vector<bool> Oracle = std::get<1>(result);
  vector<vector<int>> trace_prime =  std::get<0>(result);

  float unsat_percent = 0;
  float sat_percent = 0;
  int sat_frames = 0;
  int unsat_frames = 0;
  int total_frames = 0;
  
  string outfilename = filename +"_"+ to_string(formula_number) + ".txt";

  ofstream outfile (outfilename);
  /*outfile << "#Generation Report for " << filename << " mltl::Formula " << formula_number << endl;
  outfile << "#" <<  pattern << " Pattern" << endl;
  outfile << "#Desired Trace Length: " << trace_length <<  endl;
  outfile << "#Random Seed: " << rand_seed << endl;
  outfile << "#Generated trace written to " << filename << "_trace.csv" << endl;*/
  
  //outfile << "#\n#Results: " << endl;
  //outfile << percentage_reporter(Oracle);
  
  // outfile << "#\n#Final formula after progression:" << endl;
  //  outfile << "#" << out_formula_point->to_string() << endl;

  //outfile <<"#\n#Oracle:" << endl;
  //outfile << "#timestamp\tverdict"<<endl;
  for(int i =0; i < Oracle.size(); i++){
    outfile << i << ",";
    if (Oracle[i] == 0){
      outfile << "F" << endl;
    }
    else{
      outfile << "T" << endl;
    }
    //outfile << Oracle[i] << endl;
    // printf("%s\n",Oracle[i]?"true":"false");
  }
  outfile.close();
}
/////////////////////////////////////////// //TODO - needs updated for new patterns
void parse_args (int argc, char *argv[]) 
{
  if (argc < 4){
      cout << argc << " arguments provided, minimum 4 expected" << endl;
      print_helper ();
      exit (0);
    }
  if (strcmp(argv[1],"--help")==0){
    // cout << "--help selected" << endl;
      print_helper ();
      exit (0);
  }
  if(argc == 4){
    if(strcmp(argv[1],"-S")==0){
      pattern_sat = true;
    }
    
    else if (strcmp(argv[1],"-U")==0){
      pattern_unsat = true;
    }
    else if (strcmp(argv[1],"-M1")==0){
      pattern_median1 = true;
    }
    
    else if (strcmp(argv[1],"-M2")==0){
      pattern_median2 = true;
    }
    
    else if( strcmp(argv[1],"-R") ==0 ){
      pattern_random_trace = true;
      rand_seed = 0;
      cout << "Using random seed: " << rand_seed << endl;
      cout << "Random pattern not currently supported" << endl;
    }
    else if ( strcmp(argv[1],"-SR") == 0){
      rand_seed = 0;
      cout << "Using random seed: " << rand_seed << endl;
      pattern_random_oracle = true;
    }
  }
  else if((argc == 5) &&( strcmp(argv[1],"-R")) ==0 ){
    pattern_random_trace = true;
    //try to read in the random seed
    char* endp;
    rand_seed = strtol(argv[4],&endp,0);
    if(endp==argv[4] || *endp){
      cerr << "Failed to get random seed" << endl;
      print_helper();
      exit(0);
    }
    cout << "Random trace assignment pattern not currently supported" << endl;
  }
  else if((argc == 5) &&( strcmp(argv[1],"-SR")) ==0 ){
    pattern_random_oracle = true;
    //try to read in the random seed
    char* endp;
    rand_seed = strtol(argv[4],&endp,0);
    if(endp==argv[4] || *endp){
      cerr << "Failed to get random seed" << endl;
      print_helper();
      exit(0);
    }
  }
  else {
    print_helper ();
    exit (0);
  }
  
  formula_file = string (argv[2]);
  // cout << "input: " << input<< endl;
  // finitetracefile = string (argv[3]);
  trace_length = atoi(argv[3]);
  //cout << "trace_length "<< trace_length << endl;
  //sanity check the strings before saving
}
//////////////////////////////////////////////////////
// mltl::Formula Checker
// Recursively Checks a formula for
//intervals of the form [x -1]
//temporal operators without intervals
//intervals with components longer than the given trace length
// if pattern_median1 or pattern_median2 are true, check for intervals with widths greater than half the given trace length
void formula_checker(mltl::Formula *f){
  //base cases
  if(f->interval() != NULL && f->interval()->_right == -1){ //check for intervals of [x -1] form
    int temp = f->interval()->_left;
    f->interval()->_right = temp;
    f->interval()->_left = 0;
  }//end if
 
  if( f->oper() >= 6 && f->oper() <= 7 && f->interval() == NULL){ //check for temporal operators with undefined intervals
    f->interval()->_left = 0;
    f->interval()->_right = trace_length;
    cerr << "WARNING: Detected temporal operator with no declared bounds, assigning interval [0," << trace_length << "]" << endl;
  }
 
  if(f->interval() != NULL && f->interval()->_right > trace_length){
    cerr << "ERROR: mltl::Formula contains a bound greater than the completed trace length" << endl;
    cerr << "Max bound: " << f->max_bound() << " trace length: " << trace_length << endl;
    exit(1);
  }
  
  if(pattern_median1 || pattern_median2){//check for interval width over half the length of the total
    //   cerr << "entering median case of formula_checker" << endl;
    
    if(f->interval() == NULL){//if f has no interval, do nothing
    }
    else if(f->interval()->_right - f->interval()->_left > (trace_length / 2)){
      cerr << "WARNING: mltl::Formula contains an interval step size greater than half the total trace length. This may cause issues with obtaining median-satisfiability" << endl; 
    }
  }
  if(f->interval() !=NULL && f->interval()->_left > f->interval()->_right){
    cerr << "WARNING: Malformed interval with left bound greater than right bound detected" << endl;
  }
  //recursive case
  if(f->l_mf() != NULL){
    formula_checker(f->l_mf());
    }//end if
  
  if(f->r_mf() != NULL){
    formula_checker(f->r_mf());
  }//end if
  
}//end formula_checker

////////////////////////////////////////////////////////////
//Takes formula in from .mltl file 
vector<mltl::Formula*> formula_loader(string filename){
  vector<mltl::Formula*> formula_vector;
  mltl::Formula *f;
  ifstream mltl_file(filename);
  if(mltl_file.bad()){
    std::cerr << "Error: could not open " << filename << ".\n";
    print_helper();
  }
  if(mltl_file.is_open()){
    std::string line;
    while (mltl_file.peek() != EOF){
      std::getline(mltl_file,line);
      if(line[0] == '#' || line.length() == 0){
	//skip line if comment or blank
      }
      else{
	//trim the ';" from the end of the line;
	//	line.pop_back();
	f = mltl::Formula (line.c_str()).unique();
  //std::cerr << "Built" << endl;
	formula_checker(f);
  //std::cerr << "Checked." << endl;
	//implicit lower bound declaration detection
	if(f->interval() != NULL && f->interval()->_right == -1){
	  int temp = f->interval()->_left;
	  f->interval()->_right = temp;
	  f->interval()->_left = 0;
	}//end if
	formula_vector.push_back(f->clone());
	//TODO -- may need to add interval formalization here for temporal operators declared without interval
      } //end else
      // cout << "edge of while loop" << endl;
    }//end while loop
  }
  // cout << "formula loaded." << endl;
  mltl_file.close();
  //add check for close
  return formula_vector;
  }

/////////////////////////////////////////////////////
//Atom counter counts the number of atoms in a formula and sets the global atom count to that
//should be run only once at the beginning of the program to properly set max_atoms and num_atoms,
//so the trace can be constructed without a template file
void count_atoms(mltl::Formula *f){
  if(f->oper() != mltl::Formula::True && f->oper() != mltl::Formula::False){
     if(f->l_mf() != NULL){
      count_atoms(f->l_mf());
    }//end if
    if(f->r_mf() != NULL){
      count_atoms(f->r_mf());
    }//end if
  }//end if
 
  if (f->oper() > max_atoms){
    atom_names.push_back(mltl::Formula::get_name(f->oper()));
    max_atoms = f->oper();
    num_atoms = max_atoms - 8;
    //cout << num_atoms << endl;
  }
 
}//end atom counter



//////////////////////////////////////////////////////
//trace_generator(num_atoms,length) generates a template trace with a given number of atoms and a given length
//This is to replace the trace_loader function and move away from requiring a template trace file.
vector<vector<int>> trace_generator(int num_atoms, int length){
  vector<vector<int>> finite_trace;
  for(int i =0; i < length; i++){
    vector<int> finite_frame;
    for(int j =0; j < num_atoms; j++){
      finite_frame.push_back(0);
    }//end inner for
    finite_trace.push_back(finite_frame);
  }//end outer for
  return finite_trace;
}//end trace_generator
///////////////////////////////////////////////////////
//Recursively computes the worst-case propagation delay for a formula
int max_pd(mltl::Formula *f){
  if (f != NULL){
    if (f->interval() != NULL){
      //if this operator has a bound, add it to the worst of its subformulas
      int base = f->interval()->_right;
      return base + max(max_pd(f->l_mf()),max_pd(f->r_mf()));
    }
    else
      return max(max_pd(f->l_mf()),max_pd(f->r_mf()));
    //end else
  }//end if
  else
    return 0;
}//end max_pd

int min_pd(mltl::Formula *f){
  if (f != NULL){
    if (f->interval() != NULL){
      //if this operator has a bound, add it to the worst of its subformulas
      int base = f->interval()->_left;
      return base + min(min_pd(f->l_mf()),min_pd(f->r_mf()));
    }
    else
      return min(max_pd(f->l_mf()),min_pd(f->r_mf()));
    //end else
  }//end if
  else
  return 0;
}

//////////////////////////////////////////////////////
//Calls the z3 solver via the command line and temporary files via the SMTlibV2 standard
//returns 0 for unsat
//returns 1 for sat
//returns 2 for unknown
int  smt_SAT(mltl::Formula *f, std::vector<std::vector<int> > &trace_prime, int k){
//Compute the max propagation delay for the formula to know how far to iterate
unsigned long start = clock();
  //int maxdelay = max_pd(f);
  int maxdelay = min_pd(f) +1;

//To get the translation into z3_intputs.txt, the formula must first be written to a temp file where the translator 
//from Gokul Hariharan can translate the formula into the smt format. 
//TODO: Change the pathing options to be general and not specific to this machine
//Open the temp file
ofstream translation ("translation.txt");

//write the formula to it

//Preprocess the formula to remove the Weak operators
/*if(f->oper() == 10){
    mltl::Formula neg_a = mltl::Formula(mltl::Formula::Not,NULL,f->l_mf(),NULL);
    mltl::Formula* neg_a_point = &neg_a;
    mltl::Formula neg_b = mltl::Formula(mltl::Formula::Not,NULL,f->r_mf(),NULL);
    mltl::Formula* neg_b_point = &neg_b;
    mltl::Formula R = mltl::Formula(mltl::Formula::Release, neg_a_point, neg_b_point,f->interval());
    mltl::Formula* R_point = &R;
    mltl::Formula update = mltl::Formula(mltl::Formula::Not,NULL,R_point,NULL);
    f= &update;
  }*/
translation << f->to_string() << endl;
//close the tempfile
translation.close();
//exec the statement
int trans_exit_code = system("timeout 15s  ./Experiments/FPROGG/translator/main -f translation.txt -o z3_input.txt -t SMT");
  //There will be num_atoms atoms to be evaluated from the atom_names vector
  
  ofstream z3_inputs ("z3_input.txt", std::ios::app);
  if(z3_inputs.is_open()){
    //z3_inputs << " (declare-const t Int)" << endl;
    //z3_inputs << smt_model << endl;
    //z3_inputs << f->to_string() << endl;
    
    //Compile the satisfying trace by evaluating every atomic at every time
    //z3_inputs << "(get-model)\n";
    for(int i = 0; i <= maxdelay; i++){
    for(int j =0; j < atom_names.size(); j++){
      z3_inputs << "(eval (" << atom_names[j] << " " << i << " ))" << endl;
      //cout << "(eval (" << atom_names[i] << " t))" << endl;
    }//end atomic query loop
    
    }//end time loop
    z3_inputs.close();
  }
  else cerr << "unable to open file" << endl;

  //run z3
  int exit_code = system("timeout 15s z3 -smt2 ./z3_input.txt > ./z3_output.txt");
  if(exit_code == 124 ){
  cerr << "TIMEOUT" << endl;
  exit(1);
  }

  //parse the outputs back in
  vector<vector<int>> return_trace;
  string line;
  vector <string> results;
  ifstream z3_outputs ("z3_output.txt");
  if (z3_outputs.is_open()){
    while (getline (z3_outputs,line) ){
      //cout << line << endl;
      results.push_back(line);
    }
  }
  else cerr << "unable to open file" << endl;
  
  if (strcmp(results[0].c_str(),"sat") !=0){
    //a case that is either unsat or unknown
    if (strcmp(results[0].c_str(),"unsat") ==0){
      //a case that is unsat
      //cerr << "unsat" << endl;
      cerr << (float) (clock() - start)/CLOCKS_PER_SEC << "\t";
      return 0;
    }
     return 2;
     //cerr << "unkown" << endl;
     cerr << (float) (clock() - start)/CLOCKS_PER_SEC << "\t";
  }
  else if(strcmp(results[0].c_str(),"sat") == 0) {
    //sat case
    //iterate through the rest of the lines and get the verdicts

    int r = 1; // r tracks position in the results vector
    int j = 0; // j tracks time position in the resultant trace
    //cout << "RESULTS: " << endl;
    //cout << results << endl;
    string true_string = "true";
    string false_string = "false";
    while(r < results.size()){
      //skip over the returned model TODO: optionally debug store the returned model
      if(results[r] == "("){
        while(results[r] != ")"){
          r++;
        }//end while
        r++;
      }//end if

      vector<int> return_frame; //temporarily stores the assignments of atomics at a particular time 
      //extract the values for all the atoms at 1 time step
      for(int i = 0; i < atom_names.size(); i++){
	//if results[r] == "true"
	if (results[r].find(true_string) != string::npos) {
	  //push 1 to the ith spot in the frame
	  //cout << atom_names[i] << " true" << endl;
	  return_frame.push_back(1);
	}//end if
	
	//else if reslults[r] == "false"
	else if(results[r].find(false_string) != string::npos) {
	  //push 0 to the ith spot in the frame
	   //cout << atom_names[i] << " false" << endl;
	  return_frame.push_back(0);
	}//end elif
	
	//otherwise its unkown
	else {
	  //push 2 to the ith spot in the frame
	  //cout << results[r] << endl;
	  return_frame.push_back(2);
	}//end else
       	//increase r
	r = r + 1;
      }//end for
      	//cout << "RETURN FRAME: " << endl;
      	//cout << return_frame << endl;
	return_trace.push_back(return_frame);

	//	cout << "r = " << r << endl;
      j = j + 1;
      //push the filled in time step to the return trace
    }//end while
    
  }//end elif "sat"
  trace_prime = return_trace;
  //cout << trace_prime << endl;
  cerr << (float) (clock() - start)/CLOCKS_PER_SEC << "\t";
  return 1;
}

///////////////////////////////////////////////////////
//prog() takes in a formula f and returns a progressed formula
//recursively calls itself 
///////////////////////////////////////////////////////
mltl::Formula prog(mltl::Formula *f,  std::vector<std::vector<int> > finite_trace){
  //cout << "FINITE_TRACE" << endl;
  //cout << finite_trace << endl;
  //constructors
  mltl::Formula result;
  mltl::Interval *I;
  mltl::Formula *l_tmp;
  mltl::Formula *r_tmp;

   //if f has an interval , ensure it is in the correct order
  // cout << "f: " << f->to_string() << endl;
  // cout << "oper: " << mltl::Formula::get_name(f->oper()) << endl;
  if (f->interval() != NULL){
    assert (f->interval()->_right >= f->interval()->_left);
  }
  //Base Case: interval = 1
  if(finite_trace.size() == 1){
    switch (f->oper()){
      
    case mltl::Formula::True:
      // cout << "Returning true" << endl;
      return mltl::Formula(mltl::Formula::True, NULL,NULL,NULL);
      //unreachable
      break;
    case mltl::Formula::False:
      // cout << "Returning false" << endl;
      return mltl::Formula(mltl::Formula::False, NULL,NULL,NULL);
      //unreachable
      break;

      //Detangling: removing the (k)not
    case mltl::Formula::Not:{
      //  cout << "Not!" << endl;
      //cout << "Not! Recursing..." << endl;
      // cerr << "f:" << f->to_string() << endl;
      mltl::Formula temp;
      mltl::Formula *temp_point = &temp;
      if(f->l_mf() == NULL && f->r_mf() != NULL){
	      //	cerr << "left_sub null  and right defined, recursing" <<endl;
	      //	cerr << "right: " << f->r_mf()->to_string() << endl;
	      mltl::Formula temp = prog(f->r_mf(), finite_trace);
	      temp_point = &temp;
      }
      else if(f->l_mf() != NULL && f->r_mf() == NULL){
	    //	cerr << "right_sub null  and left defined, recursing" <<endl;
	    //	cerr << "left: " << f->l_mf()->to_string() << endl;
	      mltl::Formula temp = prog(f->l_mf(), finite_trace);
	      temp_point = &temp;
      }
      else {
	      //cerr << "prog error in not case" << endl;
      }
      
      // cerr << "l_mf(): " << f->l_mf()->to_string() << endl;
      // cerr << "r_mf(): " << f->r_mf()->to_string() << endl;

      //	cout << "not recurse returned" << endl;
      //	cout << temp_point->to_string() << endl;
	return mltl::Formula(mltl::Formula::Not,NULL, temp_point,NULL);
	//	  }
      //unreachable
     
      break;
    }
    case mltl::Formula::Or:
    case mltl::Formula::And:{
      // cout << "Or/And! Recursing..." << endl;
      // cout << mltl::Formula::get_name(f->oper()) << endl;
      mltl::Formula left_sub = prog(f->l_mf(), finite_trace);
      mltl::Formula right_sub = prog(f->r_mf(), finite_trace);
      mltl::Formula *left_sub_point = &left_sub;
      mltl::Formula *right_sub_point = &right_sub;
      return mltl::Formula(f->oper(),left_sub_point,right_sub_point,NULL);
      //return result;
      //unreachable
      break;
    }
    //case mltl::Formula::Weak_Until:
    case mltl::Formula::Until:
      // cout << "Until! Recursing..." << endl;
      if(f->interval()->_left > 0 && f->interval()->_left <= f->interval()->_right){
		    //cout << "1st branch" << endl;
        f->interval()->_left -= 1;
        f->interval()->_right -= 1;
	      //return result;
		    return mltl::Formula(mltl::Formula::Until,f->l_mf(),f->r_mf(),f->interval());
      }
      else if(f->interval()->_left == 0 && f->interval()->_left < f->interval()->_right){
		//cout << "2nd branch" << endl;
    //cout << f->to_string() << endl;
    //cout << "f->interval()->_left = " << f->interval()->_left << endl;
    //cout << "f->interval()->_right = " << f->interval()->_right << endl;
	//create the reduced right-side interval
	mltl::Interval J = mltl::Interval (0, f->interval()->_right - 1);
	mltl::Interval *K = &J;
	//check the pointers
	if(K->_left != 0 || K->_right != f->interval()->_right-1){
    cerr << "K-J Disagree! in prog() Until" <<endl;
	  cerr << "K = (" << K->_left << ", " << K->_right << ")\n";
	}
	
	mltl::Formula sub_3 = mltl::Formula (mltl::Formula::Until, f->l_mf(),f->r_mf(),K);
	mltl::Formula *sub_3_point = &sub_3;
	//cout << "progging l_mf()\t";
	mltl::Formula sub_2 = prog(f->l_mf()->unique(),finite_trace);
	mltl::Formula *sub_2_point = &sub_2;
     	mltl::Formula sub_1 = mltl::Formula (mltl::Formula::And, sub_2_point,sub_3_point,NULL);
	mltl::Formula *sub_1_point = &sub_1;
	//cout << "progging r_mf()" << endl;
	mltl::Formula sub_4 = prog(f->r_mf(),finite_trace);
	mltl::Formula *sub_4_point = &sub_4;
	//cout << "returning from prog" << endl;
       	return mltl::Formula(mltl::Formula::Or, sub_4_point,sub_1_point,NULL);
	//return result;
	
      }
      else if(f->interval()->_left == 0 && f->interval()->_right == 0){
	//*result = mltl::Formula (*prog(f->r_mf(), finite_trace));
	//	return result;
	//	cout << "3rd branch" << endl;
	return prog(f->r_mf(),finite_trace);
      }
      else{
	//unreachable
	cerr << "Unreachable state in prog() Until reached" << endl;
//	cerr << "f: " << f->to_string() << endl;
//	cerr << "f->interval()->right = " << f->interval()->_right << endl;
//	cerr << "f->interval()->leftt = " << f->interval()->_left << endl;
	exit(1);
      }
      break;
      
    //case mltl::Formula::Weak_Release:
    case mltl::Formula::Release:{
      // cout << "Release!" << endl;
      //   cerr << "f->l_mf(): " << f->l_mf() << "f->r_mf(): " << f->r_mf() << endl;
      mltl::Formula l_temp = mltl::Formula (mltl::Formula::Not,NULL,f->l_mf(),NULL);
      mltl::Formula *l_temp_point = &l_temp;
      // cerr << "l_temp built" << endl;
      mltl::Formula r_temp = mltl::Formula (mltl::Formula::Not,NULL,f->r_mf(),NULL);
      mltl::Formula *r_temp_point = &r_temp;
      // cerr << "r_temp built" << endl;
      
      mltl::Formula c_temp = mltl::Formula (mltl::Formula::Until,l_temp_point,r_temp_point,f->interval());
      mltl::Formula *c_temp_point = &c_temp;
 
      mltl::Formula progged = prog(c_temp_point,finite_trace);
      mltl::Formula *progged_point = &progged;
      // cout << "Returning from Release" << endl;
      mltl::Formula print_temp = mltl::Formula(mltl::Formula::Not,NULL,progged_point,NULL);
      mltl::Formula *print_temp_point = &print_temp;
      // cout << print_temp_point->to_string() << endl;
      return mltl::Formula (mltl::Formula::Not,NULL,progged_point,NULL);
     // unreachable
      cerr << "Unreachable state in prog() Until reached" << endl;
	    cerr << "f: " << f->to_string() << endl;
	    cerr << "f->interval()->right = " << f->interval()->_right << endl;
	    cerr << "f->interval()->leftt = " << f->interval()->_left << endl;
	    exit(1);
    }

    default:
      //Literal case located here due to how operator ids are assigned
      if(f->oper() > 12){

	      //return true iff p exists in pi[0] (is true at the start of the trace)
	      if(finite_trace[0][f->oper()-13] == 1){
	        return *mltl::Formula::TRUE();
	      }//end if
	
	    else if (finite_trace[0][f->oper()-13] == 1){	  
	      return *mltl::Formula::FALSE();
	    }//end else if
      else {
        return mltl::Formula (f->oper(), NULL,NULL,NULL);
      } //end else

	//unreachable
	cerr << "Unreachable state in literal case reached" << endl;
      }
      break;
      
      /*Globally and Finally cases are preemptively translated into U and R operators,
       so they are not addressed here*/
      
    }//end switch statement
  }//end base case

//recursive case            
else{
  //cout << "Finite trace size: " << finite_trace.size() << endl;
  //cout << "Recursing to reduce length" << endl;
  mltl::Interval J = mltl::Interval (0,1);
  mltl::Interval *K = &J;
 
  if(K->_left != 0 || K->_right != 1){
    cerr << "K-J Disagree! in prog() recursion" <<endl;
    cerr << "K = (" << K->_left << ", " << K->_right << ")\n";
  }
  vector <vector<int> > finite_zero;
  vector <vector<int> > finite_one;
  
  finite_zero.push_back(finite_trace[0]);
  int i = 1;
  for (i=1 ; i < finite_trace.size(); i++){
    finite_one.push_back(finite_trace[i]);
  }//end for
 
  mltl::Formula temp = prog(f,finite_zero);
  mltl::Formula *temp_point = &temp;
  // cout << "1st recurse returned" << endl;
  //cout << "temp: " << temp_point->to_string() << endl;
  //cout << "f = " << f->to_string() << endl;
  //mltl::Formula temp = temp_point->clone();
  //cout << "Finite Trace" << finite_trace << endl;
  return prog(temp_point,finite_one);
 } //end recursive case
  
  //unreachable
  //cerr << "unreachable state in prog() reached beyond recursive case.\n";
   return result;
}/*end prog*/

//////////////////////////////////////////////////////////
//Almost-Satisfiable Pattern function
//Line count listed in comments corresponds to the lines of the Algorithm presented in "MLTL Benchmark Generation via mltl::Formula Progression" paper

std::tuple< vector<vector<int>>,vector<bool>, mltl::Formula> Almost_SAT(mltl::Formula *f,  std::vector<std::vector<int> > in_trace,  vector<bool> Oracle, int start_time, int end_time){
  //cout << "Entering Almost-SAT" << endl;
   assert (f != NULL);
   
   vector<vector<int>> return_trace;
   vector<vector<int>> finite_trace;
   
  //create the trace_prime object to pass into SAT
  vector<vector<int>> trace_prime;
  int max_size = max_pd(f);
  if(smt_SAT(f, finite_trace, 0) == 0){
    //if the 1st frame is not satisfiable, then the pattern has failed and a trace cannot be generated
    pattern_fail = true;
    cerr << "ERROR: PATTERN FAILURE. UNABLE TO COMPLETE TRACE GENERATION" << endl;
    //return the model, empty trace, and oracle
    return_trace = finite_trace;
    return std::make_tuple(return_trace,Oracle,*f); 
  }
  int k = start_time+1;
  mltl::Formula fup;
  mltl::Formula *fup_point;
  Oracle.push_back(true);
  //Line 5
  /*cout << "Entering main loop" << endl;
  cout << "start_time = " << start_time << endl;
  cout << "end_time = " << end_time << endl;
  cout << "finite_trace.size(): " << finite_trace.size() << endl;*/
   while (start_time <= k && k < finite_trace.size()){
    mltl::Formula *fcopy = f->clone();
    //Line 6:
    //cout << "finite_trace.size()" << finite_trace.size() << "\t trace_length:" << trace_length << endl;
    if(finite_trace.size() > trace_length){
      //Line 7: return f,empty trace, oracle
      //cout << "trace length exceeds given bound. Complete" << endl;
      return_trace = finite_trace;
      return std::make_tuple(return_trace,Oracle,fup);
    }//Line 8: end if
    //Line 9: phi' = prog(phi, pi_k)
    //create the trace suffix finite trace starting at k
    vector <vector<int>> finite_k;
    int k_size = start_time+finite_trace.size();
    /*cout << "k = " << k << endl;
    cout << "k_size = " << k_size << endl;
    cout << "i-start_time = " << k - start_time << endl; */
    for(int i =k; i < k_size; i++){
      finite_k.push_back(finite_trace[i-start_time]);
    }
    //cout << "finite_trace: " << finite_trace << endl;
    //cout << "finite_k: " << finite_k << endl;
    //cout << "f: " << f->to_string() << endl;
    /*cout << "Proggers" << endl;*/
    fup = prog(fcopy, finite_k);
    fup_point = &fup;
    //cerr << "fup: " << fup_point->to_string() << endl;
    //Line 10
        if(smt_SAT(fup_point, trace_prime, k) == 0){
      //push (k,false) into oracle stream
	  Oracle.push_back(false);
	}

    //Line 12
    else{
      bool conflict_assign = false;
      //Line 14: pi = pi * pi'<-- not actually multiplication
      /*cout << "Concatenating..." << endl;
      cout << "Finite Trace Size: " << finite_trace.size() << endl;
      cout << "Trace Prime Size: " << trace_prime.size() << endl;
      cout << "k: " << k << endl;
      cout << "Finite Trace: " << finite_trace << endl;
      cout << "Trace Prime: " << trace_prime << endl;*/

      //copy(trace_prime.begin(), trace_prime.end(), finite_trace.begin()+(k-start_time));
      for(int i =0; i < trace_prime.size(); i++){
        //cout << "i = " << i << endl;
        //cout << "trace prime: " << trace_prime << endl;
        if(i< finite_trace.size()-k){
          //Assign the variable assignments, unless the element in prime is 2 and in trace is not 2
          //If the value in trace is 0 and in prime is 1, or vice-versa, then we have a bigger issue
          for(int j = 0; j < finite_trace[k].size(); j++){
            if(trace_prime[i][j] == 2 && finite_trace[k+i][j] != 2){
              //Do not overwrite known values with unknown ones
            }
            else if  (finite_trace[k+i][j] !=2 && trace_prime[i][j]!=2 && finite_trace[k+i][j] != trace_prime[i][j]){
              //cerr << "DEBUG: Error overwriting \t";

              //Oracle.push_back(false);
              conflict_assign = true;
            }
            else{
              finite_trace[k+i] = trace_prime[i];
            }//end else
          }
        }
        else{
          finite_trace.push_back(trace_prime[i]);
        }
      }//end for loop
      //cout << "updated finite trace:" << finite_trace << endl;
      //cout << "Concatenated" << endl;

    //cout << "finite_trace[" << k << "] = " << trace_prime[0] << endl;
      //Line 15: push (k,true) into oracle stream
      if(!conflict_assign){
        Oracle.push_back(true);
      }
      else{
        Oracle.push_back(false);
      } 
    }//Line 16: end else
	 //cout << "k: " << k << endl;
    k++; //increment k
   
    
  }//Line 17: end while
  //cout << "While loop exited" << endl;
  //Line 18: return the model
  return std::make_tuple(finite_trace,Oracle,fup);
}

//////////////////////////////////////////////////////////
///Almost UNSAT pattern
std::tuple< vector<vector<int>>,vector<bool>, mltl::Formula> Almost_UNSAT(mltl::Formula *f,  std::vector<std::vector<int> > finite_trace,  vector<bool> in_Oracle,int start_time, int end_time){
   vector<vector<int>> unsat_trace;
  //negate the formula
  mltl::Formula negation = mltl::Formula(mltl::Formula::Not,NULL,f,NULL);
  mltl::Formula* negation_point = &negation;

  vector<vector<int>> sat_trace;
  vector<bool> Oracle;
  //run Almost_SAT
  auto sat_result = Almost_SAT(negation_point,finite_trace,in_Oracle,start_time, end_time);
  sat_trace = std::get<0>(sat_result);
  Oracle = std::get<1>(sat_result);
  mltl::Formula sat_progged = std::get<2>(sat_result);
  //  cout << "unsat_trace: " << sat_trace << endl;
  //flip the output oracle values
  //unsat_oracle = {(k,true)|(k,false) E Oracle} U {(k,false)|(k,true) E Oracle}
  int i =0;
  for (i=0;i<Oracle.size();i++){
    if (Oracle[i]){
      Oracle[i] = false;
      //cout << "reassinging Oracle" << endl;
    }//end if
  else if (!Oracle[i]){
      Oracle[i] = true;
      // cout << "reassinging Oracle" << endl;
   }
    else
      cerr << "Error in updating Oracle for Unsat pattern at timestep "<< i << endl;
  }//end for
  
  return std::make_tuple(sat_trace,Oracle,sat_progged);
}
/////////////////////////Median SAT Patern
std::tuple< vector<vector<int>>,vector<bool>, mltl::Formula> Median_SAT(mltl::Formula *f,  std::vector<std::vector<int> > finite_trace,  vector<bool> in_Oracle,int order, int max_time){
  vector<vector<int>> sat_trace,unsat_trace,med_trace;
  vector<bool> sat_Oracle,unsat_Oracle,med_Oracle;
  mltl::Formula sat_form,unsat_form;
  mltl::Formula* sat_form_point;
  mltl::Formula* unsat_form_point;
  mltl::Formula return_form;

  //determine the sizes of the half-length traces
  int start_time = 0;
  int first_length = max_time /2;
  int second_length = max_time /2;
  if(max_time % 2 != 0){
    second_length++;
  }
  //Run the first pattern, as given by the order parameter
  if( order == 1 ){
    //  cout << "Calling Almost SAT with start time: " << start_time << " and end time: " << first_length << endl;
    tie(sat_trace,sat_Oracle,sat_form) = Almost_SAT(f, finite_trace, in_Oracle,start_time,first_length);
    if (pattern_fail){
      cerr << "PROVIDED FORMULA IS ARBITRARILY FALSE" << endl;
      return std::make_tuple(sat_trace,sat_Oracle,sat_form);
    }
    //run Almost_UNSAT
    //store the trace and oracle object
    sat_form_point = &sat_form;
    // cout << "Calling Almost UNSAT with start time: " << (first_length+1) << " and end time: " << max_time << endl;
    tie(unsat_trace,unsat_Oracle,unsat_form) = Almost_UNSAT(sat_form_point, finite_trace, in_Oracle, (first_length) ,max_time);
    // cout << "unsat_trace: " << sat_trace << endl;
    if (pattern_fail){
      cerr << "PROGRESSED FORMULA AT TIMESTEP " << (first_length+1) << " IS ARBITRARILY TRUE" << endl;
      return std::make_tuple(sat_trace,sat_Oracle,unsat_form);
    }
    //cout << "Unsat Trace:\n" << unsat_trace << endl;
    unsat_form_point = &unsat_form;
    sat_trace.insert( sat_trace.end(),unsat_trace.begin(),unsat_trace.end() );
    sat_Oracle.insert( sat_Oracle.end(),unsat_Oracle.begin(),unsat_Oracle.end() );
    med_trace = sat_trace;
    med_Oracle = sat_Oracle;
    return_form = unsat_form;
  }
  else{
    //Run almost UNSAT
    // cout << "Calling Almost UNSAT with start time: " << start_time << " and end time: " << first_length << endl;
    tie(unsat_trace,unsat_Oracle,unsat_form) = Almost_UNSAT(f, finite_trace, in_Oracle, start_time, first_length );
    
    if (pattern_fail){
      cerr << "PROVIDED FORMULA IS ARBITRARILY TRUE" << endl;
      return std::make_tuple(unsat_trace,unsat_Oracle,unsat_form);
    }
   
    unsat_form_point = &unsat_form;
    //  cout << "unsat formula: " << unsat_form_point->to_string() << endl;
    //  cout << "unsat oracle:\n";
    //Run Almost_SAT
    // cout << "Calling Almost SAT with start time: " << (first_length+1) << " and end time: " << max_time << endl;
    tie(sat_trace,sat_Oracle,sat_form) = Almost_SAT(unsat_form_point, finite_trace, in_Oracle,(first_length), max_time);
    if (pattern_fail){
          cerr << "PROGRESSED FORMULA AT TIMESTEP " << (first_length+1) << " IS ARBITRARILY FALSE" << endl;
      return std::make_tuple(sat_trace,sat_Oracle,sat_form);
    }
   
    unsat_trace.insert( unsat_trace.end(),sat_trace.begin(),sat_trace.end() );
    unsat_Oracle.insert( unsat_Oracle.end(),sat_Oracle.begin(),sat_Oracle.end() );
    med_trace = unsat_trace;
    med_Oracle = unsat_Oracle;
    return_form = sat_form;
  }
  int total = 0;
  int sat_count = 0;
  int unsat_count = 0;
  //generate the new oracle
  vector<vector<int>> frame_results;
  frame_results.push_back(finite_trace[0]);
  //  percentage_reporter(med_Oracle);
  return std::make_tuple(med_trace,med_Oracle,return_form);
}
////////////////////////////////////////////////Random Pattern
std::tuple< vector<vector<int>>,vector<bool>, mltl::Formula>  Rand_SAT(mltl::Formula *f,  std::vector<std::vector<int> > finite_trace,  vector<bool> in_Oracle){
 
  cerr << "Random pattern not yet implemented." << endl;
  vector<vector<int>> rand_trace;
  vector<bool> Oracle;
  //build the random trace
   for(int j =0; j < trace_length; j++){
    vector<int> frame;
    for(int i = 0; i < num_atoms; i++){
      frame.push_back(rand()%2);
    }
    rand_trace.push_back(frame);
    cerr << "Frame: ";
    cerr << frame << endl;
  }

  //check its satisfiability
  int total = 0;
  int sat_count = 0;
  int unsat_count = 0;
  int unkow_count = 0;
 
  
  for(int i = 0; i < rand_trace.size(); i++){ 
    //inefficient way to look only at the finite trace starting at time index i
    vector <vector<int>> finite_i;
    for(int j =i; j < rand_trace.size(); j++){
      finite_i.push_back(rand_trace[j]);
    }
    
   int result = smt_SAT(f,finite_i, finite_i.size());
    if (result == 1){
      Oracle.push_back(true);
      sat_count++;
    }
    else if (result == 0){
      Oracle.push_back(false);
      unsat_count++;
    }
    else{ //unkown case - treating as unsatisfiable, though this is something that could be made configureable
      Oracle.push_back(false);
      unkow_count++;
    }
    total++;
  }//end new oracle generation for loop
  //cout << sat_count << " frames satisfiable of " << total << " total frames" << endl;
  return std::make_tuple(rand_trace,Oracle,*f);

}// end rand trace sat

//////////////////////////////////////////////////////////
//Random Oracle Pattern - attempts to build a trace who's satisfiability shifts randomly between timesteps
std::tuple< vector<vector<int>>,vector<bool>, mltl::Formula> Semi_Rand_Oracle_SAT(mltl::Formula *f,  std::vector<std::vector<int> > finite_trace,  vector<bool> in_Oracle, int start_time, int end_time){
  //cerr << "Entering Semi Random pattern" << endl;
  mltl::Formula *g = f;
   vector<bool> Oracle;
   vector<vector<int>> temp_trace;
   //randomly choose opening pattern
   bool sat_choice;
   bool unsat_choice;
   if(rand()%2 ==0){
     //almost sat pattern
     sat_choice = true;
     //  cerr << "Almost SAT chosen" << endl;
   }
   else{
     //almost unsat pattern
     unsat_choice =true;
     //  cerr << "Almost Unsat chosen" << endl;
   }
   //cerr << "Entering main decision loop" << endl;
   for(int i = start_time; i < end_time; i++){
     // cerr << "i = " << i << endl;
    //end interval
     //  cerr << "endtime-i:" << end_time-i << endl;
    int end_interval = i + (rand()%(end_time-i));
    if((end_interval == i) || (end_interval+1 == end_time) || (end_interval - i <= 1)){
      //  cerr << "ending interval updated from " << end_interval;
      
      end_interval++;
      if (end_interval-i ==1){
	end_interval++;
      }
    }
    // cerr << "to ending interval time of " << end_interval << endl;
    //  cerr << "g: " << g->to_string() << endl;
    //execute the desired pattern
    if(sat_choice){
      
      auto result = Almost_SAT(g, finite_trace,in_Oracle,i,end_interval);
      
      vector<vector<int>> temp = get<0>(result);
      temp_trace.insert( temp_trace.end(),get<0>(result).begin(),get<0>(result).end() );
      Oracle.insert(Oracle.end(),get<1>(result).begin(),get<1>(result).end() );
      // g = &std::get<mltl::Formula>(result);
       mltl::Formula temp_form;
      mltl::Formula *temp_form_point = &std::get<mltl::Formula>(result);
      g = temp_form_point->clone();
      
    }
    else{
      
      auto result = Almost_UNSAT(g, finite_trace,in_Oracle,i,end_interval);
      vector<vector<int>> temp = get<0>(result);
      temp_trace.insert( temp_trace.end(),get<0>(result).begin(),get<0>(result).end() );
      Oracle.insert(Oracle.end(),get<1>(result).begin(),get<1>(result).end() );
      mltl::Formula temp_form;
      mltl::Formula *temp_form_point = &std::get<mltl::Formula>(result);
      g = temp_form_point->clone();
      
    }
    if (pattern_fail){
      cerr << "PROVIDED FORMULA IS ARBITRARILY TRUE/FALSE" << endl;
      return std::make_tuple(temp_trace,Oracle,*g);
    }
    
    //check if now arbitrary - or dont
    //otherwise flip pattern
    if(sat_choice){
      sat_choice = false;
    }
    else{
      sat_choice = true; 
    }
    //update current time
    i=end_interval-1;
  }//end for loop
  //cerr << "for loop ended" << endl;
   return std::make_tuple(temp_trace,Oracle,*g);
}//end semi rand oracle sat

//////////////////////////////////////////////////////////
/////////////       Main     /////////////////////////////
//////////////////////////////////////////////////////////
int main(int argc, char** argv){
  //Parse the runtime arguments
  parse_args (argc, argv);
  //load the formula and check not null
  vector<mltl::Formula* > formula_vector = formula_loader(formula_file);
  string extension = ".mltl";
  formula_file.erase(formula_file.find_last_not_of(extension)+1);
  //for each formula in the vector, run the declared pattern
  srand(rand_seed);
  unsigned long start,end;
  //cout << "Clocks: " << CLOCKS_PER_SEC << endl;
  for(int i =0; i < formula_vector.size(); i++){
    start = clock();
    //cout << "Entering loop" << endl;
    mltl::Formula *f = formula_vector[i];
    assert (f != NULL);
    //cout << f->to_string() << endl;
    //determine its number of atoms
    max_atoms = 8;
    atom_names.clear();
    count_atoms(f);
    //generate the finite trace array
    vector<vector<int>> finite_trace = trace_generator(num_atoms, trace_length);
    //create the oracle object
    vector<bool> in_Oracle;
    
    if(pattern_sat){
      //cout << "Running almost sat pattern" << endl;
      auto result = Almost_SAT(f, finite_trace, in_Oracle,0, trace_length);
      //model_printer(result);
      output_result_printer(result,i,"Almost Satisfiable",formula_file);
      output_trace_printer(result,i,formula_file);
    }
    
    else if(pattern_unsat){
      //cout << "Running almost unsat pattern" << endl;
      auto result = Almost_UNSAT(f, finite_trace, in_Oracle,0, trace_length);
      //model_printer(result);
      output_result_printer(result,i,"Almost Unsatisfiable",formula_file);
      output_trace_printer(result,i,formula_file);
    }
    
    else if(pattern_median1){
      //cout << "Running median pattern with Almost_SAT in first half of time" << endl;
      int order = 1;
      auto result = Median_SAT(f,finite_trace,in_Oracle,order, trace_length);
      //model_printer(result);
      output_result_printer(result,i,"Median Satisfiable",formula_file);
      output_trace_printer(result,i,formula_file);
    }
    
    else if(pattern_median2){
      //cerr << "Running median pattern with Almost_UNSAT in first half of time" << endl;
      int order = 2; 
      auto result = Median_SAT(f,finite_trace,in_Oracle,order, trace_length);
      //model_printer(result);
      output_result_printer(result,i,"Median Satisfiable",formula_file);
      output_trace_printer(result,i,formula_file);
    }
    
    else if(pattern_random_trace){
      //cerr << "Running random pattern with random seed: " << rand_seed << endl;
      auto result = Rand_SAT(f,finite_trace,in_Oracle);
      //model_printer(result);
      output_result_printer(result,i,"Random Trace",formula_file);
      output_trace_printer(result,i,formula_file);
    }
    else if(pattern_random_oracle){
      auto result =  Semi_Rand_Oracle_SAT(f,finite_trace,in_Oracle,0, trace_length);
      //model_printer(result);
      output_result_printer(result,i,"Semi-Random Desired Oracle",formula_file);
      output_trace_printer(result,i,formula_file);
    }
    
    else{
      //should be unreachable, if no pattern is declared execution should stop at argument parsing
      cerr << "No pattern type found" << endl;
    }
    cerr << "\n";
    cout <<"Formula "<< i << " | Benchmark Time: " << (float) (clock() - start)/CLOCKS_PER_SEC << endl;
    
  }//end for loop
  return 0;
}
