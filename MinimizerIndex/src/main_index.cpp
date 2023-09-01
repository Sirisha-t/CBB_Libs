#include "seq_index.h"
#include "loader.h"
#include "util_func.h"
#include "bio_alphabet.h"
#include "minimizer.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <ctime>
#include <cmath>
#include <omp.h>

using namespace std;

static string workspace_dir = ".";
static string target_file;
static string query_file;
static string seq_type;
static int num_blocks;
static int num_threads;
static int kmer_size;
static int win_size;
bool dump_idx = false;


int main(int argc, char** argv)
{
  boost::program_options::options_description desc("Options for indexing");
  desc.add_options()
      ("help", "print the help message")
      ("seq_type", boost::program_options::value<string>(&seq_type), "[required] type of sequence (nucl/prot)")
      ("target_file", boost::program_options::value<string>(&target_file), "[required] target sequencing reads (in FASTQ/FASTA)")
      ("query_file", boost::program_options::value<string>(&query_file), "[required] query sequencing reads (in FASTQ/FASTA)")
      ("k", boost::program_options::value<int>(&kmer_size)->default_value(15), "[optional] kmer size")
      ("w", boost::program_options::value<int>(&win_size)->default_value(5), "[optional] window size (usually 2/3 of kmer size.")
      ("work_space", boost::program_options::value<string>(&workspace_dir)->default_value(workspace_dir), "[optional] working directory for indexing file dump")
      ("num_blocks", boost::program_options::value<int>(&num_blocks)->default_value(8), "[optional] num. of index blocks (should match the number of threads to be used during search)")
      ("num_threads", boost::program_options::value<int>(&num_threads)->default_value(8), "[optional] num. of threads to use")
      ("dump_idx", boost::program_options::value<bool>(&dump_idx)->default_value(false), "[optional] option to dump index to file")
  ;

  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::parse_command_line(argc, argv, desc), vm
  );

  boost::program_options::notify(vm);
    if(vm.count("help"))  {
      cout << endl;
      cout << "Usage: ./seq_index [seq_type (\"nucl\" or \"prot\")]  [target_file (FASTQ/A)] [query_file (FASTQ/A)]" << endl << endl;
      cout << desc << endl;
      return 0;
    }
    // check options validity
    if(seq_type != "nucl" && seq_type != "prot")  {
      cout << "Error: The specified seq_type is not supported (currently we only support \"nucl\" and \"prot\"):\t" << seq_type << endl;
      cout << "Please use \'--help\' for more details." << endl;
      exit(0);
    }
    
    if(!boost::filesystem::exists(target_file))  {
      cout << "Error: target_file does not exist:\t" << target_file << endl;
      cout << "Please use \'--help\' for more details." << endl;
      exit(0);
    }
    if(!boost::filesystem::exists(query_file))  {
      cout << "Error: query_file does not exist:\t" << query_file << endl;
      cout << "Please use \'--help\' for more details." << endl;
      exit(0);
    }
    boost::filesystem::path abs_workspace = workspace_dir;
    if(!boost::filesystem::is_directory(workspace_dir))  {
      cout << workspace_dir << endl;
      cout << "Error: working space does not exist (please provide full path)." << endl;
      cout << "Please use \'--help\' for more details." << endl;
      exit(0);
    }
    if(num_blocks <= 0)  {
      cout << "Error: num_blocks should be set to a positive value." << endl;
      cout << "Please use \'--help\' for more details." << endl;
      exit(0);
    }

  //cout << "============================================================" << endl;
  //cout << "Indexing target sequences.." << endl;
  double start_time;
  double current_time;
  BioAlphabet alphabet;
  if(seq_type == "nucl")  alphabet = BioAlphabet(DNA);
  else alphabet = BioAlphabet(PROT);
  // load in sequence from file for forward sequences
  UtilFunc util;
  string target_stem = util.GetFileStem(target_file);
  //string query_stem = util.GetFileStem(query_file);

  start_time = util.MyTime();
  MMBuild target_seq(alphabet);
  //MMBuild query_seq(alphabet);

  std::vector<FastaInfo> Target;
  std::vector<FastaInfo> Query;

  if(seq_type == "nucl" || seq_type == "prot"){
    target_seq.LoadSequences(target_file, false);
    //Query = query_seq.LoadSequences(query_file, false);
  }
  else{
    cout<<"Error: Invalid input sequence type. Please input nucl or prot sequences\n\n";
    exit(0);
  }

  current_time = util.MyTime();
  util.PrintElapsed(start_time, current_time, "Loaded Sequences in ");
  cout << "============================================================" << endl;
  //target_seq.PrintAllSeqsInfo();

  /** Index target sequences **/
  start_time = util.MyTime();
  std::cout<<"Building minimizer index..\n";
  target_seq.BuildIndex(kmer_size, win_size,workspace_dir, target_stem );

  if(dump_idx)
  {
    start_time = util.MyTime();
    target_seq.DumpIndextoFile(workspace_dir, target_stem);
    current_time = util.MyTime();
    util.PrintElapsed(start_time, current_time, "Dumped index to file in ");
  }

  current_time = util.MyTime();
  util.PrintElapsed(start_time, current_time, "Built Index in ");

  start_time = util.MyTime();
  
  //Load index from file
  //TODO:: Load only one entry at a time and do mapping
  //target_seq.LoadIndex(workspace_dir, target_stem);

  current_time = util.MyTime();
  util.PrintElapsed(start_time, current_time, "Loaded index in ");
  cout << "============================================================" << endl;

  /** Map query sequence to index **/
  start_time = util.MyTime();
  std::cout<<"Mapping query to target sequences..\n";

  int nreads = target_seq.GetNumSeqs();
  //Map each query seq to the minimizer indices  
  for (int i=0; i < nreads; i++){
    target_seq.MapQuery(kmer_size, win_size, i);  
  }

  current_time = util.MyTime();
  util.PrintElapsed(start_time, current_time, "Mapped Query sequences in ");
  cout << "============================================================" << endl;

return 0;
}
