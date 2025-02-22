#include "BioAlphabet.h"

using namespace std;

BioAlphabet::BioAlphabet(const BioSequence s) {
  seq_type_ = s;
  char_map_.assign(256, -1);
  inv_char_map_.assign(256, -1);
  srand(time(NULL));
  switch (s)  {
    case PROT: InitProt();  break;
    case DNA: InitDNA();  break;
    case RNA: InitRNA();  break;
    default: std::cout << "Unknown BIOSEQUENCE: the program only supports PROT, DNA, or RNA." << std::endl; break;
  }
  return;
}

BioAlphabet::BioAlphabet(ReducedAlphabet &reduced)  {
  seq_type_ = OTHER;
  char_map_.assign(256, -1);
  inv_char_map_.assign(256, -1);
  srand(time(NULL));
  unordered_map<char, bool> nr_char;
  for(auto it = reduced.alphabet_map_.begin(); it != reduced.alphabet_map_.end(); ++ it) {
    char c = (char) (65 + it->second);
    char_map_[c] = nr_char.size(); inv_char_map_[nr_char.size()] = c; nr_char[c] = true;
  }
  alphabet_size_ = nr_char.size();
  return;
}

void BioAlphabet::InitProt(void) {
  alphabet_size_ = 25;
  num_bits_ = 5;
  char_map_[(int) '*'] = 0; inv_char_map_[0] = (int) '*';
  char_map_[(int) 'A'] = 1; inv_char_map_[1] = (int) 'A';
  char_map_[(int) 'B'] = 2; inv_char_map_[2] = (int) 'B';
  char_map_[(int) 'C'] = 3; inv_char_map_[3] = (int) 'C';
  char_map_[(int) 'D'] = 4; inv_char_map_[4] = (int) 'D';
  char_map_[(int) 'E'] = 5; inv_char_map_[5] = (int) 'E';
  char_map_[(int) 'F'] = 6; inv_char_map_[6] = (int) 'F';  
  char_map_[(int) 'G'] = 7; inv_char_map_[7] = (int) 'G';
  char_map_[(int) 'H'] = 8; inv_char_map_[8] = (int) 'H';
  char_map_[(int) 'I'] = 9; inv_char_map_[9] = (int) 'I';
  char_map_[(int) 'J'] = 10; inv_char_map_[10] = (int) 'J';
  char_map_[(int) 'K'] = 11; inv_char_map_[11] = (int) 'K';  
  char_map_[(int) 'L'] = 12; inv_char_map_[12] = (int) 'L';
  char_map_[(int) 'M'] = 13; inv_char_map_[13] = (int) 'M';
  char_map_[(int) 'N'] = 14; inv_char_map_[14] = (int) 'N';
  char_map_[(int) 'P'] = 15; inv_char_map_[15] = (int) 'P';
  char_map_[(int) 'Q'] = 16; inv_char_map_[16] = (int) 'Q';
  char_map_[(int) 'R'] = 17; inv_char_map_[17] = (int) 'R';
  char_map_[(int) 'S'] = 18; inv_char_map_[18] = (int) 'S';
  char_map_[(int) 'T'] = 19; inv_char_map_[19] = (int) 'T';
  char_map_[(int) 'V'] = 20; inv_char_map_[20] = (int) 'V';
  char_map_[(int) 'W'] = 21; inv_char_map_[21] = (int) 'W';
  char_map_[(int) 'X'] = 22; inv_char_map_[22] = (int) 'X';
  char_map_[(int) 'Y'] = 23; inv_char_map_[23] = (int) 'Y';
  char_map_[(int) 'Z'] = 24; inv_char_map_[24] = (int) 'Z';
  
  return;
}
void BioAlphabet::InitDNA(void) {
  alphabet_size_ = 5;
  num_bits_ = 3;
  char_map_[(int) 'A'] = 0; inv_char_map_[0] = (int) 'A';
  char_map_[(int) 'C'] = 1; inv_char_map_[1] = (int) 'C';
  char_map_[(int) 'G'] = 2; inv_char_map_[2] = (int) 'G';
  char_map_[(int) 'N'] = 3; inv_char_map_[3] = (int) 'N';
  char_map_[(int) 'T'] = 4; inv_char_map_[4] = (int) 'T';
    //ambiguity code
  char_map_[(int) 'R'] = 5; inv_char_map_[5] = (int) 'R';/*A or G */
  char_map_[(int) 'Y'] = 6; inv_char_map_[6] = (int) 'Y';/*C or T*/ 
  char_map_[(int) 'W'] = 7; inv_char_map_[7] = (int) 'W';/*A or T*/ 
  char_map_[(int) 'S'] = 8; inv_char_map_[8] = (int) 'S';/*C or G*/ 
  char_map_[(int) 'K'] = 9; inv_char_map_[9] = (int) 'K';/*G or T*/ 
  char_map_[(int) 'M'] = 10; inv_char_map_[10] = (int) 'M';/*A or C*/ 
  char_map_[(int) 'B'] = 11; inv_char_map_[11] = (int) 'B';/*C or G T*/ 
  char_map_[(int) 'D'] = 12; inv_char_map_[12] = (int) 'D';/*A G or T*/ 
  char_map_[(int) 'H'] = 13; inv_char_map_[13] = (int) 'H';/*A C or T*/ 
  char_map_[(int) 'V'] = 14; inv_char_map_[14] = (int) 'V';/*A C or G*/ 
  return;
}

void BioAlphabet::InitRNA(void) {
  alphabet_size_ = 5;
  num_bits_ = 3;
  char_map_[(int) 'A'] = 0; inv_char_map_[0] = (int) 'A';
  char_map_[(int) 'C'] = 1; inv_char_map_[1] = (int) 'C';
  char_map_[(int) 'G'] = 2; inv_char_map_[2] = (int) 'G';
  char_map_[(int) 'N'] = 3; inv_char_map_[3] = (int) 'N';
  char_map_[(int) 'U'] = 4; inv_char_map_[4] = (int) 'U';
  return;
}


bool BioAlphabet::isUnambiguous(char c)
{
    switch(c)
    {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return true;
        default:
            return false;
    }
}

// Returns true if c is a valid ambiguity code
bool BioAlphabet::isAmbiguous(char c)
{
    switch(c)
    {
        case 'M':
        case 'R':
        case 'W':
        case 'S':
        case 'Y':
        case 'K':
        case 'V':
        case 'H':
        case 'D':
        case 'B':
        case 'N':
            return true;
        default:
            return false;
    }
}

// Returns true if c is a valid symbol in this alphabet
bool BioAlphabet::isValidLetter(char c)
{
    return isUnambiguous(c) || isAmbiguous(c);
}

//
std::string BioAlphabet::getPossibleSymbols(char c)
{
    switch(c)
    {
        case 'A':
            return "A";
        case 'C':
            return "C";
        case 'G':
            return "G";
        case 'T':
            return "T";
        case 'M':
            return "AC";
        case 'R':
            return "AG";
        case 'W':
            return "AT";
        case 'S':
            return "CG";
        case 'Y':
            return "CT";
        case 'K':
            return "GT";
        case 'V':
            return "ACG";
        case 'H':
            return "ACT";
        case 'D':
            return "AGT";
        case 'B':
            return "CGT";
        case 'N':
            return "ACGT";
        default:
            return "";
    }
}


bool BioAlphabet::IsValid(const char c) {
  return (char_map_[(int) c] != -1);
}


char BioAlphabet::RandomChar()  {
  int pos = rand() % alphabet_size_;
  return (char) inv_char_map_[pos];
}
