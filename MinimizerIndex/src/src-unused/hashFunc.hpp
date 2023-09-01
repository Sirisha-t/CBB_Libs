/**
* =====================================================================================
*
*       Filename:  hashutil.h
*
*    Description:
*
*        Version:  1.0
*        Created:  04/18/2016 04:49:32 PM
*       Revision:  none
*       Compiler:  gcc
*
*
*         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
*                  Rob Patro (rob.patro@cs.stonybrook.edu)
*                  Rob Johnson (rob@cs.stonybrook.edu)
*   Organization:  Stony Brook University
*   Edited by: Mohamed Abuelanin (mabuelanin@gmail.com) UC Davis
*
* =====================================================================================
*/

#ifndef _HASHUTIL_H_
#define _HASHUTIL_H_

#include <sys/types.h>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include "Utils/kmer.h"
#include <algorithm>
#include <iostream>


using namespace std;

/* ---------------------
Class Hasher: Parent class
--------------------- */

class Hasher{
private:
    uint64_t seed;
public:
    explicit Hasher(uint64_t Iseed) { seed = Iseed; }

    Hasher *clone() override { return new Hasher(seed); }

    uint64_t hash(const string & Skey) { return 0; };

    ~MumurHasher(){}
};
