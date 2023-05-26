#ifndef CODE_bamRemoveDuplicates
#define CODE_bamRemoveDuplicates
#include <string>
#include "Parameters.h"

using namespace std;
int funCompareNames(const void *a, const void *b);
uint32 funStartExtendS(const uint32* const p);
uint32 funCigarExtendS(const uint32* const p, uint32* cout);
int funCompareCigarsExtendS(const uint32* const pa, const uint32* const pb);
int funCompareCoordFlagCigarSeq(const void *a, const void *b);
void bamRemoveDuplicates(const string bamFileName, const string bamFileNameOut, Parameters &P);

#endif
