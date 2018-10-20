#ifndef T_CLASSIFY_H
#define T_CLASSIFY_H
//t_classify.h
//utility function for transcript classification vs reference
//

#include "gff.h"

int classcode_rank(char c);

char getOvlCode(GffObj& m, GffObj& r, int& ovlen); //returns class code

bool singleExonTMatch(GffObj& m, GffObj& r, int& ovlen);

#endif
