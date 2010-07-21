#ifndef FLOWPATH_H_
#define FLOWPATH_H_

#include <list>

struct cp	{
	int row;
	int col;
};

typedef struct cp INTPAIR;

typedef std::list<INTPAIR> flowpath;

flowpath* get_flowpath (int t);
#endif // FLOWPATH_H_
