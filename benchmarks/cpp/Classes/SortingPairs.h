#ifndef SORTINGPAIRS_H
#define SORTINGPAIRS_H
#include <iostream>
using namespace std;

class GreaterOnTop
{
	public:
	bool operator()(pair<double,int> p1, pair<double,int> p2);
};

class LowestOnTop
{
	public:
	bool operator()(pair<double,int> p1, pair<double,int> p2);
};

#endif
