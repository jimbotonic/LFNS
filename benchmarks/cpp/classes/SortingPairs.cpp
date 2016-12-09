#include "SortingPairs.h"
using namespace std;

bool GreaterOnTop::operator()(pair<double,int> p1, pair<double,int> p2)
{
	if(p1.first<p2.first)
	{
		return true;
	}else{
		return false;
	}
}

bool LowestOnTop::operator()(pair<double,int> p1, pair<double,int> p2)
{
	if(p1.first>p2.first)
	{
		return true;
	}else{
		return false;
	}
}
