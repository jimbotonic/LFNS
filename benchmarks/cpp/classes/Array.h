//Classe pour manipuler des vecteurs 
#ifndef ARRAY_H
#define ARRAY_H
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iterator> 
#include <cmath>  
using namespace std;

class Array{

	public:
	Array();
	Array(const Array&);	//Constructuer
	Array(int taille); //Constructuer 
	Array(int taille, const vector<double>& VecteurInit); //Constructeur avec arguments
	Array(string fileName);// Constructeur a partir d'un fichier, determine la taille du vecteur a partir du nombre de lignes dans le fichier
	void showArray();
	double getComposante(int i) const {return array.at(i);};
	void setComposante(int i,double value) { array[i]=value; };
	int getSize() const {return size;};	
	double abs_max(); //returns the element of the array with the largest absolute value
	double max(); //returns the largest element of the array
	double var(); //returns the variance of the array
	double average(); //returns the average value of the array
	double total(); //returns the average value of the array
	
	//Surcharge d'operateurs
	Array operator+(const Array & b) const;
	void operator+=(const Array& b);
	Array operator*(double lambda) const;
    Array operator*(const Array & a) const;
	//Fonction amie
	friend Array operator*(double lambda, const Array& b);
	friend ostream& operator<<(ostream& os, const Array& ToPrint);

	private:

	int size;
	vector<double> array;

};

#endif
