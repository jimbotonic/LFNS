//Classe pour manipuler matrices sparses
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include <iostream>
#include "Case.h"
#include "Array.h"
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

class SparseMatrix{

	public:
	SparseMatrix(int taille, const vector<Case>& VecteurInit); //Constructeur avec arguments
	SparseMatrix(int taille, string fileName);// Constructeur a partir d'un fichier
	void showNormal(); // affiche comme vraie matrice
	void show(); // affiche an tant que sparse
	double getNumElements() const {return list.size();};	// renvoye le nombre d'elements non nuls dans la matrice sparse
	Case getCase(int i) const {return list[i];};	// renvoye la case i

//Surcharge d'operateurs
	Array operator*(const Array& v) const;
	//Fonction amie
//	friend Array operator*(double lambda, const Array& b);

	private:

	int size_matrice;
	vector<Case> list;

};

#endif
