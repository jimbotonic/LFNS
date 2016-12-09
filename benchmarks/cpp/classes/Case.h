//Classe pour qui construit les objects à l'interieur des Sparse array {int,int,double} (les deux premiers etant les positions dans le vecteur)
#ifndef CASE_H
#define CASE_H
using namespace std;

class Case{

	public:
	Case();	
	Case(int ligne, int colonne, double element); //Constructeur avec arguments

	int getLigne() const {return Ligne;};	
	int getColonne() const {return Colonne;};
	double getElement() const {return Element;};		

	private:

	int Ligne, Colonne;
	double Element;

};

#endif
