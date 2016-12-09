#include "Case.h"

// Public Methods
Case::Case() //Constructeur sans parametres: construis une case contenant que des zeros
{
	Ligne=0;
	Colonne=0;
	Element=0.0;	
}

Case::Case(int ligne, int colonne, double element) //Constructeur sans parametres: construis une case contenant que des zeros
{
	Ligne=ligne;
	Colonne=colonne;
	Element=element;	
}

