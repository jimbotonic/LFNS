#include "SparseMatrix.h"

// Public Methods

SparseMatrix::SparseMatrix(int taille, const vector<Case>& VecteurInit) //Constructeur avec parametres
{
	size_matrice=taille;

	for(unsigned int i=0; i<VecteurInit.size();i++)
	{
		list.push_back(VecteurInit[i]);
	}
	
}


SparseMatrix::SparseMatrix(int taille, string fileName)
{	
	size_matrice=taille;
	
	ifstream infile;
   	infile.open(fileName.c_str());

	int ligne, colonne;
	double element;
	while(true)
	{
		infile >> ligne >> colonne >> element;
		if(infile.eof()) break;
		Case c(ligne,colonne,element);
		list.push_back(c);
	}
	infile.close();
}
//*******************************************************************************
//************************ Surcharge Operateurs *********************************
//*******************************************************************************

Array SparseMatrix::operator*(const Array& v) const //Multiplication d'une matrice sparse par un array
{
	vector<double> product(v.getSize(),0.0); //initializes vector of zeros
	for(unsigned int i=0;i<list.size();i++)
	{
		product[list[i].getLigne()]+=list[i].getElement()*v.getComposante(list[i].getColonne());
	}

	Array Product(v.getSize(),product);
	return Product;
}


//*******************************************************************************
//************************ Methodes d'affichage *********************************
//*******************************************************************************

void SparseMatrix::showNormal() //Affiche la matrice sparse comme une vraie matrice
{
	for(int i=1; i<=size_matrice; i++)
	{
		for(int j=1; j<=size_matrice; j++)
		{
			bool temp(true);
			for(unsigned int k=0; k<list.size(); k++)
			{
				if(list[k].getLigne()==i && list[k].getColonne()==j)
				{
					cout<<list[k].getElement()<<" "; 
					temp=false;
				}
			}
			if(temp)
			{
				cout<<"0 ";
			}

		}
		cout<<endl;
	}
}

void SparseMatrix::show() //Affiche la matrice sparse comme sparse
{
	for(unsigned int k=0; k<list.size(); k++)
	{
		cout<<"("<<list[k].getLigne()<<","<<list[k].getColonne()<<") "<<list[k].getElement()<<endl;
	}
}


