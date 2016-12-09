#include "Array.h"

// Public Methods
Array::Array(){}

Array::Array(int taille) //Constructeur sans parametres: construis un array de zeros
{
	size=taille;
	for(int i(0); i<size;i++)
	{
		array.push_back(0);
	}
}

Array::Array(const Array& A)//Constructuer
{
	size=A.getSize();
	for(int i=0; i<size;i++)
	{
		array.push_back(A.getComposante(i));
	}
}	

Array::Array(int taille, const vector<double>& VecteurInit) //Constructeur avec parametres
{
	size=taille;
	for(int i=0; i<size;i++)
	{
		array.push_back(VecteurInit[i]);
	}
	
}

Array::Array(string fileName) //Constructeur depuis fichier
{
	int compteur(0);
	
	ifstream infile;
   	infile.open(fileName.c_str());

	double element;
	while(true)
	{
		infile >> element;
		if(infile.eof()) break;
		array.push_back(element);
		compteur++;
	}
	infile.close();
	size=compteur;
}

void Array::showArray() //Affiche le vecteur
{
	for(int i=0; i<size; i++)
	{
		cout<<array[i]<<" ";
	}
	cout<<endl;
}

Array Array::operator+(const Array& b) const // Surcharge operateur +
{
	vector<double> sum(size);

	for(int i=0; i<size; i++)
	{
		sum[i]=array[i]+b.getComposante(i);
	}
	Array Sum(size,sum);
	return Sum;
}

void Array::operator+=(const Array& b)// Surcharge operateur +=
{
	for(int i=0; i<size; i++)
	{
		array[i] += b.getComposante(i);
	}
}

Array Array::operator*(double lambda) const // Surcharge multiplication par un scalaire
{
	vector<double> product(size);
	for(int i=0; i<size; i++)
	{
		product[i]=lambda*array[i];
	}
	Array Product(size,product);
	return Product;
}


Array Array::operator*(const Array & a) const // Surcharge multiplication de deux array composante par composante
{
    vector<double> product(size);
    for(int i=0; i<size; i++)
    {
        product[i]=a.getComposante(i)*array[i];
    }
    Array Product(size,product);
    return Product;
}

//Fonction amie, permet d'implementer la commutativite de la multiplication
Array operator*(double lambda, const Array& b)
{
	return b*lambda;
}

//Surchage << 
ostream& operator<<(ostream& os, const Array& ToPrint)
{	
	for(int i=0; i<ToPrint.getSize(); i++)
	{
		os<<ToPrint.getComposante(i)<<" ";
	}
    return os;
}

//returns the element of the array with the largest absolute value
double Array::abs_max()
{
	double max_temp(0);
	double max (0);
	for(unsigned int i(0); i<array.size();i++)
	{
		max_temp=std::abs(array[i]);
		if(max_temp>max)
		{
			max=max_temp;
		}
	}
	return max;
}

//returns the largest element of the array
double Array::max()
{
	double max(array[0]),temp(0);
	for(unsigned int i(0); i<array.size();i++)
	{
		temp=array[i];
		if(temp>max)
		{
			max=temp;
		}
	}
	return max;
}

//returns the variance of the array
double Array::var()
{
	double first_moment(0);
	double second_moment(0);
	for(unsigned int i(0); i<array.size();i++)
	{
		first_moment+=array[i];
		second_moment+=pow(array[i],2);
	}
	second_moment=second_moment/array.size();
	first_moment=first_moment/array.size();
	
	return second_moment-pow(first_moment,2);
}

//returns the average of the array
double Array::average()
{
	double first_moment(0);
	for(unsigned int i(0); i<array.size();i++)
	{
		first_moment+=array[i];
	}
	first_moment=first_moment/array.size();
	
	return first_moment;
}

//returns the sum of the elements of the array
double Array::total()
{
	double sum(0);
	for(unsigned int i(0); i<array.size();i++)
	{
		sum+=array[i];
	}
	return sum;
}
