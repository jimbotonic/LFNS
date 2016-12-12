#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "../classes/Case.h"
#include "../classes/SparseMatrix.h"
#include "../classes/Array.h"
using namespace std;


//************************* FUNCTION DEFINITIONS **************************************
// Function to evaluate the RHS
Array rhs(const Array & y, const Array & P, const SparseMatrix & K, const SparseMatrix & G); 

// Runge-Kutta step
bool runge_kutta(Array & y, Array & ydot, double deltat, const Array & P, const SparseMatrix & B, const SparseMatrix & G, double precision);
//*************************************************************************************

int main()
{
	int node_num;//nombre de noeuds de la boucle
    double deltat;// pas temporel
    double precision;// precision
    string B_fn, G_fn, P_fn; //Enter the file names for B, G, and P
    string Theta_fn, Res_fn; //Enter the file names for Powers and Result file
    
    cin>> node_num >> deltat >> precision >> B_fn >> G_fn >> P_fn >> Res_fn ;

	vector< double > ydot(node_num,0); //vecteurs de conditions initiales
	vector< double > y(node_num,0); //vecteurs de conditions initiales

	Array Y(node_num,y);
	Array Ydot(node_num,ydot);

    SparseMatrix B(node_num,B_fn),G(node_num,G_fn);
    Array P(P_fn);
	
	auto start = chrono::steady_clock::now();
	
    ofstream Results;//Creates file to store the resuts, the file name contains the information about the connectivity on the first link
    Results.open(Res_fn);
	Results.precision(12);
	
	int step(0),real_step(0),step_max(1000000);//
	bool stop(0);
    while(step<step_max)
	{
        if(step%1000==0)
		{
			cout<<"iteration "<< step<<"  max(abs(y_dot))="<<Ydot.abs_max()<<endl;
		}

        //Results << step*deltat<<" "<< Y<<" "<< Ydot <<endl;
        stop=runge_kutta(Y,Ydot,deltat,P,B,G,precision);//Runge-Kutta step
        if (stop){
			real_step=step;;
			step=step_max;}
        else{step++;};

	}
	cout<<"REAL STEP "<<real_step<<endl;
	// Results << real_step*deltat<<" "<< Y << " " << Ydot <<endl;
	Results << Y << endl;
	Results.close();
	
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout <<"Duration="<< (chrono::duration <double, milli> (diff).count())/1000 << "seconds" << endl;
	return 0;
}


//************************* RUNGE-KUTTA STEP ****************************
//************************* RUNGE-KUTTA STEP ****************************
//************************* RUNGE-KUTTA STEP ****************************
bool runge_kutta(Array & y, Array & ydot, double deltat, const Array & P, const SparseMatrix & B, const SparseMatrix & G, double precision)
{
	int size(y.getSize());
	Array k1(size),k2(size),k3(size),k4(size);
	Array ydot_old(ydot);
	
	k1=deltat*rhs(y,P,B,G);
	k2=deltat*rhs(y+(1.0/2.0)*k1,P,B,G);
	k3=deltat*rhs(y+(1.0/2.0)*k2,P,B,G);
	k4=deltat*rhs(y+k3,P,B,G);
	
	y=y+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
	ydot=rhs(y,P,B,G);
	if ((ydot+(-1.0)*ydot_old).abs_max() < precision && ydot.var() < precision){ return true; }
	// if (ydot.abs_max() < precision){ return true; }
	else{ return false; }
}

//********************** EVALUATION RHS *********************************
//********************** EVALUATION RHS *********************************
//********************** EVALUATION RHS *********************************
Array rhs(const Array & y, const Array & P, const SparseMatrix & B, const SparseMatrix & G)
{
	int num_nodes(y.getSize());
    vector<double> cosy(num_nodes,0),siny(num_nodes,0),interaction(num_nodes,0),v(num_nodes,1.0);

    for(int j=0;j<num_nodes;j++)
    {
        cosy[j]=cos(y.getComposante(j));
        siny[j]=sin(y.getComposante(j));
    }

	Array SINY(num_nodes,siny),COSY(num_nodes,cosy),V(num_nodes,v);
	Array b_s(B*SINY),b_c(B*COSY),g_s(G*SINY),g_c(G*COSY);

	for(int j=0;j<num_nodes;j++)
    {
        interaction[j]=cosy[j]*(b_s.getComposante(j)-g_c.getComposante(j))-siny[j]*(g_s.getComposante(j)+b_c.getComposante(j));
    }
	Array INT_1(num_nodes,interaction),INT_2(G*V);
	
    return P+INT_1+INT_2;
}
