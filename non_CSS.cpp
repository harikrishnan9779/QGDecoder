#include <armadillo>
using namespace std;

#include "QGDecoder/Auxillaries.h"
#include "QGDecoder/Standard_codes.h"
#include "QGDecoder/Additive.h"
			
int main()
{
unsigned int d=5; //code distance
std::vector<std::string> stabs,Lz,Lx;
//optimal_codes(d,stabs,Lz,Lx); //optimal codes of d=3,5,7,9,11 are available
XZZX_code(d,stabs,Lz,Lx); // XZZX code on d*d square lattice

stab_to_graph S(stabs,Lz,Lx,d);
unsigned int N=S.N; // No. of physical qubits
unsigned int T=S.t; // BDD target weight S.t = (d-1)/2

unsigned int w=N; // Maximum weight of Pauli error combination 
double p=0.1; // single qubit channel probability p 

for(int i=0;i<pow(10,3);i++)
	{

	// Generate N qubit, weight w Pauli error with single qubit channel probability p 
	// flags 'D' - Depolarization, 'X', 'Y' or 'Z' - Bit flip, Bit-phase flip and phase flip channels 

	std::string E=get_pauli_error_vector(N,w,p,'D');

	// Obtaining stabilizer syndrome and graph syndromes

	arma::urowvec stab_syn;
	arma::umat graph_syn;
	get_stabilizer_syndrome(E,S,stab_syn);
	get_graph_syndrome(stab_syn,S,graph_syn);

	// Bounded distance decoding

	std::string C = decode(graph_syn,T,S);	

	// Checking logical error		
	
	bool Log_err=logical_error(E,C,Lz,Lx);

	std::cout<<"Error:"<<E<<" Correction:"<<C<<" Weight:"<<Pauli_error(E).w<<" Logical error:"<<Log_err<<std::endl;
	}
return(0);
}
