#include <armadillo>
using namespace std;

#include "QGDecoder/Auxillaries.h"
#include "QGDecoder/Standard_codes.h"
#include "QGDecoder/CSS.h"
			
int main()
{
unsigned int d=7; //code distance
std::vector<std::string> stabs,Lz,Lx;
triangular_color_code(d,stabs,Lz,Lx); 
//rotated_surface_code(d,stabs,Lz,Lx); 

stab_to_graph_CSS S(stabs,Lz,Lx,d);
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
	
	arma::urowvec stab_syn,graph_syn_Z,graph_syn_X;
	get_stabilizer_syndrome_CSS(E,S,stab_syn);
	get_graph_syndrome_CSS(stab_syn,S,graph_syn_Z,graph_syn_X);	

	// Bounded distance decoding

	string C = decode_CSS(graph_syn_Z,graph_syn_X,T,S);
	
	// Checking logical error		
		
	bool Log_err=logical_error(E,C,Lz,Lx);

	std::cout<<"Error:"<<E<<" Correction:"<<C<<" Weight:"<<Pauli_error(E).w<<" Logical error:"<<Log_err<<std::endl;
	}
return(0);
}
