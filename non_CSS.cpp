#include <armadillo>
using namespace std;

#include "Auxillaries.h"
#include "Standard_codes.h"
#include "Additive.h"
			
int main()
{
std::vector<std::string> stabs,Lz,Lx;
stabs={"XZZXI","IXZZX","XIXZZ","ZXIXZ"};
Lz={"ZZZZZ"}; Lx={"XXXXX"};
unsigned int d=3;//[[5,1,3]]

unsigned int t=(d-1)/2;
unsigned int tmax=t;
double p=0.09;
stab_to_graph S(stabs,Lz,Lx);
for(int i=0;i<pow(10,3);i++)
	{
	std::string E=get_pauli_error_vector(S.N,S.N,p,'D');
	arma::urowvec stab_syn;
	arma::umat graph_syn;
	get_stabilizer_syndrome(E,S,stab_syn);
	get_graph_syndrome(stab_syn,S,graph_syn);
	std::string C = decode(graph_syn,t,tmax,S);	
	bool Log_err=logical_error(E,C,Lz,Lx);
	std::cout<<"Error:"<<E<<" Correction:"<<C<<" Weight:"<<Pauli_error(E).w<<" Logical error:"<<Log_err<<std::endl;
	}
return(0);
}
