void triangular_color_code(unsigned int d,std::vector<std::string> &stabs,std::vector<std::string> &Lz,std::vector<std::string> &Lx)
	{
	stabs={};Lz={};Lx={};
	if(d%2==0)
		{std::cout<<"d must be odd!!"<<std::endl; exit(0);}
	unsigned int N = (3*d*d+1)/4;
	std::vector<arma::uvec> R,G,B,Plaquettes;
	arma::uvec qubits = arma::regspace<arma::uvec> (0,N-1),Lx_ind,Lz_ind;
	unsigned int j = 1;	
	while(qubits.n_elem>0)
		if(j%2==1)
			{
			arma::uvec dummy;
			for(unsigned int k=0;k<j;k++)
				dummy.insert_rows(dummy.n_elem,arma::urowvec({qubits(k)}));
			R.push_back(dummy);
			for(unsigned int k=0;k<j;k++)
				qubits.shed_row(0);
			if(qubits.n_elem==0)
				continue;
			dummy = arma::uvec({});
			for(unsigned int k=0;k<j;k++)
				dummy.insert_rows(dummy.n_elem,arma::urowvec({qubits(k)}));
			G.push_back(dummy);			
			for(unsigned int k=0;k<j;k++)
				qubits.shed_row(0);				
			j++;
			}	
		else
			{
			arma::uvec dummy;
			for(unsigned int k=0;k<j;k++)
				dummy.insert_rows(dummy.n_elem,arma::urowvec({qubits(k)}));
			B.push_back(dummy);	
			for(unsigned int k=0;k<j;k++)
				qubits.shed_row(0);				
			j++;
			}
	for(int i=0;i<R.size();i++)  // Lx
		Lx_ind.insert_rows(Lx_ind.n_elem,arma::urowvec({R[i](0)}));
	for(int i=0;i<G.size();i++)
		Lx_ind.insert_rows(Lx_ind.n_elem,arma::urowvec({G[i](0)}));
	for(int i=0;i<R.size();i++)  // Lz
		Lz_ind.insert_rows(Lz_ind.n_elem,arma::urowvec({R[i](R[i].n_elem-1)}));
	for(int i=0;i<B.size();i++)
		Lz_ind.insert_rows(Lz_ind.n_elem,arma::urowvec({B[i](B[i].n_elem-1)}));
	Lx_ind = sort(Lx_ind);
	Lz_ind = sort(Lz_ind);
	std::string S(N,'I');
	for(unsigned int i:Lz_ind){S[i]='Z';}
	Lz.push_back(S);
	S = std::string(N,'I');
	for(unsigned int i:Lx_ind){S[i]='X';}
	Lx.push_back(S);
	for(int i=0;i<R.size()-1;i++) // Red plaquettes
		{
		unsigned int s = R[i].size();
		j = 1;
		while(j<=s)
			{
			if(j==s)
				{
				arma::uvec dummy = {R[i](R[i].size()-1),G[i](G[i].size()-1),B[i](B[i].size()-2),B[i](B[i].size()-1)};
				Plaquettes.push_back(dummy);
				j++;
				}			
			else
				{
				arma::uvec dummy = {R[i](j-1),R[i](j),G[i](j-1),G[i](j),B[i](j-1),B[i](j)};
				Plaquettes.push_back(dummy);
				j+=2;
				}
			}
		}							
	for(int i=0;i<G.size();i++)  // Green plaquettes
		{
		unsigned int s = G[i].size();
		j = 1;
		while(j<=s)
			{
			if(j==1)
				{
				arma::uvec dummy = {G[i](0),B[i](0),R[i+1](0),R[i+1](1)};
				Plaquettes.push_back(dummy);
				j++;
				}			
			else
				{
				arma::uvec dummy = {G[i](j-1),G[i](j),B[i](j-1),B[i](j),R[i+1](j),R[i+1](j+1)};
				Plaquettes.push_back(dummy);
				j+=2;
				}
			}
		}
	for(int i=0;i<B.size();i++)  // Blue plaquettes
		{
		if(i==(B.size()-1))
			{
			j = 1;
			while(j<=B[i].size())
				{
				arma::uvec dummy = {B[i](j-1),B[i](j),R[i+1](j),R[i+1](j+1)};
				Plaquettes.push_back(dummy);
				j+=2;
				}
			}	
		else
			{
			j = 1;
			while(j<=B[i].size())
				{
				arma::uvec dummy = {B[i](j-1),B[i](j),R[i+1](j),R[i+1](j+1),G[i+1](j),G[i+1](j+1)};
				Plaquettes.push_back(dummy);
				j+=2;
				}
			}
		}			
	for(int i=0;i<Plaquettes.size();i++)
		{
		std::string Sz(N,'I'),Sx(N,'I');
		for(unsigned int j:Plaquettes[i])
			{
			Sz[j]='Z';
			Sx[j]='X';
			}
		stabs.push_back(Sz);
		stabs.push_back(Sx);
		}
	return;
	} 

void rotated_surface_code(unsigned int d,std::vector<std::string> &stabs,std::vector<std::string> &Lz,std::vector<std::string> &Lx)
	{
	stabs={};Lz={};Lx={};
	if(d%2==0)
		{std::cout<<"d must be odd!!"<<std::endl; exit(0);}
	unsigned int N = d*d;
	std::vector<arma::uvec> Z_Plaquettes,X_Plaquettes;
	arma::uvec Lx_ind,Lz_ind;
	for(unsigned int i=0;i<(d-1)/2;i++) // Top row 2 local Z
		{arma::uvec dummy={2*i+1,2*i+2};Z_Plaquettes.push_back(dummy);}
	for(unsigned int i=0;i<(d-1);i++)   // Bulk 4 local Z
		for(unsigned int j=0;j<(d-1)/2;j++)
			{arma::uvec dummy={i*d+2*j+i%2,i*d+2*j+1+i%2,(i+1)*d+2*j+i%2,(i+1)*d+2*j+1+i%2};Z_Plaquettes.push_back(dummy);}
	for(unsigned int i=0;i<(d-1)/2;i++) // Bottom row 2 local Z
		{arma::uvec dummy={d*(d-1)+2*i,d*(d-1)+2*i+1};Z_Plaquettes.push_back(dummy);}
	for(unsigned int i=0;i<(d-1)/2;i++) //Left column 2 local X
		{arma::uvec dummy={2*i*d,(2*i+1)*d};X_Plaquettes.push_back(dummy);}
	for(unsigned int i=0;i<(d-1);i++)   // Bulk 4 local X
		for(unsigned int j=0;j<(d-1)/2;j++)
			{arma::uvec dummy={i*d+2*j+(i+1)%2,i*d+2*j+1+(i+1)%2,(i+1)*d+2*j+(i+1)%2,(i+1)*d+2*j+1+(i+1)%2};X_Plaquettes.push_back(dummy);}
	for(unsigned int i=0;i<(d-1)/2;i++) // Right column 2 local X
		{arma::uvec dummy={2*(i+1)*d-1,2*(i+1)*d-1+d};X_Plaquettes.push_back(dummy);}
	for(unsigned int i=0;i<d;i++)       // Logical operators
		{Lz_ind.insert_rows(Lz_ind.n_elem,arma::urowvec({(i+1)*d-i-1}));
		Lx_ind.insert_rows(Lx_ind.n_elem,arma::urowvec({i*(d+1)}));}
	std::string S(N,'I');
	for(unsigned int i:Lz_ind){S[i]='Z';}
	Lz.push_back(S);
	S = std::string(N,'I');
	for(unsigned int i:Lx_ind){S[i]='X';}
	Lx.push_back(S);
	for(int i=0;i<X_Plaquettes.size();i++)
		{
		std::string S(N,'I');
		for(unsigned int j:X_Plaquettes[i])
			S[j]='X';
		stabs.push_back(S);
		}
	for(int i=0;i<Z_Plaquettes.size();i++)
		{
		std::string S(N,'I');
		for(unsigned int j:Z_Plaquettes[i])
			S[j]='Z';
		stabs.push_back(S);
		}		
	return;
	}

void XZZX_code(unsigned int d,std::vector<std::string> &stabs,std::vector<std::string> &Lz,std::vector<std::string> &Lx)
	{
	rotated_surface_code(d,stabs,Lz,Lx);
	for(int i=0;i<d;i++)
		for(int j=0;j<d;j++)
			if((i%2==0 and j%2==0) or (i%2==1 and j%2==1))
				{
				for(int s=0;s<stabs.size();s++)
					if(stabs[s][i*d+j]=='X')
						stabs[s][i*d+j]='Z';
					else if(stabs[s][i*d+j]=='Z')
						stabs[s][i*d+j]='X';
				for(int s=0;s<Lz.size();s++)
					if(Lz[s][i*d+j]=='Z')
						Lz[s][i*d+j]='X';
				for(int s=0;s<Lx.size();s++)
					if(Lx[s][i*d+j]=='X')
						Lx[s][i*d+j]='Z';
				}
	}

void optimal_codes(unsigned int d,std::vector<std::string> &stabs,std::vector<std::string> &Lz,std::vector<std::string> &Lx)
	{
	if(d==3) //[[5,1,3]]
		{
		stabs={"XZZXI","IXZZX","XIXZZ","ZXIXZ"};
		Lz={"YYIXI"};
		Lx={"IYIZZ"};
		}
	else if(d==5) //[[11,1,5]] https://www.codetables.de/
		{		
		stabs={"XIIIIXZZIXX","ZIIIIZXYYIX","IXIIIXZYZIY","IZIIIZZZYYI","IIXIIXZXXZI","IIZIIZYIYZZ","IIIXIXYIXXY","IIIZIZIZXZX","IIIIXXXIZZX","IIIIZZIYZYZ"};
		Lz={"IZIIXYZYIII"};
		Lx={"IXXIIIZIXIX"};
		}
	else if(d==7) //[[17,1,7]] https://www.codetables.de/	
		{
stabs={"XIIIIIIIXZYYIIYYZ","ZIIIIIIIZYXXIIXXY","IXIIIIIIZZYZYIXZI","IZIIIIIIYYXYXIZYI","IIXIIIIIIZZYZYIXZ","IIZIIIIIIYYXYXIZY","IIIXIIIIZYYYYZZXZ","IIIZIIIIYXXXXYYZY","IIIIXIIIZXZZYYYYZ","IIIIZIIIYZYYXXXXY","IIIIIXIIZXIYZYZZI","IIIIIZIIYZIXYXYYI","IIIIIIXIIZXIYZYZZ","IIIIIIZIIYZIXYXYY","IIIIIIIXZYYIIYYZX","IIIIIIIZYXXIIXXYZ"};
		Lz={"IXIIZYYXIIIZIIIIX"};
		Lx={"YIYIIIIXIIIYZZXII"};
		}
	else if(d==9) //[[25,1,9]] https://www.codetables.de/
		{
stabs={"XIIYYIIIIIIIZYXIIXZYIIYXZ","ZIIXXIIIIIIIYXZIIZYXIIXZY","IXIYZIIIIIIIYXZIIZYXIIXZY","IZIXYIIIIIIIXZYIIYXZIIZYX","IIXZYIIIIIIIYXZIIZYXIIXZY","IIZYXIIIIIIIXZYIIYXZIIZYX","IIIIIXIIYYIIYXZIIYXZIIYXZ","IIIIIZIIXXIIXZYIIXZYIIXZY","IIIIIIXIYZIIXZYIIXZYIIXZY","IIIIIIZIXYIIZYXIIZYXIIZYX","IIIIIIIXZYIIXZYIIXZYIIXZY","IIIIIIIZYXIIZYXIIZYXIIZYX","IIIIIIIIIIXIYZXIIIIIIIIII","IIIIIIIIIIZIXYZIIIIIIIIII","IIIIIIIIIIIXXXXIIIIIIIIII","IIIIIIIIIIIZZZZIIIIIIIIII","IIIIIIIIIIIIIIIXIYZXIIIII","IIIIIIIIIIIIIIIZIXYZIIIII","IIIIIIIIIIIIIIIIXXXXIIIII","IIIIIIIIIIIIIIIIZZZZIIIII","IIIIIIIIIIIIIIIIIIIIXIYZX","IIIIIIIIIIIIIIIIIIIIZIXYZ","IIIIIIIIIIIIIIIIIIIIIXXXX","IIIIIIIIIIIIIIIIIIIIIZZZZ"};
		Lz={"XIIYYIYZXIIIIIIXIZIZIIIII"};
		Lx={"IIIIIYIYYIIIIIIIYIZXIYXIZ"};
		d=9;
		}
	else if(d==11)	//[[29,1,11]] https://www.codetables.de/
		{
stabs={"XIIIIIIIIIIIIIXYZZXZZXXZZXZZY","ZIIIIIIIIIIIIIZXYYZYYZZYYZYYX","IXIIIIIIIIIIIIYYZYXIYXZIYXIYI","IZIIIIIIIIIIIIXXYXZIXZYIXZIXI","IIXIIIIIIIIIIIIYYZYXIYXZIYXIY","IIZIIIIIIIIIIIIXXYXZIXZYIXZIX","IIIXIIIIIIIIIIYZZZXZIYIIYYZIZ","IIIZIIIIIIIIIIXYYYZYIXIIXXYIY","IIIIXIIIIIIIIIZZXXIZXZXYYXIXX","IIIIZIIIIIIIIIYYZZIYZYZXXZIZZ","IIIIIXIIIIIIIIXXIYIZIIYYXZYZZ","IIIIIZIIIIIIIIZZIXIYIIXXZYXYY","IIIIIIXIIIIIIIZIZYXYXZZIIYXIY","IIIIIIZIIIIIIIYIYXZXZYYIIXZIX","IIIIIIIXIIIIIIYIXYIIZZXYXYZIZ","IIIIIIIZIIIIIIXIZXIIYYZXZXYIY","IIIIIIIIXIIIIIZZYZXYYIIZIYIXX","IIIIIIIIZIIIIIYYXYZXXIIYIXIZZ","IIIIIIIIIXIIIIXXIXYYXZXZIXXZZ","IIIIIIIIIZIIIIZZIZXXZYZYIZZYY","IIIIIIIIIIXIIIZIZYYIIYIZXZZZY","IIIIIIIIIIZIIIYIYXXIIXIYZYYYX","IIIIIIIIIIIXIIYIXYIZXYIXYZYYI","IIIIIIIIIIIZIIXIZXIYZXIZXYXXI","IIIIIIIIIIIIXIIYIXYIZXYIXYZYY","IIIIIIIIIIIIZIIXIZXIYZXIZXYXX","IIIIIIIIIIIIIXYZZXZZXXZZXZZYX","IIIIIIIIIIIIIZXYYZYYZZYYZYYXZ"};
		Lz={"IIIZIIIIIXIXYYZIYIIIIXIXIIXIY"};
		Lx={"IIIIZIIXYIIXIIYIIXIZIIIYIYXIZ"};
		}
	else
		{
		std::cout<<"Code unavailable!!"<<std::endl; exit(0);
		}
	}			
