void triangular_color_code(unsigned int d,std::vector<std::string> &stabs,std::vector<std::string> &Lz,std::vector<std::string> &Lx)
	{
	stabs={};Lz={};Lx={};
	if(d%2==0)
		{throw std::runtime_error("Only odd d values are allowed!");}	
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
		{throw std::runtime_error("Only odd d values are allowed!");}
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
	
void cluster_CWS(unsigned int d,unsigned int &N,std::vector<std::string> &stabs,std::vector<std::string> &Lz,std::vector<std::string> &Lx,int l=-1)
	{
	if(d%2==0)
		{throw std::runtime_error("Only odd d values are allowed!");}
	arma::umat Adj;
	arma::uvec A,Ap;
	if(d==3)
		{
		l = 7;
		N = l; 
		A = arma::regspace<arma::uvec> (0,6);
		Adj=arma::umat(N,N,arma::fill::zeros);
		for(int k=0;k<l;k++)
			{
			Adj(k,mod_index(k+1,l))=1;
			Adj(k,mod_index(k-1,l))=1;
			}
		}	
	else if(d==5)
		{
		unsigned int l = 5;
		N = l*l; // =25 min |A|=21
		A = arma::regspace<arma::uvec> (0,20);
		Adj=arma::umat(N,N,arma::fill::zeros);
		for(int j=0;j<l;j++)
			{
			for(int k=0;k<l;k++)
				{
				Adj(j*l+k,j*l+mod_index(k-1,l))=1;
				Adj(j*l+k,j*l+mod_index(k+1,l))=1;
				Adj(j*l+k,k+l*mod_index(j+1,l))=1;
				Adj(j*l+k,k+l*mod_index(j-1,l))=1;
				}
			}						
		}	 
	else if(d==7)
		{
		unsigned int l = 4;
		N = l*l*l; // =64 min |A|=43
		A = arma::regspace<arma::uvec> (0,42);
		Adj=arma::umat(N,N,arma::fill::zeros);
		for(int i=0;i<l;i++)
			{
			for(int j=0;j<l;j++)
				{
				for(int k=0;k<l;k++)
					{
					Adj(i*l*l+j*l+k,i*l*l+j*l+mod_index(k-1,l))=1;
					Adj(i*l*l+j*l+k,i*l*l+j*l+mod_index(k+1,l))=1;
					Adj(i*l*l+j*l+k,i*l*l+mod_index(j-1,l)*l+k)=1;
					Adj(i*l*l+j*l+k,i*l*l+mod_index(j+1,l)*l+k)=1;
					Adj(i*l*l+j*l+k,mod_index(i-1,l)*l*l+j*l+k)=1;
					Adj(i*l*l+j*l+k,mod_index(i+1,l)*l*l+j*l+k)=1;
					}
				}
			}									
		}	 
	else if(d==9)
		{
		unsigned int l = 3;
		N=l*l*l*l; // =81 min |A|=73
		A = arma::regspace<arma::uvec> (0,72);
		Adj=arma::umat(N,N,arma::fill::zeros);
		for(int h=0;h<l;h++)
			{
			for(int i=0;i<l;i++)
				{
				for(int j=0;j<l;j++)
					{
					for(int k=0;k<l;k++)
						{
						Adj(h*l*l*l+i*l*l+j*l+k,h*l*l*l+i*l*l+j*l+mod_index(k-1,l))=1;
						Adj(h*l*l*l+i*l*l+j*l+k,h*l*l*l+i*l*l+j*l+mod_index(k+1,l))=1;
						Adj(h*l*l*l+i*l*l+j*l+k,h*l*l*l+i*l*l+mod_index(j-1,l)*l+k)=1;
						Adj(h*l*l*l+i*l*l+j*l+k,h*l*l*l+i*l*l+mod_index(j+1,l)*l+k)=1;
						Adj(h*l*l*l+i*l*l+j*l+k,h*l*l*l+mod_index(i-1,l)*l*l+j*l+k)=1;
						Adj(h*l*l*l+i*l*l+j*l+k,h*l*l*l+mod_index(i+1,l)*l*l+j*l+k)=1;
						Adj(h*l*l*l+i*l*l+j*l+k,mod_index(h-1,l)*l*l*l+i*l*l+j*l+k)=1;
						Adj(h*l*l*l+i*l*l+j*l+k,mod_index(h+1,l)*l*l*l+i*l*l+j*l+k)=1;
						}
					}
				}
			}										
		}
	arma::umat A_G(2*N,N,arma::fill::zeros),A_S(2*N,N,arma::fill::zeros);
	A_G(arma::span(0,N-1),arma::span::all)=Adj;
	A_G(arma::span(N,2*N-1),arma::span::all)=arma::eye<arma::umat> (N,N);	
	int count=0;
	for(unsigned int i=0;i<A.n_elem-1;i++)
		{A_S.col(i) = mod2(arma::uvec(A_G.col(A(i)) + A_G.col(A(i+1))));count++;}
	for(unsigned int i=0;i<N;i++)
		if(!any(A==i))
			{A_S.col(count) = A_G.col(i);count++;}
	std::string S(N,'I');
	for(unsigned int i:A)
		{A_S(i,N-1) = 1; S[i]='Z';}
	Lz.push_back(S);
	S=std::string(N,'I');
	for(unsigned int i=0;i<N;i++)
		if(i==A(0))
			S[i]='X';
		else if(Adj(i,A(0))==1)
			S[i]='Z';
	Lx.push_back(S);

	for(unsigned int j=0;j<N-1;j++)
		{
		S =std::string(N,'I');
		for(unsigned int i=0;i<N;i++)
			if(A_S(i,j)==0 and A_S(N+i,j)==1)
				S[i]='X';
			else if(A_S(i,j)==1 and A_S(N+i,j)==1)
				S[i]='Y';
			else if(A_S(i,j)==1 and A_S(N+i,j)==0)
				S[i]='Z';
		stabs.push_back(S);
		}
	return;
	}		
