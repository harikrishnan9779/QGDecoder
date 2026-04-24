struct stab_to_graph_CSS
	{
	std::vector<std::string> stabs,Lz,Lx;
	unsigned int N,t,N_stab,N_log,Nx_Z,Nx_X,N_leftZ,N_leftX;
	arma::uvec stabs_orderZ,qubits_rev_orderZ,stabs_orderX,qubits_rev_orderX;
	arma::umat Xs,Zs,AdjZ,JZ,AdjX,JX;
	bool print_graph_status=false;
	stab_to_graph_CSS(std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,unsigned int);
	};

stab_to_graph_CSS::stab_to_graph_CSS(std::vector<std::string> stabilizer_gen,std::vector<std::string> logicalZ_op,std::vector<std::string> logicalX_op,unsigned int d)
	{
	stabs = stabilizer_gen;
	Lz = logicalZ_op;	
	Lx = logicalX_op;	
	N_stab=stabs.size(), N_log=Lz.size();
	N=N_stab+N_log;
	t = (d-1)/2;
	arma::umat A_SZ(2*N,N,arma::fill::zeros),A_SX(2*N,N,arma::fill::zeros);
	for(unsigned int j=0;j<N_stab;j++) //stabilizers
		{
		if(stabs[j].length()!=N)
			{throw std::runtime_error("All stabilizer generators must be same length N!");}
		for(unsigned int i=0;i<N;i++) //i th row qubit
			{
			if(stabs[j].at(i)!='I' and stabs[j].at(i)!='X' and stabs[j].at(i)!='Y' and stabs[j].at(i)!='Z')
				{throw std::runtime_error("Illegal entry in the stabilizer generator list. Only I,X,Y,Z are allowed!");}
			if(stabs[j].at(i)=='Z')
				A_SZ(i,j)=1;
			else if(stabs[j].at(i)=='X')
				{A_SZ(i+N,j)=1;}
			else if(stabs[j].at(i)=='Y')
				{throw std::runtime_error("For CSS code no Y stabilizers allowed!");}
			}
		}
	Zs=A_SZ(arma::span(0,N-1),arma::span(0,N_stab-1));	
	Xs=A_SZ(arma::span(N,2*N-1),arma::span(0,N_stab-1));	
	A_SX=A_SZ;
	for(unsigned int j=0;j<N_log;j++) //logical operators
		{
		if(Lz[j].length()!=N)
			{throw std::runtime_error("All logical operators must be length N!");}
		for(unsigned int i=0;i<N;i++) //i th row qubit
			{
			if(Lz[j].at(i)!='I' and Lz[j].at(i)!='Z')
				{throw std::runtime_error("Lz should only have I and Z!");}
			else if(Lz[j].at(i)=='Z')
				A_SZ(i,N_stab+j)=1;
			if(Lx[j].at(i)!='I' and Lx[j].at(i)!='X')
				{throw std::runtime_error("Lx should only have I and X!");}
			else if(Lx[j].at(i)=='X')
				A_SX(i+N,N_stab+j)=1;
			}
		}
	stabs_orderZ = arma::regspace<arma::uvec> (0,N-1); //Reordering to [0 | Z]
	stabs_orderX = arma::regspace<arma::uvec> (0,N-1); //              [X | 0] form
	unsigned int swap_index=0;             
	for(unsigned int j=0;j<N;j++)		
		{
		if(accu(A_SZ(arma::span(0,N-1),j))==0 and accu(A_SZ(arma::span(N,2*N-1),j))>0)
			{stabs_orderZ.swap_rows(swap_index,j);swap_index++;}
		else if(accu(A_SZ(arma::span(0,N-1),j))>0 and accu(A_SZ(arma::span(N,2*N-1),j))>0)
			{throw std::runtime_error("Y type stabilizer detected in CSS code A_SZ!");}
		}	
	swap_index=0;             
	for(unsigned int j=0;j<N;j++)		
		{
		if(accu(A_SX(arma::span(0,N-1),j))==0 and accu(A_SX(arma::span(N,2*N-1),j))>0)
			{stabs_orderX.swap_rows(swap_index,j);swap_index++;}
		else if(accu(A_SX(arma::span(0,N-1),j))>0 and accu(A_SX(arma::span(N,2*N-1),j))>0)
			{throw std::runtime_error("Y type stabilizer detected in CSS code A_SX!");}
		}				
	////////////////////////////////     GRAPH  FOR   Z   ERRORS        /////////////////////////////////
	arma::umat Z=A_SZ(arma::span(0,N-1),arma::span::all);
	arma::umat X=A_SZ(arma::span(N,2*N-1),arma::span::all);
	if(any(vectorise(mod2(arma::umat(Z.t()*X+X.t()*Z)))))
		{throw std::runtime_error("Stabilizer generators do not mutually commute! (Lz in list)");}	
	Nx_Z=0;
	for(unsigned int j=0;j<N;j++)		
		if(accu(X.col(j))>0)
			Nx_Z++;			
	arma::umat X1 = GF2_gauss(X);
	N_leftZ = GF2_rank(X1);	
	arma::uvec qubits_order=arma::regspace<arma::uvec>(0,N-1);
	invertible_block(X,N_leftZ,Nx_Z,qubits_order,stabs_orderZ); 
	X = X(qubits_order,stabs_orderZ);
	Z = Z(qubits_order,stabs_orderZ);
	if(N_leftZ==0)//X is 0 rank
		{throw std::runtime_error("X is zero rank! Need atleast one X stabilizer!");}
	else if(N_leftZ==N)//X is full rank
		{throw std::runtime_error("X is full rank! Need atleast one Z stabilizer!");}
	else if(accu(Z(arma::span::all,arma::span(0,N_leftZ-1)))>0 or accu(X(arma::span::all,arma::span(N_leftZ,N-1)))>0)	
		{throw std::runtime_error("Z1 or X2 is not zero! Not a CSS code!");}								
	qubits_rev_orderZ = arma::uvec(qubits_order.n_elem);
	qubits_rev_orderZ(qubits_order) = arma::regspace <arma::uvec> (0, qubits_order.n_elem - 1);			
	//Constructing the graph	
	arma::umat X1L = X(arma::span(0,N_leftZ-1),arma::span(0,N_leftZ-1));
	arma::umat X1R = X(arma::span(N_leftZ,N-1),arma::span(0,N_leftZ-1));
	arma::umat Z2L = Z(arma::span(0,N_leftZ-1),arma::span(N_leftZ,N-1));
	arma::umat Z2R = Z(arma::span(N_leftZ,N-1),arma::span(N_leftZ,N-1));
	arma::umat X1L_inv = GF2_inverse(X1L);
	arma::umat Z2R_inv = GF2_inverse(Z2R);
	arma::umat B = mod2(arma::umat(X1R*X1L_inv));
	AdjZ = B.t();	
	JZ = arma::umat(N,N,arma::fill::zeros);	
	JZ(arma::span(0,N_leftZ-1),arma::span(0,N_leftZ-1)) = X1L_inv; 	
	JZ(arma::span(N_leftZ,N-1),arma::span(N_leftZ,N-1)) = Z2R_inv;					
	arma::umat H=arma::eye<arma::umat> (2*N,2*N);	//Hadamard 
	for(unsigned int i=N_leftZ;i<N;i++)
			{H(i,i)=0;H(N+i,N+i)=0;H(i,N+i)=1;H(N+i,i)=1;}	
	arma::umat A_G(2*N,N,arma::fill::zeros);
	A_G(arma::span(0,N_leftZ-1),arma::span(N_leftZ,N-1))=AdjZ;
	A_G(arma::span(N_leftZ,N-1),arma::span(0,N_leftZ-1))=AdjZ.t();
	A_G(arma::span(N,2*N-1),arma::span::all)=arma::eye<arma::umat> (N,N);		
	Z = A_SZ(arma::span(0,N-1),arma::span::all); Z = Z(qubits_order,stabs_orderZ); A_SZ(arma::span(0,N-1),arma::span::all) = Z; 
	X = A_SZ(arma::span(N,2*N-1),arma::span::all); X = X(qubits_order,stabs_orderZ); A_SZ(arma::span(N,2*N-1),arma::span::all) = X;
	if(!all(vectorise(A_G==mod2(arma::umat(H*A_SZ*JZ)))))
		{throw std::runtime_error("A_S and A_G do not map each other!");}	
	if(print_graph_status)
		std::cout<<"Equivalent graph successfully generated for Z"<<std::endl;
	////////////////////////////////     GRAPH  FOR   X   ERRORS        /////////////////////////////////
	Z=A_SX(arma::span(0,N-1),arma::span::all);
	X=A_SX(arma::span(N,2*N-1),arma::span::all);
	if(any(vectorise(mod2(arma::umat(Z.t()*X+X.t()*Z)))))
		{throw std::runtime_error("Stabilizer generators do not mutually commute! (Lx in list)");}
	Nx_X=0;
	for(unsigned int j=0;j<N;j++)		
		if(accu(X.col(j))>0)
			Nx_X++;			
	X1 = GF2_gauss(X);
	N_leftX = GF2_rank(X1);	
	qubits_order=arma::regspace<arma::uvec>(0,N-1);
	invertible_block(X,N_leftX,Nx_X,qubits_order,stabs_orderX); 
	X = X(qubits_order,stabs_orderX);
	Z = Z(qubits_order,stabs_orderX);
	if(N_leftX==0)//X is 0 rank
		{throw std::runtime_error("X is zero rank! Need atleast one X stabilizer!");}
	else if(N_leftX==N)//X is full rank
		{throw std::runtime_error("X is full rank! Need atleast one Z stabilizer!");}
	else if(accu(Z(arma::span::all,arma::span(0,N_leftX-1)))>0 or accu(X(arma::span::all,arma::span(N_leftX,N-1)))>0)	
		{throw std::runtime_error("Z1 or X2 is not zero! Not a CSS code!");}						
	qubits_rev_orderX = arma::uvec(qubits_order.n_elem);
	qubits_rev_orderX(qubits_order) = arma::regspace <arma::uvec> (0, qubits_order.n_elem - 1);
	//Constructing the graph	
	X1L = X(arma::span(0,N_leftX-1),arma::span(0,N_leftX-1));
	X1R = X(arma::span(N_leftX,N-1),arma::span(0,N_leftX-1));
	Z2L = Z(arma::span(0,N_leftX-1),arma::span(N_leftX,N-1));
	Z2R = Z(arma::span(N_leftX,N-1),arma::span(N_leftX,N-1));
	X1L_inv = GF2_inverse(X1L);
	Z2R_inv = GF2_inverse(Z2R);
	B = mod2(arma::umat(X1R*X1L_inv));
	AdjX = B;	
	JX = arma::umat(N,N,arma::fill::zeros);	
	JX(arma::span(0,N_leftX-1),arma::span(0,N_leftX-1)) = X1L_inv; 	
	JX(arma::span(N_leftX,N-1),arma::span(N_leftX,N-1)) = Z2R_inv;			
	H=arma::eye<arma::umat> (2*N,2*N);	//Hadamard 
	for(unsigned int i=N_leftX;i<N;i++)
			{H(i,i)=0;H(N+i,N+i)=0;H(i,N+i)=1;H(N+i,i)=1;}	
	A_G = arma::umat(2*N,N,arma::fill::zeros);
	A_G(arma::span(0,N_leftX-1),arma::span(N_leftX,N-1))=AdjX.t();
	A_G(arma::span(N_leftX,N-1),arma::span(0,N_leftX-1))=AdjX;
	A_G(arma::span(N,2*N-1),arma::span::all)=arma::eye<arma::umat> (N,N);		
	Z = A_SX(arma::span(0,N-1),arma::span::all); Z = Z(qubits_order,stabs_orderX); A_SX(arma::span(0,N-1),arma::span::all) = Z; 
	X = A_SX(arma::span(N,2*N-1),arma::span::all); X = X(qubits_order,stabs_orderX); A_SX(arma::span(N,2*N-1),arma::span::all) = X;
	if(!all(vectorise(A_G==mod2(arma::umat(H*A_SX*JX)))))
		{throw std::runtime_error("A_S and A_G do not map each other!");}	
	if(print_graph_status)
		std::cout<<"Equivalent graph successfully generated for X"<<std::endl;
	return;
	}  


template <typename MatType>
void return_layers_CSS(arma::urowvec graph_syn,MatType &G,unsigned int t,std::vector<arma::uvec> &layers,unsigned int &n_layers,unsigned int &n_qubits)
	{
	arma::umat GT=G.t(); 
	unsigned int Nt=G.n_rows;
	unsigned int Nc=G.n_cols;
	n_layers=0;
	n_qubits=0;
	if(sum(graph_syn)==0)
		return;
	arma::uvec lt=find(graph_syn>0),t_qubits_covered,lc,c_qubits_covered;
	while(c_qubits_covered.n_elem<Nc)
		{
		lc={};
		for(unsigned int i:lt)
			{
			arma::uvec neighbors=nbr(G,i);
			for(unsigned int j:neighbors)
				{
				if(any(c_qubits_covered==j))
					continue;
				else
					c_qubits_covered.insert_rows(c_qubits_covered.n_elem,arma::urowvec({j}));
				if(any(lc==j))
					continue;
				else
					lc.insert_rows(lc.n_elem,arma::urowvec({j}));
				}	
			}			
		arma::uvec nei_density(lc.n_elem,arma::fill::zeros);
		for(unsigned int i=0;i<lc.n_elem;i++)
			{			
			arma::uvec neighbors=nbr(GT,lc(i));
			for(unsigned int j:neighbors)
				if(any(lt==j))
					nei_density(i)=nei_density(i)+1;
			}
		lc=lc(sort_index(nei_density,"descend"));
		if(lc.n_elem==0)	
			break;
		layers.push_back(lc);

		n_layers++;
		n_qubits+=lc.n_elem;
		if(n_layers==t)
			break;

		t_qubits_covered.insert_rows(t_qubits_covered.n_elem,lt);
		lt={};
		for(unsigned int i:lc)
			{
			arma::uvec neighbors=nbr(GT,i);
			for(unsigned int j:neighbors)
				{
				if(any(t_qubits_covered==j))
					continue;
				else
					lt.insert_rows(lt.n_elem,arma::urowvec({j}));
				}
			}				
		}	
	return;
	}
	
std::string decode_CSS_subroutine(const unsigned int t,const unsigned int max_wt,arma::urowvec graph_syn,arma::umat Adj,char error='X',bool dynamic_check=false)
	{
	unsigned int Nc=Adj.n_cols; //Control 
	unsigned int Nt=Adj.n_rows; //Target
	unsigned int N=Nc+Nt; 
	std::vector<arma::uvec> feed_forward;
	unsigned int N_layers,N_q_layers,max_layers;
	return_layers_CSS(graph_syn,Adj,max_wt,feed_forward,N_layers,N_q_layers);
	max_layers=min(arma::uvec({N_layers,max_wt}));

	//Decoding begins
	arma::urowvec mu_best(Nc,arma::fill::zeros); 
	arma::urowvec nu_best = arma::urowvec(graph_syn);
	unsigned int w_best = sum(graph_syn);
	if(w_best<=t)
		{goto correction_identified;}
	
	for(int cur_l=0;cur_l<max_layers;cur_l++)
		{
		arma::uvec curr_layer = feed_forward[cur_l];
		arma::uvec past_layer = {};
		for(int q=0;q<cur_l;q++)
			past_layer.insert_rows(past_layer.n_rows,feed_forward[q]);		
		int past_k_max = min(arma::uvec({past_layer.n_elem,max_wt}));
		for(int past_k=0;past_k<=past_k_max;past_k++) //Past layer all combination loop
			{
			int curr_k_max = min(arma::uvec({curr_layer.n_elem,max_wt-past_k}));
			for(int curr_k=0;curr_k<=curr_k_max;curr_k++) //Current layer all combination loop
				{
				if((curr_k+past_k)==0){continue;}
				std::vector<arma::uvec> past_perm_list=generate_excitations(past_layer.n_elem, past_k);
				std::vector<arma::uvec> curr_perm_list=generate_excitations(curr_layer.n_elem, curr_k);
				for(arma::uvec past_layer_ones:past_perm_list)
					{
					for(arma::uvec curr_layer_ones:curr_perm_list)
						{
						arma::urowvec mu_i(Nc,arma::fill::zeros);
						for(unsigned int q:past_layer_ones)
							mu_i(past_layer(q))=1;
						for(unsigned int q:curr_layer_ones)
							mu_i(curr_layer(q))=1;
						arma::urowvec nu_i = mod2(arma::urowvec(mu_i*Adj.t() + graph_syn));   //Adding to the first
						unsigned int w_i = sum(mu_i)+sum(nu_i);
						if(w_i<w_best)
							{
							mu_best = mu_i;
							nu_best = nu_i;
							w_best = w_i;
							}
						if(w_best<=t)
							{goto correction_identified;}
							
						if(dynamic_check) //Adding to the best
							{
							arma::urowvec mu_curr = mod2(arma::urowvec(mu_i+mu_best));
							arma::urowvec nu_curr = mod2(arma::urowvec(mu_i*Adj.t() + nu_best));
							unsigned int w_curr = sum(mu_curr)+sum(nu_curr);
							if(w_curr<w_best)
								{
								mu_best = mu_curr;
								nu_best = nu_curr;
								w_best = w_curr;
								}
							if(w_best<=t)
								{goto correction_identified;}
							}
						
						}
					}	
				}
			}				
		}
	correction_identified:
	std::string C(N,'I');
	if(error=='X')
		{
		arma::urowvec T=mu_best;
		T.insert_cols(T.n_elem,nu_best);
		for(int i=0;i<N;i++)
			if(T(i)==1)
				C[i]='X';
		}
	else if(error=='Z')
		{
		arma::urowvec T=nu_best;
		T.insert_cols(T.n_elem,mu_best);
		for(int i=0;i<N;i++)
			if(T(i)==1)
				C[i]='Z';
		}
	return(C);
	}

std::string decode_CSS(arma::urowvec Z_syn, arma::urowvec X_syn,const unsigned int T,stab_to_graph_CSS &S)
	{
	std::string C(S.N,'I');
	if(accu(Z_syn)>0 and accu(X_syn)==0) //Only Z errors
		{
		std::string CG = decode_CSS_subroutine(S.t,T,Z_syn,S.AdjZ,'Z');
		for(int i = 0; i < S.N; i++)
	    	C[i] = CG[S.qubits_rev_orderZ(i)];
		}
	else if(accu(Z_syn)==0 and accu(X_syn)>0) //Only X errors
		{
		std::string CG = decode_CSS_subroutine(S.t,T,X_syn,S.AdjX,'X');
		for(int i = 0; i < S.N; i++)
	    	C[i] = CG[S.qubits_rev_orderX(i)];
		}
	else if(accu(Z_syn)>0 and accu(X_syn)>0) //Both X,Z errors
		{
		std::string CG_Z = decode_CSS_subroutine(S.t,T,Z_syn,S.AdjZ,'Z');
		std::string CG_X = decode_CSS_subroutine(S.t,T,X_syn,S.AdjX,'X');

		std::string CZ(S.N,'I');
		for(int i = 0; i < S.N; i++)
	    	CZ[i] = CG_Z[S.qubits_rev_orderZ(i)];

		std::string CX(S.N,'I');
		for(int i = 0; i < S.N; i++)
	    	CX[i] = CG_X[S.qubits_rev_orderX(i)];

		for(int i=0;i<S.N;i++)
			{
			if(CZ[i]=='Z' and CX[i]=='I')
				C[i]='Z';
			else if(CZ[i]=='I' and CX[i]=='X')
				C[i]='X';
			else if(CZ[i]=='Z' and CX[i]=='X')
				C[i]='Y';
			}
		}
	return C;
	}
	
void get_stabilizer_syndrome_CSS(std::string &E,stab_to_graph_CSS &S,arma::urowvec &stab_syn)
	{
	Pauli_error Er(E);
	arma::urowvec stab_synZ = mod2(arma::urowvec(Er.nu*S.Xs));
	arma::urowvec stab_synX = mod2(arma::urowvec(Er.mu*S.Zs));
	stab_syn = mod2(arma::urowvec(stab_synZ+stab_synX));	
	return;
	}
	
void get_graph_syndrome_CSS(arma::urowvec stab_syn,stab_to_graph_CSS &S,arma::urowvec &graph_syn_Z,arma::urowvec &graph_syn_X)
	{
	stab_syn.insert_cols(stab_syn.n_elem,arma::urowvec(S.Lz.size(),arma::fill::zeros)); //Appending with trivial logical syndrome
	arma::urowvec stab_syn_Z_graph = stab_syn.cols(S.stabs_orderZ);
	arma::urowvec stab_syn_X_graph = stab_syn.cols(S.stabs_orderX);
	graph_syn_Z = mod2(arma::urowvec(stab_syn_Z_graph*S.JZ));
	graph_syn_X = mod2(arma::urowvec(stab_syn_X_graph*S.JX));
	graph_syn_Z = graph_syn_Z(arma::span(0,S.N_leftZ-1));	
	graph_syn_X = graph_syn_X(arma::span(S.N_leftX,S.N-1));	
	return;
	}			
