struct stab_to_graph
	{
	std::vector<std::string> stabs,Lz,Lx;
	unsigned int N,t,N_stab,N_log,Nx,N_left;
	arma::uvec stabs_order,qubits_rev_order;
	arma::umat Xs,Zs,Adj,J;
	arma::uvec H,P;
	bool print_graph_status=false;
	stab_to_graph(std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,unsigned int,bool);
	};
		
stab_to_graph::stab_to_graph(std::vector<std::string> stabilizer_gen,std::vector<std::string> logicalZ_op,std::vector<std::string> logicalX_op,unsigned int d,bool use_Lz=true)
	{
	stabs = stabilizer_gen;
	Lz = logicalZ_op;	
	Lx = logicalX_op;	
	N_stab=stabs.size(), N_log=Lz.size();
	N = N_stab+N_log;
	t = (d-1)/2;
	arma::umat A_S(2*N,N,arma::fill::zeros);
	for(unsigned int j=0;j<N_stab;j++) //stabilizers
		{
		if(stabs[j].length()!=N)
			{throw std::runtime_error("All stabilizer generators must be same length N!");}
		for(unsigned int i=0;i<N;i++) //i th row qubit
			{
			if(stabs[j].at(i)!='I' and stabs[j].at(i)!='X' and stabs[j].at(i)!='Y' and stabs[j].at(i)!='Z')
				{throw std::runtime_error("Illegal entry in the stabilizer generator list. Only I,X,Y,Z are allowed!");}
			if(stabs[j].at(i)=='Z')
				A_S(i,j)=1;
			else if(stabs[j].at(i)=='X')
				A_S(i+N,j)=1;
			else if(stabs[j].at(i)=='Y')
				{A_S(i,j)=1;A_S(i+N,j)=1;}
			}
		}
	Zs=A_S(arma::span(0,N-1),arma::span(0,N_stab-1));	
	Xs=A_S(arma::span(N,2*N-1),arma::span(0,N_stab-1));	
	std::vector<std::string> log_op;
	if(use_Lz)
		log_op=Lz;
	else
		log_op=Lx;	
	for(unsigned int j=0;j<N_log;j++) //logical operators
		{
		if(log_op[j].length()!=N)
			{throw std::runtime_error("All logical operators must be length N!");}
		for(unsigned int i=0;i<N;i++) //i th row qubit
			{
			if(log_op[j].at(i)!='I' and log_op[j].at(i)!='X' and log_op[j].at(i)!='Y' and log_op[j].at(i)!='Z')
				{throw std::runtime_error("Illegal entry in the stabilizer generator list. Only I,X,Y,Z are allowed!");}
			if(log_op[j].at(i)=='Z')
				A_S(i,N_stab+j)=1;
			else if(log_op[j].at(i)=='X')
				A_S(i+N,N_stab+j)=1;
			else if(log_op[j].at(i)=='Y')
				{A_S(i,N_stab+j)=1;A_S(i+N,N_stab+j)=1;}
			}
		}
	stabs_order = arma::regspace<arma::uvec> (0,N-1); //Reordering to [Z1 | Z2]
	                                      //              [X1 | 0 ] form		
	Nx=0;             
	for(unsigned int j=0;j<N;j++)		
		if(accu(A_S(arma::span(N,2*N-1),j))>0)
			{stabs_order.swap_rows(Nx,j);Nx++;}
	arma::umat Z=A_S(arma::span(0,N-1),arma::span::all);
	arma::umat X=A_S(arma::span(N,2*N-1),arma::span::all);
	if(any(vectorise(mod2(arma::umat(Z.t()*X+X.t()*Z)))))
		{throw std::runtime_error("Stabilizer generators do not commute!");}
	else if(accu(X)==0)//X is 0 rank
		{throw std::runtime_error("X is zero rank! Need atleast one X stabilizer!");}
	else if(accu(Z)==0)//Z is full rank
		{throw std::runtime_error("Z is zero rank! Need atleast one Z stabilizer!");}	
	arma::umat X1 = GF2_gauss(X);
	N_left = GF2_rank(X1);	
	arma::uvec qubits_order=arma::regspace<arma::uvec>(0,N-1);
	invertible_block(X,N_left,Nx,qubits_order,stabs_order); 
	X = X(qubits_order,stabs_order);
	Z = Z(qubits_order,stabs_order);
	arma::umat J_ini = arma::eye<arma::umat> (N,N); //Recombinations to make X2=0
	if(N_left<N and accu(X(arma::span::all,arma::span(N_left,N-1)))>0)
		{
		arma::umat A = X(arma::span::all,arma::span(0,N_left-1));
		for(unsigned int i=N_left;i<N;i++)
			{
			arma::uvec b = X.col(i);
			arma::uvec x = GF2_solve(A,b);
			J_ini(arma::span(0,N_left-1),i) = x; 
			}
		}

	Z=A_S(arma::span(0,N-1),arma::span::all);
	X=A_S(arma::span(N,2*N-1),arma::span::all);
	X = X(qubits_order,stabs_order);
	Z = Z(qubits_order,stabs_order);
	X = mod2(arma::umat(X*J_ini)); 			
	Z = mod2(arma::umat(Z*J_ini)); 			
	qubits_rev_order = arma::uvec(qubits_order.n_elem);
	qubits_rev_order(qubits_order) = arma::regspace <arma::uvec> (0, qubits_order.n_elem - 1);
	////////////////////    Constructing the graph     ////////////////////////
	if(N_left==N)//X is full rank
		{
		arma::umat X_inv = GF2_inverse(X);
		arma::umat C = mod2(arma::umat(Z*X_inv));
		for (unsigned int i=0;i<N;i++)
			if(C(i,i)==1)
				{
				P.insert_rows(P.n_elem,arma::urowvec({i}));
				C(i,i)=0;
				}
		Adj = C;
		J = X_inv; 	
		arma::umat P_gate=arma::eye<arma::umat> (2*N,2*N); //Phase gate
		for(unsigned int i:P)
			{P_gate(i,N+i)=1;}	
		arma::umat A_G(2*N,N,arma::fill::zeros);
		A_G(arma::span(0,N-1),arma::span::all)=Adj;
		A_G(arma::span(N,2*N-1),arma::span::all)=arma::eye<arma::umat> (N,N);		
		Z = A_S(arma::span(0,N-1),arma::span::all); Z = Z(qubits_order,stabs_order); A_S(arma::span(0,N-1),arma::span::all) = Z; 
		X = A_S(arma::span(N,2*N-1),arma::span::all); X = X(qubits_order,stabs_order); A_S(arma::span(N,2*N-1),arma::span::all) = X;
		if(!all(vectorise(A_G==mod2(arma::umat(P_gate*A_S*J)))))
			{throw std::runtime_error("A_S and A_G do not map each other!");}
		if(print_graph_status)
			std::cout<<"Equivalent graph successfully generated"<<std::endl;
		return;
		}
	//X is non full rank
	arma::umat X1L = X(arma::span(0,N_left-1),arma::span(0,N_left-1));
	arma::umat X1R = X(arma::span(N_left,N-1),arma::span(0,N_left-1));
	arma::umat Z1L = Z(arma::span(0,N_left-1),arma::span(0,N_left-1));
	arma::umat Z1R = Z(arma::span(N_left,N-1),arma::span(0,N_left-1));
	arma::umat Z2L = Z(arma::span(0,N_left-1),arma::span(N_left,N-1));
	arma::umat Z2R = Z(arma::span(N_left,N-1),arma::span(N_left,N-1));
	arma::umat X1L_inv = GF2_inverse(X1L);
	arma::umat Z2R_inv = GF2_inverse(Z2R);
	arma::umat B = mod2(arma::umat(X1R * X1L_inv));
	arma::umat C = mod2(arma::umat((Z1L+X1L_inv.t()*X1R.t()*Z1R)*X1L_inv));
	for (unsigned int i=0;i<N;i++)
		{
		if(i<N_left and C(i,i)==1)
			{
			P.insert_rows(P.n_elem,arma::urowvec({i}));
			C(i,i)=0;
			}
		else if(i>=N_left)
			{
			H.insert_rows(H.n_elem,arma::urowvec({i}));
			}
		}
	Adj = arma::umat(N,N,arma::fill::zeros);	
	J = arma::umat(N,N,arma::fill::zeros);	
	Adj(arma::span(0,N_left-1),arma::span(0,N_left-1)) = C; 	
	Adj(arma::span(N_left,N-1),arma::span(0,N_left-1)) = B; 
	Adj(arma::span(0,N_left-1),arma::span(N_left,N-1)) = B.t();
	J(arma::span(0,N_left-1),arma::span(0,N_left-1)) = X1L_inv; 	
	J(arma::span(N_left,N-1),arma::span(0,N_left-1)) = mod2(arma::umat(Z2R_inv*Z1R*X1L_inv)); 
	J(arma::span(N_left,N-1),arma::span(N_left,N-1)) = Z2R_inv;		
	J = mod2(arma::umat(J_ini*J));
	
	arma::umat H_gate=arma::eye<arma::umat> (2*N,2*N); //Hadamard 
	arma::umat P_gate=arma::eye<arma::umat> (2*N,2*N); //Phase gate
	for(unsigned int i:H)
		{H_gate(i,i)=0;H_gate(N+i,N+i)=0;H_gate(i,N+i)=1;H_gate(N+i,i)=1;}	
	for(unsigned int i:P)
		{P_gate(i,N+i)=1;}	
	arma::umat A_G(2*N,N,arma::fill::zeros);
	A_G(arma::span(0,N-1),arma::span::all)=Adj;
	A_G(arma::span(N,2*N-1),arma::span::all)=arma::eye<arma::umat> (N,N);		
	Z = A_S(arma::span(0,N-1),arma::span::all); Z = Z(qubits_order,stabs_order); A_S(arma::span(0,N-1),arma::span::all) = Z; 
	X = A_S(arma::span(N,2*N-1),arma::span::all); X = X(qubits_order,stabs_order); A_S(arma::span(N,2*N-1),arma::span::all) = X;	
	if(!all(vectorise(A_G==mod2(arma::umat(P_gate*H_gate*A_S*J)))))
		{throw std::runtime_error("A_S and A_G do not map each other!");}
	if(print_graph_status)
		std::cout<<"Equivalent graph successfully generated"<<std::endl;
	return;
	}  


template <typename MatType>
void return_layers(arma::urowvec graph_syn,MatType &G,unsigned int t,std::vector<arma::uvec> &layers,unsigned int &n_layers,unsigned int &n_qubits) //returns feed forward graph upto depth t
	{
	unsigned int N=G.n_cols;
	if(sum(graph_syn)==0)
		{
		n_layers=0;
		n_qubits=0;
		return;
		}
	arma::uvec cur_l={},cur_l_nei={},qubits_covered={},nei_density={};
	for(unsigned int i=0;i<N;i++)
		if(graph_syn(i)==1)
			{
			cur_l.insert_rows(cur_l.n_elem,arma::urowvec({i}));
			qubits_covered.insert_rows(qubits_covered.n_elem,arma::urowvec({i}));
			}
	for(unsigned int i=0;i<cur_l.n_elem;i++)
		{
		arma::uvec neighbors = nbr(G,cur_l(i));
		for(unsigned int j:neighbors)
			if(!any(qubits_covered==j))
				{
				cur_l_nei.insert_rows(cur_l_nei.n_elem,arma::urowvec({j}));
				qubits_covered.insert_rows(qubits_covered.n_elem,arma::urowvec({j}));
				}
		}
	nei_density=arma::uvec(cur_l.n_elem+cur_l_nei.n_elem,arma::fill::zeros);
	for(unsigned int i=0;i<cur_l.n_elem;i++)
		{
		arma::uvec neighbors = nbr(G,cur_l(i));
		for(unsigned int j:neighbors)
			if(any(cur_l==j))
				nei_density(i)=nei_density(i)+1;
		}
	for(unsigned int i=0;i<cur_l_nei.n_elem;i++)
		{
		arma::uvec neighbors = nbr(G,cur_l_nei(i));
		for(unsigned int j:neighbors)
			if(any(cur_l==j))
				nei_density(cur_l.n_elem+i)=nei_density(cur_l.n_elem+i)+1;
		}	
	cur_l.insert_rows(cur_l.n_elem,cur_l_nei);			
	cur_l=cur_l(sort_index(nei_density,"descend"));
	layers.push_back(cur_l);
	n_layers=1;
	n_qubits=cur_l.n_elem;

	while(qubits_covered.n_elem<N)
		{
		cur_l = arma::uvec({});
		for(unsigned int i=0;i<cur_l_nei.n_elem;i++)
			{
			arma::uvec neighbors = nbr(G,cur_l_nei(i));
			for(unsigned int j:neighbors)
				if(!any(qubits_covered==j))
					{
					cur_l.insert_rows(cur_l.n_elem,arma::urowvec({j}));
					qubits_covered.insert_rows(qubits_covered.n_elem,arma::urowvec({j}));
					}
			}
		nei_density=arma::uvec(cur_l.n_elem,arma::fill::zeros);
		for(unsigned int i=0;i<cur_l.n_elem;i++)
			{
			arma::uvec neighbors = nbr(G,cur_l(i));
			for(unsigned int j:neighbors)
				if(any(cur_l_nei==j))
					nei_density(i)=nei_density(i)+1;
			}
		cur_l=cur_l(sort_index(nei_density,"descend"));
	
		cur_l_nei = arma::uvec({});
		for(unsigned int i=0;i<cur_l.n_elem;i++)
			{
			arma::uvec neighbors = nbr(G,cur_l(i));
			for(unsigned int j:neighbors)
				if(!any(qubits_covered==j))
					{
					cur_l_nei.insert_rows(cur_l_nei.n_elem,arma::urowvec({j}));
					qubits_covered.insert_rows(qubits_covered.n_elem,arma::urowvec({j}));
					}
			}
		nei_density=arma::uvec(cur_l_nei.n_elem,arma::fill::zeros);
		for(unsigned int i=0;i<cur_l_nei.n_elem;i++)
			{
			arma::uvec neighbors = nbr(G,cur_l_nei(i));
			for(unsigned int j:neighbors)
				if(any(cur_l==j))
					nei_density(i)=nei_density(i)+1;
			}
		cur_l_nei=cur_l_nei(sort_index(nei_density,"descend"));
		cur_l.insert_rows(cur_l.n_elem,cur_l_nei);			
		layers.push_back(cur_l);
		n_layers+=1;
		n_qubits+=cur_l.n_elem;
		if(n_layers==t)
			break;
		}	
	return;
	}

Pauli_error decode(unsigned int t,unsigned int T,arma::urowvec &graph_syn,arma::umat &Adj)
	{
	unsigned int N=graph_syn.n_cols; 
	arma::urowvec mu_best(N,arma::fill::zeros); 
	arma::urowvec nu_best = arma::urowvec(graph_syn);
	Pauli_error C_best(mu_best,nu_best);
	if(C_best.w<=t)
		{return C_best;}

	std::vector<arma::uvec> feed_forward;
	unsigned int N_layers,N_q_layers,max_layers;
	return_layers(graph_syn,Adj,T,feed_forward,N_layers,N_q_layers);
	max_layers=min(arma::uvec({N_layers,T}));
	
	//Weight reduction begins
	for(int cur_l=0;cur_l<max_layers;cur_l++)
		{
		arma::uvec curr_layer = feed_forward[cur_l];
		arma::uvec past_layer = {};
		for(int q=0;q<cur_l;q++)
			past_layer.insert_rows(past_layer.n_rows,feed_forward[q]);		
		int past_k_max = min(arma::uvec({past_layer.n_elem,T})); //past_layer.n_elem
		for(int past_k=0;past_k<=past_k_max;past_k++) //Past layer all combination loop
			{
			int curr_k_max = min(arma::uvec({curr_layer.n_elem,T-past_k})); //curr_layer.n_elem
			for(int curr_k=0;curr_k<=curr_k_max;curr_k++) //Current layer all combination loop
				{
				if((curr_k+past_k)==0){continue;}
				std::vector<arma::uvec> past_perm_list=generate_excitations(past_layer.n_elem, past_k);
				std::vector<arma::uvec> curr_perm_list=generate_excitations(curr_layer.n_elem, curr_k);
				for(arma::uvec past_layer_ones:past_perm_list)
					{
					for(arma::uvec curr_layer_ones:curr_perm_list)
						{
						arma::urowvec mu_i(N,arma::fill::zeros);
						for(unsigned int q:past_layer_ones)
							mu_i(past_layer(q))=1;
						for(unsigned int q:curr_layer_ones)
							mu_i(curr_layer(q))=1;
						arma::urowvec nu_i = mod2(arma::urowvec(mu_i*Adj + graph_syn));
						Pauli_error C_i(mu_i,nu_i);
						if(C_i.w<C_best.w)
							{C_best = C_i;}
						if(C_best.w<=t)
							{return C_best;}
						}
					}	
				}
			}				
		}
	return(C_best);
	}

std::string decode(arma::umat graph_syndromes,const unsigned int T,stab_to_graph &S)
	{
	Pauli_error C_list[graph_syndromes.n_rows]; //Correction for every graph syndrome 
	arma::uvec weights(graph_syndromes.n_rows);
	for(int i=0;i<graph_syndromes.n_rows;i++) //Going through each graph syndrome
		{
		arma::urowvec graph_syn=graph_syndromes.row(i);
		C_list[i]=decode(S.t,T,graph_syn,S.Adj);
		weights(i)=C_list[i].w;
		}
	std::string CG = C_list[weights.index_min()].E;
	for(int i=0;i<S.N;i++) // Applying H and P
		{
		if(any(S.H==i))
			{
			if(CG[i]=='Z'){CG[i]='X';}
			else if(CG[i]=='X'){CG[i]='Z';}
			}
		else if(any(S.P==i))
			{
			if(CG[i]=='X'){CG[i]='Y';}
			else if(CG[i]=='Y'){CG[i]='X';}
			}
		}	
	std::string C(S.N,'I');
	for(int i = 0; i < S.N; i++) //Reordering qubits
    	C[i] = CG[S.qubits_rev_order(i)];
	return(C);
	}

void get_stabilizer_syndrome(std::string &E,stab_to_graph &S,arma::urowvec &stab_syn)
	{
	Pauli_error Er(E);
	arma::urowvec stab_synZ = mod2(arma::urowvec(Er.nu*S.Xs));
	arma::urowvec stab_synX = mod2(arma::urowvec(Er.mu*S.Zs));
	stab_syn = mod2(arma::urowvec(stab_synZ+stab_synX));	
	return;
	}

void get_graph_syndrome(arma::urowvec &stab_syn,stab_to_graph &S,arma::umat &graph_syn)
	{
	arma::umat tot_syn(pow(2,S.N_log),stab_syn.n_elem);
	tot_syn.each_row() = stab_syn;
	tot_syn.insert_cols(tot_syn.n_cols,arma::umat(pow(2,S.N_log),S.N_log,arma::fill::zeros));
	for(int j=0;j<pow(2,S.N_log);j++)
		{
		arma::uvec lz_vals = dec_to_base_n(j, 2, S.N_log);
		for(int k=0;k<S.N_log;k++)
			tot_syn(j,S.N_stab+k)=lz_vals(k);
		}
	tot_syn = tot_syn.cols(S.stabs_order);
	graph_syn = mod2(arma::umat(tot_syn*S.J));
	return;
	}	
