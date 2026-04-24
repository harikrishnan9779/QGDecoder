template <typename T>	
void prnvec(const std::vector<T>& vecObj)
	{
	cout<<"\n\n"<<vecObj.size()<<" element(s)\n";
	for(int i=0; i<vecObj.size(); i++)
		{cout<<vecObj[i]<<"  ";}
	return;
	}	


template <typename MatType>
MatType mod2(MatType M)
	{
	if constexpr (arma::is_Col<MatType>::value || arma::is_Row<MatType>::value)
		{
		for(int i=0;i<M.n_elem;i++)
			{
			if(M(i)%2==0)
				M(i)=0;
			else
				M(i)=1;	
			}		
		}	
	else if constexpr (arma::is_Mat<MatType>::value)
		{
		for(int i=0;i<M.n_rows;i++)
			{
			for(int j=0;j<M.n_cols;j++)
				{
				if(M(i,j)%2==0)
					M(i,j)=0;
				else
					M(i,j)=1;	
				}	
			} 
		}
	return M;
	}	


std::vector<arma::uvec> generate_excitations(int N, int k)
	{
    std::vector<arma::uvec> result;
    std::vector<int> mask(N, 0);
    std::fill(mask.begin(), mask.begin() + k, 1);
    do 
    {
        arma::uvec v;
        for(unsigned int i = 0; i < N; i++)
        	{if (mask[i]) v.insert_rows(v.n_elem,arma::uvec({i}));}
        result.push_back(v);
    } while(std::prev_permutation(mask.begin(), mask.end()));
    return result;
	}
	
arma::uvec dec_to_base_n(int N, int b, int m) // Number N in base b bit length m 
	{
	arma::uvec c(m);
	for(int i=0;i<m;i++)
		{
		c(i) = N % b;
		N = N/ b;
		}
	return c;
	}		

unsigned int mod_index(int i,unsigned int m)
	{
	if(i>=0)
		return(i%m);
	else
		return(m-(abs(i)%m));	
	}

template <typename MatType>
arma::uvec nbr(MatType& G,unsigned int j)
	{
	arma::urowvec R=arma::urowvec(G.row(j));
	arma::uvec indices=arma::find(R>0);
	return(indices);
	}

		
struct Pauli_error
	{
	std::string E;
	arma::urowvec mu,nu;
	unsigned int N,w;
	Pauli_error(){};
	Pauli_error(std::string);
	Pauli_error(arma::urowvec,arma::urowvec);
	};

Pauli_error::Pauli_error(std::string e)
	{
	E=e;
	N=E.size();
	w=0;
	mu=arma::urowvec(N,arma::fill::zeros);
	nu=arma::urowvec(N,arma::fill::zeros);
	for(int i=0;i<N;i++)
		{
		if(E[i]=='X')
			{mu(i)=1;w++;}
		else if(E[i]=='Y')
			{mu(i)=1;nu(i)=1;w++;}
		else if(E[i]=='Z')
			{nu(i)=1;w++;}
		}			 
	return;
	}

Pauli_error::Pauli_error(arma::urowvec mu,arma::urowvec nu)
	{
	unsigned int N=mu.n_elem;
	if(nu.n_elem!=N){throw std::runtime_error("Pauli_error mu and nu must be of same length!");}
	w=0;
	E=std::string(N,'I');
	for(int i=0;i<N;i++)
		{
		if(mu(i)==1 and nu(i)==0)
			{E[i]='X';w++;}
		else if(mu(i)==1 and nu(i)==1)
			{E[i]='Y';w++;}
		else if(mu(i)==0 and nu(i)==1)
			{E[i]='Z';w++;}
		}			 
	return;
	}	


std::string get_pauli_error_vector(int N,int w,double p,char error_type='D',int seed_random=-1) //D - depolarization, X -	Bit fip, Z - Phase flip, Y - Bit-Phase flip 	
	{	
	if(error_type!='X' and error_type!='Y' and error_type!='Z' and error_type!='D')
		throw std::runtime_error("Not a valid error model choose X,Y,Z or D only");																				  			
	if(seed_random==-1)
		arma::arma_rng::set_seed_random();
	else
		arma::arma_rng::set_seed(seed_random);	
	arma::ivec Indices;
	while(Indices.n_elem<w)
		{
		int a=randi(arma::distr_param(0,N-1));
		if(any(Indices==a))
			continue;
		Indices.insert_rows(Indices.n_elem,arma::irowvec({a}));	
		}
	std::string E(N,'I');
	if(error_type=='D')
		for(int i:Indices)
			{
			double number = arma::randu();	
			if(number<(1.0-p))
				{continue;}
			else if((number>=(1.0-p)) and (number<(1.0-2.0*p/3.0)))	
				{E[i]='X';}
			else if((number>=(1.0-2.0*p/3.0)) and (number<(1.0-1.0*p/3.0)))	
				{E[i]='Y';}
			else
				{E[i]='Z';}
			}
	else
		for(int i:Indices)
			{
			double number = arma::randu();	
			if(number<(1.0-p))
				{continue;}
			else
				{E[i]=error_type;}
			}			
	return(E);
	}


arma::umat GF2_gauss(arma::umat output)
	{
	unsigned int m = output.n_rows, n = output.n_cols;
	if (m == 0 or n==0)
		{throw std::runtime_error("Empty matix in GF2_gauss");output=arma::umat(m,n,arma::fill::zeros);}

	int row=-1;
	for (int j=0;j<n;j++)
		{
		int pivot=-1;
		for (int i=row+1;i<m;i++)
			{
			if(output(i,j)==1)
				{
				if (pivot==-1)
					{
					pivot=i;
					row+=1;
					output.swap_rows(pivot, row);
					}
				else
					output.row(i)=mod2(arma::urowvec(output.row(i)+output.row(row)));
				}
			}
		}		                       
    return (output);
	}

unsigned int GF2_rank(arma::umat A)
	{
	unsigned int rank=0;
	for (int i=0;i<A.n_rows;i++)
		for (int j=0;j<A.n_cols;j++)
		    if (A(i,j)==1)
		        {rank+=1;break;}
	return (rank);
	}

bool GF2_is_invertible(const arma::umat& A)
	{
	arma::umat B = GF2_gauss(A);
	if(B.n_cols==GF2_rank(B))
		return(1);
	return(0);
	}
	
arma::umat GF2_inverse(arma::umat A)
	{
	unsigned int m = A.n_rows, n = A.n_cols;
	arma::umat output;
	if (m!=n)
		{throw std::runtime_error("Not a square matrix. Does not have inverse!");output=arma::uvec(m,arma::fill::zeros);}
	A.insert_cols(A.n_cols,arma::eye<arma::umat>(m,m));
	A = GF2_gauss(A);
	if (GF2_rank(A(arma::span::all,arma::span(0,m-1))) < m)
		{throw std::runtime_error("Singular matrix!");}
	A = reverse(A).t();
	arma::umat dummy = arma::reverse(A(arma::span(0,m-1),arma::span::all));
	dummy.insert_rows(m,arma::reverse(A(arma::span(m,2*m-1),arma::span::all)));
	A = dummy.t();
	output = GF2_gauss(A);
	output = output(arma::span::all,arma::span(m,2*m-1));
	output = reverse(output);output = reverse(output,1);
	return (output);
	}
	

void invertible_block(arma::umat &X,unsigned int n,unsigned int last_non_zer_col,arma::uvec& row_perm,arma::uvec& col_perm) 
	{
	unsigned int N = X.n_rows;
	unsigned int M = X.n_cols;
	if (n == M) //Already full rank
		{return;}
	if (n > N || n > M)
		throw std::runtime_error("n cannot exceed matrix dimensions");
	bool inv_block;
	int shift_ind_count=0;
	for (unsigned int k = 0; k < n; k++)
		{
		arma::uvec temp_col_perm = col_perm;
		arma::uvec temp_row_perm = row_perm;
		arma::umat temp_X = X(temp_row_perm,temp_col_perm);
		arma::umat submat = temp_X(arma::span(0,k),arma::span(0,k));
		inv_block = GF2_is_invertible(submat);
		while(!inv_block) //If submat not invertible
			{
			arma::uvec froz_ind_left = {};
			if(k>0)
				froz_ind_left = temp_col_perm(arma::span(0,k-1));
			arma::uvec shift_ind = temp_col_perm(arma::span(k,last_non_zer_col-1));
			arma::uvec froz_ind_right = {};
			if(last_non_zer_col<=(M-1))
				froz_ind_right = temp_col_perm(arma::span(last_non_zer_col,M-1));
			for (unsigned int i=k;i<N;i++)
				if(X(temp_row_perm(i),shift_ind(0)))
					{
					temp_row_perm.swap_rows(i,k);					
					temp_X = X(temp_row_perm,temp_col_perm); 
					submat = temp_X(arma::span(0,k),arma::span(0,k));
					if(GF2_is_invertible(submat))
						{inv_block=true;break;}
					else	
						temp_row_perm.swap_rows(i,k);
					} 
			if(!inv_block)
				{
				if(shift_ind.n_elem>1 and (shift_ind_count > shift_ind.n_elem)) //Shifting columns after k  to left by one if there is a reduction
					{
					shift_ind = shift(shift_ind, -1);
					temp_col_perm = arma::join_cols(froz_ind_left, arma::join_cols(shift_ind, froz_ind_right));
					}
				else  // Last chance try swapping last row with all remaining rows
					for (unsigned int i=k+1;i<N;i++) 
							{
							temp_row_perm.swap_rows(i,k);					
							temp_X = X(temp_row_perm,temp_col_perm); 
							submat = temp_X(arma::span(0,k),arma::span(0,k));
							if(GF2_is_invertible(submat))
								{inv_block=true;break;}
							else	
								temp_row_perm.swap_rows(i,k);
							}
				}
			shift_ind_count = shift_ind.n_elem;
			}
		col_perm = temp_col_perm;
		row_perm = temp_row_perm;			
		}	
	return;
	}


arma::uvec GF2_solve(const arma::umat& A, const arma::uvec& b) //Solve A x = b over GF(2)
	{
	unsigned int N = A.n_rows;
	unsigned int n = A.n_cols;
	if ((unsigned int)b.n_elem != N)
		throw std::runtime_error("Dimension mismatch in GF2_solve");
	arma::umat Aug(N, n+1, arma::fill::zeros);// Augmented matrix [A | b]
	Aug.cols(0, n-1) = A;
	Aug.col(n) = b;
	unsigned int row = 0;
	for (unsigned int col = 0; col < n && row < N; col++)// --- Forward elimination ---
		{
		unsigned int pivot = -1;
		for (unsigned int r = row; r < N; r++) // Find pivot
			if (Aug(r, col) == 1)
				{pivot = r;break;}
		if (pivot == -1)
			throw std::runtime_error("Matrix not full rank or inconsistent");
		if (pivot != row) // Swap into position
			{Aug.swap_rows(pivot, row);}
		for (unsigned int r = row + 1; r < N; r++)// Eliminate below
			if (Aug(r, col) == 1)
				Aug.row(r)=mod2(arma::urowvec(Aug.row(r)+Aug.row(row)));
		row++;
		}
	arma::uvec x = Aug(arma::span(0,n-1),n);
	for (int i = n-1; i >= 0; i--)// --- Back substitution ---
		{
		unsigned int pivot_row = -1;
		for (int r = N-1; r >= 0; r--)// Find pivot row for column i
			if (Aug(r, i) == 1)
				{
				pivot_row = r;
				break;
				}
		if (pivot_row == -1)
			throw std::runtime_error("Unexpected: missing pivot in GF2_solve");
		unsigned int val = Aug(pivot_row, n);
		for (unsigned int j = i+1; j < n; j++)// subtract known contributions (XOR)
			val = (val + (Aug(pivot_row,j) & x(j)))%2;
		x(i) = val;
		}
    return x;
	}


bool logical_error(std::string E,std::string C,std::vector<std::string> Lz,std::vector<std::string> Lx)
	{
	unsigned int N=E.size();
	std::string T(N,'I');
	for(int i=0;i<N;i++)
		{
		if(E[i]=='I'){T[i]=C[i];}
		else if(C[i]=='I'){T[i]=E[i];}
		else if((E[i]=='X' and C[i]=='Y') or (E[i]=='Y' and C[i]=='X'))
			{T[i]='Z';}
		else if((E[i]=='X' and C[i]=='Z') or (E[i]=='Z' and C[i]=='X'))
			{T[i]='Y';}
		else if((E[i]=='Y' and C[i]=='Z') or (E[i]=='Z' and C[i]=='Y'))
			{T[i]='X';}
		}
	arma::uvec countZ(Lz.size(),arma::fill::zeros);
	arma::uvec countX(Lz.size(),arma::fill::zeros);
	for(int j=0;j<Lz.size();j++)
		{
		for(int i=0;i<N;i++)
			{
			if(T[i]=='X' and (Lz[j][i]=='Y' or Lz[j][i]=='Z'))
				countZ(j)=countZ(j)+1;
			if(T[i]=='Y' and (Lz[j][i]=='X' or Lz[j][i]=='Z'))
				countZ(j)=countZ(j)+1;
			if(T[i]=='Z' and (Lz[j][i]=='X' or Lz[j][i]=='Y'))
				countZ(j)=countZ(j)+1;
			if(T[i]=='X' and (Lx[j][i]=='Y' or Lx[j][i]=='Z'))
				countX(j)=countX(j)+1;
			if(T[i]=='Y' and (Lx[j][i]=='X' or Lx[j][i]=='Z'))
				countX(j)=countX(j)+1;
			if(T[i]=='Z' and (Lx[j][i]=='X' or Lx[j][i]=='Y'))
				countX(j)=countX(j)+1;
			}		
		}		  
	countZ=mod2(countZ);
	countX=mod2(countX);
	if(any(countZ)==1 or any(countX)==1)
		return(1);
	return(0);
	}		
