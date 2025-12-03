//STAN model with correlated random effects
//function to compute kronecker product
functions {
  matrix kronecker_prod(matrix A, matrix B) {
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
}

  real matMean(matrix A){
    real sums;
    sums = 0;
    for (i in 1:rows(A)){
      for (j in 1:cols(A)){
        sums = sums + A[i,j];
      }
    }
    return(sums/(cols(A)*rows(A)));
  }
}

data {
	int<lower=0> I;                  // number of subjects
	int<lower=0> J;                  // number of visits per subject
	int<lower=0> IJ;                 // total number of subjects
	int<lower=0> D;                  // grid length
	int<lower=0> p;                  // number of parameters

	int Kt;                          // number of spline basis functions

	matrix [IJ, D] Y;                // outcome matrix  
	matrix [I, p] X;                //fixed effects design matrix
  //matrix [IJ, I] Z;                //random effects design matrix
  array[IJ] int<lower=1, upper=I> subj_id;
	
	matrix[D,3*Kt] THETA;                 // B-spline evaluation matrix
	cov_matrix[Kt] PenMatInv;        //inverse of the penalty matrix
	//matrix[3,3] Ident2;             //2*2 identity matrix
	
	//matrix[IJ, D] Y_true;
	//matrix[p, D] beta_true;
}

//transpose theta bc we use it a lot
transformed data {
  matrix [3*Kt, D] THETA_t = THETA';
}

parameters {
  //matrices centered
	//matrix[3*Kt, I] B_Z;               //matrix of fixed effect spline coefficients
	//matrix[3*Kt, p] B_W;               //matrix of random effect spline coefficients
	
	//covariance for random effects
	cholesky_factor_corr[3] L_Omega_z;
	vector<lower=0>[3] tau_z;
	
	//covariance for fixed effects
	array[p] cholesky_factor_corr[3] L_Omega_W;
  array[p] vector<lower=0>[3] tau_W;
	
	//covariances
	//cov_matrix [3] Omega_z; //covariance matrix for the B_Z component
	//array holding the cov matrices for each 
	//array [p] cov_matrix [3] Omega_W;
	
	//uncentered spline coefficinets
	matrix[3*Kt, p] B_W_raw;
	matrix[3*Kt, I] B_Z_raw;
	
	//level 1 standard deviation
	real<lower=0> lev1_sigma;
}

transformed parameters {
  //definitions
  //covariance matrices
  matrix[3, 3] Omega_z;
  array[p] matrix[3, 3] Omega_W;
  
  matrix[3*Kt, 3*Kt] B_Z_cov_mat;
  matrix[3*Kt, 3*Kt] L_BZ;   //Cholesky 
  
  array[p] matrix[3*Kt, 3*Kt] B_W_cov_mat;
  array[p] matrix[3*Kt, 3*Kt] L_BW;   //cholesky
  
  //actual spline coefficients
  matrix[3*Kt, p] B_W;  
  matrix[3*Kt, I] B_Z;
  
  //beta matrix
  matrix[p, D] BETA;
  
  
  //Omega Z
  Omega_z = diag_pre_multiply(tau_z, L_Omega_z)
            * diag_pre_multiply(tau_z, L_Omega_z)';
            
  //omega W for each predictor
  for (k in 1:p) {
    Omega_W[k] = diag_pre_multiply(tau_W[k], L_Omega_W[k])
                 * diag_pre_multiply(tau_W[k], L_Omega_W[k])';
  }
  
  B_Z_cov_mat = kronecker_prod(Omega_z, PenMatInv);
  L_BZ = cholesky_decompose(B_Z_cov_mat);
  
  for(k in 1:p){
    B_W_cov_mat[k] = kronecker_prod(Omega_W[k], PenMatInv);
    L_BW[k] = cholesky_decompose(B_W_cov_mat[k]);
  }
  
  //noncentered portion
  for (k in 1:p) {
    B_W[, k] = L_BW[k] * B_W_raw[, k];
  }
  
  for (i in 1:I) {
    vector[3 * Kt] mu_i = B_W * X[i]';
    B_Z[, i] = mu_i + L_BZ * B_Z_raw[, i];
  }
  
  BETA = (THETA*B_W)';
  
  //matrices
  //cov_matrix[3*Kt] B_Z_cov_mat;
  //array[p] cov_matrix[3*Kt] B_W_cov_mat;
  
  //constructing the covariance structures
  //B_Z_cov_mat = kronecker_prod(Omega_z, PenMatInv);
  
  
  //beta matrix
  //matrix[p, D] BETA;
  //BETA = (THETA*B_W)';
}

model {
  //priors
  L_Omega_z ~ lkj_corr_cholesky(1.0);
  tau_z     ~ normal(0, 0.5);
  
  for (k in 1:p) {
    L_Omega_W[k] ~ lkj_corr_cholesky(1.0);
    //tau_W[k]     ~ normal(0, 1);
    tau_W[k] ~ normal(0, 0.5);
  }
  
  to_vector(B_W_raw) ~ normal(0, 1);
  to_vector(B_Z_raw) ~ normal(0, 1);
  
  //idk what prior for a sd?  half normal?
  //lev1_sigma ~ normal(0, 1);
  lev1_sigma ~ student_t(3, 0, 1);
  
  //likelihood
  for (n in 1:IJ) {
    int i = subj_id[n];                   
    row_vector[D] mu_n = (B_Z[, i]') * THETA_t;  
    Y[n] ~ normal(mu_n, lev1_sigma);
  }
}
