function L = Neg_log_lik_eigen_with_nugget(param, p, input1, input2, R01, R02, output_mat,output)
    n1 = length(input1);
    n2 = length(input2);
    num_obs = n1*n2;
    %X = ones(1,num_obs);
    X = ones(num_obs,1);
    q_X = size(X,2);
    X_list = cell(1,q_X);
    for i = 1:q_X
        X_list{i} = reshape(X(:,i),[n1,n2]);
    end
    beta = exp(param(1:p));
    nu = exp(param(p+1));
    R1 = Matern_5_2(R01, beta(1));
    R2 = Matern_5_2(R02, beta(2));
    [eigen_R1_V,eigen_R1_D] = eig(R1);% diagonal matrix D of eigenvalues 
    %and matrix V whose columns are the corresponding right eigenvectors, 
    %so that A*V = V*D.
    [eigen_R2_V,eigen_R2_D] = eig(R2);
    U_x = [];
    for i = 1:q_X
        U_x = eigen_R1_V' * X_list{1} * eigen_R2_V;
        U_x = reshape(U_x,[1,num_obs]);
    end
    Lambda_tilde_inv = (1./(kron(diag(eigen_R2_D),diag(eigen_R1_D))+nu))';
    Lambda_tilde_inv_U_x = Lambda_tilde_inv.*U_x;
    X_R_tilde_inv_X_inv = inv(U_x * Lambda_tilde_inv_U_x');
   
    output_tilde = eigen_R1_V' * output_mat * eigen_R2_V;
    output_tilde = reshape(output_tilde,[1,num_obs]);
   
    theta_hat = X_R_tilde_inv_X_inv * (Lambda_tilde_inv_U_x * output_tilde');
   
    output_mat_normalized = output'-X.*theta_hat;
    output_mat_normalized = reshape(output_mat_normalized,[n1,n2]);
    
    output_normalize_tilde = eigen_R1_V' * output_mat_normalized * eigen_R2_V;
    output_normalize_tilde = reshape(output_normalize_tilde,[1,num_obs]);
    
    S_2 = sum(output_normalize_tilde .* Lambda_tilde_inv .* output_normalize_tilde);
    
    L = -(1/2*sum(log(Lambda_tilde_inv))-num_obs/2*log(S_2));
end