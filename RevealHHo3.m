% It is built for linear reconstruction based on cov and partial cov.
% The idea of this function is to regress y, tau and C based on the
% knowledge of Symmetric SC, Estimated Jacobian Matrix and Stable Points.

% Firstly, we will build Matrix Mat_cont for regression, that is, to comput
% Matrix_cont * y = S_vec, where Matrix_cont, y and S_vec stand for the
% Symmetric Constraints, Fixed-points related Heterogeneity and Vectorized
% Symmetric SC.
% After getting y, we can reveal timescale Tau and 
% Asymmetric Structural Connectivity C.

% Input:
% SC: prior knowledge from Structural Connectivity (Symmetric);
% W: Reconstructed Jacobian Matrix(Have same size as S).

% Output:
% H: Reconstructed hierarchical heterogeneity,
% Tau_recon: Reconstructed timescales,
% C_recon: Reconstructed connectivity.

% Inter Parameters:
% S_vec: reshaped S, easier for estimation;
% W_trans: Transpose of W;
% Mat_cont: Temporary Matrix for defining R;
% R: Concatenated Matrix for revealing H;
% S_index: non-zero element indeices;
% y_st: 1/(gamma*G*J*(1-S_star)*h_i), where h_i is the slope of H_i at x_i;

function [y_st,C_recon] = RevealHHo3(SC, W)
N = size(W,1);
W_off = W-diag(diag(W));
y_st = zeros(N,1);
for i = 1:N
    y_st(i) = W_off(i,:)/SC(i,:);
end
C_recon = 1./y_st.*W_off;

% % Use lsqr for iterative solution
% maxit = 1000; % You can adjust this value
% tol = 1e-6; % You can adjust this value
% y_st = lsqr(SC_vec,R1, tol, maxit);

end
