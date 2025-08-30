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

function [y_st, C_recon] = RevealHHo2(SC, W)
N = size(W,1);
SC_vec = reshape(SC,N^2,1);
W_off = W-diag(diag(W));
R = reshape(W_off/2+W_off'/2,N^2,1);
% Eliminating Diagonal Elements
SC_index = (SC_vec~=0);
SC_vec = SC_vec(SC_index);
R1 = R(SC_index,:); 

% Homogeneous h
y_st = SC_vec\R1;
C_recon = W ./y_st;
end
