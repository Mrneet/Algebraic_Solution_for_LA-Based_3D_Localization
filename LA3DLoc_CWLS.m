function [u] = LA3DLoc_CWLS(thetaM,a,bs,Q)

[~,M] = size(bs);

ss = diag(bs'*bs);
ct = cos(thetaM);
st = sin(thetaM);
as = diag(a'*bs);

b1 = -ct.^2.*ss + as.^2;
A1 = zeros(M, 9);
for i = 1:M
    A1(i,:) = [2*(as(i)*a(:,i) - ct(i)^2*bs(:,i))', ct(i).^2 - a(:,i)'.^2, ...
              -2*a(1,i)*a(2,i), -2*a(1,i)*a(3,i), -2*a(2,i)*a(3,i)];
end

C = eye(M);
E = eye(3);
E = padarray(E, [6, 6], 0, 'post');
P = [0; 0; 0; -1; -1; -1; 0; 0; 0];

for iter = 1:2
    W = inv(C*Q*C');
    % Calculate eigenvector matrix U and eigenvalue diagonal matrix V
    [U, V] = eig((A1'*W*A1)\E);

    c = (P'*U)';
    g = U\inv(A1'*W*A1)*P;
    e = (b1'*W*A1*U)';
    f = U\inv(A1'*W*A1)*A1'*W*b1;

    result = zeros(1, 19);
    % Two polynomials are multiplied to be zero, use the conv function
    for i = 1:9
        tmpa(1,i) = -0.25*V(i,i)*c(i,1)*g(i,1);
        tmpb(1,i) = 0.5*V(i,i)*c(i,1)*f(i,1) - 0.5*c(i,1)*g(i,1) - 0.5*V(i,i)*e(i,1)*g(i,1);
        tmpc(1,i) = c(i,1)*f(i,1) + V(i,i)*e(i,1)*f(i,1);
    end 
    for i=1:9
    % Calculate the coefficients of \(\prod_{j \neq i} (1 + {\gamma }_j \lambda)^{2}\)
    prod_coeffs = [1];  % Initialize as 1
        for j=1:9
            if j ~= i
                prod_coeffs = conv(prod_coeffs, [V(j,j)^2, 2*V(j,j), 1]);
            end
        end
    % Calculate the coefficients of \(\mathcal {A}_i \lambda^{2} + \mathcal {B}_i \lambda + \mathcal {C}_i\)
        poly_coeffs = [tmpa(1,i), tmpb(1,i), tmpc(1,i)];
    % Calculate the convolution of two polynomials
        tmpconv = conv(prod_coeffs, poly_coeffs);
    % Accumulate to the result vector
        result = result + tmpconv;
    end

    lambda = roots(result);

    % delete complex roots
    lambda = lambda(imag(lambda)==0);
    
    if ~isempty(lambda)
        for i = 1:length(lambda)
            psi(:,i) = (A1'*W*A1 + lambda(i,1)*E)\(A1'*W*b1 - 0.5*lambda(i,1)*P);
        end
    else
        lambda = zeros(8,1);
        for i = 1:length(lambda)
            psi(:,i) = (A1'*W*A1 + lambda(i,1)*E)\(A1'*W*b1 - 0.5*lambda(i,1)*P);
        end
    end
    
    objective_values = zeros(length(lambda),1);
    % Calculate the objective function value for each solution
    for il = 1:length(lambda)
        psi_i = psi(:, il);
        objective_values(il) = (b1 - A1*psi_i)' * W * (b1 - A1*psi_i);
    end
    % Find the minimum value and its corresponding index
    [min_value, min_index] = min(objective_values);

    % Optimal solution
    optimal_psi = psi(:, min_index);
    optimal_lamda = lambda(min_index, 1);
     
    u = optimal_psi(1:3);
    r_2 = sum((u-bs).^2,1)';
    C = 2*diag(r_2.*st.*ct);
end