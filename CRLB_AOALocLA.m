function [CRB] = CRLB_AOALocLA(bs, u, a, Q)
% The Cramer-Rao Bound of AOA localization by LA
% Input:
%   bs: sensor positions, NxM
%   u: source position, Nx1
%   a: direction vector of LA
% Output:
%   CRB: CRB matrix, NxN

[N,M] = size(bs);
r = sqrt(sum((u-bs).^2,1))';

for i = 1:M
    theta(i,1) = acos(a(:,i)'*(u-bs(:,i))/norm(u-bs(:,i)));
    A(i,:) = a(:,i)'*((u-bs(:,i))*(u-bs(:,i))'-r(i)^2*eye(N)) / (r(i)^3*abs(sin(theta(i))));
end

CRB = inv(A'/Q*A);
