function [uCF,u] = LA3DLoc_WLS(thetaM,a,bs,Q)

[~,M] = size(bs);

ss = diag(bs'*bs);
ct = cos(thetaM);
st = sin(thetaM);
as = diag(a'*bs);

b1 = -ct.^2.*ss + as.^2;
for i = 1:M
    A1(i,:) = [2*(as(i)*a(:,i) - ct(i)^2*bs(:,i))', ct(i).^2 - a(:,i)'.^2, ...
              -2*a(1,i)*a(2,i), -2*a(1,i)*a(3,i), -2*a(2,i)*a(3,i)];
end

C = eye(M);
for iter = 1:2
    W = inv(C*Q*C');
    psi = (A1'*W*A1)\A1'*W*b1;

    u = psi(1:3);
    r_2 = sum((u-bs).^2,1)';
    C = 2*diag(r_2.*st.*ct);
end

H = [eye(3); 2*diag(u); u(2), u(1), 0; u(3), 0, u(1); 0, u(3), u(2)];
psi2 = [u; u.*u; u(1)*u(2); u(1)*u(3); u(2)*u(3)];
b2 = b1 - A1*psi2;
A2 = -A1*H;
du = (A2'*W*A2)\A2'*W*b2;

uCF = u - du;