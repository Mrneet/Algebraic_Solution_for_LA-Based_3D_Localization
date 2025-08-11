clear; close; clc;
rng('default');

u = [424;519;375];

scen = 'general';
switch scen
    case 'general'
        M = 12;
        aTmp = ((0:1/M:0.99))'*2*pi;
        bs = [300, 187,  -67, -270, -270,  -67,  187,  232, -170, 129, -106,   40;
              0,   235,  292,  130, -130, -292, -235, -150, -171,  55,  126, -101;
              94,  122, -112,  124,   39, -121,  -66,   49, -116, 233, -203,  107];
        alpha = aTmp;
        beta = [0; 0.813; 0.228; 2.460; 2.270; 0.438; 1.500; 1.232; 0.145; 0.305; 2.985; 0.108];
        nsepwr_dB = -40:5:15;
        plotLabel = 1:3;
        ylimMSE = [0.05,1e8];
        ylimBias = [0.1,1e4];
end

[N,~] = size(bs);
a = [cos(alpha).*cos(beta),sin(alpha).*cos(beta),sin(beta)]';

% true angles
for m = 1:M
    theta(m,1) = acos(a(:,m)'*(u-bs(:,m))/norm(u-bs(:,m)));
end

nsepwr_deg = 10.^(nsepwr_dB/10);
nsepwr = nsepwr_deg*(pi/180)^2;

NumEnsembles = 1000;
rng('default');
nseTmp = randn(M,NumEnsembles);
nse = nseTmp - mean(nseTmp,2);

rng('default');
atmp = rand(M,1)+1/9;
Ra = diag(roundn(atmp/mean(atmp),-2));

for n = 1:length(nsepwr)
    disp([num2str(n),'/',num2str(length(nsepwr))]);
    Q = nsepwr(n)*Ra;
    CRB = CRLB_AOALocLA(bs, u, a, Q);
    crlb(n) = trace(CRB);
   
    % WLS, CWLS, CWLS-C
    for m = 1:NumEnsembles
        thetaM = theta + sqrtm(Q)*nse(:,m);
        ia = 1;
        
        % WLS
        [uCF,u1] = LA3DLoc_WLS(thetaM,a,bs,Q);
        uEst(:,ia,m) = u1; ia = ia + 1;

        % CWLS
        [u2] = LA3DLoc_CWLS(thetaM,a,bs,Q);
        uEst(:,ia,m) = u2; ia = ia + 1;

        % CWLS-C
        [uCWLS] = LA3DLoc_CWLS_C(thetaM,a,bs,Q);
        uEst(:,ia,m) = uCWLS;

   end
    mse(:,n) = mean(sum((uEst-u).^2,1),3);
    bias(:,n) = mean(sqrt(sum((uEst-u).^2,1)),3);
end

symbs = ['x','^','o','s','*','v','+','p','d','h','>','<'];
names = {'WLS','CWLS','CWLS-C'};

figure;
for m = 1:ia-1
    semilogy(nsepwr_dB,(mse(m,:)),symbs(m),'Linewidth',1.5,'DisplayName',names{m}); hold on; grid on;
end
semilogy(nsepwr_dB,(mse(m+1,:)),symbs(m+1),'Linewidth',1.5,'LineStyle','--','DisplayName',names{m});
semilogy(nsepwr_dB,(crlb),'Linewidth',1.5,'DisplayName','CRLB');
legend('Show','Location','Northwest','fontsize',11,'NumColumns',2);
xlabel('10log(\sigma^2(rad^2))','fontsize',12); ylabel('MSE(u) (m^2)','fontsize',12);
xlim([-40,15])
ylim([0.01,1e6])

figure;
for i = 1:ia
    semilogy(nsepwr_dB,(bias(i,:)),symbs(i),'Linewidth',1.5,'DisplayName',names{i}); hold on; grid on;
end
legend('Show','Location','Northwest','fontsize',11,'NumColumns',2);
xlabel('10log(\sigma^2(rad^2))','fontsize',12); ylabel('Bias(u) (m)','fontsize',12);
xlim([-40,15])
ylim([0.1,1e3])