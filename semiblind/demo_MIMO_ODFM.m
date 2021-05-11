% clear all
% clc;
% close all;

%%
tic
N_t  = 2;          % number of TX
N_r  = 3;          % number of RX
M    = 5;          % number of multipaths
L    = 4;          % channel order (L <= M)
K       = 64;      % ODFM sub-carriers
N_pilot = 4;       % Data symbols
N_data  = 48;      % Data symbols
sigma_x = 5;       % Signal Power
W       = dftmtx(K);
F_L     = W(:,1:L);

%% Generate MIMO
[fading,delay,DOA] = gen_specular_para(M,N_t);

% AOA = pi * rand(M,N_r);
% 
% fading = rand(M,N_t);
% 
% delay  = rand(M,N_t);
% 
% DOA = rand(M,N_t);
% 
% AOA = rand(M,N_r);
[H, h_true]        = gen_chan_specular(fading,delay,DOA,N_r,L,N_t);

%% Derivative
h_derivative  = repmat(transpose(h_true),M*2,1);

% w.r.t. Fading
Hderiv_fading  = [];
M_fading  = cell(1,N_t);
for jj = 1 : N_t
    f_jj = 1./fading(:,jj);
    M_fading{1,jj} = repmat(f_jj,[1 L]);
end
D1            = blkdiag(M_fading{:});
D_fading      = repmat(D1, [1 N_r]);
Hderiv_fading = D_fading .* h_derivative;

% w.r.t. DOA
M_doa = cell(1,N_t);
for jj = 1 : N_t
    f_jj = -1i*pi*cos(DOA(:,1));
    M_doa{1,jj} = repmat(f_jj,[1 L]);
end
D2    = blkdiag(M_doa{:});
D_doa = repmat(D2,[1 N_r]);
RR    = repmat([0:N_r-1]',[1,L*N_t])';
rr    = RR(:)';
RR_doa = repmat(rr,M*2,1);
Hderiv_doa = RR_doa .* D_doa .* h_derivative;
% 
% % w.r.t. AOA
% M_aoa = cell(1,N_r);
% for ii = 1 : N_r
%     f_ii = -1i*pi*cos(AOA(:,1));
%     M_aoa{1,ii} = repmat(f_ii,[1 L]);
% end
% D4    = blkdiag(M_aoa{:});
% D_aoa = repmat(D4,[1 N_t]);
% RR2    = repmat([0:N_t-1]',[1,L*N_r])';
% rr2    = RR2(:)';
% RR_aoa = repmat(rr2,M*2,1);
% Hderiv_aoa = RR_aoa .* D_aoa .* h_derivative;

Hderiv_aoa = spec_chan_derive_AOA(fading,delay, DOA,AOA, N_r, L,M, N_t)
% w.r.t. delay
M_sinc = cell(1,N_t);
for jj = 1 : N_t
    f_jj = 1./sinc(delay(:,jj));
    dev  = [];
    for m = 1 : M
        dev_m = sin(delay(m,jj))/delay(m,jj)^2 - cos(delay(m,jj))/delay(m,jj);
        dev   = [dev; dev_m];
    end
    f_jj = f_jj .* dev;
    M_sinc{1,jj} = repmat(f_jj,[1 L]);
end
D3     = blkdiag(M_sinc{:});
D_sinc = repmat(D3,[1 N_r]);
Hderiv_sinc = D_sinc .* h_derivative;

G = [Hderiv_fading', Hderiv_doa', Hderiv_sinc' Hderiv_aoa];
P = G * G';

%% OFDM
% Channel Matrix
LAMBDA  = [];
for jj = 1 : N_t
    lambda_j =[];
    for r = 1 : N_r
        h_rj       = transpose(H(r,:,jj));
        lambda_rj  = diag(F_L*h_rj);
        lambda_j   = [lambda_j; lambda_rj];
    end
    LAMBDA = [LAMBDA lambda_j];
end
% Signal
% Use Zadoff-Chu sequence to generate signal
U = [1 3 5 7 9];
ZC_p = [];
for jj = 1 : N_t
    for k = 1 : K
        ZC(k,jj)  = 1 * exp( (-1i * pi * U(jj) * (k-1)^2 ) / K );
    end
end
x = []; % vector form
X = []; % Matrix form
for jj = 1 : N_t
    X_jj = 1/K * W * diag(ZC(:,jj)) * F_L ;
    X    = [X X_jj];
end
XXX = kron(eye(N_r),X);

%% CRB
% Partial derivative of LAMBDA w.r.t. h_i
partial_LAMBDA  = cell(1,N_r*N_t*L);
idx = 1;
while  idx <= N_r*N_t*L
    for ll = 1 : L
        partial_LAMBDA_ll = [];
        for jj = 1 : N_t
            lambda_jj =[];
            for r = 1 : N_r
                lambda_rj_ll = diag(F_L(:,ll));
                lambda_jj    = [lambda_jj; lambda_rj_ll];
            end
            partial_LAMBDA_ll = [partial_LAMBDA_ll lambda_jj];
        end
        partial_LAMBDA{1,idx} = partial_LAMBDA_ll;
        idx = idx + 1;
    end
end

% CRB Computation
SNR = -20:5:20;
CRB_OP      = zeros(length(SNR),1);
CRB_OP_Spec = zeros(length(SNR),1);
CRB_SB      = zeros(length(SNR),1);
CRB_SB_Spec = zeros(length(SNR),1);

for snr_i = 1 : length(SNR)
    sigma_v2 = 10^(-SNR(snr_i)/10);
    
    % For Data Symbol
    Cyy      = sigma_x * LAMBDA * LAMBDA'  + sigma_v2 * eye(N_r*K);
    Cyy_inv  = pinv(Cyy);
    
    I_D = zeros(N_r*N_t*L);
    for ii = 1 : N_r*N_t*L
        partial_Cyy_hii = sigma_x * LAMBDA * partial_LAMBDA{1,ii}';
        for jj = ii : N_r*N_t*L
            partial_Cyy_hjj = sigma_x * LAMBDA * partial_LAMBDA{1,jj}';
            % Slepian-Bangs Formula
            I_D(ii,jj) = trace(Cyy_inv * partial_Cyy_hii * Cyy_inv * partial_Cyy_hjj);
            I_D(ii,jj) = I_D(ii,jj)';
        end
    end
    
    % Pilot
    I_OP       = N_pilot * XXX'*XXX/sigma_v2;
    I_OP_Spec  = P * I_OP * P';
    I_SB       = N_data*I_D + I_OP;
    I_SB_Spec  = P * I_SB * P';
    
    
    CRB_OP_i           = pinv(I_OP);
    CRB_OP_Spec_i      = pinv(I_OP_Spec);
    CRB_SB_i           = pinv(I_SB);
    CRB_SB_Spec_i      = pinv(I_SB_Spec);
    
    CRB_OP(snr_i)      = abs(trace(CRB_OP_i));
    CRB_OP_Spec(snr_i) = abs(trace(CRB_OP_Spec_i));
    CRB_SB(snr_i)      = abs(trace(CRB_SB_i));
    CRB_SB_Spec(snr_i) = abs(trace(CRB_SB_Spec_i));
    
end

fig = figure;
semilogy(SNR,CRB_OP,'-r+');
hold on; semilogy(SNR,CRB_OP_Spec,'--r+');
hold on; semilogy(SNR,CRB_SB,'-bs');
hold on; semilogy(SNR,CRB_SB_Spec,'--bs');

 legend('CRB_{pilot usual}','CRB_{pilot specular}','CRB_{semiblind usual}','CRB_{semiblind specular}')
axis([-20 20 1e-6 1e10 ])
xlabel('SNR'); ylabel('trace(CRB)')
set(gca,'FontSize', 13,'FontName','Times New Roman');
grid on;


