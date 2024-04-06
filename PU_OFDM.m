dmin = 2; drange = 2;
phimin = 30;    phirange = 30;
psimin = 30;    psirange = 30;
fbmin = 80e6;   fbrange = 40e6;
rho = 0.5; delta_t = 5e-9;
snr_dB = 0:1:30;
EbN0 = 10.^(snr_dB./10);
M = 64;
bitPerSymbol = log2(M); % 每符号比特数
N = 2048; % 子载波数
N_data_symbol = N/2-1;
Nsym = 100; % 每次仿真发的码元数
N_iteration = 1000*ones(size(snr_dB));
Es = 1; % 每个码元能量归一化为1
Eb = Es/bitPerSymbol;   % 每比特能量
% up = 1.6;    % 截断上限   
% pLM = 0.3;  % 线性映射的PAM信号幅度
% pELM = 0.4;
% qLM = 0.2;    % 极性余量
% qELM = 0.2;
% qNM = 0.25;
% L = 5;  % 非线性映射PAM信号阶数
up = 2.3;    % 截断上限   
pLM = 0.5;  % 线性映射的PAM信号幅度
pELM = 0.65;
qLM = 0.3;    % 极性余量
qELM = 0.35;
qNM = 0.3;
L = 4;  % 非线性映射PAM信号阶数
pn = (1:L-1)*up/L;    % 非线性映射的PAM信号幅度
deltapmax = up-pn-qNM;    % 每一阶非线性映射压缩的幅值
deltaP = cumsum(deltapmax); % 非线性映射压缩的总量
fs = 500e6; % 码元采样速率
Lmax = 10;    % LED拖尾长度
CP = Lmax;   % 循环前缀长度
leng = 5; width = 5;  height = 3;
Arx = 1e-4; Aroom = 2*(leng*width + leng*height + height*width);
t_bracket = 2/3e8* leng*width*height/Aroom;
tau = -t_bracket / log(rho);
f_0 = 1/(2*pi*tau);
eta_diff = Arx/Aroom *rho/(1-rho);
phi_half = 60;
m = -log(2)/log(cosd(phi_half));

d = rand(1,max(N_iteration))*drange + dmin;
phi = rand(1,max(N_iteration))*phirange + phimin;
psi = rand(1,max(N_iteration))*psirange + psimin;
fb = rand(1,max(N_iteration))*fbrange + fbmin;
ber_CE = zeros(size(snr_dB));
ber_U = zeros(size(snr_dB));    ber_LM = zeros(size(snr_dB));   ber_ELM = zeros(size(snr_dB));  ber_NM = zeros(size(snr_dB));
ber_theLM = zeros(size(EbN0));  ber_theNM = zeros(size(EbN0));  ber_theU = zeros(size(EbN0));
for i = 1:length(snr_dB)
    N0 = Eb./EbN0(i);   % 噪声平均功率    
    errorU = 0; errorLM = 0; errorELM = 0; errorNM = 0; errorCE = 0;
    ber_iterU = zeros(1,N_iteration(i)); ber_iterLM = ber_iterU;   ber_iterNM = ber_iterU;
    parfor k = 1:N_iteration(i)
        eta_LoS = (m+1)*Arx/(2*pi*d(k).^2) * cosd(phi(k)).^m * cosd(psi(k));
        f = (-N/2:N/2-1)*fs/N;
        H_diff = eta_diff * exp(-1j*2*pi*f*delta_t) ./ (1+1j*f/f_0);
        h_LED = exp(-2*pi*fb(k)/fs*(0:N-1));
        H_LED = fft(h_LED,N);
        H = ifftshift(eta_LoS + H_diff) .* H_LED;
        h_all = ifft(H);
        h = real(h_all(1:Lmax));
        h = h / sqrt(sum(h.^2));  % 归一化，保证功率不变
        H = fft(h,N);   % 频域特性，用于频域均衡
        fsC = 1.5*fs;   % CEO-OFDM需3个block，符号速率*1.5以保证信息速率相同
        h_LEDC = exp(-2*pi*fb(k)/fsC*(0:N-1));
        H_LEDC = fft(h_LEDC,N);
        HC = ifftshift(eta_LoS + H_diff) .* H_LEDC;
        h_allC = ifft(HC);
        LmaxC = round(1.5*Lmax);    % LED拖尾长度
        CPC = LmaxC;   % 循环前缀长度
        hC = real(h_allC(1:LmaxC));
        hC = hC./sqrt(sum(hC.^2));  % 归一化，保证功率不变
        HC = fft(hC,N);   % 频域特性，用于频域均衡
        % 产生QAM信号
        bits = round(rand(bitPerSymbol, N_data_symbol*Nsym));
        X_QAM = reshape(qammod(bits,M,'UnitAveragePower',true,'InputType','bit'), N_data_symbol,Nsym);
        % 赫米特对称
        X = [zeros(1,Nsym); X_QAM; zeros(1,Nsym); conj(flipud(X_QAM))];
        % IFFT
        x_ori = ifft(X,'symmetric')*sqrt(N);
        xp_u = zeros(N,Nsym);  xn_u = zeros(N,Nsym);
        xp_u(x_ori > 0) = x_ori(x_ori > 0);
        xn_u(x_ori < 0) = -x_ori(x_ori < 0);
        x_u = zeros(N,2*Nsym);
        x_u(:, 1:2:end) = xp_u;    x_u(:, 2:2:end) = xn_u;
        % PU调制
        xp_lm = xp_u;   xn_lm = xn_u;
        index_p1 = xp_lm > up; index_n1 = xn_lm > up; % 需进行压缩的index
        xp_lm(index_p1) = xp_lm(index_p1)-up+pLM+qLM;
        xn_lm(index_p1) = pLM;
        xn_lm(index_n1) = xn_lm(index_n1)-up+pLM+qLM;
        xp_lm(index_n1) = pLM;
        x_lm = zeros(N,2*Nsym);
        x_lm(:, 1:2:end) = xp_lm;    x_lm(:, 2:2:end) = xn_lm;
        
        xp_elm = xp_u;   xn_elm = xn_u;
        index_p1 = (xp_elm <= up-qELM)&(xp_elm>0); index_n1 = (xn_elm <= up-qELM)&(xn_elm>0); % 需时间分集的index
        index_p2 = xp_elm > up-qELM;  index_n2 = xn_elm > up-qELM; % 需进行压缩的index
        xp_elm(index_p2) = xp_elm(index_p2)-(up-qELM)+pELM+qELM;
        xn_elm(index_p1) = xp_elm(index_p1);
        xp_elm(index_p1) = xp_elm(index_p1)+qELM;
        xn_elm(index_n2) = xn_elm(index_n2)-(up-qELM)+pELM+qELM;
        xp_elm(index_n1) = xn_elm(index_n1);
        xn_elm(index_n1) = xn_elm(index_n1)+qELM;
        x_elm = zeros(N,2*Nsym);
        x_elm(:, 1:2:end) = xp_elm;    x_elm(:, 2:2:end) = xn_elm;
        
        xp_nm = xp_u;   xn_nm = xn_u;
        index_p1 = xp_nm > up; index_n1 = xn_nm > up; % 需进行压缩的index
        lnm = zeros(N,Nsym); xa_abs_nm = zeros(N,Nsym);
        for j = 1:L-1
            xa_abs_pos = xp_nm - deltaP(j);
            index_p2 = xa_abs_pos < up & lnm == 0 & index_p1;   % 第j阶压缩后不clipping的index
            lnm(index_p2) = j;
            xa_abs_nm(index_p2) = xa_abs_pos(index_p2);
            xa_abs_neg = xn_nm - deltaP(j);
            index_n2 = xa_abs_neg < up & lnm == 0 & index_n1;
            lnm(index_n2) = j;
            xa_abs_nm(index_n2) = xa_abs_neg(index_n2);
        end
        index_p3 = index_p1 & xa_abs_nm == 0;  index_n3 = index_n1 & xa_abs_nm == 0;    % 最后一阶压缩后仍clipping的index
        lnm(index_p3) = L-1;    xa_abs_nm(index_p3) = xa_abs_pos(index_p3);
        lnm(index_n3) = L-1;    xa_abs_nm(index_n3) = xa_abs_neg(index_n3);
        lnm(lnm == 0) = 1;  pnm = pn(lnm);  % 为确保lnm可作为索引，将零值置1
        xp_nm(index_p1) = xa_abs_nm(index_p1);
        xn_nm(index_p1) = pnm(index_p1);
        xn_nm(index_n1) = xa_abs_nm(index_n1);
        xp_nm(index_n1) = pnm(index_n1);
        x_nm = zeros(N,2*Nsym);
        x_nm(:, 1:2:end) = xp_nm;    x_nm(:, 2:2:end) = xn_nm;
        % 截断
        x_u(x_u > up) = up;
        x_lm(x_lm > up) = up;
        x_elm(x_elm > up) = up;
        x_nm(x_nm > up) = up;
        % 加循环前缀
        xu_cp = [x_u(N-CP+1:N, :); x_u];
        xlm_cp = [x_lm(N-CP+1:N, :); x_lm];
        xelm_cp = [x_elm(N-CP+1:N, :); x_elm];
        xnm_cp = [x_nm(N-CP+1:N, :); x_nm];
        % 信道
        xu_tx = conv2(xu_cp,h.'); xu_tx = xu_tx(1:end-Lmax+1,:);
        xu_rx = xu_tx + sqrt(N0)*randn(size(xu_tx));   % 加噪
        xlm_tx = conv2(xlm_cp,h.'); xlm_tx = xlm_tx(1:end-Lmax+1,:);
        xlm_rx = xlm_tx + sqrt(N0)*randn(size(xlm_tx));   % 加噪
        xelm_tx = conv2(xelm_cp,h.'); xelm_tx = xelm_tx(1:end-Lmax+1,:);
        xelm_rx = xelm_tx + sqrt(N0)*randn(size(xelm_tx));   % 加噪
        xnm_tx = conv2(xnm_cp,h.'); xnm_tx = xnm_tx(1:end-Lmax+1,:);
        xnm_rx = xnm_tx + sqrt(N0)*randn(size(xnm_tx));   % 加噪
        % 去循环前缀
        xu_rx = xu_rx(CP+1:end, :);
        xlm_rx = xlm_rx(CP+1:end, :);
        xelm_rx = xelm_rx(CP+1:end, :);
        xnm_rx = xnm_rx(CP+1:end, :);
        % FFT
        XU_rx = fft(xu_rx)/sqrt(N);
        XLM_rx = fft(xlm_rx)/sqrt(N);
        XELM_rx = fft(xelm_rx)/sqrt(N);
        XNM_rx = fft(xnm_rx)/sqrt(N);
        % 频域均衡
        XU_rx = XU_rx ./ H.';
        XLM_rx = XLM_rx ./ H.';
        XELM_rx = XELM_rx ./ H.';
        XNM_rx = XNM_rx ./ H.';
        x_eqU = ifft(XU_rx)*sqrt(N);
        x_eqLM = ifft(XLM_rx)*sqrt(N);
        x_eqELM = ifft(XELM_rx)*sqrt(N);
        x_eqNM = ifft(XNM_rx)*sqrt(N);
        
        x_pU = x_eqU(:, 1:2:end);   x_nU = x_eqU(:, 2:2:end);
        x_rxU = zeros(N,Nsym);
        x_rxU(x_pU >= x_nU) = x_pU(x_pU >= x_nU);
        x_rxU(x_pU < x_nU) = -x_nU(x_pU < x_nU);
        X_rxU = fft(x_rxU)/sqrt(N);
        
        x_pLM = x_eqLM(:, 1:2:end); x_nLM = x_eqLM(:, 2:2:end);
        index = x_nLM > x_pLM;  % 同时隙中负>正的索引
        x_aLM = zeros(N,Nsym); x_bLM = zeros(N,Nsym);
        x_aLM(~index) = x_pLM(~index);  x_aLM(index) = x_nLM(index);
        x_bLM(~index) = x_nLM(~index);  x_bLM(index) = x_pLM(index);
        x_rxLM = zeros(N,Nsym);
        index2 = x_bLM >= pLM/2;  % 存在压缩的索引
        x_rxLM(index2) = x_aLM(index2) +up-pLM-qLM;
        x_rxLM(~index2) = x_aLM(~index2);
        x_rxLM(index) = -x_rxLM(index);
        X_rxLM = fft(x_rxLM)/sqrt(N);
        
        x_pELM = x_eqELM(:, 1:2:end);   x_nELM = x_eqELM(:, 2:2:end);
        index = x_nELM > x_pELM;  % 同时隙中负>正的索引
        x_aELM = zeros(N,Nsym); x_bELM = zeros(N,Nsym);
        x_aELM(~index) = x_pELM(~index);  x_aELM(index) = x_nELM(index);
        x_bELM(~index) = x_nELM(~index);  x_bELM(index) = x_pELM(index);
        x_dELM = x_aELM -qELM -x_bELM;
        x_rxELM = zeros(N,Nsym);
        index2 = x_dELM >= pELM/2;  % 存在压缩的索引
        x_rxELM(index2) = x_aELM(index2) +up-pELM-2*qELM;
        x_rxELM(~index2) = (x_aELM(~index2) -qELM +x_bELM(~index2))/2;
        x_rxELM(index) = -x_rxELM(index);
        X_rxELM = fft(x_rxELM)/sqrt(N);
        
        x_pNM = x_eqNM(:, 1:2:end); x_nNM = x_eqNM(:, 2:2:end);
        index = x_nNM > x_pNM;  % 同时隙中负>正的索引
        x_aNM = zeros(N,Nsym); x_bNM = zeros(N,Nsym);
        x_aNM(~index) = x_pNM(~index);  x_aNM(index) = x_nNM(index);
        x_bNM(~index) = x_nNM(~index);  x_bNM(index) = x_pNM(index);
        lNM = round(x_bNM*L/up);    index2 = lNM > 0;   % 需解限幅的索引
        lNM(lNM > L-1) = L-1;   lNM(lNM <= 0) = 1;  deltapNM = deltaP(lNM); % 为确保lNM可作为索引，将大于L-1的限幅为L-1，非负的限幅为1
        x_rxNM = x_aNM;
        x_rxNM(index2) = x_aNM(index2) + deltapNM(index2);
        x_rxNM(index) = -x_rxNM(index);
        X_rxNM = fft(x_rxNM)/sqrt(N);
        % 解调
        X_deU = X_rxU(2:N_data_symbol+1, :);
        bits_U = reshape(qamdemod(X_deU,M,'UnitAveragePower',true,'OutputType','bit'), bitPerSymbol,N_data_symbol*Nsym);
        X_deLM = X_rxLM(2:N_data_symbol+1, :);
        bits_LM = reshape(qamdemod(X_deLM,M,'UnitAveragePower',true,'OutputType','bit'), bitPerSymbol,N_data_symbol*Nsym);
        X_deELM = X_rxELM(2:N_data_symbol+1, :);
        bits_ELM = reshape(qamdemod(X_deELM,M,'UnitAveragePower',true,'OutputType','bit'), bitPerSymbol,N_data_symbol*Nsym);
        X_deNM = X_rxNM(2:N_data_symbol+1, :);
        bits_NM = reshape(qamdemod(X_deNM,M,'UnitAveragePower',true,'OutputType','bit'), bitPerSymbol,N_data_symbol*Nsym);
        % CEO-OFDM
        x_ce_3 = zeros(N,Nsym);
        x_ce_3(xp_u > up) = xp_u(xp_u > up);    x_ce_3(xn_u > up) = xn_u(xn_u > up);
        x_ce_3 = x_ce_3 - up;   x_ce_3(x_ce_3 < 0) = 0;
        x_ce = zeros(N,3*Nsym);
        x_ce(:, 1:3:end) = xp_u;    x_ce(:, 2:3:end) = xn_u;    x_ce(:, 3:3:end) = x_ce_3;
        x_ce(x_ce > up) = up;
        xce_cp = [x_ce(N-CPC+1:N, :); x_ce];
        xce_tx = conv2(xce_cp,hC.');  xce_tx = xce_tx(1:end-LmaxC+1, :);
        xce_rx = xce_tx + sqrt(N0)*randn(size(xce_tx)); xce_rx = xce_rx(CPC+1:end, :);
        XCE_rx = fft(xce_rx)/sqrt(N);
        XCE_rx = XCE_rx ./ HC.';    x_eqCE = ifft(XCE_rx)*sqrt(N);
        x_pCE = x_eqCE(:, 1:3:end); x_nCE = x_eqCE(:, 2:3:end); x_3CE = x_eqCE(:, 3:3:end);
        index = x_nCE > x_pCE;
        x_rxCE = -x_nCE .* index + x_pCE .* ~index + (-1).^index .* x_3CE;
        X_rxCE = fft(x_rxCE)/sqrt(N);   
        X_deCE = X_rxCE(2:N_data_symbol+1, :);
        bits_CE = reshape(qamdemod(X_deCE,M,'UnitAveragePower',true,'OutputType','bit'), bitPerSymbol,N_data_symbol*Nsym);
        errorCE = errorCE + sum(bits_CE ~= bits, 'all');
        % 计算误码率
        errorU = errorU + sum(bits_U ~= bits, 'all');
        errorLM = errorLM + sum(bits_LM ~= bits, 'all');
        errorELM = errorELM + sum(bits_ELM ~= bits, 'all');
        errorNM = errorNM + sum(bits_NM ~= bits, 'all');
    end
    ber_U(i) = errorU / (bitPerSymbol*N_data_symbol*Nsym*N_iteration(i));
    ber_LM(i) = errorLM / (bitPerSymbol*N_data_symbol*Nsym*N_iteration(i));
    ber_ELM(i) = errorELM / (bitPerSymbol*N_data_symbol*Nsym*N_iteration(i));
    ber_NM(i) = errorNM / (bitPerSymbol*N_data_symbol*Nsym*N_iteration(i));
    ber_CE(i) = errorCE / (bitPerSymbol*N_data_symbol*Nsym*N_iteration(i));
end