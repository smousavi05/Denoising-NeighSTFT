 function [denoised] = neigBlock(data);
 % neighboring thresholding function 
 % Mostafa Mousavi 2015
 
%% STFT
% time_win = 1000;
time_win = 500;
factor_redund = 1;
stftCoef = STFT(data.x, time_win, factor_redund, 1/data.dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate noise level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha     = 0.92;                            % Decision-directed a priori SNR 
dB_xi_min = -25;                             % Minimum a priori SNR in dB
imcra.IS                 = 1 %10;            % Number of frames used to initialize noise
imcra.w                  = 1;                % How many adjacent bins (before and after bin k) are used to smooth spectrogram in frecuency (S_f(k,l))
imcra.alpha_s            = 0.9;              % Adaptation rate for the spectrogram smoothing in time
imcra.U                  = 8;                % U spectrogram values are stored to perform minimum tracking
imcra.V                  = 15;               % Each V frames minimum tracking will take place
imcra.Bmin               = 1.66;             % Bias of minimum noise estimate
imcra.Gamma0             = 4.6;              % A priori signal to noise ratio treshold for rough speech absence I (other conditions must also be met)
imcra.Gamma1             = 3;                % A priori signal to noise ratio treshold for speech absence q computation (other conditions must also be met)
imcra.zeta0              = 1.67;             % Threshold for "something similar to minimun estimated a priori signal to noise ratio but with spectrogram instead of signal [eq 18]"
imcra.alpha_d            = 0.85;             % Adaptation rate for speech probability dependent time smoothing paramater
imcra.beta               = 1.5 %1.47;       % Bias compensation when speech is absent for noise variance estimation
imcra.normalized_hanning = 1;                % Allows hanning window to be normalized
imcra.K = size(stftCoef, 1);                 % Get number of frequency bins
imcra.l = 0;                                 % Set counter for number of frames processed
imcra.j = 0;                                 % Set counter for the buffer update 
% To do the smoothing in frequency we use a matrix of indices on each frame
% and a corresponding matrix of row-stacked hanning windows
if imcra.normalized_hanning
    imcra.smth_win = repmat((hanning(2*imcra.w+1)/sum(hanning(2*imcra.w+1)))',imcra.K,1);
else
    imcra.smth_win = repmat(hanning(2*imcra.w+1)',K,1);
end
% Create coresponding matrix of indices
imcra.smth_mat = (repmat(1:imcra.K,2*imcra.w+1,1) + repmat(-imcra.w:imcra.w,imcra.K,1)')';
% Ignore indices out of bounds
idx              = find(imcra.smth_mat<=0 | imcra.smth_mat > imcra.K);
imcra.smth_mat(idx) = 1;
imcra.smth_win(idx) = 0;
% Smoothed spectrograms
imcra.S             = [];    % Smoothed Spectrogram first iteration
imcra.tilde_S       = imcra.S;  % Smoothed Spectrogram second iteration
imcra.Smin          = imcra.S;  % Smoothed Spectrogram minimum first iteration
imcra.tilde_Smin    = imcra.S;  % Smoothed Spectrogram minimum first second iteration
imcra.Smin_sw       = imcra.S;  % Smoothed Spectrogram minimum running minimum
imcra.tilde_Smin_sw = imcra.S;  % Second smoothed Spectrogram minimum running minimum
% Store buffers
imcra.Storing       = [];    % Smoothed Spectrogram minimum first iteration store buffer
imcra.tilde_Storing = [];    % Smoothed Spectrogram minimum second iteration store buffer
% Other parameters
imcra.ov_Lambda_D   = [];    % Biased noise variance estimate
imcra.Lambda_D      = [];    % Unbiased noise variance estimate 
imcra.p             = 1;     % A posteriori speech presence probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Noise Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y = stftCoef; 
    [K,L]  = size(Y);
    
    % 1) Compute Wiener posterior
    % This will hold the Wiener estimated clean signal
    hat_X_W        = zeros(size(Y));
    % This will hold the residual estimation uncertainty, in other words
    % the variance of the Wiener posterior
    Lambda         = zeros(size(Y));
    % Initialize noise power
    imcra.Lambda_D = abs(Y(:,1)).^2;
    % Initialize Gain and a posteriori SNR
    GH1            = ones(size(Y(:,1)));
    Gamma          = GH1;
    
    % Loop over frames
    for l=1:L 
        
        % SNR ESTIMATION (II)
        % A posteriori SNR [2, eq.3]
        new_Gamma = (abs(Y(:,l)).^2)./imcra.Lambda_D;                      
        % Decision-directed a priori SNR, with flooring [2, eq.32]
        xi        = alpha*(GH1.^2).*Gamma + (1-alpha)*max(new_Gamma-1,0);  
        xi        = max(xi,10^(dB_xi_min/20));                             
        % Update Gamma
        Gamma     = new_Gamma;

        % WIENER Posterior
        % Mean (Wiener filter)
        hat_X_W(:,l) = xi./(1+xi).*Y(:,l);
        % Variance (residual MSE)
        Lambda(:,l)  = xi./(1+xi).*imcra.Lambda_D;
        % Get the gain as well
        GH1          = xi./(1+xi);
        
        % SNR ESTIMATION (I), yes it is done in this order
        % IMCRA estimation of noise variance
        imcra = IMCRA(imcra,Y(:,l),Gamma,xi) 
    end
    
%     sigma = mean(imcra.Lambda_D);
%     stftCoef = stftCoef./sqrt(sigma);
    
sigma = (imcra.Lambda_D);
sigma = sqrt(sigma);
mean(sigma);

for q = 1:length(sigma)
    stftCoef(q,:) = (stftCoef(q,:))./sigma(q);
end

 
%% blocking 
% Maximum block length and width
[nR, nC] = size(stftCoef);

 W = 16;
 L = 8;

nb_L = floor(nC / L);
nb_W = floor(nR / W);

rr=0; cc=0;
wr = (nR - (nb_W.*W));
if wr > 0
rr = W - wr;
end 

lr = (nC - (nb_L.*L));
if lr > 0
cc = L - lr;
end

pad_coef = zeros(nR+rr, nC+cc);
pad_coef(1:nR,1:nC) = stftCoef;

[nR, nC] = size(pad_coef);
nb_L = floor(nC / L);
nb_W = floor(nR / W);

% nb_W = floor(((nR-1)/2) / W);
m = 0;
for i = 1: nb_L
    for j = 1: nb_W 
        m = m + 1;
        b1 = (j-1)*W; b2 = (i-1)*L; 
        H{m} = pad_coef(b1+1:b1+W,b2+1:b2+L);
    end
end 

%% Threshold coefficients 
for i = 1: length(H)

    [LB{i}, thres{i},T_H{i},a_m{i}] = Thresh(H{i});
end

%% Generating the thresholded coeff. matrix

% N = size(stftCoef,2);
t_decom = zeros(nR, nC);
a_map = zeros(nR, nC);
tlMap = [0,0];

m = 0;
for i = 1: nb_L
    for j = 1: nb_W 
        m = m + 1;
        b1 = (j-1)*W; b2 = (i-1)*L; 
        t_decom(b1+1:b1+W,b2+1:b2+L)= T_H{m};
    end
end 


dnCoef = t_decom(1:nR-rr,1:nC-cc);



denoised = inverseSTFT(dnCoef, time_win, factor_redund, 1/data.dt, length(data.x));
denoised =real(denoised);
% plot(denoised)


function [L, thr,t_subb a_map] = Thresh(subb)
[nRow,nCol] = size(subb); 

% Estimating optimal block size and threshold level 
[L,thr] = optimalTL(subb); 
L = round(L);

% Thresholding
 subb_real = real(subb);  
 ext_coef = padarray(subb_real,[L,L],0);
 t_subb = zeros(nRow,nCol); 
 a_map = zeros(nRow,nCol);
for i = L+1:L+nRow
    for j = L+1:L+nCol
        ux = ext_coef(i-L:i+L,j-L:j+L);
%         ux = ux./sig;
        S = sum(ux(:).^2);
        factor = max(1-thr^2/S,0);
%            a = (1 - thr*L*L*(sig)^2./S);
% %          a = (1 - thr*L*L./S); %%% CAI(1999)page 904
%          factor = a * (a > 0);
        t_subb(i-L,j-L) = factor*subb(i-L,j-L);
        a_map(i-L,j-L) = factor;
    end
end
    
 
function [LB,th] = optimalTL(subb)
% Estimatting the optimal threshold and block size using the SURE principle

% Input:
% subb: noisy coefficients

% Outputs;
% LB: Optimal block size
% th: Optimal threshold value


% Computting different block sizes based on the restirictions on Cai and Zhuo(2005) page 10 
% Lmax = ceil(size(subb,2)^0.75); 
% 
% p = floor(log2(Lmax));
% for i = 1:p
%     L_set(i) = 2^i;
% end

L_set = [3 5 7 9 11 13 15 17]; 
k = 0;
for L = L_set 
    k = k + 1;
    [risk(k),thr(k)] = SURE(subb,L);
end

% Obtain the optimal block size and the corresponding threshold
[guess,ibest] = min(risk);
th = thr(ibest); LB = L_set(ibest);



function [risk,th] = SURE(subb,L)
subb = real(subb);
[nRow,nCol] = size(subb);
Ns = nRow*nCol; % subband length


% if L == 1
% 
%     sb = subb(:).^2;
%     if min(sb) > 0 
%        Thres = [sb; 0];
%     end
%   else 

ext_sub_coef = padarray(subb,[L,L],0);

% Compute all S's and S2's of all neighbouring coefficients
S = zeros(Ns,1); m = 0;
for j = L+1:L+nCol
    for i = L+1:L+nRow            
        ux = ext_sub_coef(i-L:i+L,j-L:j+L);
        m = m + 1; S(m) = sum(ux(:).^2);      
    end
end
S2 = S.^2; tc2 = subb(:).^2;


Thres = L+1:0.1:(L+1)*3;
%% selectting the threshold coresponding to minimum SURE
risk = zeros(1,length(Thres)); m = 0;
for th = Thres
        
    th2 = th*th; temp1 = 0; temp2 = 0; temp3 = 0;
    for k = 1:Ns
        if S(k) > th2
           temp1 = temp1 + 1/S(k)-2*tc2(k)/S2(k);
           temp2 = temp2 + tc2(k)/S2(k);
       else
           temp3 = temp3 + tc2(k)-2;
       end
   end
    m = m + 1; risk(m) = Ns-2*th2*temp1+th2*th2*temp2+temp3;
    
end

[risk,ibest] = min(risk);
th = Thres(ibest);



% function cf = IMCRA(cf,Y_l,Gamma,xi)
%
% Implementation of Improved Minima Controlled Recursive Averaging(IMCRA) 
% after
%
%  [1] Israel Cohen, Noise Spectrum estimation in Adverse Environments:
%  Improved Minima Controlled Recursive Averaging. IEEE. Trans. Acoust.
%  Speech Signal Process. VOL. 11, NO. 5, Sep 2003.
%
% Input: cf     IMCRA configuration from previous iteration or obtained from
%               init_IMCRA.m
%
%        Y_l    Current frame of the STFT of size [K, 1]
%
%        Gamma  A posteriori signal to noise ratio    
%
%           xi  A priori signal to noise ratio 
%
% Output:   
%
%         cf.Lambda_D is the noise variance estimate 
%         cf.p        is the a posteriori speech probability
      
%
% Ramón F. Astudillo

function cf = IMCRA(cf,Y_l,Gamma,xi)

% Increase frame counter
cf.l = cf.l + 1;

% If in first frame, initialize buffers
if cf.l == 1    
    % Smoothed spectrograms
    cf.S             = sum(cf.smth_win.*abs(Y_l(cf.smth_mat)).^2,2); % Smoothed Spectrogram first iteration [3,eq.14]
    cf.tilde_S       = cf.S;                                         % Smoothed Spectrogram second iteration
    cf.Smin          = cf.S;                                         % Smoothed Spectrogram minimum first iteration
    cf.tilde_Smin    = cf.S;                                         % Smoothed Spectrogram minimum first second iteration
    cf.Smin_sw       = cf.S;                                         % Smoothed Spectrogram minimum running minimum
    cf.tilde_Smin_sw = cf.S;                                         % Second smoothed Spectrogram minimum running minimum
    % Store buffers
    cf.Storing       = [];                                           % Smoothed Spectrogram minimum first iteration store buffer
    cf.tilde_Storing = [];                                           % Smoothed Spectrogram minimum second iteration store buffer
    % Other parameters
    cf.ov_Lambda_D   = abs(Y_l).^2;                                  % Biased noise variance estimate
    cf.Lambda_D      = cf.ov_Lambda_D;                               % Unbiased noise variance estimate
    cf.p             = ones(size(Y_l));                              % A posteriori signal presence probability
end
    
% If in initialization segment, update noise stats only
% Note: Keep in mind that IS might be zero
if cf.l <= cf.IS
    % Compute minimum statistics smoothed spectrograms
    Sf         = sum(cf.smth_win.*abs(Y_l(cf.smth_mat)).^2,2);             % [3,eq.14]
    % Time smoothing
    cf.S       = cf.alpha_s*cf.S + (1-cf.alpha_s)*Sf;                      % [3,eq.15]
    % Update running minimum
    cf.Smin    = min(cf.Smin,cf.S);                                        % [3,eq.16]
    cf.Smin_sw = min(cf.Smin_sw,cf.S);
    % Compute smoothed spectrogram
    cf.Lambda_D = cf.alpha_d*cf.Lambda_D + (1-cf.alpha_d)*abs(Y_l).^2;     % [3,eq.8]
    % Set a posteriori speech probability to zero
    cf.p           = zeros(size(Y_l));
    
else
    
    % FIRST MINIMA CONTROLLED VAD
    % This provides a rough VAD to eliminate relatively strong signal
    % components towards the second power spectrum estimation
    Sf         = sum(cf.smth_win.*(abs(Y_l(cf.smth_mat)).^2),2);           % [3,eq.14]
    % Time smoothing
    cf.S       = cf.alpha_s*cf.S+(1-cf.alpha_s)*Sf;                        % [3,eq.15]
    % update running minimum
    cf.Smin    = min(cf.Smin,cf.S);                                        % [3,eq.16]
    cf.Smin_sw = min(cf.Smin_sw,cf.S);
    % Indicator function for VAD
    Gamma_min     = (abs(Y_l).^2)./(cf.Bmin*cf.Smin);                      % [3,eq.18]
    zeta          = cf.S./(cf.Bmin*cf.Smin);                               
    I             = zeros(size(Y_l));
    I((Gamma_min < cf.Gamma0 ) & (zeta < cf.zeta0)) = 1;                   % [3,eq.21]
    
    % SECOND MINIMA CONTROLLED VAD
    % This provides the signal probability needed to compute the final
    % noise estimation. The hard VAD index I, computed in the previous
    % estimation, is here used to exclude strong signal components.
    idx              = find(sum(I(cf.smth_mat),2) == 0);                   % [3,eq.26]
    warning off
    cf.tilde_Sf      = sum(cf.smth_win.*I(cf.smth_mat).*(abs(Y_l(cf.smth_mat)).^2),2)./sum(cf.smth_win.*I(cf.smth_mat),2);
    warning on
    cf.tilde_Sf(idx) = cf.tilde_S(idx);
    % Time smoothing
    cf.tilde_S       = cf.alpha_s*cf.tilde_S+(1-cf.alpha_s)*cf.tilde_Sf;   % [3,eq.27]
    % Update running minimum
    cf.tilde_Smin    = min(cf.tilde_Smin,cf.tilde_S);                      % [3,eq.26]
    cf.tilde_Smin_sw = min(cf.tilde_Smin_sw,cf.tilde_S);                   % [3,eq.27]
    
    % A PRIORI SIGNAL ABSENCE
    tilde_Gamma_min  = (abs(Y_l).^2)./(cf.Bmin*cf.tilde_Smin);             % [3,eq.28]
    tilde_zeta       = cf.S./(cf.Bmin*cf.tilde_Smin);                      
    
    % Signal absence
    q                = zeros(size(Y_l));
    idx_1            = find((tilde_Gamma_min <= 1) & (tilde_zeta < cf.zeta0));          % [3,eq.29]
    q(idx_1)         = ones(size(q(idx_1)));
    idx_2            = find((1 < tilde_Gamma_min) & (tilde_Gamma_min < cf.Gamma1) & (tilde_zeta < cf.zeta0));
    q(idx_2)         = (repmat(cf.Gamma1,size(q(idx_2)))-tilde_Gamma_min(idx_2))./(repmat(cf.Gamma1,size(q(idx_2)))-ones(size(q(idx_2))));
    
    % A POSTERIORI SIGNAL PROBABILITY
    nu               = Gamma.*xi./(1+xi);
    cf.p             = zeros(size(Y_l));
    cf.p(q < 1)      = (1+(q(q < 1)./(1-q(q < 1))).*(1+xi(q < 1)).*exp(-nu(q < 1))).^(-1); % [3,eq.7]
    
    % PROBABILITY DRIVEN RECURSIVE SMOOTHING
    % Smoothing parameter
    tilde_alpha_d    = cf.alpha_d+(1-cf.alpha_d)*cf.p;                                     % [3,eq.11]
    % UPDATE NOISE SPECTRUM ESTIMATE
    cf.ov_Lambda_D   = tilde_alpha_d.*cf.ov_Lambda_D + (1-tilde_alpha_d).*abs(Y_l).^2;     % [3,eq.10]
    % Bias correction
    cf.Lambda_D      = cf.beta*cf.ov_Lambda_D;                                             % [3,eq.12]
end

% UPDATE MINIMUM TRACKING
cf.j = cf.j+1;
if cf.j == cf.V
    % Minimum traking for the first estimation
    % store Smin_sw
    if size(cf.Storing,2) <= cf.U
        cf.Storing = [cf.Storing cf.Smin_sw];
    else
        cf.Storing = [cf.Storing(:,2:end) cf.Smin_sw];
    end
    % Set Smin to minimum
    cf.Smin = min(cf.Storing,[],2);
    % Let Smin_sw = S
    cf.Smin_sw = cf.S;
    % Minimum traking for the second estimation
    % store Smin_sw
    if size(cf.tilde_Storing,2) <= cf.U
        cf.tilde_Storing = [cf.tilde_Storing cf.tilde_Smin_sw];
    else
        cf.tilde_Storing = [cf.tilde_Storing(:,2:end) cf.tilde_Smin_sw];
    end
    % Set Smin to minimum
    cf.tilde_Smin    = min(cf.tilde_Storing,[],2);
    % Let Smin_sw = tilde_S
    cf.tilde_Smin_sw = cf.tilde_S;
    % reset counter
    cf.j=0;
end



function f_rec = inverseSTFT(STFTcoef, time_win, factor_redund, f_sampling, length_f)
%
% Inverse windowed Fourier transform. 
%
% Input:
% - STFTcoef: Spectrogram. Column: frequency axis from -pi to pi. Row: time
%   axis. (Output of STFT). 
% - time_win: window size in time (in millisecond).
% - factor_redund: logarithmic redundancy factor. The actual redundancy
%   factor is 2^factor_redund. When factor_redund=1, it is the minimum
%   twice redundancy. 
% - f_sampling: the signal sampling frequency in Hz.
% - length_f: length of the signal. 
%
% Output:
% - f_rec: reconstructed signal. 

% Guoshen Yu
% Version 1, Sept 15, 2006

% Window size
size_win = round(time_win/1000 * f_sampling);

% Odd size for MakeHanning
if mod(size_win, 2) == 0
    size_win = size_win + 1;
end
halfsize_win =  (size_win - 1) / 2;

Nb_win = floor(length_f / size_win * 2);

% Reconstruction
f_rec = zeros(1, length_f);

shift_k = round(halfsize_win / 2^(factor_redund-1));

% Loop over windows 
for k = 1 : 2^(factor_redund-1)
    for j = 1 : Nb_win - 1
        f_win_rec = ifft(STFTcoef(:, (k-1)+2^(factor_redund-1)*j));
        f_rec(shift_k*(k-1)+(j-1)*halfsize_win+1 : shift_k*(k-1)+(j-1)*halfsize_win+size_win) =  f_rec(shift_k*(k-1)+(j-1)*halfsize_win+1 : shift_k*(k-1)+(j-1)*halfsize_win+size_win) +  (f_win_rec');
    end
end

f_rec = f_rec / 2^(factor_redund-1);
