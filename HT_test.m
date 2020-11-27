function [tau_n]=HT_test(x)
% hilbert-based tau estimator for frequency ramp
% gives the index of the rocof step location in tau_n
% tau_n = NaN if no step is found

% Double differentiate input signal
dx = gradient(x);
dx2 = gradient(dx);

% take the analytical signal of dx2
% it is amplitude-modulated by Rf
z = hilbert(dx2);
psi = unwrap(angle(z));
f_i = gradient(psi); 

    NSamples = length(f_i);
    br = 0.10;%80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples)+1; % number of samples to be ignored
    brmask = [brn+1:NSamples-brn];

    az = abs(z);
% compensate the f_i signal with az
    f_icomp = az.*f_i/median(az);
    figure(1); hold off; plot(f_i(brmask));
    
%calculate ri signal
    ri2 = gradient(f_icomp);  
%     dmag = gradient(az);
%     figure(2); hold on; plot((dmag(brmask)));
    
    threshold = 1.5e-6;
    %detection signal d
    d = zeros(1,NSamples);
    ri2br = gradient(ri2(brmask));
    d(brmask(2:end-1)) = abs(ri2br(2:end-1));
    
    % --- figure ---
    figure(3); plot(d); title('d');hold on; plot(threshold*ones(1,NSamples),'k--');
    % ----
    
    % find tau
    [dmax,tau_n] = max(d);
    if dmax<threshold
        tau_n = NaN;
    end
    