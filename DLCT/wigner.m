function varargout=wigner(varargin)
%WIGNER Wigner distribution
%
%   WD=WIGNER(X,Fs) computes the Wigner distribution of X. By default,
%   the number of frequency points used to compute the discrete Fourier
%   transform is set to 512 (NFFT=512). The returned matrix is of
%   the form WD(F,T).
%
%   WD=WIGNER(X,Fs,XFLAG) gives the option to use the original input signal
%   in the Wigner computation. This is achieved by setting XFLAG to 'NHT'.
%   If XFLAG is not specified, the Hilbert transform of X is used. Type
%   'help hilbert' for more details on Hilbert transform.
%
%   WD=WIGNER(X,Fs,XFLAG,NFFT) specifies the number of frequency points
%   used to calculate the discrete Fourier transform. If NFFT is not
%   specified, the default NFFT is used.
%
%   [WD,F,T]=WIGNER(...) returns a vector of frequencies F and a
%   vector of times T at which the distribution is computed.
%
%   WIGNER(...) with no output arguments plots the distribution along with
%   its frequency marginal and X.
%
%   Example: Display the Wigner distribution of a two-tone signal.
%       [x,fs]=synthetic_test_signal('2tone');
%       wigner(x,fs)
%
%   Reference page in browser: <a href="wigner_doc.html">doc wigner</a>
%   see also <a href="wigner.html">wigner</a>

%   Copyright (c)  2010 Arash Mahboobin, Patrick Loughlin,
%   Applied Signals and Systems Analysis Lab
%   University of Pittsburgh
%   Time-Frequency Signal Analysis, BIOE 3528

%   ACKNOWLEDGEMENT: Portions of this code are extracted from the
%   Time-Frequency Toolbox (TFTB) developed at Centre National de la
%   Recherche Scientifique (CNSR), France.

error(nargchk(1,4,nargin)); %validate number of input arguments
error(nargoutchk(0,3,nargout)); %validate number of output arguments
[msg,x,fs,xflag,nfft]=input_chk(varargin);
error(msg)

[xr,xc]=size(x); %get size of signal
if xc>xr
    x=x'; %transpose signal
    xr=xc;
end
if strcmp(xflag,'ht')
    x=hilbert(x); %Hilbert transform
    x=x./norm(x); %normalize to unit energy
else
    if isreal(x)
        msg1='Warning: ';
        msg2='Results may not be correct for frequencies greater fs/4';
        disp([msg1 msg2])
    end
end
Ex=sum(abs(x).^2); %energy of signal
t=1:xr; %time vector with ts=1 [s]
[tr,tc]=size(t); %get size of time
wd=zeros(nfft,tc); %initialize matrix to store data. Note wd(f,t)
for n=1:tc
tn=t(n); %time sample
max_tau=min([tn-1,xr-tn,round(nfft/2)-1]); %maximum points allowed to shift
tau=-max_tau:max_tau; %shift range
m=rem(nfft+tau,nfft)+1;
wd(m,n)=x(tn+tau).*conj(x(tn-tau));
end
wd=fft(wd);
wd=real(wd);
wd=Ex*wd/sum(sum(wd)); %normalize to signal energy
[fdim,tdim]=size(wd);
F=(0:1:fdim-1)/2/fdim*fs; %frequency vector
T=(0:tdim-1)/fs; %time vector with ts=1/fs [s]
switch nargout
    case 0
        fmarg=sum(wd,2); %frequency marginal
        L1='Time (sec)';
        L2='Freq (Hz)'; %set up some labels for plots
        L3='Wigner distribution of input signal';
        tfd_plot(x,wd,fmarg,T,F,L1,L2,L3)
        colormap('jet')
    case 1
        varargout={wd};
    case 2
        varargout={wd,F};
    case 3
        varargout={wd,F,T};
end

function [msg,x,fs,xflag,nfft]=input_chk(P)
msg=[];
x=P{1}; %signal
if sum(abs(x).^2) ~= 1 %check if signal is normalized
    x=x./norm(x); %normalize signal
end
if (length(P) > 1) && ~isempty(P{2})
    fs=P{2}; %sampling frequency
else
    msg='You must specify a sampling frequency';
end
if (length(P) > 2) && ~isempty(P{3})
    xflag=lower(P{3});
else
    xflag='ht'; %default setting: use Hilbert transform
end
if (length(P) > 3) && ~isempty(P{4})
    nfft=P{4};
else
    nfft=512; %default value
end