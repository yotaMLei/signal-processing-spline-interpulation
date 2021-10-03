%% Question 1
% Submitted by:
% Yotam Leibovitch
% Asaf Bar-El
clear all;
clc;
%% 1.2
N = 200;% number of coefficients.
T = -50:0.1:50;% time vector, N supposed to be greater then T
a = 10*rand([1 N]);% coefficients
tau = 10;

% calculate x(t)
x = zeros(size(T));
for t = 1:length(T)
    for n = -N/2:(N/2-1)
       g_t_n = (1/tau)*exp(-(T(t)-n)/tau)*heaviside(T(t)-n);% g(t-n)
       x(t) = x(t)+tau*a(n+N/2+1)*g_t_n;
    end
end

X = fftshift(fft(x));% frequency domain
f = linspace(0,length(T),length(X));

figure(1);
subplot(1,2,1);plot(T,x);axis ([-50 50 20 70]);title('1.2:  x(t)');
subplot(1,2,2);plot(f,abs(X));axis tight;title('X(f)');

%% 1.3
% Nyquist reconstruction

s = sinc(T);
c = conv(x,s,'same')*(T(2)-T(1));% sample using the sampling filter s.

g = s;
x_hat = conv(c,g,'same')*(T(2)-T(1));% reconstruct using the reconstruction filter g.

X_hat = fftshift(fft(x_hat));% frequency domain.
f_hat = linspace(0,length(T),length(X_hat));

figure(2);
subplot(1,2,1);plot(T,x_hat);axis ([-50 50 20 70]);title('1.3:  $$\hat{x}$$(t) using sinc','Interpreter','Latex');
subplot(1,2,2);plot(f_hat,abs(X_hat));axis tight;title('$$\hat{X}$$(f)','Interpreter','Latex');

%% 1.5

%rect=@(x,a) ones(1,numel(x)).*(abs(x)<a/2); % a is the width of the pulse

s_n = zeros(size(T));
g_n = zeros(size(T));

n = -N/2:(N/2-1);

S_n = zeros(length(n), length(T));
G_n = zeros(length(n), length(T));

% calculating the digital correction filter H in frequency domain.
for i = 1:length(n)
    for t = 1:length(T)
        s_n(t) = sinc(T(t))*exp(2*pi*1i*n(i)*T(t));% frequency shifted sinc
        g_n(t) = (1/tau)*exp(-T(t)/tau)*heaviside(T(t))*exp(2*pi*1i*n(i)*T(t));% frequency shifted g
    end
    S_n(i,:) = fft(s_n);% the sampling filter in frequency domain shifted by n.
    G_n(i,:) = fft(g_n);% the reconstruction filter in frequency domain shifted by n.
end

S_G = S_n.*G_n;% S(f-n)*G(f-n)
S_G_sum = sum(S_G,1);% sum over all n
H = 1./S_G_sum;

% calculating the sampling and reconstruction filters s and g in time domain.
s = sinc(T);
g = zeros(1,length(T));
for t = 1:length(T)
    g(t) = (1/tau)*exp(-T(t)/tau)*heaviside(T(t));
end

% calculating the sampling and reconstruction filters S and G in frequency domain.
S = fftshift(fft(s));
G = fftshift(fft(g));

% sampling, correcting and reconstructing x in frequency domain.
X_hat = X.*S.*H.*G;
f_hat = linspace(0,length(T),length(X_hat));

% reconstructed x in time domain.
x_hat = ifft(X);
t_hat = linspace(T(1),T(end),length(x_hat));

figure(3);
subplot(1,2,1);plot(t_hat,abs(x_hat));axis ([-50 50 20 70]);title('1.5:  $$\hat{x}$$(t) using g(t)','Interpreter','Latex');
subplot(1,2,2);plot(f_hat,abs(X_hat));axis tight;title('$$\hat{X}$$(f)','Interpreter','Latex');

