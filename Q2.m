%% Question 2
% Submitted by:
% Yotam Leibovitch
% Asaf Bar-El
clear all;
clc;
%% 2.2a
N = 100;% number of coefficients
P = 10;
d = P*rand([1 N]);% coefficients

n = 3; %cubic spline
M = 20;
t = -N-2:(1/M):N+2;% time vector 

% delta between time samples ((t(end)-t(1))/length(t));
step = t(2)-t(1);

% number of time samples in a time interval of length 1
unit_step = floor(1/step); 

x = SplineExpansion(d,t,n);% calculating the spline expansion.
figure(1);
plot(t,x);axis ([-1 N+2 -1 P]);title('2.2a: x(t)');

% calculating shifted spline expansions of order n-1
s = SplineExpansion(d,t,n-1);
k=0.5;
s_shift_r = [zeros(1,ceil(unit_step*k)) s(1:end - unit_step*k)];
s_shift_l = [s(1+unit_step*k:end) zeros(1,ceil(unit_step*k))];

% calculating the derivative of x using shifted splines.
dx_dt = s_shift_l - s_shift_r;

figure(2);
subplot(3,1,1);
plot(t,dx_dt);axis ([-1 N+2 -1 P]);title('2.2b:  dx(t)/dt');

% calculating the derivative of x using differences.
x_diff = diff(x)/(t(2)-t(1));
subplot(3,1,2);
plot(t(2:length(x_diff)+1),x_diff);axis ([-1 N+2 -1 P]);title('dx(t)/dt using differences');
x_diff = [zeros(1) x_diff];

subplot(3,1,3);
plot(t,dx_dt-x_diff); axis ([-1 N+2 -1 P]);title('dx(t)/dt - x_{diff}');
%% 
% 2.4.a
%create a cubic spline and sample it
clear all;
clc;

N = 100;% number of coefficients
P = 10;
d = P*rand([1 N]);% coefficients

n = 3; %cubic spline
M = 20;
t = -N-2:1/M:N+2;% time vector

% calculating the cubic spline expansion
x = SplineExpansion(d,t,n);

figure(3);
subplot(2,1,1);
plot(t,x);axis ([-1 N+2 -1 P]);title('2.4a:  x(t)');

%sample x at the integers
sample_indexes = find(abs(t-floor(t)) <= min(t-floor(t)));

first_index = find(abs(t) <= min(abs(t))); %index for t = 1
last_index = find(abs(t-(N)) <= min(abs(t-(N+2)))); %index for t = N+1

a = find(sample_indexes >= first_index);
sample_indexes = sample_indexes(a);
a = find(sample_indexes <= last_index);
sample_indexes = sample_indexes(a); %all the indexes of t=-1,0,1,..,N+2

c = x(sample_indexes);% sample x

%reconstruct x using the samples
x_rec = interpCubic(c,t);

subplot(2,1,2);
plot(t,x_rec);axis ([-1 N+2 -1 P]);title('x(t) reconstructed');

%% 
% 2.5.b

%insert the h filter from 2.4 a delta function to get the series h(k) of
%order 3
clear all;
clc;

N = 10;
M = 20;

t = -N-2:1/M:N+2; 

delta = [zeros(1,floor(N/2)) 1 zeros(1,floor(N/2))];% delta function

%filter delta to get h(k)
c = delta;
d = filter(1,[1 2-sqrt(3)],c);% filter using the causal filter
d = filter(1,[(2+sqrt(3))/6 1/6],wrev(d));% filter the reverse signal
                                          % using the non causal filter.
d = wrev(d);% reverse back.
h_k =d;

%insert h(k) to spline expansion of order 3
eta_3 = SplineExpansion(h_k,t,3);

%shift the time so the result will be centered around zero
figure(4);
subplot(2,1,1);
plot(t-find(h_k == max(h_k)),eta_3);axis ([-5 5 -1 1]);title('2.5b: \eta^3(t)');

subplot(2,1,2);
plot(t,sinc(t));axis ([-5 5 -1 1]);title('sinc(t)');


%% 2.6
clear all;
clc;

I = im2double(imread('lena.png'));% load the original image
figure(5);
subplot(1,3,1);
imshow(I)
title('2.6b: Original Image');

N = length(I);
M = 20;
t = -N-2:(1/M):N+2;% time vector

% apply the downsample
M = 2;% downsample factor
I_down = I((1:M:end),(1:M:end));% downsampling by factor M
subplot(1,3,2);
imshow(I_down)
title('Downsampled Image');

% interpolation (for upsample)
I_up = zeros(size(I));
I_rows = zeros(N/M,N);

first_index = find(abs(t-1) <= min(abs(t-1))); % index for t = 1
last_index = find(abs(t-(N/M)) <= min(abs(t-(N/M)))); % index for t = N/M

% interpolate each row in the downsampled image
for i = 1:N/M
    d = interpSquare(I_down(i,:),t);% interpolating using spline of order 2
    I_rows(i,:) = d(round(linspace(first_index,last_index,N)));% Upsample
end
% interpolate each column after we interpolated the rows
for i = 1:N
    d = interpSquare(I_rows(:,i),t);% interpolating using spline of order 2
    I_up(:,i) = d(round(linspace(first_index,last_index,N)));% Upsample
end
subplot(1,3,3);
imshow(I_up);
title('Upsampled Image');

%% 2.7
clear all;
clc;

I = im2double(imread('lena.png'));% load the original image
figure(6);
imshow(I)
title('2.7b:  Original Image');

N = length(I);
M = 20;
t = -N-2:(1/M):N+2;

% Upsample with spline of order 3
M = 2;% upsample factor
I_up = zeros(N*M, N*M);
I_rows = zeros(N, N*M);

first_index = find(abs(t-1) <= min(abs(t-1))); % index for t = 1
last_index = find(abs(t-N) <= min(abs(t-N))); % index for t = N

% interpolate each row in the image
for i = 1:N
    d = interpCubic(I(i,:),t);% interpolating using spline of order 3
    I_rows(i,:) = d(round(linspace(first_index,last_index,N*M)));% Upsample the row
end
% interpolate each column after we interpolated the rows
for i = 1:(N*M)
    d = interpCubic(I_rows(:,i),t);% interpolating using spline of order 3
    I_up(:,i) = d(round(linspace(first_index,last_index,N*M)));% Upsample the column
end
figure(7);
imshow(I_up);
title('2.7b: Upsampled Image');

%Downsample using spline of order 1
I_down = zeros(N, N);
I_rows = zeros(N*M, N);
t = -N*M-2:(1/M):N*M+2;% time voctor for the upsampled image
first_index = find(abs(t-1) <= min(abs(t-1))); %index for t = 1
last_index = find(abs(t-N*M) <= min(abs(t-N*M))); %index for t = N*M

% interpolate each row in the upsampled image
for i = 1:(N*M)
    d = SplineExpansion(I_up(i,:),t,1);% interpolating using spline of order 1
    I_rows(i,:) = d(round(linspace(first_index,last_index,N)));% Downsample the row
end
% interpolate each column after we interpolated the rows
for i = 1:N
    d = SplineExpansion(I_rows(:,i),t,1);% interpolating using spline of order 1
    I_down(:,i) = d(round(linspace(first_index,last_index,N)));% Downsample the column
end
figure(8);
imshow(I_down);
title('2.7b: Downsampled Image');