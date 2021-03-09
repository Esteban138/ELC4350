%% 7.5 a
% Trying f = 200, 300, 450, 550, 600, 800, 2200 Hz
f = [200, 300, 450, 550, 600, 800, 2220]; 
Ts=1/1000; 
time=10.0; % f req , sampling inte rval , time
t=Ts : Ts : time ; % def ine a time vector
for x = 1:7
 w=sin (2* pi* f(x) *t ) ; % de f ine the s inuso id
N=2^10; % s i z e of a n a l y s i s window
ssf=(-N/2:N/2-1)/(Ts*N) ; % f requency vector
fw=fft (w( 1 :N)) ; % do DFT/FFT
fws=fftshift( fw ) ; % s h i f t i t for plot t ing
figure
plot (ssf , abs( fws ))
title(['Frequency =  ',num2str(f(x)),'Hz'])
end



%% 7.5 b
% Testing different Ts values
f = 100;
Ts= [(1/500), (1/250), (1/50)];
time=50.0; % f req , sampling inte rval , time
for x = 1:3
t=Ts(x) : Ts(x) : time ; % def ine a time vector
 w=sin (2* pi* f *t ); % de f ine the s inuso id
N=2^10; % s i z e of a n a l y s i s window
ssf=(-N/2:N/2-1)/(Ts(x)*N) ; % f requency vector
fw=fft (w(1:N)) ; % do DFT/FFT
fws=fftshift(fw) ; % s h i f t i t for plot t ing
figure
plot (ssf , abs( fws ))
title(['Trying period Ts =  ',num2str(Ts(x)),'s'])
end
%% 7.5 c
% Testing different values of N
N_values = [11, 14, 8, 4, 2, 20];
f = 100;
Ts= 1/1000;
time=2000.0; % f req , sampling inte rval , time
for x = 1:6
t=Ts:Ts:time ; % def ine a time vector
 w=sin (2* pi* f *t ); % de f ine the s inuso id
N=2^N_values(x); % s i z e of a n a l y s i s window
ssf=(-N/2:N/2-1)/(Ts*N) ; % f requency vector
fw=fft (w(1:N)) ; % do DFT/FFT
fws=fftshift(fw) ; % s h i f t i t for plot t ing
figure
plot (ssf ,abs( fws ))
title(['Trying value of N=  ',num2str(N_values(x))])
end


%% Exercise 7.6
f = 100;
Ts= 1/1000;
time=10.0; % f req , sampling inte rval , time
t=Ts:Ts:time ; % def ine a time vector
w=sin (2* pi* f *t ).^2; % de f ine the s inuso id
N=2^10; % s i z e of a n a l y s i s window
ssf=(-N/2:N/2-1)/(Ts*N) ; % f requency vector
fw=fft (w(1:N)) ; % do DFT/FFT
fws=fftshift(fw) ; % s h i f t i t for plot t ing
figure
plot (ssf ,abs( fws ))
title('Spectrum sin ^2')

w=sin (2* pi* f *t ).^3;
fw=fft (w(1:N)) ; % do DFT/FFT
fws=fftshift(fw) ; % s h i f t i t for plot t ing
figure
plot (ssf ,abs( fws ))
title('Spectrum sin ^3')

k = 0.0001;
w=sin (2* pi* f *t ).^k;
fw=fft (w(1:N)) ; % do DFT/FFT
fws=fftshift(fw) ; % s h i f t i t for plot t ing
figure
plot (ssf ,abs( fws ))
title('Spectrum sin ^k, where k = 0.0001')

% The bigger the k value the smaller the amplitude value at each
% peak and the smaller the k value it convolves into a single point

%% Exercise 7.7
f = 100;
Ts= 1/1000;
time=20.0; % f req , sampling inte rval , time
t=Ts:Ts:time ; % def ine a time vector
w=sinc(2* pi* f *t ); % de f ine the s inuso id
N=2^10; % s i z e of a n a l y s i s window
ssf=(-N/2:N/2-1)/(Ts*N) ; % f requency vector
fw=fft (w(1:N)) ; % do DFT/FFT
fws=fftshift(fw) ; % s h i f t i t for plot t ing
figure
plot (ssf ,abs( fws ))
title('Spectrum sinc ')

% Testing sinc ^2
w=sinc(2* pi* f *t ).^2; % de f ine the s inuso id
N=2^10; % s i z e of a n a l y s i s window
ssf=(-N/2:N/2-1)/(Ts*N) ; % f requency vector
fw=fft (w(1:N)) ; % do DFT/FFT
fws=fftshift(fw) ; % s h i f t i t for plot t ing
figure
plot (ssf ,abs( fws ))
title('Spectrum sinc ^2')

%% Exercise 7.8
f = 100;
Ts= 1/1000;
time=10.0; % f req , sampling inte rval , time
t=Ts:Ts:time ; % def ine a time vector
w=sin(t) + x*exp(-t); % de f ine the s inuso id
N=2^10; % s i z e of a n a l y s i s window
ssf=(-N/2:N/2-1)/(Ts*N) ; % f requency vector
fw=fft (w(1:N)) ; % do DFT/FFT
fws=fftshift(fw) ; % s h i f t i t for plot t ing
figure
plot (ssf ,abs( fws ))
title('Spectrum for w(t) = sin(t) + je^ -t ')

% Using specsin2.m is probably the way to go because of the fftshift 
% organizing the frequencies nicely.

%% Exercise 7.9
f = 100;
Ts= 1/1000;
phi = [0, 0.2, 0.4, 0.8, 1.5, 3.14];
time=10.0; % f req , sampling inte rval , time
t=Ts:Ts:time ; % def ine a time vector
for x = 1:6
    w=sin((2*pi*f*t) + phi(x)); % de f ine the s inuso id
    N=2^10; % s i z e of a n a l y s i s window
    ssf=(-N/2:N/2-1)/(Ts*N) ; % f requency vector
    fw=fft (w(1:N)) ; % do DFT/FFT
    fws=fftshift(fw) ; % s h i f t i t for plot t ing
    figure
    plot (ssf ,unwrap(angle( fws )))
    title(['Phi value of ',num2str(phi(x)), ' for sin'])
end


% Finding phase output of sin.^2
for x = 1:6
    w=sin((2*pi*f*t) + phi(x)).^2; % de f ine the s inuso id
    N=2^10; % s i z e of a n a l y s i s window
    ssf=(-N/2:N/2-1)/(Ts*N) ; % f requency vector
    fw=fft (w(1:N)) ; % do DFT/FFT
    fws=fftshift(fw) ; % s h i f t i t for plot t ing
    figure
    plot (ssf ,unwrap(angle( fws )))
    title(['Phi value of ',num2str(phi(x)), ' for sin ^2'])
end
%% Exercise 7.10
% Finding the value of N for 0.1 s
% s = (1/Ts) * 2 ^ N 
filename='gong.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;                           % sample interval and # of samples
N=2^12; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);  % establish time base for plotting
figure 
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('First 0.1s');

% Exercise 7.10 Taking Middle Sample
% Taking value in the middle of the sound
filename='gong.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;  
N=2^12;
x=x(N/2:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ;  % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Middle 0.1s');
%% Exercise 7.11

filename='gong.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;                           % sample interval and # of samples
N=2^15; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), semilogy( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Semilogy');


%% Exercise 7.13
filename='medieval.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;                           % sample interval and # of samples
N=2^18; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);  % establish time base for plotting
figure 
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('First 0.1s');

% Exercise 7.10 Taking Middle Sample
% Taking value in the middle of the sound
filename='gong.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;  
N=2^12;
x=x(N/2:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ;  % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Middle 0.1s');

%% Exercise 7.15
clc 
clear all
a=[0.9] ; lena=length ( a)-1; % autor egr es s ive c o e f f i c i e n t s
b= [2] ; lenb=length (b ) ; % moving average c o e f f i c i e n t s
d=randn ( 1 ,20) ;
h=impz (b , a );
if lena>=lenb  % dimpulse needs lena>=lenb % impulse response of f i l t e r
yfilt=filter(h , 1 , d) % f i l t e r x [ k ] with h [ k ]
end
IIR=filter(b , a , d)



% Creating FIR
FIR = conv(h,d);
disp('IIR Results');
disp(IIR)
disp('FIR Results');
disp(FIR);
disp('Results are nearly identical');
