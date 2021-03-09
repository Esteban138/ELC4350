%% Exercise 7.12, N = 7
filename='gong2.wav' ;                % name of wave file goes here
N_values = [2^7, 2^9, 2^12, 2^14, 2^16, 2^18];
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;                           % sample interval and # of samples
N=2^7; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Testing N = 7');

%% Exercise 7.12, N = 9
filename='gong2.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;                           % sample interval and # of samples
N=2^9; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Testing N = 9')
%% Exercise 7.12, N = 14
filename='gong2.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;                           % sample interval and # of samples
N=2^14; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Testing N = 14')
%% Exercise 7.12, N = 18
filename='gong2.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr ;                           % sample interval and # of samples
N=2^18; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Testing N = 18')
%% Exercise 7.12, half the time
filename='gong2.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr *2 ;                           % sample interval and # of samples
N=2^15; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Testing N = 15, half time period')
%% Exercise 7.12, third the time
filename='gong2.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr *3 ;                           % sample interval and # of samples
N=2^15; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Testing N = 15, one third of time period')
%% Exercise 7.12, fourth the time
filename='gong2.wav' ;                % name of wave file goes here
[ x , sr ]=audioread( filename ) ;     % read in wavefile
Ts=1/ sr *4 ;                           % sample interval and # of samples
N=2^15; x=x(1:N)';                   % length for analysis
sound(x , 1 / Ts )                   % play sound , if sound card installed
time=Ts * ( 0 : length(x)-1);        % establish time base for plotting
figure
subplot ( 2 , 1 , 1 ) , plot ( time , x )          % and plot top figure
magx=abs ( fft(x) ) ;                % take FFT magnitude
ssf =(0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot( 2 , 1 , 2 ), plot( ssf, magx ( 1:N/2) ) % plot mag spectrum
title('Testing N = 15, fourth time period')