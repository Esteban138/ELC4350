%% Exercise 11.3
alphabet = [2, 3, 5];
variance = [1, 5, 35/3];
for iter = 1:3
    figure()                     % used to plot figure eyediag3
    N=1000; m=pam(N,alphabet(iter),variance(iter));          % random +/-1 signal of length N
    M=20; mup=zeros(1,N*M); mup(1:M:N*M)=m;  % oversampling by factor of M
    ps=ones(1,M);                            % square pulse width M
    x=filter(ps,1,mup);            % convolve pulse shape with mup
    neye=5;
    c=floor(length(x)/(neye*M))
    xp=x(N*M-neye*M*c+1:N*M);      % dont plot transients at start
    q=reshape(xp,neye*M,c);        % plot in clusters of size 5*Mt=(1:198)/50+1;
    subplot(3,1,1), plot(q)
    title(['Eye diagram for rectangular pulse shape using ',num2str(alphabet(iter))])

    N=1000; m=pam(N,alphabet(iter),variance(iter));          % random +/-1 signal of length N
    M=20; mup=zeros(1,N*M); mup(1:M:N*M)=m;  % oversampling by factor of M
    ps=hamming(M);                           % square pulse width M
    x=filter(ps,1,mup);            % convolve pulse shape with mup
    %x=x+0.15*randn(size(x));
    neye=5;
    c=floor(length(x)/(neye*M))
    xp=x(N*M-neye*M*c+1:N*M);      % dont plot transients at start
    q=reshape(xp,neye*M,c);        % plot in clusters of size 5*Mt=(1:198)/50+1;
    subplot(3,1,2), plot(q)
    title(['Eye diagram for hamming pulse shape using ',num2str(alphabet(iter))])

    N=1000; m=pam(N,alphabet(iter),variance(iter));          % random +/-1 signal of length N
    M=20; mup=zeros(1,N*M); mup(1:M:N*M)=m;  % oversampling by factor of M
    L=10; ps=srrc(L,0,M,50);
    ps=ps/max(ps);         % sinc pulse shape L symbols wide
    x=filter(ps,1,mup);    % convolve pulse shape with mup
    %x=x+0.15*randn(size(x));
    neye=5;
    c=floor(length(x)/(neye*M))
    xp=x(N*M-neye*M*(c-3)+1:N*M);  % dont plot transients at start
    q=reshape(xp,neye*M,c-3);      % plot in clusters of size 5*Mt=(1:198)/50+1;
    subplot(3,1,3), plot(q)
    axis([0,100,alphabet(iter) *-1,alphabet(iter)])
    title(['Eye diagram for sinc pulse shape using ',num2str(alphabet(iter))])
end

%% Exercise 11.4
random_values = [0.075,0.15, 0.35 , 0.45, 0.55];
for iter = 1:5
    figure                      % used to plot figure eyediag3
    N=1000; m=pam(N,2,1);          % random +/-1 signal of length N
    M=20; mup=zeros(1,N*M); mup(1:M:N*M)=m;  % oversampling by factor of M
    ps=ones(1,M);                            % square pulse width M
    x=filter(ps,1,mup);            % convolve pulse shape with mup
    x=x+random_values(iter)*randn(size(x));
    neye=5;
    c=floor(length(x)/(neye*M))
    xp=x(N*M-neye*M*c+1:N*M);      % dont plot transients at start
    q=reshape(xp,neye*M,c);        % plot in clusters of size 5*Mt=(1:198)/50+1;
    subplot(3,1,1), plot(q)
    title(['Eye diagram for rectangular pulse shape using ',num2str(random_values(iter))])

    N=1000; m=pam(N,2,1);          % random +/-1 signal of length N
    M=20; mup=zeros(1,N*M); mup(1:M:N*M)=m;  % oversampling by factor of M
    ps=hamming(M);                           % square pulse width M
    x=filter(ps,1,mup);            % convolve pulse shape with mup
    x=x+random_values(iter)*randn(size(x));
    neye=5;
    c=floor(length(x)/(neye*M))
    xp=x(N*M-neye*M*c+1:N*M);      % dont plot transients at start
    q=reshape(xp,neye*M,c);        % plot in clusters of size 5*Mt=(1:198)/50+1;
    subplot(3,1,2), plot(q)
    title(['Eye diagram for hamming pulse shape using ',num2str(random_values(iter))])

    N=1000; m=pam(N,2,1);          % random +/-1 signal of length N
    M=20; mup=zeros(1,N*M); mup(1:M:N*M)=m;  % oversampling by factor of M
    L=10; ps=srrc(L,0,M,50);
    ps=ps/max(ps);         % sinc pulse shape L symbols wide
    x=filter(ps,1,mup);    % convolve pulse shape with mup
    x=x+random_values(iter)*randn(size(x));
    neye=5;
    c=floor(length(x)/(neye*M))
    xp=x(N*M-neye*M*(c-3)+1:N*M);  % dont plot transients at start
    q=reshape(xp,neye*M,c-3);      % plot in clusters of size 5*Mt=(1:198)/50+1;
    subplot(3,1,3), plot(q)
    axis([0,100,random_values(iter) *-1,random_values(iter)])
    title(['Eye diagram for sinc pulse shape using ',num2str(random_values(iter))])
end

% ANSWER
% It appears that the largest value of v is around 0.5, anything after that
% and it appears that the "eye" is no longer visible.

%% Exercise 12.1 a
% clockrecDD.m: clock recovery minimizing 4-PAM cluster variance
% to minimize J(tau) = (Q(x(kT/M+tau))-x(kT/M+tau))^2

% prepare transmitted signal
n=10000;                         % number of data points
m=2;                             % oversampling factor
beta=0.3;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap=srrc(l,beta,m,toffset);  % srrc pulse shape with timing offset
s=pam(n,4,5);                    % random data sequence with var=5
sup=zeros(1,n*m);                % upsample the data by placing...
sup(1:m:n*m)=s;                  % ... m-1 zeros between each data point
hh=conv(pulshap,chan);           % ... and pulse shape
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,m,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

% clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu_values = [0.005, 0.01, 0.02, 0.05]
mu=2;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-2*l*m          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% plot results
figure();
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

%ANSWER
% It appears that mu seems the affect the amount of noise present in the
% signal, around 0.2 noise can start to appear visible. However if mu goes
% as large as 2 it no longer converges to the correct value of the offset. 

%% Exercise 12.1 b Testing 2-PAM
% prepare transmitted signal
n=10000;                         % number of data points
m=2;                             % oversampling factor
beta=0.3;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap=srrc(l,beta,m,toffset);  % srrc pulse shape with timing offset
s=pam(n,2,1);                    % random data sequence with var=5
sup=zeros(1,n*m);                % upsample the data by placing...
sup(1:m:n*m)=s;                  % ... m-1 zeros between each data point
hh=conv(pulshap,chan);           % ... and pulse shape
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,m,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

% clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.01;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-2*l*m          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% plot results
figure();
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')
title('Testing 2-PAM');
%% Exercise 12.1 b Testing 6-PAM
% prepare transmitted signal
n=10000;                         % number of data points
m=2;                             % oversampling factor
beta=0.3;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap=srrc(l,beta,m,toffset);  % srrc pulse shape with timing offset
s=pam(n,6,11);                    % random data sequence with var=5
sup=zeros(1,n*m);                % upsample the data by placing...
sup(1:m:n*m)=s;                  % ... m-1 zeros between each data point
hh=conv(pulshap,chan);           % ... and pulse shape
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,m,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

% clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.01;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-2*l*m          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% plot results
figure();
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')
title('Testing 6-PAM');
 % Answer
 % It appears that the more levels are used to represent a signal more
 % noise appears in the signal as well.
 
 %% Exercise 12.2 
 % prepare transmitted signal
n=10000;                         % number of data points
m=2;                             % oversampling factor
beta=0.3;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
x = (m+toffset)/2:beta:(m+toffset)/2;
%pulshap = 1*rectpuls(x, 1);
pulshap = rect_pulse_maker(1, beta, m, toffset);
s=pam(n,4,5);                    % random data sequence with var=5
sup=zeros(1,n*m);                % upsample the data by placing...
sup(1:m:n*m)=s;                  % ... m-1 zeros between each data point
hh=conv(pulshap,chan);           % ... and pulse shape
r=conv(hh,sup);                  % ... to get received signal
%matchfilt=srrc(l,beta,m,0);      % matched filter = srrc pulse shape
matchfilt = rect_pulse_maker(1, beta, m, 0);
x=conv(r,matchfilt);             % convolve signal with matched filter

% clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.01;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-2*l*m          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% ANSWER
% It appears that a rectuangular pulse has very poor spectral properites
% and because it varies so high from the maxima you end with a lot of
% interference. A raised consine can control the excess bandwidth but that
% is hard to do in a rectangular pulse leading to a noisy signal that does
% not converge properly.

%% Exercise 12.3

% clockrecDD.m: clock recovery minimizing 4-PAM cluster variance
% to minimize J(tau) = (Q(x(kT/M+tau))-x(kT/M+tau))^2

% prepare transmitted signal
n=10000;                         % number of data points
m=2;                             % oversampling factor
beta=0.3;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap=srrc(l,beta,m,toffset);  % srrc pulse shape with timing offset
s=pam(n,4,5);                    % random data sequence with var=5
sup=zeros(1,n*m);                % upsample the data by placing...
sup(1:m:n*m)=s;                  % ... m-1 zeros between each data point
hh=conv(pulshap,chan);           % ... and pulse shape
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,m,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);        % convolve signal with matched filter
x=x+0.15*randn(size(x));
% clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.01;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-2*l*m          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% plot results
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')
% ANSWER
% It appears that it stills converges to the correct value however it is a
% lot more noisy then before.

%% Exercise 13.1 
% LSequalizer.m find a LS equalizer f for the channel b
 b=[0.5 1 -0.6];                    % define channel
    figure
    plot(freqz(b));
    title('Plotting channel ');
    [hb, wb] = freqz(b);
    m=1000; s=sign(randn(1,m));        % binary source of length m
    r=filter(b,1,s);                   % output of channel
    n=3;                               % length of equalizer - 1
for iter = 1:4
    delta=iter-1;                           % use delay <=n*length(b)
    p=length(r)-delta;
    R=toeplitz(r(n+1:p),r(n+1:-1:1));  % build matrix R
    S=s(n+1-delta:p-delta)';           % and vector S
    f=inv(R'*R)*R'*S;                   % calculate equalizer f
    Jmin=S'*S-S'*R*inv(R'*R)*R'*S      % Jmin for this f and delta
    y=filter(f,1,r);                   % equalizer is a filter
    dec=sign(y);                       % quantize and find errors
    err=0.5*sum(abs(dec(delta+1:m)-s(1:m-delta)))
    figure
    plot(freqz(b,f))
    title(['Trying equalizer with delta value ', num2str(iter)])
    
    [h,w] = freqz(b,f);
    figure
    plot(abs(h.*hb))
    title(['Plotting magnitude of frequency response ', num2str(iter)]);
end
% ANSWER
% Even though the magnitude of teh channel and the equalizers seems to be
% opposite they never really cancel out and reach unity. 

%% Exercise 13.2 a 
% LSequalizer.m find a LS equalizer f for the channel b
b=[0.5 1 -0.6];                    % define channel
m=1000; s=sign(randn(1,m));        % binary source of length m
sd = 0.3;
r=filter(b,1,s)+sd*randn(size(s));                   % output of channel
n=3;                               % length of equalizer - 1
delta=2;                           % use delay <=n*length(b)
p=length(r)-delta;
R=toeplitz(r(n+1:p),r(n+1:-1:1));  % build matrix R
S=s(n+1-delta:p-delta)';           % and vector S
f=inv(R'*R)*R'*S                   % calculate equalizer f
Jmin=S'*S-S'*R*inv(R'*R)*R'*S      % Jmin for this f and delta
y=filter(f,1,r);                   % equalizer is a filter
dec=sign(y);                       % quantize and find errors
err=0.5*sum(abs(dec(delta+1:m)-s(1:m-delta)))
% ANSWER
% It appears that the value where error first start appearing is around 0.3

%% Exercise 13.2b
x = []
for sd = 0.1:0.05:1
    b=[0.5 1 -0.6];                    % define channel
    m=1000; s=sign(randn(1,m));        % binary source of length m
    r=filter(b,1,s)+sd*randn(size(s));                   % output of channel
    n=3;                               % length of equalizer - 1
    delta=2;                           % use delay <=n*length(b)
    p=length(r)-delta;
    R=toeplitz(r(n+1:p),r(n+1:-1:1));  % build matrix R
    S=s(n+1-delta:p-delta)';           % and vector S
    f=inv(R'*R)*R'*S                   % calculate equalizer f
    Jmin=S'*S-S'*R*inv(R'*R)*R'*S      % Jmin for this f and delta
    y=filter(f,1,r);                   % equalizer is a filter
    dec=sign(y);                       % quantize and find errors
    err=0.5*sum(abs(dec(delta+1:m)-s(1:m-delta)))
    x = [x, Jmin]

end

figure
plot(0.1:0.05:1, x)
title('Jmin vs sd plot')

%% Exercise 13.2 c
% LSequalizer.m find a LS equalizer f for the channel b
b=[0.5 1 -0.6];                    % define channel
m=1000; s=sign(randn(1,m));        % binary source of length m
sd = 0.19;
r=filter(b,1,s)+sd*randn(size(s));                   % output of channel
n=3;                               % length of equalizer - 1
delta=1;                           % use delay <=n*length(b)
p=length(r)-delta;
R=toeplitz(r(n+1:p),r(n+1:-1:1));  % build matrix R
S=s(n+1-delta:p-delta)';           % and vector S
f=inv(R'*R)*R'*S                   % calculate equalizer f
Jmin=S'*S-S'*R*inv(R'*R)*R'*S      % Jmin for this f and delta
y=filter(f,1,r);                   % equalizer is a filter
dec=sign(y);                       % quantize and find errors
err=0.5*sum(abs(dec(delta+1:m)-s(1:m-delta)))


%ANSWER
% It appears that when sd is around 0.2 error appears

%% Exercise 13.2 d
% ANSWER
% The first one with the larger delay seems to perform better, because the
% value of sd can be larger than the one with the smaller delay. This means
% it can take larger noise values and is therefore a more robust system.
