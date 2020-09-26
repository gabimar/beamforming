clear all

%coef

Na=25;%number of antennas
Ns=3;%total signals recived by atennas array Ns<=Na-2
Nu=1;%useful signals
Ni=2;%interfering signals
N=1024;%samples number
P=floor(Na/2)+1;%antennas number per subarray P=(Na/2)+1
L=Na-P+1;%correlation matrix number
f=10^9;%carrier frequency
fe=48000;%sampling frequency
lambda=(3*(10^8))/f; %wavelength
beta=(2*pi)/lambda; %propagation constant
dist=lambda/2;%distance between antennas

%emission signal matrix
A=[1 0.9 0.9];%amplitude vector
fs=[3000 6000 9000];%frequency vector
theta=[70 150 30];%DOA vector
s=zeros(Ns,N-1);
for i=1:Ns
    s(i,:)=A(i)*exp(1i*2*pi*fs(i)/fe*(1:N-1));
end
e1=zeros(Na,Ns);
for i=1:Na
    for j=1:Ns
        e1(i,j)=exp((-1i)*beta*dist*(i-1)*cos(theta(j)/180*pi));%???(vector directie)
    end
end
y1=e1*s;%emission signals
y2=awgn(y1,10);%reception signals(added noise SNR=10-important!!!)

%CAPON algorithm
R=zeros(Na,Na);
R=(1/(N-1))*(y2*y2');
e2=zeros(Na,181);
for i=1:Na
    for j=1:181
        e2(i,j)=exp((-1i)*beta*dist*(i-1)*cos((j)/180*pi));%???
    end
end
f1=zeros(1,181);
for j=1:181
    f1(j)=1/((((e2(:,j))')*inv(R)*e2(:,j)));
end

%MUSIC algorithm
R=(1/(N-1))*(y2*y2');
[v,d]=eig(R);%???
d2=diag(d);
[d3,poz]=sort(d2,'descend');
poz2=poz((Ns+1):Na);
q=zeros(Na,Na-Ns);

for i=1:Na-Ns
    q(:,i)=v(:,poz2(i));
end
e2=zeros(Na,181);
for i=1:Na
    for j=1:181
        e2(i,j)=exp((-1i)*beta*dist*(i-1)*cos((j)/180*pi));
    end
end
f2=zeros(1,181);
for j=1:181
    f2(j)=1/((((e2(:,j))')*q*(q')*e2(:,j)));
end

%Smooth MUSIC algorithm
Re=zeros(P,P);%P-antennas number per subarray!!!
y3=zeros(P,N);
for i=1:L
    y3=y2(i:(i+P-1),:);
    Re=Re+(1/(N-1))*(y3*y3');
end
Re=(1/L)*Re;
[v2,d4]=eig(Re);
d5=diag(d4);
[d6,poz3]=sort(d5,'descend');
poz4=poz3((Ns+1):P);
q2=zeros(P,P-Ns);
for i=1:P-Ns
    q2(:,i)=v2(:,poz4(i));
end
e3=zeros(P,181);
for i=1:P;
    for j=1:181
        e3(i,j)=exp((-1i)*beta*dist*(i-1)*cos((j)/180*pi));
    end
end
f3=zeros(1,181);
for j=1:181
    f3(j)=1/((((e3(:,j))')*q2*(q2')*e3(:,j)));
end
    %sortarea DOA(direction of arrival)
[peaks,thetad]=findpeaks(real(f3));
[peaks2,locs]=sort(peaks,'descend');
%locs2=sort(locs(1:Ns),'descend');
thetaf=thetad(locs(1:Ns));%by me

%Beamforming
W=zeros(1,Na);
W((Ns+1):Na)=1;
delta=zeros(Ns);
%   matricea delta
for i=1:Ns
    for j=1:Ns
        delta(i,j)=exp((-1i)*beta*dist*(j-1)*cos(thetaf(i)/180*pi));
    end
end
    %free terms
V=zeros(Ns,1);
V(1:Nu)=Na/Nu;
V2=zeros(Ns,1);
for i=1:Ns
    for j=(Ns+1):Na
        V2(i)=V(i)-exp((-1i)*beta*dist*(j-1)*cos(thetaf(i)/180*pi));
    end
end
    %deltaN matrix and coef W determination
deltaN=zeros(Ns);
deltaN=delta;
for j=1:Ns
    deltaN(:,j)=V2;
    W(j)=(det(deltaN))/(det(delta));
    deltaN(:,j)=delta(:,j);
end
C=zeros(1,181);
e4=zeros(Na,181);
for i=1:Na
    for j=1:181
        e4(i,j)=exp((-1i)*beta*dist*(i-1)*cos(j/180*pi));
    end
end
C=W*e4;
y4=sum(y2); %unfiltered signal with beamforming
yf=W*y2; %filtered signal using beamforming

%samples transformation in frequency domain for ploting TF for y4 and yf signals
N2=200; %number of samples for y4 and yf graphs and Fourier Transformation of them
xfft=fe/N*(1:N-1);
yfft1=abs(fft(y4(1:N-1)));
yfft2=abs(fft(yf(1:N-1)));
%frequency finder of each signal with TF and signal sorting
[ampld_y4,fd_y4]=findpeaks(yfft1);
[amplf_y4,locf]=sort(ampld_y4,'descend');
ff_y4=fd_y4(locf(1:Ns));
ff_y4f=fe/N*(ff_y4);

% [ampld_yf,fd_yf]=findpeaks(yfft2);
% [amplf_yf,locf2]=sort(ampld_yf,'descend');
% ff_yf=fd_yf(locf2(1:Nu));
% ff_yff=fe/N*(ff_yf);
% yfft3=zeros(1,1023);
%     for i=1:Nu
%         yfft3=yfft3+[zeros(1,ff_yf(i)-2) abs(fft(yf((ff_yf(i)-1):(ff_yf(i)+1))) zeros(1,N-1-ff_yf(i))];
%             end
%[peaksfftd1,locsfftd1]=findpeaks(abs(fft(yf))); %??
%[peaksfftd2,locsfftd2]=sort(peaksfftd1,'descend');

%Comparation between the three methods
figure
subplot(3,1,1);
set(gcf,'color','w');
plot(10*log10(real((f1(1:180)))));
title(['Capon', ' (Na=',num2str(Na),', Ns=',num2str(Ns),', Nu=',num2str(Nu),', Ni=',num2str(Ni),')'],'fontweight','bold','fontsize',14)
xlabel('Direction of Arrival[grades]','fontsize',14);
ylabel('Capon Spectrum[dB]','fontsize',14);
subplot(3,1,2);
set(gcf,'color','w');
plot(10*log10((real(f2(1:180)))));
title(['MUSIC' ' (Na=',num2str(Na),', Ns=',num2str(Ns),', Nu=',num2str(Nu),', Ni=',num2str(Ni),')'],'fontweight','bold','fontsize',14)
xlabel('Direction of Arrival[grades]','fontsize',14);
ylabel('MUSIC Spectrum[dB]','fontsize',14);
subplot(3,1,3);
plot(10*log10((real(f3(1:180)))));
title(['Smooth MUSIC' ' (Na=',num2str(Na),', Ns=',num2str(Ns),', Nu=',num2str(Nu),', Ni=',num2str(Ni),', P=',num2str(P),')'],'fontweight','bold','fontsize',14);
xlabel('Direction of Arrival[degrees]','fontsize',14);
ylabel('Sm. MUSIC Spectrum[dB]','fontsize',14);

%Antenna array radiation pattern graph
figure
set(gcf,'color','w');
plot(abs(C(1:180)));
title(['Radiation patern C' ' (maximum direction=',num2str(thetaf(1:Nu)),', minimum direction=',num2str(thetaf(Nu+1:Ns)),')'],'fontweight','bold','fontsize',14)
xlabel('Angle of incidence [degrees]','fontsize',14);
ylabel('Array transfer factor','fontsize',14);

%yf, y4, TF(yf) and TF(y4) graphs
figure
subplot(2,2,1);
set(gcf,'color','w');
plot(real(y4(1:200)));
title('y4 signal recived without Beamforming','fontweight','bold','fontsize',14);
xlabel('Samples number','fontsize',14);
ylabel('Amplitude','fontsize',14);
subplot(2,2,2);
set(gcf,'color','w');
plot(real(yf(1:200)));
title('yf signal recived with Beamforming','fontweight','bold','fontsize',14);
xlabel('Samples number','fontsize',14);
ylabel('Amplitude','fontsize',14);
subplot(2,2,3);
set(gcf,'color','w');
fft_y4=plot(xfft-fe/N,yfft1);
title(['Fourier transfor of y4',' (fu=',num2str(ff_y4f(1:Nu)),' fi=',num2str(ff_y4f(Nu+1:Ns)),')'],'fontweight','bold','fontsize',14);
xlabel('Frequency[Hz]','fontsize',14);
ylabel('Amplitude','fontsize',14);
subplot(2,2,4);
set(gcf,'color','w');
fft_yf=plot(xfft-fe/N,yfft2);
title('Fourier transform of yf','fontweight','bold','fontsize',14);
xlabel('Frequency[Hz]','fontsize',14);
ylabel('Amplitude','fontsize',14);
