clc;
close all;
clear all;
%% setting part
j=10000;
Q=zeros(1,j);
A=zeros(1,j);
sigma=1;
K=normrnd(0,sigma,1,j);
B=zeros(1,j);
C1=zeros(1,j);
C2=zeros(1,j);
En=zeros(1,j);
P=zeros(1,j);
P(1)=1;
En(1)=0;
TH=zeros(1,j);
delta=zeros(1,j+1);
R1=2;
M1=2;
M2=2;
E=zeros(1,j);
I=zeros(1,j);
U1=zeros(1,j);
U2=zeros(1,j);
Ia=zeros(1,j);
Imax=0;
a=zeros(1,j);
b=zeros(1,j);
n=zeros(1,j);
n(1)=1;


W=zeros(1,j+1);
col=j;
lin=1;
S=unifrnd(0,1,lin,col);
S1=unifrnd(0,1,lin,col);
S2=unifrnd(0,1,lin,col);
p=0.8;
if S1(1)<=p
   S(1)=1;
else
   S(1)=0;
end
A(1)=K(1);
phi1=0.8;
phi2=0.7;
V=1;
Z=1;
W(1)=1;
m=100;
if (S2(1)>phi1)
        B(1)=0;
    else
        B(1)=1;
end
W(1)=1;
for t=2:j
%     if (t>3950)&&(t<4000)
%         W(t)=100;
%     elseif (t>6950)&&(t<7000)
%         W(t)=100;
%     else
%         W(t)=1;
%     end
W(t)=1;
    determine=mod(t,m);
    if (determine==0)
        delta(t)=1;
    else
        delta(t)=0;
    end
    if (S2(t)>phi1)
        B(t)=0;
    else
        B(t)=1;
    end
end
Wa=mean(W);
if (P(1)>=1)
    C1(1)=1;
else
    C1(1)=0;
end
if ((P(1)-U1(1))>=M1)
    C2(1)=1;
else
    C2(1)=0;
end
theta=(2*Wa*R1/m)*(((1+M1)/(M1*p*phi1))-1);
beta=theta+(2*Wa*R1*M2*(1+M1)/(m*phi1))+2*Wa*R1*(1-M2)/m;

a(1)=V*En(1)-Z*TH(1)*S(1)*C1(1)-0.5*S(1)*C1(1)*theta*Q(1)^2-W(2)*R1*delta(2)*S(1)*C1(1)*Q(1)^2;
b(1)=V*M1*En(1)+0.5*theta*n(1)*C2(1)*sigma^2-0.5*beta*C2(1)*n(1)*sigma^2+(1-M2)*delta(2)*R1*W(2)*n(1)*sigma^2;

if b(1)>=0
    U1(1)=0;
    U2(1)=0;
else
    U1(1)=0;
    U2(1)=1;
end
E(1)=Q(1)^2+M2*A(1)^2;
I(1)=W(1)*(Q(1)^2+M2*A(1)^2);
Ia(1)=I(1);
if U2(1)==1
    n(2)=1;
else
    n(2)=2;
end
P(2)=P(1)+B(1)-U1(1)-M1*U2(1);
En(2)=max(En(1)-B(1)+U1(1)+M1*U2(1),0);
Q(2)=(1-S(1)*C1(1)*U1(1))*Q(1)+U2(1)*C2(1)*A(1);
TH(2)=max(TH(1)+phi2-U1(1),0);

% Iteration of our algoritm
for t=2:j
    p=0.8;
    if S1(t)<=p
       S(t)=1;
    else
       S(t)=0;
    end
    
if (P(t)>=1)
    C1(t)=1;
else
    C1(t)=0;
end
    A(t)=(1-U2(t-1)*C2(t-1))*A(t-1)+K(t);
    I(t)=W(t)*(Q(t)^2+M2*A(t)^2);
    Ia(t)=((t-1)*Ia(t-1)+I(t))/t;
    E(t)=Q(t)^2+M2*A(t)^2;
    a(t)=V*En(t)-Z*TH(t)*S(t)*C1(t)-0.5*S(t)*C1(t)*theta*Q(t)^2-W(t+1)*R1*delta(t+1)*S(t)*C1(t)*Q(t)^2;
    



    if a(t)<0
       U1(t)=1;
    else
       U1(t)=0;
    end
    
if ((P(t)-U1(t))>=M1)
    C2(t)=1;
else
    C2(t)=0;
end
b(t)=V*M1*En(t)+0.5*theta*n(t)*C2(t)*sigma^2-0.5*beta*C2(t)*n(t)*sigma^2+(1-M2)*delta(t+1)*R1*W(t+1)*C2(t)*n(t)*sigma^2;
    if b(t)<0
       U2(t)=1;
    else
       U2(t)=0;
    end

    if U2(t)==1
        n(t+1)=1;
    else
        n(t+1)=n(t)+1;
    end
    P(t+1)=P(t)+B(t)-U1(t)-M1*U2(t);
    En(t+1)=max(En(t)-B(t)+U1(t)+M1*U2(t),0);
    Q(t+1)=(1-S(t)*C1(t)*U1(t))*Q(t)+U2(t)*C2(t)*A(t);
    TH(t+1)=max(TH(t)+phi2-U1(t),0);
end
Imax=max(I,[],2);
figure
plot(1:1:j,I);
xlim([0 10000]);
%title('Value of UOI(Lyapunov)')
xlabel('Time Slot'); 
ylabel('Value of UOI'); 
pause;
plot(1:1:j, Ia);
xlim([0 10000]);
%title('Average UOI(Lyapunov)')
xlabel('Time Slot'); 
ylabel('Average of UOI'); 
pause;
plot(1:1:j+1,En);
xlim([0 10000]);
%title('Length of Virtual Queue: Energy')
xlabel('Time Slot'); 
ylabel('Length of Virtual Queue: Energy E(t)'); 
pause;
plot(1:1:j+1,TH);
xlim([0 10000]);
%title('Length of Virtual Queue: Feedback')
xlabel('Time Slot'); 
ylabel('Length of Virtual Queue: Throughput TH_{t}'); 
pause;
yyaxis left
plot(1:1:j,E);
xlim([0 10000]);
%title('compare of error and weight of urgency \delta')
xlabel('Time Slot'); 
ylabel('{Q_{t}}^2+M*{A_{t}}^2');
yyaxis right
plot(1:1:j+1,delta);
xlim([0 10000]);
ylabel('Query Time delta');
legend({'Square Error','Query Time'},'Location','northeast');
pause;
yyaxis left
plot(1:1:j+1,En);
xlim([0 10000]);
%title('compare of energy queue and weight of urgency W')
xlabel('Time Slot'); 
ylabel('Energy Queue En(t)');
yyaxis right
plot(1:1:j+1,delta);
xlim([0 10000]);
ylabel('Query Time');
legend({'Energy Virtue Queue','Query Time'},'Location','northeast');
pause;
yyaxis left
plot(1:1:j+1,TH);
xlim([0 10000]);
%title('compare of throughput queue and weight of urgency W')
xlabel('Time Slot'); 
ylabel('Throughput Queue TH_{t}');
yyaxis right
plot(1:1:j+1,delta);
xlim([0 10000]);
ylabel('Query Time');
legend({'Throughput Virtue Queue','Query Time'},'Location','northeast');
pause;

