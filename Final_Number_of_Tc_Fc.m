clc;
close all;
clear all;
%% setting part
x=10;
e=10;
j=10000;
Q=zeros(x,j);
A=zeros(x,j);
n=zeros(x,j);
sigma=1;
K=normrnd(0,sigma,x,j);
G=zeros(x,j);
H=zeros(x,j);
I=zeros(x,j);
U1=zeros(1,j);
U2=zeros(1,j);
Ia=zeros(x,1);
B=zeros(x,j);
C1=zeros(x,j);
C2=zeros(x,j);
m=10;
R1=2;
M1=2;
M2=2;
En=zeros(x,j);
P=zeros(x,j);

TH=zeros(x,j);
delta=zeros(1,j+1);
NZT=zeros(10,x);

NZF=zeros(10,x);
phi2=zeros(x,1);
phi1=zeros(x,1);
p=1;
Z=1000;
V=1;
theta=zeros(x,1);
beta=zeros(x,1);


col=j;
lin=x;
W=unifrnd(0,1,lin,col+1);
W1=unifrnd(0,1,lin,col+1);
for i=1:x
    for t=1:j+1
%         if W1(i,t)<=0.01
%             W(i,t)=100;
%         else
            W(i,t)=1;
%         end
determine=mod(t,m);
    if (determine==0)
        delta(i,t)=1;
    else
        delta(i,t)=0;
    end
    end
end

Wa=mean(W,2);
S=unifrnd(0,1,lin,col);
S1=unifrnd(0,1,lin,col);
S2=unifrnd(0,1,lin,col);
for d=1:e
   if d==1
     phi2(d,1)=0.1;
   else
     phi2(d,1)=phi2(d-1,1)+0.1;
   end

for i=1:x
   
if S1(i,1)<=p
   S(i,1)=1;
else
   S(i,1)=0;
end
if (S2(i,1)>phi1(i,1))
        B(i,1)=0;
    else
        B(i,1)=1;
end
if (P(i,1)>=1)
    C1(i,1)=1;
else
    C1(i,1)=0;
end
if ((P(i,1)-U1(1))>=M1)
    C2(i,1)=1;
else
    C2(i,1)=0;
end
A(i,1)=K(i,1);
if i==1
   phi1(i,1)=0.1;
else
   phi1(i,1)=phi1(i-1,1)+0.1;
end
theta(i,1)=(2*Wa(i,1)*R1/m)*(((1+M1)/(M1*p*phi1(i,1)))-1);
beta(i,1)=theta(i,1)+(2*Wa(i,1)*R1*M2*(1+M1)/(m*phi1(i,1)))+2*Wa(i,1)*R1*(1-M2)/m;
a=V*En(i,1)-Z*TH(i,1)*S(i,1)*C1(i,1)-0.5*S(i,1)*C1(i,1)*theta(i,1)*Q(i,1)^2-W(i,2)*R1*delta(i,2)*S(i,1)*C1(i,1)*Q(i,1)^2;
b=V*M1*En(i,1)+0.5*theta(i,1)*n(i,1)*C2(i,1)*sigma^2-0.5*beta(i,1)*C2(i,1)*n(i,1)*sigma^2+(1-M2)*delta(i,2)*R1*W(i,2)*n(i,1)*sigma^2;

if a>=0
     U1(1)=0;
else
     U1(1)=1;
end
if b>=0
    U2(1)=0;
else
    U2(1)=1;
end
if U2(1)==1
    n(i,2)=1;
else
    n(i,2)=2;
end
I(i,1)=delta(i,1)*W(i,1)*(Q(i,1)^2+M2*A(i,1)^2);
Ia(i,1)=I(i,1);
P(i,2)=P(i,1)+B(i,1)-U1(1)-M1*U2(1);
En(i,2)=max(En(i,1)-B(i,1)+U1(1)+M1*U2(1),0);
Q(i,2)=(1-S(i,1)*C1(i,1)*U1(1))*Q(i,1)+U2(1)*C2(i,1)*A(i,1);
TH(i,2)=max(TH(i,1)+phi2(d,1)-U1(1),0);
% Iteration of our algoritm
for t=2:j
    if S1(i,t)<=p
       S(i,t)=1;
    else
       S(i,t)=0;
    end 
    
    if (S2(i,t)>phi1(i,1))
        B(i,t)=0;
    else
        B(i,t)=1;
    end
if (P(i,t)>=1)
    C1(i,t)=1;
else
    C1(i,t)=0;
end
    A(i,t)=(1-U2(t-1)*C2(i,t-1))*A(i,t-1)+K(i,t);
       
    I(i,t)=delta(i,t)*W(i,t)*(Q(i,t)^2+M2*A(i,t)^2);
    Ia(i,1)=((t-1)*Ia(i,1)+I(i,t))/t;
   a=V*En(i,t)-Z*TH(i,t)*S(i,t)*C1(i,t)-0.5*S(i,t)*C1(i,t)*theta(i,1)*Q(i,t)^2-W(i,t+1)*R1*delta(i,t+1)*S(i,t)*C1(i,t)*Q(i,t)^2;
     if a>=0
        U1(t)=0;
    else
        U1(t)=1;
    end
    
if ((P(i,t)-U1(t))>=M1)
    C2(i,t)=1;
else
    C2(i,t)=0;
end
    
    b=V*M1*En(i,t)+0.5*theta(i,1)*n(i,t)*C2(i,t)*sigma^2-0.5*beta(i,1)*C2(i,t)*n(i,t)*sigma^2+(1-M2)*delta(i,t+1)*R1*W(i,t+1)*n(i,t)*sigma^2;
    if b>=0
        U2(t)=0;
    else
        U2(t)=1;
    end 
    if U2(t)==1
        n(i,t+1)=1;
    else
        n(i,t+1)=n(i,t)+1;
    end
   P(i,t+1)=P(i,t)+B(i,t)-U1(t)-M1*U2(t);
En(i,t+1)=max(En(i,t)-B(i,t)+U1(t)+M1*U2(t),0);
Q(i,t+1)=(1-S(i,t)*C1(i,t)*U1(t))*Q(i,t)+U2(t)*C2(i,t)*A(i,t);
TH(i,t+1)=max(TH(i,t)+phi2(d,1)-U1(t),0);
end
NZT(d,i)=sum(U1~=0);
NZF(i,d)=sum(U2~=0);

end
end
%Imax=max(I,[],2);
figure
plot(phi1, NZT);
ylim([0 10000]);
xlim([0.1 1]);
%title('number of transmission when v=1')
xlabel('Energy Constraint '); 
ylabel('Number of Transimission'); 
legend({'\phi_{2}=0.1','\phi_{2}=0.2','\phi_{2}=0.3','\phi_{2}=0.4','\phi_{2}=0.5','\phi_{2}=0.6','\phi_{2}=0.7','\phi_{2}=0.8','\phi_{2}=0.9','\phi_{2}=1.0',},'Location','northwest');
pause;
plot(phi1, NZF);
ylim([0 10000]);
xlim([0.1 1]);
%title('number of sense when z=1')
xlabel('Throughput Constraint'); 
ylabel('Number of Sense'); 
legend({'\phi_{2}=0.1','\phi_{2}=0.2','\phi_{2}=0.3','\phi_{2}=0.4','\phi_{2}=0.5','\phi_{2}=0.6','\phi_{2}=0.7','\phi_{2}=0.8','\phi_{2}=0.9','\phi_{2}=1.0',},'Location','northwest');