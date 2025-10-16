clc;
close all;
clear all;
%% setting part
x=100;
j=1000000;
Q=zeros(x,j);
A=zeros(x,j);
sigma=1;
K=normrnd(0,sigma,x,j);
G=zeros(x,j);
H=zeros(x,j);
I=zeros(x,j);
n=zeros(x,j);
U1=zeros(x,j);
U2=zeros(x,j);
Ia=zeros(x,1);
Ia1=zeros(x,1);
Ia2=zeros(x,1);
Ia3=zeros(x,1);
Ia4=zeros(x,1);
Ia5=zeros(x,1);
Tc=zeros(x,1);
theta=zeros(x,1);
beta=zeros(x,1);
phi1=zeros(x,1);
B=zeros(x,j);
C1=zeros(x,j);
C2=zeros(x,j);
En=zeros(x,j);
P=zeros(x,j);

TH=zeros(x,j);
delta=zeros(1,j+1);
R1=2;
M1=2;
M2=2;
p=0.8;
m=10;
V=1;
Z=1;
col=j;
lin=x;
W=unifrnd(0,1,lin,col+1);
S=unifrnd(0,1,lin,col);
S1=unifrnd(0,1,lin,col);
S2=unifrnd(0,1,lin,col);
% W1=unifrnd(0,1,lin,col+1);
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

for d=1:5
    if d==1
        phi2=0.2;
     elseif d==2
        phi2=0.4;
     elseif d==3
        phi2=0.6;
     elseif d==4
        phi2=0.8;
     elseif d==5
        phi2=1.0;
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
if ((P(i,1)-U1(i,1))>=M1)
    C2(i,1)=1;
else
    C2(i,1)=0;
end
A(i,1)=K(i,1);
if i==1
   phi1(i,1)=0.01;
else
   phi1(i,1)=phi1(i-1,1)+0.01;
end
n(i,1)=1;
theta(i,1)=(2*Wa(i,1)*R1/m)*(((1+M1)/(M1*p*phi1(i,1)))-1);
beta(i,1)=theta(i,1)+(2*Wa(i,1)*R1*M2*(1+M1)/(m*phi1(i,1)))+2*Wa(i,1)*R1*(1-M2)/m;
a=V*En(i,1)-Z*TH(i,1)*S(i,1)*C1(i,1)-0.5*S(i,1)*C1(i,1)*theta(i,1)*Q(i,1)^2-W(i,2)*R1*delta(i,2)*S(i,1)*C1(i,1)*Q(i,1)^2;
b=V*M1*En(i,1)+0.5*theta(i,1)*n(i,1)*C2(i,1)*sigma^2-0.5*beta(i,1)*C2(i,1)*n(i,1)*sigma^2+(1-M2)*delta(i,2)*R1*W(i,2)*n(i,1)*sigma^2;
if a>=0
     U1(i,1)=0;
else
     U1(i,1)=1;
end
if b>=0
    U2(i,1)=0;
else
    U2(i,1)=1;
end
if U2(i,1)==1
    n(i,2)=1;
else
    n(i,2)=2;
end
I(i,1)=W(i,1)*(Q(i,1)^2+M2*A(i,1)^2);
Ia(i,1)=I(i,1);
P(i,2)=P(i,1)+B(i,1)-U1(i,1)-M1*U2(i,1);
En(i,2)=max(En(i,1)-B(i,1)+U1(i,1)+M1*U2(i,1),0);
Q(i,2)=(1-S(i,1)*C1(i,1)*U1(i,1))*Q(i,1)+U2(i,1)*C2(i,1)*A(i,1);
TH(i,2)=max(TH(i,1)+phi2-U1(i,1),0);
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
    A(i,t)=(1-U2(i,t-1)*C2(i,t-1))*A(i,t-1)+K(i,t);
       
    I(i,t)=W(i,t)*(Q(i,t)^2+M2*A(i,t)^2);
    Ia(i,1)=((t-1)*Ia(i,1)+I(i,t))/t;
   a=V*En(i,t)-Z*TH(i,t)*S(i,t)*C1(i,t)-0.5*S(i,t)*C1(i,t)*theta(i,1)*Q(i,t)^2-W(i,t+1)*R1*delta(i,t+1)*S(i,t)*C1(i,t)*Q(i,t)^2;
    if a>=0
        U1(i,t)=0;
    else
        U1(i,t)=1;
    end
if ((P(i,t)-U1(i,t))>=M1)
    C2(i,t)=1;
else
    C2(i,t)=0;
end
    
    b=V*M1*En(i,t)+0.5*theta(i,1)*n(i,t)*C2(i,t)*sigma^2-0.5*beta(i,1)*C2(i,t)*n(i,t)*sigma^2+(1-M2)*delta(i,t+1)*R1*W(i,t+1)*n(i,t)*sigma^2;
    
    if b>=0
        U2(i,t)=0;
    else
        U2(i,t)=1;
    end 
    
    if U2(i,t)==1
        n(i,t+1)=1;
    else
        n(i,t+1)=n(i,t)+1;
    end
    P(i,t+1)=P(i,t)+B(i,t)-U1(i,t)-M1*U2(i,t);
En(i,t+1)=max(En(i,t)-B(i,t)+U1(i,t)+M1*U2(i,t),0);
Q(i,t+1)=(1-S(i,t)*C1(i,t)*U1(i,t))*Q(i,t)+U2(i,t)*C2(i,t)*A(i,t);
TH(i,t+1)=max(TH(i,t)+phi2-U1(i,t),0);
end
end
    if d==1
        Ia1=Ia;
    elseif d==2
       Ia2=Ia;
    elseif d==3
        Ia3=Ia;
    elseif d==4
       Ia4=Ia;
    elseif d==5
        Ia5=Ia;
    end
end

figure
plot(0.01:0.01:1, Ia1,0.01:0.01:1, Ia2,0.01:0.01:1, Ia3,0.01:0.01:1, Ia4, 0.01:0.01:1, Ia5);
ylim([1 200]);
xlim([0 1]);
%title('average of UOI(Lyapunov)')
xlabel('Energy constraint'); 
ylabel('Average UoI'); 
legend({'\phi_{2}=0.2','\phi_{2}=0.4','\phi_{2}=0.6','\phi_{2}=0.8','\phi_{2}=1.0',},'Location','northeast');

pause;