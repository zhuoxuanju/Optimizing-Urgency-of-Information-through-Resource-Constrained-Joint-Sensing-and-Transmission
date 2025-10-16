clc;
close all;
clear all;
%% setting part
x=100;
j=1000000;
Q=zeros(x,j);
A=zeros(x,j);
Q1=zeros(x,j);
A1=zeros(x,j);
Q2=zeros(x,j);
A2=zeros(x,j);
sigma=1;
K=normrnd(0,sigma,x,j);

I=zeros(x,j);
I1=zeros(x,j);
I2=zeros(x,j);
n=zeros(x,j);
U1=zeros(x,j);
U2=zeros(x,j);
U11=zeros(x,j);
U21=zeros(x,j);
U12=zeros(x,j);
U22=zeros(x,j);
Ia=zeros(x,1);
Ia1=zeros(x,1);
Ia2=zeros(x,1);

theta=zeros(x,1);
beta=zeros(x,1);
phi1=zeros(x,1);
pr1=zeros(x,1);
pr2=zeros(x,1);
B=zeros(x,j);
C1=zeros(x,j);
C2=zeros(x,j);
C11=zeros(x,j);
C21=zeros(x,j);
a=zeros(x,j);
b=zeros(x,j);
En=zeros(x,j);
P=zeros(x,j);
En1=zeros(x,j);
P1=zeros(x,j);
TH=zeros(x,j);
TH1=zeros(x,j);
delta=zeros(1,j+1);
p=0.8;
R1=2;
M1=2;
M2=2;
phi2=0.5;
m=10;
V=1;
Z=1;
col=j;
lin=x;
W=unifrnd(0,1,lin,col+1);
S=unifrnd(0,1,lin,col);
S1=unifrnd(0,1,lin,col);
S2=unifrnd(0,1,lin,col);
S3=unifrnd(0,1,lin,col);
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
   pr1(i,1)=M1*phi1(i,1)/(1+M1);
  pr2(i,1)=phi1(i,1)/(1+M1);
else
   phi1(i,1)=phi1(i-1,1)+0.01;
   pr1(i,1)=M1*phi1(i,1)/(1+M1);
  pr2(i,1)=phi1(i,1)/(1+M1);
end

n(i,1)=1;
theta(i,1)=(2*Wa(i,1)*R1/m)*(((1+M1)/(M1*p*phi1(i,1)))-1);
beta(i,1)=theta(i,1)+(2*Wa(i,1)*R1*M2*(1+M1)/(m*phi1(i,1)))+2*Wa(i,1)*R1*(1-M2)/m;
a(i,1)=V*En(i,1)-Z*TH(i,1)*S(i,1)*C1(i,1)-0.5*S(i,1)*C1(i,1)*theta(i,1)*Q(i,1)^2-W(i,2)*R1*delta(i,2)*S(i,1)*C1(i,1)*Q(i,1)^2;
b(i,1)=V*M1*En(i,1)+0.5*theta(i,1)*n(i,1)*C2(i,1)*sigma^2-0.5*beta(i,1)*C2(i,1)*n(i,1)*sigma^2+(1-M2)*delta(i,2)*R1*W(i,2)*n(i,1)*sigma^2;
if a(i,1)>=0
     U1(i,1)=0;
else
     U1(i,1)=1;
end
if b(i,1)>=0
    U2(i,1)=0;
else
    U2(i,1)=1;
end
if U2(i,1)==1
    n(i,2)=1;
else
    n(i,2)=2;
end
I(i,1)=delta(i,1)*W(i,1)*(Q(i,1)^2+M2*A(i,1)^2);
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
       
    I(i,t)=delta(i,t)*W(i,t)*(Q(i,t)^2+M2*A(i,t)^2);
    Ia(i,1)=((t-1)*Ia(i,1)+I(i,t))/t;
   a(i,t)=V*En(i,t)-Z*TH(i,t)*S(i,t)*C1(i,t)-0.5*S(i,t)*C1(i,t)*theta(i,1)*Q(i,t)^2-W(i,t+1)*R1*delta(i,t+1)*S(i,t)*C1(i,t)*Q(i,t)^2;
    if a(i,t)>=0
        U1(i,t)=0;
    else
        U1(i,t)=1;
    end
if ((P(i,t)-U1(i,t))>=M1)
    C2(i,t)=1;
else
    C2(i,t)=0;
end
    
    b(i,t)=V*M1*En(i,t)+0.5*theta(i,1)*n(i,t)*C2(i,t)*sigma^2-0.5*beta(i,1)*C2(i,t)*n(i,t)*sigma^2+(1-M2)*delta(i,t+1)*R1*W(i,t+1)*n(i,t)*sigma^2;
    
    if b(i,t)>=0
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
    
for i=1:x
A1(i,1)=K(i,1);
Q1(i,1)=0;
I1(i,1)=delta(i,1)*(Q1(i,1)^2+M2*A1(i,1)^2);
A2(i,1)=K(i,1);
Q2(i,1)=0;
I2(i,1)=delta(i,1)*(Q2(i,1)^2+M2*A2(i,1)^2);
n1=0;
n2=0;
Ia1(i,1)=I1(i,1);
U11(i,1)=0;
U21(i,1)=0;
U12(i,1)=0;
U22(i,1)=0;

    
%Q1(i,2)=max((1-S(i,1)*U21(i,1))*Q1(i,1)+U22(i,1)*A1(i,1),0);
Q1(i,2)=(1-S(i,1)*U11(i,1))*Q1(i,1)+U21(i,1)*A1(i,1);
Q2(i,2)=(1-S(i,1)*U12(i,1))*Q2(i,1)+U22(i,1)*A2(i,1);

for t=2:j
    A1(i,t)=(1-U21(i,t-1))*A1(i,t-1)+K(i,t);
   A2(i,t)=(1-U22(i,t-1))*A2(i,t-1)+K(i,t);
    
    
    I1(i,t)=delta(i,t)*(Q1(i,t)^2+M2*A1(i,t)^2);
    I2(i,t)=delta(i,t)*(Q2(i,t)^2+M2*A2(i,t)^2);
    Ia1(i,1)=((t-1)*Ia1(i,1)+I1(i,t))/t;
   Ia2(i,1)=((t-1)*Ia2(i,1)+I2(i,t))/t;
    
        if (S2(i,t)>=pr1(i,1))
            U12(i,t)=0;
        else
            U12(i,t)=1;
        end
         if (S3(i,t)>=pr2(i,1))
             U22(i,t)=0;
             
         else
            
            U22(i,t)=1;
         end
        
if (P1(i,t)>=1)&&(n1<M1)
    U11(i,t)=1;
    n1=n1+1;
else
    U11(i,t)=0;
end
if (P1(i,t)-U11(i,t)>=M1)&&(n1==3)
    U21(i,t)=1;
    n1=0;
else
    U21(i,t)=0;
end
P1(i,t+1)=P1(i,t)+B(i,t)-U11(i,t)-M1*U21(i,t);
         Q1(i,t+1)=(1-S(i,t)*U11(i,t))*Q1(i,t)+U21(i,t)*A1(i,t);
         Q2(i,t+1)=(1-S(i,t)*U12(i,t))*Q2(i,t)+U22(i,t)*A2(i,t);
        %Q1(i,t+1)=max((1-S(i,t)*U21(i,t))*Q1(i,t)+U22(i,t)*A1(i,t),0);
end
end



figure


plot(0.01:0.01:1, Ia,0.01:0.01:1, Ia1);
%ylim([0 200]);
xlim([0.1 1]);
%ylim([-100 2000]);
%title('average of error')
xlabel('Frequency Constraint of Transmission'); 
ylabel('{Q_{t}}^2+M2*{A_{t}}^2'); 
legend({'Lyapunov','Round Robin'},'Location','northeast');
pause;