clc;
close all;
clear all;
x=100;
j=1000;
Tc=0.3;
Fc=0.3;
p=0.8;
Z=30;
V=30;
q=Tc*p;
Q=zeros(x,j);
n=zeros(x,j);
A=zeros(x,j);
Q1=zeros(x,j);
Q2=zeros(x,j);
A1=zeros(x,j);
A2=zeros(x,j);
G=zeros(x,j);
H=zeros(x,j);
Ia=zeros(x,1);
Wa=zeros(x,1);
I=zeros(x,j);
Ia1=zeros(x,1);
Ia2=zeros(x,1);
theta=zeros(x,1);
beta=zeros(x,1);
I1=zeros(x,j);
I2=zeros(x,j);
U1=zeros(x,j);
U2=zeros(x,j);
U21=zeros(x,j);
U22=zeros(x,j);
sigma=1;
K=normrnd(0,sigma,x,j);
col=j;
lin=x;
S=zeros(x,j);
S1=unifrnd(0,1,lin,col);
W=unifrnd(0,1,lin,col+1);
W1=unifrnd(0,1,lin,col+1);
for i=1:x
    for t=1:j+1
        if W1(i,t)<=0.01
            W(i,t)=1000;
        else
            W(i,t)=1;
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
A(i,1)=K(i,1);
n(i,1)=1;
theta(i,1)=2*Wa(i,1)*(1-p*Tc)/(p*Tc);
beta(i,1)=2*Wa(i,1)*(1-p*Tc)/(p*Tc)+2*Wa(i,1)/Fc;
a=V*H(i,1)-0.5*p*theta(i,1)*Q(i,1)^2-W(i,2)*p*Q(i,1)^2;
b=Z*G(i,1)+0.5*theta(i,1)*n(i,1)*sigma^2-0.5*beta(i,1)*n(i,1)*sigma^2;
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
I(i,1)=(Q(i,1)^2+A(i,1)^2);
Ia(i,1)=I(i,1);
    Q(i,2)=max((1-S(i,1)*U1(i,1))*Q(i,1)+U2(i,1)*A(i,1),0);
    H(i,2)=max(H(i,1)-Tc+U1(i,1),0);
    G(i,2)=max(G(i,1)-Fc+U2(i,1),0);
% Iteration of our algoritm
for t=2:j
    if S1(i,t)<=p
       S(i,t)=1;
    else
       S(i,t)=0;
    end 
    A(i,t)=(1-U2(i,t-1))*A(i,t-1)+K(i,t);
       
    I(i,t)=(Q(i,t)^2+A(i,t)^2);
    Ia(i,1)=((t-1)*Ia(i,1)+I(i,t))/t;
    
    a=V*H(i,t)-0.5*p*theta(i,1)*Q(i,t)^2-W(i,t+1)*p*Q(i,t)^2;
    b=Z*G(i,t)+0.5*theta(i,1)*n(i,t)*sigma^2-0.5*beta(i,1)*n(i,t)*sigma^2;
    if a>=0
        U1(i,t)=0;
    else
        U1(i,t)=1;
    end
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
    Q(i,t+1)=max((1-S(i,t)*U1(i,t))*Q(i,t)+U2(i,t)*A(i,t),0);
    H(i,t+1)=max(H(i,t)-Tc+U1(i,t),0);
    G(i,t+1)=max(G(i,t)-Fc+U2(i,t),0);
end
end

for i=1:x
A1(i,1)=K(i,1);
Q1(i,1)=0;
I1(i,1)=Q1(i,1)^2+A1(i,1)^2;
n1=0;
n2=0;
Ia1(i,1)=I1(i,1);

U21(i,1)=0;

    if A1(i,1)>0
            U22(i,1)=1;
    else
            U22(i,1)=0;
    end
Q1(i,2)=max((1-S(i,1)*U21(i,1))*Q1(i,1)+U22(i,1)*A1(i,1),0);
for t=2:j
    A1(i,t)=(1-U22(i,t-1))*A1(i,t-1)+K(i,t);
   
    
    
    I1(i,t)=Q1(i,t)^2+A1(i,t)^2;
    
    Ia1(i,1)=((t-1)*Ia1(i,1)+I1(i,t))/t;
   
    
        if (n2/t)<Fc
            U22(i,t)=1;
        else
            U22(i,t)=0;
        end
         if (n1/t)<q
            U21(i,t)=1;
        else
            U21(i,t)=0;
         end
         if U21(i,t)==1
             n1=n1+1;
         end 
         if U22(i,t)==1
             n2=n2+1;
         end 
         Q1(i,t+1)=max((1-S(i,t)*U21(i,t))*Q1(i,t)+U22(i,t)*A1(i,t),0);
end
end
figure
plot(1:1:100, Ia,1:1:100, Ia1);
%ylim([-100 2000]);
%title('average of error')
xlabel('Number of Rounds'); 
ylabel('Average of Error'); 
legend({'Lyapnov','greedy'},'Location','northeast');
pause;