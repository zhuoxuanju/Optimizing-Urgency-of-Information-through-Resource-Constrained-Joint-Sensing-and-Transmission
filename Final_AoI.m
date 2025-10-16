clc;
close all;
clear all;
%% setting part
x=100;
j=1500000;
Q=zeros(x,j);
A=zeros(x,j);
sigma=1;
K=normrnd(0,sigma,x,j);
G=zeros(x,j);
G1=zeros(x,j);
H=zeros(x,j);
H1=zeros(x,j);
I=zeros(x,j);
I1=zeros(x,j);
I2=zeros(x,j);
Q1=zeros(x,j);
Q2=zeros(x,j);
A1=zeros(x,j);
A2=zeros(x,j);
n=zeros(x,j);
U1=zeros(x,j);
U2=zeros(x,j);
U21=zeros(x,j);
U22=zeros(x,j);
Ia=zeros(x,1);
Ia1=zeros(x,1);
Ia2=zeros(x,1);
R=2;
M=2;
Tc=zeros(x,1);
theta=zeros(x,1);
beta=zeros(x,1);
Fc=0.5;
p=0.6;
Z=20;
V=20;
col=j;
lin=x;
W=unifrnd(0,1,lin,col+1);
W1=unifrnd(0,1,lin,col+1);
for i=1:x
    for t=1:j+1
        if W1(i,t)<=0.01
            W(i,t)=100;
        else
            W(i,t)=1;
        end
    end
end
Wa=mean(W,2);
S=unifrnd(0,1,lin,col);
S1=unifrnd(0,1,lin,col);

for i=1:x
   
if S1(i,1)<=p
   S(i,1)=1;
else
   S(i,1)=0;
end
A(i,1)=1;
if i==1
   Tc(i,1)=0.01;
else
   Tc(i,1)=Tc(i-1,1)+0.01;
end
n(i,1)=1;
theta(i,1)=2*R*Wa(i,1)*(1-p*Tc(i,1))/(p*Tc(i,1));
beta(i,1)=theta(i,1)+2*((1/Fc)-1)*R*M*Wa(i,1)+2*R*Wa(i,1);
a=V*H(i,1)-0.5*p*theta(i,1)*Q(i,1)^2-W(i,2)*R*p*Q(i,1)^2;
b=Z*G(i,1)+0.5*theta(i,1)*n(i,1)*sigma^2-0.5*beta(i,1)*n(i,1)*sigma^2+(1-M)*R*W(i,2)*n(i,1)*sigma^2;
a1=V*H1(i,1)-0.5*p*theta(i,1)*Q1(i,1)^2-R*W(i,2)*Q1(i,1);
b1=Z*G1(i,1)-0.5*beta(i,1)*A1(i,1)^2-beta(i,1)*A1(i,1)+0.5*theta(i,1)*A1(i,1)^2+theta(i,1)*Q1(i,1)*A1(i,1)+R*W(i,2)*A1(i,1)-R*M*W(i,2)*A1(i,1);
c1=-theta(i,1)*p*Q1(i,1)*A1(i,1);
if a>=0
     U1(i,1)=0;
else
     U1(i,1)=1;
end
if b>=0
    U2(i,1)=0;
    U22(i,1)=0;
else
    U2(i,1)=1;
    U22(i,1)=1;
end
if U2(i,1)==1
    n(i,2)=1;
else
    n(i,2)=2;
end
I(i,1)=W(i,1)*(Q(i,1)^2+M*A(i,1)^2);
I1(i,1)=W(i,1)*(Q1(i,1)^2+M*A1(i,1)^2);
I2(i,1)=W(i,1)*(Q1(i,1)+M*A1(i,1));
Ia(i,1)=I(i,1);
Ia1(i,1)=I1(i,1);
Ia2(i,1)=I2(i,1);
    %Q(i,2)=max((1-S(i,1)*U1(i,1))*Q(i,1)+U2(i,1)*A(i,1),0);
    Q(i,2)=(1-S(i,1)*U1(i,1))*Q(i,1)+U2(i,1)*A(i,1);
    H(i,2)=max(H(i,1)-Tc(i,1)+U1(i,1),0);
    G(i,2)=max(G(i,1)-Fc+U2(i,1),0);
    Q1(i,2)=(1-S(i,1)*U21(i,1))*Q1(i,1)+U22(i,1)*A1(i,1);
    H1(i,2)=max(H1(i,1)-Tc(i,1)+U21(i,1),0);
    G1(i,2)=max(G1(i,1)-Fc+U22(i,1),0);
% Iteration of our algoritm
for t=2:j
    if S1(i,t)<=p
       S(i,t)=1;
    else
       S(i,t)=0;
    end 
    A(i,t)=(1-U2(i,t-1))*A(i,t-1)+K(i,t);
       A1(i,t)=(1-U22(i,t-1))*A1(i,t-1)+1;
    I(i,t)=W(i,t)*(Q(i,t)^2+M*A(i,t)^2);
    Ia(i,1)=((t-1)*Ia(i,1)+I(i,t))/t;
    I1(i,t)=W(i,t)*(Q1(i,t)^2+M*A1(i,t)^2);
    Ia1(i,1)=((t-1)*Ia1(i,1)+I1(i,t))/t;
    I2(i,t)=W(i,t)*(Q1(i,t)+M*A1(i,t));
    Ia2(i,1)=((t-1)*Ia2(i,1)+I2(i,t))/t;
    a=V*H(i,t)-0.5*p*theta(i,1)*Q(i,t)^2-W(i,t+1)*R*p*Q(i,t)^2;
    b=Z*G(i,t)+0.5*theta(i,1)*n(i,t)*sigma^2-0.5*beta(i,1)*n(i,t)*sigma^2+(1-M)*R*W(i,t+1)*n(i,t)*sigma^2;
    a1=V*H1(i,t)-0.5*p*theta(i,1)*Q1(i,t)^2-R*W(i,t+1)*Q1(i,t);
    b1=Z*G1(i,t)-0.5*beta(i,1)*A1(i,t)^2-beta(i,1)*A1(i,t)+0.5*theta(i,1)*A1(i,t)^2+theta(i,1)*Q1(i,t)*A1(i,t)+R*W(i,t+1)*A1(i,t)-R*M*W(i,t+1)*A1(i,t);
    c1=-theta(i,1)*p*Q1(i,t)*A1(i,t);
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
    
    if (a1>=0)&&(b1>=0)&&(c1>=0)
        U21(i,t)=0;
        U22(i,t)=0;
    end
    if (a1>=0)&&(b1>=0)&&(c1<0)
        if (-c1<=a1+b1)
            U21(i,t)=0;
            U22(i,t)=0;
        end
        if(-c1>a1+b1)
            if Q1(i,t)~=0
            U21(i,t)=1;
            else
            U21(i,t)=0;
            end 
            U22(i,t)=1;
        end
    end
    if (a1>=0)&&(b1<0)&&(c1>=0)
        U21(i,t)=0;
        U22(i,t)=1;
    end
    if (a1>=0)&&(b1<0)&&(c1<0)
        if (-c1<=a1)
            U21(i,t)=0;
            U22(i,t)=1;
        end
        if (-c1>a1)
            if Q1(i,t)~=0
            U21(i,t)=1;
            else
            U21(i,t)=0;
            end 
            U22(i,t)=1;
            U22(i,t)=1;
        end             
    end
    if (a1<0)&&(b1>=0)&&(c1>=0)
        if Q1(i,t)~=0
            U21(i,t)=1;
            else
            U21(i,t)=0;
            end 
        U22(i,t)=0;
    end
    if (a1<0)&&(b1>=0)&&(c1<0)
        if (-c1<=b1)
            if Q1(i,t)~=0
            U21(i,t)=1;
            else
            U21(i,t)=0;
            end 
            U22(i,t)=0;
        end  
        if (-c1>b1)
            if Q1(i,t)~=0
            U21(i,t)=1;
            else
            U21(i,t)=0;
            end 
            U22(i,t)=1;
        end  
    end
    if (a1<0)&&(b1<0)&&(c1>=0)
       if (-a1>=-b1)
           if(-b1>=c1)
               if Q1(i,t)~=0
            U21(i,t)=1;
            else
            U21(i,t)=0;
            end 
               U22(i,t)=1;
           end
           if(-b1<c1)
               if Q1(i,t)~=0
                  U21(i,t)=1;
               else
                  U21(i,t)=0;
               end 
               U22(i,t)=0;
           end
       end
       if (-a1<-b1)
           if(-a1>=c1)
               if Q1(i,t)~=0
                  U21(i,t)=1;
               else
                  U21(i,t)=0;
               end 
               U22(i,t)=1;
           end
           if(-a1<c1)
               U21(i,t)=0;
               U22(i,t)=1;
           end
       end
    end
    if (a1<0)&&(b1<0)&&(c1<0)
        if Q1(i,t)~=0
            U21(i,t)=1;
            else
            U21(i,t)=0;
        end 
        U22(i,t)=1;
    end       
    
    
    %Q(i,t+1)=max((1-S(i,t)*U1(i,t))*Q(i,t)+U2(i,t)*A(i,t),0);
    Q(i,t+1)=(1-S(i,t)*U1(i,t))*Q(i,t)+U2(i,t)*A(i,t);
    H(i,t+1)=max(H(i,t)-Tc(i,1)+U1(i,t),0);
    G(i,t+1)=max(G(i,t)-Fc+U2(i,t),0);
    Q1(i,t+1)=(1-S(i,t)*U21(i,t))*Q1(i,t)+U22(i,t)*A1(i,t);
    H1(i,t+1)=max(H1(i,t)-Tc(i,1)+U21(i,t),0);
    G1(i,t+1)=max(G1(i,t)-Fc+U22(i,t),0);
end
end
    




figure


plot(0.01:0.01:1, Ia,0.01:0.01:1, Ia1,0.01:0.01:1, Ia2);
ylim([0 100]);
xlim([0.01 1]);
%ylim([-100 2000]);
%title('average of error')
xlabel('Frequency Constraint of Transmission'); 
ylabel('Average of UoI or Average Weighted Age Cost'); 
legend({'Average UoI of UoI Optimal','Average UoI of AoI Optimal', 'Average Weighted Age Cost of AoI Optimal'},'Location','northeast');
pause;