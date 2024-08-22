function dv=HighPackingModel_GA(t,a,OP)
global R b br R0 ka k1 k2 k3 k4 k5 k6 k7 k8 km1 km2  km3 km4 km5 km6 kts ktf num m1 m2 
dv = zeros(num*7,1);
i = 0;
%% before drurg
if t<0
    dv = zeros(num*7,1);
else
%% after drug
while i < num 
P = OP(i+1,:);
R0=P(1);
R=P(2);
ka=P(3);
ktf=P(4);
kts=P(5);
m1=P(6);
m2=P(7);    
k1 =P(8);
k2 =P(9);
k3 =P(10);
k4 =P(11);     
k5 =P(12);
k6 =P(13);
k7 =P(14);
k8 =P(15);
km1 =P(16);
km2 =P(17);
km3 =P(18);
km4 =P(19);
km5 =P(20);
km6 =P(21);
dv(1+7*i) = m1*ka*(R0+R) - k1*(a(1+7*i)/(km1+a(1+7*i)))*a(5+7*i) - k2*(a(1+7*i)/(km2+a(1+7*i)))*a(6+7*i) -ktf*(a(1+7*i)-a(2+7*i));
dv(2+7*i) = ktf*(a(1+7*i) - a(2+7*i));
dv(3+7*i) = m2*ka*(R0+R) - k1*(a(3+7*i)/(km1+a(3+7*i)))*a(5+7*i) - k2*(a(3+7*i)/(km2+a(3+7*i)))*a(6+7*i) -kts*(a(3+7*i)-a(4+7*i));
dv(4+7*i) = kts*(a(3+7*i) - a(4+7*i));
bg(i+1)=a(1+7*i)+a(3+7*i);
dv(5+7*i) = k3*bg(i+1) - k4*(a(5+7*i)/(a(5+7*i) + km3));
cadiffout =0;
cadiffin = 0;
cadiffout1 =0;
cadiffin1 = 0;
j = 0;
%calculating term for total outward and inward diffusion of Cacyt 
while j < num
    cadiffout1=cadiffout1-br(i+1,j+1)*(a(6+7*i)-a(6+7*j));
    cadiffin1=cadiffin1  +br(j+1,i+1)*(a(6+7*j)-a(6+7*i));
    if a(6+7*i)>a(6+7*j)
            cadiffout=cadiffout-b(i+1,j+1)*(a(6+7*i)-a(6+7*j));  
    elseif a(6+7*j)>a(6+7*i)
            cadiffin=cadiffin  +b(j+1,i+1)*(a(6+7*j)-a(6+7*i));
     end
     j=j+1;
end
dv(6+7*i) = k5*a(5+7*i)+k6*a(5+7*i)*a(6+7*i)*(a(7+7*i)/(a(7+7*i)+km4))-k7*(a(6+7*i))/(a(6+7*i)+km5)-k8*(a(6+7*i)/(a(6+7*i)+km6))...
    + cadiffin+ cadiffout+ cadiffin1+ cadiffout1 ;
dv(7+7*i) = -k6*a(5+7*i)*a(6+7*i)*(a(7+7*i)/(a(7+7*i)+km4)) + k8*(a(6+7*i)/(a(6+7*i)+km6));
i= i+1;
end
end
% t;
end