function dv=Single_Cell_Model_GA(t,a,P)
% global R b br R0 ka k1 k2 k3 k4 k5 k6 k7 k8 km1 km2  km3 km4 km5 km6 kts ktf num m1 m2 
dv=zeros(7,1);
R0= P(1);
R=P(2);
ka=P(3);
ktf= P(4);
kts=P(5);
m1=P(6);
m2= P(7);
k1 =P(8);
k2 =P(9);
k3 =P(10);
k4 = P(11);
k5 =P(12);
k6 =P(13);
k7 = P(14);
k8 =P(15);
km1 =P(16);
km2 =P(17);
km3 =P(18);
km4 =P(19);
km5 =P(20);
km6 =P(21);

if t<0
    dv = zeros(7,1);
else
    dv(1) = m1*ka*(R0+R) - k1*(a(1)/(km1+a(1)))*a(5) - k2*(a(1)/(km2+a(1)))*a(6) -ktf*(a(1)-a(2));
    dv(2) = ktf*(a(1) - a(2));
    dv(3) = m2*ka*(R0+R) - k1*(a(3)/(km1+a(3)))*a(5) - k2*(a(3)/(km2+a(3)))*a(6) -kts*(a(3)-a(4));
    dv(4) = kts*(a(3) - a(4));
    bg=a(1)+a(3);
    dv(5) = k3*bg - k4*(a(5)/(a(5) + km3));
    dv(6) = k5*a(5)+k6*a(5)*a(6)*(a(7)/(a(7)+km4))-k7*(a(6))/(a(6)+km5)-k8*(a(6)/(a(6)+km6));
    dv(7) = -k6*a(5)*a(6)*(a(7)/(a(7)+km4)) + k8*(a(6)/(a(6)+km6));
end
end