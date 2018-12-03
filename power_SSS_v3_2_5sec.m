function SSS = power_SSS_v3_2_5sec(x)

syms q1 q2 q3 q4 q5 real
q = [q1 q2 q3 q4 q5];
l2 = 8;l3 = 10;l4 = 4;l5 = 4;
a2 = 15;a3 = 15;d1 =3;d2 = 4;d5 = 8;
A01 = [cos(q1),0,sin(q1),0;
       sin(q1),0,-cos(q1),0;
       0,1,0,d1;
       0,0,0,1];
P0 = [0;0;0];
Z0 = [0;0;1];
P1 = A01(1:3,4);
Z1 = A01(1:3,3);

A12 = [cos(q2),-sin(q2),0,a2*cos(q2);
       sin(q2),cos(q2),0,a2*sin(q2);
       0,0,1,d2;
       0,0,0,1];
 
A02 = A01*A12;
P2 = A02(1:3,4);
Z2 = A02(1:3,3);

A23 = [cos(q3),-sin(q3),0,a3*cos(q3);
       sin(q3),cos(q3),0,a3*sin(q3);
       0,0,1,0;
       0,0,0,1];
A03 = A02*A23;
P3 = A03(1:3,4);
Z3 = A03(1:3,3);

A34 = [cos(q4),0,sin(q4),0;
       sin(q4),0,-cos(q4),0;
       0,1,0,0;
       0,0,0,1];
A04 = A03*A34;
P4 = A04(1:3,4);
Z4 = A04(1:3,3);

A45 = [cos(q5),-sin(q5),0,0;
       sin(q5),cos(q5),0,0;
       0,0,1,d5;
       0,0,0,1];
A05 = A04*A45;
P5 = A05(1:3,4);
Z5 = A05(1:3,3);
 J = [cross(Z0,(P5-P0)),cross(Z1,(P5-P1)),cross(Z2,(P5-P2)),cross(Z3,(P5-P3)),cross(Z4,(P5-P4));
      Z0,               Z1,               Z2,               Z3,               Z4];

% simplify(J)
%% Jacobians for torque calculations(dynamics)
%%%%%%%%%%%&&&&&&&&&&&&&&&&&
Ic1 = [0.0011303*10^7,0,0;0,0.0011303*10^7,0;0,0,0.0016747*10^7]; %(g.cm^2)
Ic2 = [0.00048527*10^7,0,0.00014106*10^7;0,0.0011004*10^7,0;0,0,0.0013838*10^7];
Ic3 = [973,-27,8.27;-27,2911,0;8.27,0,3673];
Ic4 = [178.6,0,35.8;0,283.5,0;0,0,238.6];
Ic5 = [159.6,0,33.8;0,211.5,0;0,0,144.6];
Pc1 = P1;
Pc2 = P2 + [l2*cos(q2);l2*sin(q2);0];
Pc3 = P3 + [l3*cos(q3);l3*sin(q3);0];
Pc4 = P4 + [l4*cos(q4);l4*sin(q4);0];
Pc5 = P5 + [l5*cos(q5);l5*sin(q5);0];

Pc = [Pc1,Pc2,Pc3,Pc4,Pc5];
% Jv = [Jv1,Jv2,Jv3,Jv4,Jv5];
Jv1 = [diff(Pc1,q1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)];
Jv2 = [diff(Pc2,q1),diff(Pc2,q2),zeros(3,1),zeros(3,1),zeros(3,1)];
Jv3 = [diff(Pc3,q1),diff(Pc3,q2),diff(Pc3,q3),zeros(3,1),zeros(3,1)];
Jv4 = [diff(Pc4,q1),diff(Pc4,q2),diff(Pc4,q3),diff(Pc4,q4),zeros(3,1)];
Jv5 = [diff(Pc5,q1),diff(Pc5,q2),diff(Pc5,q3),diff(Pc5,q4),diff(Pc5,q5)];
m1 = 170;m2 = 250;m3 = 100;m4 = 50;m5 = 70;
Jw1 = [Z1,zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)];
Jw2 = [Z1,Z2,zeros(3,1),zeros(3,1),zeros(3,1)];
Jw3 = [Z1,Z2,Z3,zeros(3,1),zeros(3,1)];
Jw4 = [Z1,Z2,Z3,Z4,zeros(3,1)];
Jw5 = [Z1,Z2,Z3,Z4,zeros(3,1)];
M = m1*Jv1'*Jv1 + Jw1'*Ic1*Jw1 + m2*Jv2'*Jv2 + Jw2'*Ic2*Jw2 + m3*Jv3'*Jv3 + Jw3'*Ic3*Jw3 +m4*Jv4'*Jv4 + Jw4'*Ic4*Jw4 +m5*Jv5'*Jv5 + Jw5'*Ic5*Jw5;
size(M);
% return;
%% For centrifugal calculation 
 W = sym(zeros(5,5));
b = sym(zeros(5,5));
% C = zeros(5,5);
k = 1;
for i = 1:5
    for j = 1:5
        k == j;
        W(i,j) = diff(M(i,j),q(k));
        W(i,k) = diff(M(i,k),q(j));
        W(j,k) = diff(M(j,k),q(i));
        b(i,j) = 0.5*(W(i,j)+W(i,k)-W(j,k));
    end
end
    


for i = 1:5
    for j = 1:5
        k == j;
    C(i,j)= b(i,j);
    end
end

%% for coriollis calculation

b1 = sym(zeros(5,5,5));
W1 = sym(zeros(5,5,5));
W2 = sym(zeros(5,5,5));
W3 = sym(zeros(5,5,5));
for i = 1:5
    for j = 1:4
       for k = j+1:5;
         W1(i,j,k) = diff(M(i,j),q(k));
        W2(i,k,j) = diff(M(i,k),q(j));
        W3(j,k,i) = diff(M(j,k),q(i));
        b1(i,j,k) = 0.5*(W1(i,j,k)+W2(i,k,j)-W3(j,k,i));   
       end
    end
end
% global B;
B = sym(zeros(5,10));

for i =1:5
B(i,1) = 2*b1(i,1,2);
B(i,2) = 2*b1(i,1,3);
B(i,3) = 2*b1(i,1,4);
B(i,4) = 2*b1(i,1,5);
B(i,5) = 2*b1(i,2,3);
B(i,6) = 2*b1(i,2,4);
B(i,7) = 2*b1(i,2,5);
B(i,8) = 2*b1(i,3,4);
B(i,9) = 2*b1(i,3,5);
B(i,10) = 2*b1(i,4,5);
end

thetaf = [1.707 , 1.2217 , 1.2217 , 1.2217, 1.2217]';
thetat=[0.7071 1.7071 1.7071 1.7071 1.7071]';
Asss = [2.5^5 2.5^4 2.5^3;
    5*2.5^4   4*2.5^3 3*2.5^2;
    20*2.5^3 12*2.5^2 6*2.5];
xsss1 = [ thetaf(1,1)-thetat(1,1) - x(1)*2.5^7 - x(2)*2.5^6;
    -7*x(1)*2.5^6 - 6*x(2)*2.5^5;
    42*x(1)*2.5^5 - 30*x(2)*2.5^4];
xsss2 = [ thetaf(2,1)-thetat(2,1) - x(3)*2.5^7 - x(4)*2.5^6;
    -7*x(3)*2.5^6 - 6*x(4)*2.5^5;
    42*x(3)*2.5^5 - 30*x(4)*2.5^4];
xsss3= [ thetaf(3,1)-thetat(3,1) - x(5)*2.5^7 - x(6)*2.5^6;
    -7*x(5)*2.5^6 - 6*x(6)*2.5^5;
    42*x(5)*2.5^5 - 30*x(6)*2.5^4];

xsss4 = [ thetaf(4,1)-thetat(4,1) - x(7)*2.5^7 - x(8)*2.5^6;
    -7*x(7)*2.5^6 - 6*x(8)*2.5^5;
    42*x(7)*2.5^5 - 30*x(8)*2.5^4];
xsss5 = [ thetaf(5,1)-thetat(5,1) - x(9)*2.5^7 - x(10)*2.5^6;
    -7*x(9)*2.5^6 - 6*x(10)*2.5^5;
    42*x(9)*2.5^5 - 30*x(10)*2.5^4];
%Bzsss = [csss;dsss;esss];
Zsss1 = Asss\xsss1;
Zsss2 = Asss\xsss2;
Zsss3 = Asss\xsss3;
Zsss4 = Asss\xsss4;
Zsss5 = Asss\xsss5;




%%%Joint 1%%
c1 = Zsss1(1,1);
d1= Zsss1(2,1);
e1 = Zsss1(3,1);
%%%Joint 2%%
c2 = Zsss2(1,1);
d2= Zsss2(2,1);
e2 = Zsss2(3,1);

%%%Joint 3%%
c3 = Zsss3(1,1);
d3= Zsss3(2,1);
e3 = Zsss3(3,1);

%%%Joint 4%%
c4 = Zsss4(1,1);
d4= Zsss4(2,1);
e4 = Zsss4(3,1);

%%%Joint 5%%
c5 = Zsss5(1,1);
d5= Zsss5(2,1);
e5 = Zsss5(3,1);

thetat=[0.7071 1.7071 1.7071 1.7071 1.7071];

%% joint Angle%%
t = [0:0.1:2.5];
Q1=[x(1)*t.^7+x(2)*t.^6+c1*t.^5+d1*t.^4+e1*t.^3+thetat(1,1)];
Q2=[x(3)*t.^7+x(4)*t.^6+c2*t.^5+d2*t.^4+e2*t.^3+thetat(1,2)];
Q3=[x(5)*t.^7+x(6)*t.^6+c3*t.^5+d3*t.^4+e3*t.^3+thetat(1,3)];
Q4=[x(7)*t.^7+x(8)*t.^6+c4*t.^5+d4*t.^4+e4*t.^3+thetat(1,4)];
Q5=[x(9)*t.^7+x(10)*t.^6+c5*t.^5+d5*t.^4+e5*t.^3+thetat(1,5)];
% size(Q1)
%%jotnt angle veloctty%%
Qd1=[7*x(1)*t.^6+6*x(2)*t.^5+5*c1*t.^4+4*d1*t.^3+3*e1*t.^2];
Qd2=[7*x(3)*t.^6+6*x(4)*t.^5+5*c2*t.^4+4*d2*t.^3+3*e2*t.^2];
Qd3=[7*x(5)*t.^6+6*x(6)*t.^5+5*c3*t.^4+4*d3*t.^3+3*e3*t.^2];
Qd4=[7*x(7)*t.^6+6*x(8)*t.^5+5*c4*t.^4+4*d4*t.^3+3*e4*t.^2];
Qd5=[7*x(9)*t.^6+6*x(10)*t.^5+5*c5*t.^4+4*d5*t.^3+3*e5*t.^2];
% size(Qd1)
%%jotnt angle acceleratton%%
Qdd1=42*x(1)*t.^5+30*x(2)*t.^4+20*c1*t.^3+12*d1*t.^2+6*e1*t;
Qdd2=42*x(3)*t.^5+30*x(4)*t.^4+20*c2*t.^3+12*d2*t.^2+6*e2*t;
Qdd3=42*x(5)*t.^5+30*x(6)*t.^4+20*c3*t.^3+12*d3*t.^2+6*e3*t;
Qdd4=42*x(7)*t.^5+30*x(8)*t.^4+20*c4*t.^3+12*d4*t.^2+6*e4*t;
Qdd5=42*x(9)*t.^5+30*x(10)*t.^4+20*c5*t.^3+12*d5*t.^2+6*e5*t;
 
P = 0;
% tm = zeros(5,1);
for i=1:1:length(t)
     m=double(subs(M,{q1,q2,q3,q4,q5},{Q1(1,i),Q2(1,i),Q3(1,i),Q4(1,i),Q5(1,i)}));%%Mass Matrtx%%
    Tm=m*[Qdd1(1,i);Qdd2(1,i);Qdd3(1,i);Qdd4(1,i);Qdd5(1,i)];
%     Tm = tm + Tm;
%     Tm = abs(Tm);
  cor= double(subs(C,{q1,q2,q3,q4,q5},{Q1(1,i),Q2(1,i),Q3(1,i),Q4(1,i),Q5(1,i)}));
  Tc=cor*[Qd1(1,i)^2;Qd2(1,i)^2;Qd3(1,i)^2;Qd4(1,i)^2;Qd5(1,i)^2];
%   Tc = abs(Tc);
bor= double(subs(B,{q1,q2,q3,q4,q5},{Q1(1,i),Q2(1,i),Q3(1,i),Q4(1,i),Q5(1,i)}));
p=[Qd1(1,i)*Qd2(1,i);Qd1(1,i)*Qd3(1,i);Qd1(1,i)*Qd4(1,i);Qd1(1,i)*Qd5(1,i);Qd2(1,i)*Qd3(1,i);Qd2(1,i)*Qd4(1,i);Qd2(1,i)*Qd5(1,i);Qd3(1,i)*Qd4(1,i);Qd3(1,i)*Qd5(1,i);Qd4(1,i)*Qd5(1,i)];
Tb= bor*[p(1); p(2);p(3);p(4);p(5);p(6);p(7);p(8);p(9);p(10)]; 

Tt = Tm+Tc+Tb;

Power= 1e-07*(Tt')*[(Qd1(1,i));(Qd2(1,i));(Qd3(1,i));(Qd4(1,i));(Qd5(1,i))];
P = vpa(Power,4) + P;
end
SSS = double(P);
