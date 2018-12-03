% syms x1 x2
% eqn1 = x1*t^7 + x2*t^6 + (900*x1 - 15*x2 + 0.00327)*t^5 + (75*x2 -9500*x1-0.0409)*t^4 + (24375*x1 - 125*x2 + 0.1365)*t^3 + 0.071==6.28;
% eqn2 = x1 + x2 + (900*x1 - 15*x2 + 0.00327) + (75*x2 -9500*x1-0.0409) + (24375*x1 - 125*x2 + 0.1365)+ 0.071== 0;
% sol = solve([eqn1, eqn2], [x1,x2])
% x1Sol = sol.x1
% x2Sol = sol.x2
load("x_2_5")
% x =[0.0000    0.0097    0.0000    0.0000    0.0000    0.0122    0.0000    0.0109    0.0000    0.0129];
x1 = x(1)          
x2 = x(2)
x3 =x(3)
x4 =  x(4)      

x5 = x(5)
x6 = x(6)
x7 = x(7)
x8 = x(8)
x9 =    x(9)
x10 = x(10)    
q1 = zeros(6,1);q2 = zeros(6,1);q3 = zeros(6,1);q4 = zeros(6,1);q5 = zeros(6,1);
 thetat=[0.7071 1.7071 1.7071 1.7071 1.7071];
 i = 1;
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
for tsss = 0:0.1:2.5
  q1(i,1) = x1*tsss.^7 + x2*tsss.^6 + Zsss1(1,1)*tsss.^5 + Zsss1(2,1)*tsss.^4 + Zsss1(3,1)*tsss.^3 +thetat(1,1);
  q2(i,1) = x3*tsss.^7 + x4*tsss.^6 + Zsss2(1,1)*tsss.^5 + Zsss2(2,1)*tsss.^4 + Zsss2(3,1)*tsss.^3 +thetat(2,1);
  q3(i,1) = x5*tsss.^7 + x6*tsss.^6 + Zsss3(1,1)*tsss.^5 + Zsss3(2,1)*tsss.^4 + Zsss3(3,1)*tsss.^3 + thetat(3,1);
  q4(i,1) = x7*tsss.^7 + x8*tsss.^6 + Zsss4(1,1)*tsss.^5 + Zsss4(2,1)*tsss.^4 + Zsss4(3,1)*tsss.^3 + thetat(4,1);
  q5(i,1) = x9*tsss.^7 + x10*tsss.^6 + Zsss5(1,1)*tsss.^5 + Zsss5(2,1)*tsss.^4 + Zsss5(3,1)*tsss.^3 + thetat(5,1);
  i=i+1;
end
tsss = [0:0.1:2.5]'
figure(1)
plot(tsss,q1)
hold on
plot(tsss,q2)
hold on
plot(tsss,q3)
hold on
plot(tsss,q4)
hold on
plot(tsss,q5)
xlabel('Time(seconds)')
ylabel('Joint angle(Radians)')
title("Joint Position t=2.5sec")
 legend('joint-1','joint-2','joint-3','joint-4','joint-5')
set(gcf,'color','w')
ylim([-3.14 3.14])
 QQQ=[q1 q2 q3 q4 q5];
 figure()
 %%Velocity Plot%%
 i=1;
 for tsss = 0:0.1:2.5
Qd1(i,1)=[7*x(1)*tsss.^6+6*x(2)*tsss.^5+5*Zsss1(1,1)*tsss.^4+4*Zsss1(2,1)*tsss.^3+3*Zsss1(3,1)*tsss.^2];
Qd2(i,1)=[7*x(3)*tsss.^6+6*x(4)*tsss.^5+5*Zsss2(1,1)*tsss.^4+4*Zsss2(2,1)*tsss.^3+3*Zsss2(3,1)*tsss.^2];
Qd3(i,1)=[7*x(5)*tsss.^6+6*x(6)*tsss.^5+5*Zsss3(1,1)*tsss.^4+4*Zsss3(2,1)*tsss.^3+3*Zsss3(3,1)*tsss.^2];
Qd4(i,1)=[7*x(7)*tsss.^6+6*x(8)*tsss.^5+5*Zsss4(1,1)*tsss.^4+4*Zsss4(2,1)*tsss.^3+3*Zsss4(3,1)*tsss.^2];
Qd5(i,1)=[7*x(9)*tsss.^6+6*x(10)*tsss.^5+5*Zsss5(1,1)*tsss.^4+4*Zsss5(2,1)*tsss.^3+3*Zsss5(3,1)*tsss.^2];
  i=i+1;
end
tsss = [0:0.1:2.5]';
figure(2)
plot(tsss,Qd1)
hold on
plot(tsss,Qd2)
hold on
plot(tsss,Qd3)
hold on
plot(tsss,Qd4)
hold on
plot(tsss,Qd5)
xlabel('Time(seconds)')
ylabel('Joint velocity(Radians/sec)')
title("Joint Velocity t=2.5sec")
 legend('joint-1','joint-2','joint-3','joint-4','joint-5')
set(gcf,'color','w')
%ylim([-3.14 3.14])
 QQQ_V=[q1 q2 q3 q4 q5];
 
 
 %%%Acceleration plot%%%
 i=1;
 for tsss = 0:0.1:2.5
  Qdd1(i,1)=42*x(1)*tsss.^5+30*x(2)*tsss.^4+20*Zsss1(1,1)*tsss.^3+12*Zsss1(2,1)*tsss.^2+6*Zsss1(3,1)*tsss;
 Qdd2(i,1)=42*x(3)*tsss.^5+30*x(4)*tsss.^4+20*Zsss2(1,1)*tsss.^3+12*Zsss2(2,1)*tsss.^2+6*Zsss2(3,1)*tsss;
 Qdd3(i,1)=42*x(5)*tsss.^5+30*x(6)*tsss.^4+20*Zsss3(1,1)*tsss.^3+12*Zsss3(2,1)*tsss.^2+6*Zsss3(3,1)*tsss;
 Qdd4(i,1)=42*x(7)*tsss.^5+30*x(8)*tsss.^4+20*Zsss4(1,1)*tsss.^3+12*Zsss4(2,1)*tsss.^2+6*Zsss4(3,1)*tsss;
  Qdd5(i,1)=42*x(9)*tsss.^5+30*x(10)*tsss.^4+20*Zsss5(1,1)*tsss.^3+12*Zsss5(2,1)*tsss.^2+6*Zsss5(3,1)*tsss;
  i=i+1;
end
tsss = [0:0.1:2.5]';
figure(3)
plot(tsss,Qdd1)
hold on
plot(tsss,Qdd2)
hold on
plot(tsss,Qdd3)
hold on
plot(tsss,Qdd4)
hold on
plot(tsss,Qdd5)
xlabel('Time(seconds)')
ylabel('Joint acceleration(Radians/sec^2)')
 legend('joint-1','joint-2','joint-3','joint-4','joint-5')
 title("Joint Acceleration t=2.5sec")
set(gcf,'color','w')
%ylim([-3.14 3.14])
 QQQ_A=[q1 q2 q3 q4 q5];
 