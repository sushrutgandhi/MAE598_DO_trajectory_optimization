
x0 = [5*1e-5 2*1e-4 9*1e-5 7*1e-4 5*1e-5 1.5*1e-4 6*1e-5 1*1e-4 6*1e-5 1*1e-4 ];
A = [];% inequalities
b = [];% inequalities
Aeq = [];
beq = [];
lb=[1e-05 1-05 1e-05 1e-05 1e-05 1e-05 1e-05 1e-05 1e-05 1e-05];
ub=[10 10 10 10 10 10 10 10 10 10];
% nonlcon = @nonlinear;
options = optimoptions(@fmincon,'Algorithm','sqp');
nonlcon = @nonlineq_v2;
x = fmincon(@(x)power_SSS_v3_2_5sec(x),x0,A,b,Aeq,beq,lb,ub,nonlcon,options)





function [q_SSS,qeu_SSS] = nonlineq_v2(x)
 thetaf = [1.707 , 1.2217 , 1.2217 , 1.2217, 1.2217]';%%%final position%%
 thetat=[0.7071 1.7071 1.7071 1.7071 1.7071]';%%%initial position%%
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

 time = [0:0.1:2.5];
 q_SSS = [x(1)*time.^7+x(2)*time.^6+c1*time.^5+d1*time.^4+e1*time.^3+thetat(1,1) - 3.14;-3.14 - (x(1)*time.^7+x(2)*time.^6+c1*time.^5+d1*time.^4+e1*time.^3+thetat(1,1));
     x(3)*time.^7+x(4)*time.^6+c2*time.^5+d2*time.^4+e2*time.^3+thetat(2,1) - 3.14;-3.14 - (x(3)*time.^7+x(4)*time.^6+c2*time.^5+d2*time.^4+e2*time.^3+thetat(2,1));
          x(5)*time.^7+x(6)*time.^6+c3*time.^5+d3*time.^4+e3*time.^3+thetat(3,1) - 1.57;-1.57 - (x(5)*time.^7+x(6)*time.^6+c3*time.^5+d3*time.^4+e3*time.^3+thetat(3,1));
               x(7)*time.^7+x(8)*time.^6+c4*time.^5+d4*time.^4+e4*time.^3+thetat(4,1) - 1.57;-1.57 - (x(7)*time.^7+x(8)*time.^6+c4*time.^5+d4*time.^4+e4*time.^3+thetat(4,1));
                    x(9)*time.^7+x(10)*time.^6+c5*time.^5+d5*time.^4+e5*time.^3+thetat(5,1) - 0.707;-0.707 - (x(9)*time.^7+x(10)*time.^6+c5*time.^5+d5*time.^4+e5*time.^3+thetat(5,1));
           7*x(1)*time.^6+6*x(2)*time.^5+5*c1*time.^4+4*d1*time.^3+3*e1*time.^2-1.48;-1.48-(7*x(1)*time.^6+6*x(2)*time.^5+5*c1*time.^4+4*d1*time.^3+3*e1*time.^2);
       7*x(3)*time.^6+6*x(4)*time.^5+5*c2*time.^4+4*d2*time.^3+3*e2*time.^2-1.48;-1.48-(7*x(3)*time.^6+6*x(4)*time.^5+5*c2*time.^4+4*d2*time.^3+3*e2*time.^2);
       7*x(5)*time.^6+6*x(6)*time.^5+5*c3*time.^4+4*d3*time.^3+3*e3*time.^2-1.74;-1.74-(7*x(5)*time.^6+6*x(6)*time.^5+5*c3*time.^4+4*d3*time.^3+3*e3*time.^2);
       7*x(7)*time.^6+6*x(8)*time.^5+5*c4*time.^4+4*d4*time.^3+3*e4*time.^2-1.31;-1.31-(7*x(7)*time.^6+6*x(8)*time.^5+5*c4*time.^4+4*d4*time.^3+3*e4*time.^2);
       7*x(9)*time.^6+6*x(10)*time.^5+5*c5*time.^4+4*d5*time.^3+3*e5*time.^2-2.26;-2.26-(7*x(9)*time.^6+6*x(10)*time.^5+5*c5*time.^4+4*d5*time.^3+3*e5*time.^2)
       42*x(1)*time.^5+30*x(2)*time.^4+20*c1*time.^3+12*d1*time.^2+6*e1*time-3;-2-(42*x(1)*time.^5+30*x(2)*time.^4+20*c1*time.^3+12*d1*time.^2+6*e1*time);
       42*x(3)*time.^5+30*x(4)*time.^4+20*c2*time.^3+12*d2*time.^2+6*e2*time-3;-2-(42*x(3)*time.^5+30*x(4)*time.^4+20*c2*time.^3+12*d2*time.^2+6*e2*time);
       42*x(5)*time.^5+30*x(6)*time.^4+20*c3*time.^3+12*d3*time.^2+6*e3*time-3;-2-( 42*x(5)*time.^5+30*x(6)*time.^4+20*c3*time.^3+12*d3*time.^2+6*e3*time);
       42*x(7)*time.^5+30*x(8)*time.^4+20*c4*time.^3+12*d4*time.^2+6*e4*time-3;-2-(42*x(7)*time.^5+30*x(8)*time.^4+20*c4*time.^3+12*d4*time.^2+6*e4*time);
       42*x(9)*time.^5+30*x(10)*time.^4+20*c5*time.^3+12*d5*time.^2+6*e5*time-3;-2-(42*x(9)*time.^5+30*x(10)*time.^4+20*c5*time.^3+12*d5*time.^2+6*e5*time)];
 qeu_SSS = [];
end


