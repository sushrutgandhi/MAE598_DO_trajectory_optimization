 L = Link('revolute', 'd', 0, 'a', 0, 'alpha', pi/2);
 L(2) = Link('revolute', 'd', 0.0397, 'a', 0.14611, 'alpha', 0);
 L(3) = Link('revolute', 'd', 0, 'a', 0.14551, 'alpha', 0);
 L(4) = Link('revolute', 'd', 0, 'a', 0, 'alpha', pi/2);
 L(5) = Link('revolute', 'd', 0.046, 'a', 0, 'alpha', 0);
 
 R = SerialLink(L)
R.plot(QQQ,'fps',2,'movie','P')

T = R.fkine(QQQ);
T1 = transl(T);
Tx = T1(:,1);Ty = T1(:,2);Tz = T1(:,3);
figure(2)
scatter3(Tx,Ty,Tz);
hold on
line(Tx,Ty,Tz)
title("End effector for 2.5 sec")

