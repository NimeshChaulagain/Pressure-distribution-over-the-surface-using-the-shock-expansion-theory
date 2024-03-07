%%This code is used for computing the pressure above a surface in a
%%super/hypersonic flow using shock expansion theory.
%M0 is the initial mach number
%T0 is the deflection angle at the nose, in accordance with the PM theory the deflection angle must be less than the maximum deflection angle.
%B0 is the nose wave angle
% nu refers to the PM function 
%Mn refers to the normal component of the flow across the shockwave
%%

hold all
M0 = 8; %incoming flow mach number
T0 = 30; %nose deflection angle in degree
B0 = 39; %nose wave angle

Tr0 = convangle(T0); %nose deflection angle in radian
Br0 = convangle(B0);
%calculating properties past the oblique shock
Mn0 = M0*sin(Br0);
Mn1 = ((1+0.2*Mn0^2)/(1.4*Mn0^2-0.2))^0.5;
M1 = Mn1 /sin(Br0-Tr0);

dT = 1:1:60;

[mach,nu1,mu] = flowprandtlmeyer(1.4,M1,'mach'); 

nu2 = nu1 + dT;
for j = 1:60
    
[M2(j),nu2(j),mu] = flowprandtlmeyer(1.4,nu2(j),'nu'); 

end

num = 1+0.2*M1^2;
denum = 1 + 0.2*M2.^2;
 P2 = ((num./denum).^3.5 );
 x = 0:1:60;
 Pn(1) = 1;
 for i = 2:61
     Pn(i) = P2(i-1);   
 end
 plot(x, Pn, 'k');
 xlabel('Defection angle');
 ylabel('Pressure');

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 M01 = 6; %incoming flow mach number
T01 = 30; %nose deflection angle in degree
B01 = 41; %nose wave angle

Tr01 = convangle(T01); %nose deflection angle in radian
Br01 = convangle(B01);
%calculating properties past the oblique shock
Mn01 = M01*sin(Br01);
Mn11 = ((1+0.2*Mn01^2)/(1.4*Mn01^2-0.2))^0.5;
M11 = Mn11 /sin(Br01-Tr01);

dT1 = 1:1:60;

[mach,nu11,mu] = flowprandtlmeyer(1.4,M11,'mach'); 

nu21 = nu11 + dT1;
for j = 1:60
    
[M21(j),nu21(j),mu] = flowprandtlmeyer(1.4,nu21(j),'nu'); 

end

num1 = 1+0.2*M11^2;
denum1 = 1 + 0.2*M21.^2;
 P21 = ((num1./denum1).^3.5 );
 x1 = 0:1:60;
 Pn1(1) = 1;
 for i = 2:61
     Pn1(i) = P21(i-1);   
 end
 
 

 plot(x1, Pn1, 'b');
 xlabel('Defection angle');
 ylabel('Pressure');


 M02 = 5; %incoming flow mach number
T02 = 30; %nose deflection angle in degree
B02 = 42; %nose wave angle

Tr02 = convangle(T02); %nose deflection angle in radian
Br02 = convangle(B02);
%calculating properties past the oblique shock
Mn02 = M02*sin(Br02);
Mn12 = ((1+0.2*Mn02^2)/(1.4*Mn02^2-0.2))^0.5;
M12 = Mn12 /sin(Br02-Tr02);

dT2 = 1:1:60;

[mach,nu12,mu] = flowprandtlmeyer(1.4,M12,'mach'); 

nu22 = nu12 + dT2;
for j = 1:60
    
[M22(j),nu22(j),mu] = flowprandtlmeyer(1.4,nu22(j),'nu'); 

end

num2 = 1+0.2*M12^2;
denum2 = 1 + 0.2*M22.^2;
 P22 = ((num2./denum2).^3.5 );
 x2 = 0:1:60;
 Pn2(1) = 1;
 for i = 2:61
     Pn2(i) = P22(i-1);   
 end
 
 

 plot(x2, Pn2, 'r');
 xlabel('Defection angle');
 ylabel('Pressure ratio (Pi/Pn)');
 legend('Mach 8','Mach 6','Mach 5');

 
function degtorad = convangle(x)

degtorad = (pi/180)*x;

end