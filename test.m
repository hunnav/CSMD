f = @(a,b,c,d,e,f,g) (a-10).^2+5*(b-12).^2+c.^4+3*(d-11).^2+10*e.^6+7*f.^2+g.^4-4*f.*g-10*f-8*g;

g1 = @(a,b,c,d,e,f,g) 127-2*a^2-3*b^4-c-4*d^2-5*e;
g2 = @(a,b,c,d,e,f,g) 282-7*a-3*b-10*c^2-d+e;
g3 = @(a,b,c,d,e,f,g) 196-23*a-b^2-6*f^2+8*g;
g4 = @(a,b,c,d,e,f,g) -4*a^2-b^2+3*a*b-2*c^2-5*f+11*g;

x=[2.3331998,1.942837, -0.479546,  4.387811, -0.632617,1.039775,1.600996];
A = -642.2563; B = 1.0525; C = 1.0416;

object = f(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
constraint1 = g1(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
constraint2 = g2(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
constraint3 = g3(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
constraint4 = g4(x(1),x(2),x(3),x(4),x(5),x(6),x(7));

object(:,2) = log((object)+A)/log(B)
constraint1(:,2) = log((10^(floor(1+log10(abs(object(:,2)))))*max(-constraint1,0))+C^1)/log(C)-1
constraint2(:,2) = log((10^(floor(1+log10(abs(object(:,2)))))*max(-constraint2,0))+C^1)/log(C)-1
constraint3(:,2) = log((10^(floor(1+log10(abs(object(:,2)))))*max(-constraint3,0))+C^1)/log(C)-1
constraint4(:,2) = log((10^(floor(1+log10(abs(object(:,2)))))*max(-constraint4,0))+C^1)/log(C)-1


% S = struct('add', add, 'acqui', acqui, 'Hypopt',Hypopt,'print',print,'prob',prob)
% addpoints(animatedline,S.print.x,S.print.y); drawnow;