tf = 1;
steps = 1000;
dt = tf/steps;

solinit = bvpinit(linspace(0,tf,steps), [0 0 0 0 0 0 ...
                                         0 0 0 0 0 0 ...
                                         0 0 0 0 0 0 ...
                                         0 0 0 0 0 0]); %can be chosen arbitrary
options = bvpset('Stats','on','RelTol',1e-1);
sol = bvp4c(@BVP_ode, @BVP_bc, solinit, options);
tau = sol.x;
y = sol.y;
FW4Plotfunction
figure(4)
plot3(y(1,:),y(3,:),y(5,:),'-m',linewidth=3)
grid on

function dydt = BVP_ode(t,y)
g=9.81;
m=1;
%I = [Ixx Iyy Izz]
I = [1 1 1];
u = [0 0 0 0];

u(1) = -y(18)/m;
u(2) = -y(20)/I(1);
u(3) = -y(22)/I(2);
u(4) = -y(24)/I(3);

%y(1:12) = qdot   y(13:24) = lamndadot y(25) = tf
dydt = [
y(2);
-g*y(9);
y(4);
g*y(7);
y(6);
g-u(1)/m;
y(8);
u(2)/I(1);
y(10)
u(3)/I(2);
y(12);
u(4)/I(3);
0;
-(y(2)+y(13));
0;
-(y(4)+y(15));
-(g);
-(y(6)+y(17));
-(g*y(16));
-(y(8)+y(19));
g*y(14);
-(y(10)+y(21));
0;
-(y(12)+y(23))];
end

%y(1:12) = qdot   y(13:24) = lamndadot y(25) = tf
function res = BVP_bc(ya,yb)
xtf = 10;
ytf = -1;
ztf = 200*griewank([xtf ytf]);

x0 = 0;
y0 = 0;
z0 = 200*griewank([x0 y0]);

res = [ ya(1) - x0;
ya(2) - 0;
ya(3) - y0;
ya(4) - 0;
ya(5) - z0;
ya(6) - 0;
ya(7) - 0;
ya(8) - 0;
ya(9) - 0;
ya(10) - 0;
ya(11) - 0;
ya(12) - 0;
yb(1) - xtf;
yb(2) - 0;
yb(3) - ytf;
yb(4) - 0;
yb(5) - ztf;
yb(6) - 0;
yb(7) - 0;
yb(8) - 0;
yb(9) - 0;
yb(10) - 0;
yb(11) - 0;
yb(12) - 0];
end