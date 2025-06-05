function [t,x]=Runge_Kutta4(fun,tb,te,x0,N)

%四阶龙格-库塔方法求解一阶微分方程数值解

%fun 微分方程

%tb t的取值范围的左端点

%te t的取值范围的左端点

%y0 y的迭代初始值

%N 步长

%如果函数的输入参数没有N时，步长数N取默认值200
% 
% if nargin<4
% 
% N = 200;
% 
% end



h = (te-tb)/N;%步长

t = tb+[0:N]'*h;

%x = zeros(30,N+1);

x(:,1) = x0(:); %行列变换
%x0 = [5 5 5*pi/180 0 0 0 0 0 0 0 0 0 0 0 0 3.5*10^6-5 3.5*10^6-5 8*pi*10^9/180-5*pi/180 0 0 0 0.35 0.35 pi/180 0 0 0 0 0 0];

for k=1:N

K1 = feval(fun,t(k),x(:,k));

K2 = feval(fun,t(k)+h/2,x(:,k)+h/2*K1);

K3 = feval(fun,t(k)+h/2,x(:,k)+h/2*K2);

K4 = feval(fun,t(k),x(:,k)+h*K3);

x(:,k+1) = x(:,k) + h/6*(K1+2*K2+2*K3+K4);

end

end