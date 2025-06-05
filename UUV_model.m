function dx = case1(t,x)
eta=x(1:3);
v=x(4:6);            
etam=x(7:9);          
vm=x(10:12);
zetam=x(13:15);
beta=x(16:18);      
alphag=x(19:21);     
rho=x(22:24);        
rho1=x(22);
rho2=x(23);
rho3=x(24);
eta0g=x(25:27);
K1=diag([5 5 5]);            K2=diag([10 10 10]);            K3=diag([5 5 5]);
C1=diag([1 1 1]);            C2=diag([1 1 1]);               C3=diag([1 1 1]); 
K01=15*diag([1 1 1]);        C01=diag([1 1 1]);   
K02=15*diag([1 1 1]);        C02=diag([1 1 1]);     
K03=15*diag([1 1 1]);        C03=diag([1 1 1]);
M=[25.8 0 0;0 24.6612 1.0948;0 1.0948 2.76]; 
C=[0 0 -24.6612*x(5)-1.0948*x(6);0 0 25.8*x(4);24.6612*x(5)+1.0948*x(6) -25.8*x(4) 0];
D=[
    0.7225+1.3274*abs(x(4))+5.8664*(x(4))^2 0 0;
    0 0.8612+36.2823*abs(x(5))+8.05*abs(x(6)) -0.1052-5.0437*abs(x(5))-0.13*abs(x(6));
    0 0.845*abs(x(5))+3.45*abs(x(6))-0.1079 1.9-0.08*abs(x(5))+0.75*abs(x(6))
    ];
R=[cos(x(3)) -sin(x(3)) 0;sin(x(3)) cos(x(3)) 0;0 0 1];
  d=[sin(0.3*t) sin(0.4*t) sin(0.3*t)]';
l=0.3;   l0=0.3;
T=20;       w=0.03;     TE=6;       Tm=15; 
if t>0 && t<T
 ze=exp(w*(T+TE-t))-1;
 ze_dot=-w*exp(w*(T+TE-t));
else
    ze=exp(w*(t-T));
    ze_dot=w*exp(w*(t-T));
end
etam_dot=-(K01+C01*abs(ze_dot)/ze)*(etam-eta)+R*vm;        
zetam_dot=-(K03+C03*abs(ze_dot)/ze)*R'*(etam-eta); 
etad=[0.2*t 0 0]';
etad_dot=[0.2 0 0]';
p=[x(1) x(2)]';  
x0=10; 
 y0=-1;
zhangai=[x0 y0]'; 
pv=norm(p-zhangai);  
R0g=2.6;
if pv<=R0g       
      k=2;
      alpha00=0.1;
      p0=alpha00*exp(-k*(pv-R0g)); 
else
    p0=0;
end
eta0=[p0 p0 0]';
eta0g_dot=(eta0-eta0g)/l0; 
if t>0 && t<T
    ze=exp(w*(T+TE-t))-1;            
    ze_dot=-w*exp(w*(T+TE-t));        
    ze_ddot=w^2*exp(w*(T+TE-t));      
else
    ze=exp(w*(t-T));
    ze_dot=w*exp(w*(t-T));
    ze_ddot=w^2*exp(w*(t-T)); 
end
krho=1;  crho=2;   delta1=1;  delta2=1;
erhok1=21;   erhok2=20;  erhok3=15*pi/180;
zetag=krho-crho*abs(ze_dot)/ze;
ksita=zetag^2+(pi/Tm)^2;
kbeta=zetag^2+(pi/(T-Tm))^2;
grhokx=(1/(2*zetag)-zetag/(2*kbeta))*exp(-zetag*(T-Tm))+1/(2*zetag)-zetag/(2*kbeta);
grhoks=(zetag/ksita-1/zetag)*exp(-zetag*T)+(1/ksita-1/kbeta)*exp(-zetag*(T-Tm))+1/zetag-zetag/kbeta;
grhok=grhokx/grhoks;
if t>=0 &&t<Tm
    xrhok1=(krho+crho*abs(ze_dot)/ze)*erhok1*(1-grhok*(1-cos(pi*t/Tm)));
    xrhok2=(krho+crho*abs(ze_dot)/ze)*erhok2*(1-grhok*(1-cos(pi*t/Tm)));
    xrhok3=(krho+crho*abs(ze_dot)/ze)*erhok3*(1-grhok*(1-cos(pi*t/Tm)));
elseif t>=Tm && t<T
    xrhok1=(krho+crho*abs(ze_dot)/ze)*erhok1*(1/2-grhok)*(1+cos((pi*t-pi*Tm)/(T-Tm)));
    xrhok2=(krho+crho*abs(ze_dot)/ze)*erhok2*(1/2-grhok)*(1+cos((pi*t-pi*Tm)/(T-Tm)));
    xrhok3=(krho+crho*abs(ze_dot)/ze)*erhok3*(1/2-grhok)*(1+cos((pi*t-pi*Tm)/(T-Tm)));
else
    xrhok1=0;
    xrhok2=0;
    xrhok3=0;
end
rho = [rho1 rho2 rho3]';
omega1=1;      omega2=2;        omega3=0.25;
epsilon=1;
if pv<=R0g 
     rhoT1=1+omega1*(p0+epsilon * p0^2); 
    rhoT2=1+omega2*p0;
    rhoT3=1+omega3*p0;
else
    rhoT1=1;
    rhoT2=1;
    rhoT3=1;
end
rhoT = [rhoT1 rhoT2 rhoT3]';
xrhok = [xrhok1 xrhok2 xrhok3]';
rho_dot=-(krho+crho*abs(ze_dot)/ze).*(rho-rhoT)+xrhok;
e1L=delta1*rho1;                        e1u=delta2*rho1;
e2L=delta1*rho2;                        e2u=delta2*rho2;
e3L=delta1*rho3;                        e3u=delta2*rho3;
e=eta-etad-eta0; 
fai1=1/(2*(e(1)+e1L))-1/(2*(e(1)-e1u));
fai2=1/(2*(e(2)+e2L))-1/(2*(e(2)-e2u));
fai3=1/(2*(e(3)+e3L))-1/(2*(e(3)-e3u));
fai=diag([fai1 fai2 fai3]);
f1=rho_dot(1)/rho1*fai1;   
f2=rho_dot(2)/rho2*fai2;
f3=rho_dot(3)/rho3*fai3;
f=diag([f1 f2 f3]);
z11=0.5*log((e(1)+e1L)/(e1u-e(1)))-0.5*log(delta1/delta2);
z12=0.5*log((e(2)+e2L)/(e2u-e(2)))-0.5*log(delta1/delta2);
z13=0.5*log((e(3)+e3L)/(e3u-e(3)))-0.5*log(delta1/delta2);
Z1=[z11 z12 z13]';   
alpha=R'*(-inv(fai)*(K1+C1*abs(ze_dot)/ze)*Z1+etad_dot+eta0g_dot+inv(fai)*f*e);
z2=vm-alpha-beta; 

alphag_dot=(alpha-alphag)/l; 
tauc=M*(-(K3+C3*abs(ze_dot)/ze)*z2-zetam+alphag_dot+K2*beta-R'*fai*Z1);  
 tauM=[20 20 15]';
if abs(tauc(1))>=tauM(1) 
    tau(1,:) = sign(tauc(1))*tauM(1);
else
    tau(1,:) = tauc(1);
end
if abs(tauc(2))>=tauM(2) 
    tau(2,:) = sign(tauc(2))*tauM(2);
else
    tau(2,:) = tauc(2);
end
if abs(tauc(3))>=tauM(3) 
    tau(3,:) = sign(tauc(3))*tauM(3);
else
    tau(3,:) = tauc(3);
end
vm_dot=-(K02+C02*abs(ze_dot)/ze)*R'*(etam-eta)+zetam+inv(M)*tau; 
W=tauc-tau;
r1=0.1;  r2=1;
if norm(beta) <= r1
    sbeta = 0;
elseif norm(beta) > r1 && norm(beta) < r2
    sbeta = 1 - cos((pi / 2) * sin((pi / 2) * ( norm(beta)^2-r1^2) / (r2^2 - r1^2)));
else
    sbeta = 1;
end
beta_dot=-(K2+C2*abs(ze_dot)/ze)*beta-inv(M)*W*sbeta;
eta_dot=R*v;  
zeta=inv(M)*(-C*v-D*v+d);
v_dot=inv(M)*tau+zeta;     
dx=[eta_dot;v_dot;etam_dot;vm_dot;zetam_dot;beta_dot;alphag_dot;rho_dot;eta0g_dot];
end