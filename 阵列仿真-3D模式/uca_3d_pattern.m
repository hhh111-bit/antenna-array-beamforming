% =====================================================================
% uca_3d_pattern.m  —  均匀圆阵(UCA) 三维方向图
% ---------------------------------------------------------------------
% 3D 阵因子 (含 sinθ 仰角因子, 阵列在 xy 平面):
%   AF(θ,φ) = (1/N) Σ_n exp{ j·2πR/λ·[ sinθ·cos(φ-φ_n) - sinθ_0·cos(φ_0-φ_n) ] }
%   φ_n = 2πn/N (阵元方位);  (θ_0,φ_0)=波束指向(仰角,方位)
% 符号: θ=仰角(从+z), φ=方位。 θ_0=0 => 波束指向天顶(+z)。
% =====================================================================
clear; clc; close all;

N       = 12;           % 阵元数
R       = 1.0;          % 半径 / λ
theta_0 = 0;            % 波束指向 仰角 (度)
phi_0   = 0;            % 波束指向 方位 (度)

% ---- 角度网格 (上半球) ----
th_v = deg2rad(0:1:180);        % 仰角 0..180
ph_v = deg2rad(0:1:360);        % 方位 0..360
[TH, PH] = meshgrid(th_v, ph_v);

th0 = deg2rad(theta_0);  
ph0 = deg2rad(phi_0);
phn = 2*pi*(0:N-1)/N;                          % 阵元方位
    
% ---- 3D 阵因子 (对阵元求和, 矢量化在网格上) ----
AF = zeros(size(TH));
for n = 1:N
    AF = AF + exp(1j*2*pi*R*( sin(TH).*cos(PH-phn(n)) - sin(th0).*cos(ph0-phn(n)) ));
end
P = abs(AF)/ N;                                 % 归一化幅度, 峰值=1

% ---- 转直角坐标 (半径=方向图幅度) ----
X = P.*sin(TH).*cos(PH);
Y = P.*sin(TH).*sin(PH);
Z = P.*cos(TH);

% (1) 几何俯视图
figure(1);
polarplot(phn, R*ones(1,N), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 9);
title(sprintf('阵列几何 (俯视)  N=%d, R=%.1f\\lambda', N, R));
rlim([0 R]);

% (2) 线性幅度 3D 曲面
figure(2);
surf(X,Y,Z,P,'EdgeColor','none'); hold on;
% 画指向箭头
xa=sin(th0)*cos(ph0); ya=sin(th0)*sin(ph0); za=cos(th0);
plot3([0 xa],[0 ya],[0 za],'r-','LineWidth',2);

axis equal; 
colormap(turbo); 
colorbar; 
shading interp; 
%camlight; 
%lighting gouraud;
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('UCA 3D 方向图 (线性)\nN=%d, R=%.1f\\lambda, 指向(\\theta_0=%d°,\\phi_0=%d°)', ...
      N,R,theta_0,phi_0));

figure(3);
% (3) dB 曲面
% ----  (下限 -30dB) ----
PdB = 20*log10(P+eps);  floorDB = -30;
Pd  = max(PdB, floorDB) - floorDB;             % 平移到 [0, |floor|]
Xd=Pd.*sin(TH).*cos(PH); Yd=Pd.*sin(TH).*sin(PH); Zd=Pd.*cos(TH);
surf(Xd,Yd,Zd,PdB,'EdgeColor','none');
axis equal; colormap(turbo); c=colorbar; c.Label.String='dB';
shading interp;
clim([floorDB 0]); xlabel('x'); ylabel('y'); zlabel('z');
title('UCA 3D 方向图 (dB, 下限-30)');




