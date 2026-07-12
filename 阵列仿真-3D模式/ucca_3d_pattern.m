% =====================================================================
% ucca_3d_pattern.m  —  同心圆环阵 三维方向图
% ---------------------------------------------------------------------
%  3D 阵因子 (多环 + 可选中心元):
%  AF(θ,φ) = (1/Q)[ c + Σ_m Σ_n exp{ j·2πR_m/λ·( sinθ·cos(φ-φ_mn) - sinθ_0·cos(φ_0-φ_mn) ) } ]
%  φ_mn = 2πn/N_m + β_m ;  c=中心元(1或0);  Q=总阵元数
% 每行 rings = [N_m, R_m/λ, β_m(度)]。改这张表即可换任意圆环阵。
% =====================================================================
clear; clc; close all;

% ---- 阵列配置 (可改) ----
% 例: 中心 + 内环6 + 外环12 
center = 1;                      % 1=有中心元, 0=无
rings  = [ 8,   1.0,  0;         % [阵元数, 半径/λ, 偏置角°]
           16,  2.0, 0];
theta_0 = 0;                     % 指向仰角(度)
phi_0   = 0;                     % 指向方位(度)

% ---- 角度网格 (上半球) ----
[TH,PH] = meshgrid(deg2rad(0:1:180), deg2rad(0:1:360));
th0 = deg2rad(theta_0); ph0 = deg2rad(phi_0);

% ---- 3D 阵因子 ----
AF = center * ones(size(TH));
for m = 1:size(rings,1)  % 圈数
    Nm = rings(m,1); Rm = rings(m,2); bm = deg2rad(rings(m,3));
    for n = 0:Nm-1
        phmn = 2*pi*n/Nm + bm;
        AF = AF + exp(1j*2*pi*Rm*( sin(TH).*cos(PH-phmn) - sin(th0).*cos(ph0-phmn) ));
    end
end
Q = center + sum(rings(:,1));
P = abs(AF)/Q;                   % 归一化, 峰值=1

X=P.*sin(TH).*cos(PH); 
Y=P.*sin(TH).*sin(PH); 
Z=P.*cos(TH);


% (1) 几何俯视图
figure(1);
%中心阵元黑色三角
if center, polarplot(0,0,'k^','MarkerFaceColor','k','MarkerSize',9); hold on; end
col='brgmc';
for m=1:size(rings,1)
    Nm=rings(m,1); Rm=rings(m,2); bm=deg2rad(rings(m,3));
    ph=2*pi*(0:Nm-1)/Nm + bm;
    polarplot(ph, Rm*ones(1,Nm),'o','Color',col(m),'MarkerFaceColor',col(m)); hold on;
end
title('阵列几何 (俯视)'); rlim([0 max(rings(:,2))*1.1]);


figure(2);
% (2) 3D 方向图 (线性)
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

xlabel('x'); ylabel('y'); zlabel('z'); view(135,25);
title(sprintf('同心圆环阵 3D 方向图\n总元数 Q=%d, 指向(\\theta_0=%d°,\\phi_0=%d°)', ...
      Q, theta_0, phi_0));

figure(3);
% (3) dB 曲面 (下限 -30dB)
PdB = 20*log10(P+eps);  floorDB = -30;
Pd  = max(PdB, floorDB) - floorDB;             % 平移到 [0, |floor|]
Xd=Pd.*sin(TH).*cos(PH); Yd=Pd.*sin(TH).*sin(PH); Zd=Pd.*cos(TH);

surf(Xd,Yd,Zd,PdB,'EdgeColor','none');
axis equal;
colormap(turbo);
c=colorbar; c.Label.String='dB';
shading interp;
%camlight;
%lighting gouraud;
clim([floorDB 0]); xlabel('x'); ylabel('y'); zlabel('z'); view(135,25);
title(sprintf('同心圆环阵 3D 方向图 (dB, 下限-30)\n总元数 Q=%d, 指向(\\theta_0=%d°,\\phi_0=%d°)', ...
      Q, theta_0, phi_0));
