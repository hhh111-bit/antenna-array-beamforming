% =====================================================================
% uva_3d_pattern.m  —  均匀立体阵(UVA, 长方体体阵) 三维方向图
% ---------------------------------------------------------------------
% 标准阵因子 (sin 函数形式, 阵元沿 x,y,z 三向均匀排布):
%   ψx = kdx·sinθ·cosφ + αx,   αx = -kdx·sinθ_0·cosφ_0
%   ψy = kdy·sinθ·sinφ + αy,   αy = -kdy·sinθ_0·sinφ_0
%   ψz = kdz·cosθ       + αz,   αz = -kdz·cosθ_0
%   AF = | AFx·AFy·AFz |,  AFi = sin(Ni·ψi/2)/(Ni·sin(ψi/2))
%   (θ_0,φ_0) = 波束指向 (仰角, 方位)。三向方向图相乘。
% 改 Nx,Ny,Nz / d / theta_0,phi_0 即可扫描任意指向。
% =====================================================================
clear; clc; close all;

Nx = 6;  Ny = 6;  Nz = 6;                                % 三向阵元数
dx_lambda = 0.5;  dy_lambda = 0.5;  dz_lambda = 0.5;     % 间距 / λ
theta_0 = 0;   phi_0 = 0;                              % 波束指向 (仰角,方位,度)

kdx = 2*pi*dx_lambda;  kdy = 2*pi*dy_lambda;  kdz = 2*pi*dz_lambda;
th0 = deg2rad(theta_0);  ph0 = deg2rad(phi_0);

% ---- 渐进馈电相位 (使主瓣指向 θ_0,φ_0) ----
ax_ = -kdx*sin(th0)*cos(ph0);
ay_ = -kdy*sin(th0)*sin(ph0);
az_ = -kdz*cos(th0);

% ---- 球面网格 (整球) ----
th_v = deg2rad(0:1:180);
ph_v = deg2rad(0:1:360);
[TH, PH] = meshgrid(th_v, ph_v);

% ---- 标准阵因子 (三向 sin 之积) ----
psi_x = kdx*sin(TH).*cos(PH) + ax_;
psi_y = kdy*sin(TH).*sin(PH) + ay_;
psi_z = kdz*cos(TH)          + az_;
AF_x = sin(Nx*psi_x/2)./(Nx*sin(psi_x/2));  AF_x(abs(sin(psi_x/2))<1e-9)=1;
AF_y = sin(Ny*psi_y/2)./(Ny*sin(psi_y/2));  AF_y(abs(sin(psi_y/2))<1e-9)=1;
AF_z = sin(Nz*psi_z/2)./(Nz*sin(psi_z/2));  AF_z(abs(sin(psi_z/2))<1e-9)=1;
P = abs(AF_x .* AF_y .* AF_z);
P = P / max(P(:));               % 归一化, 峰值=1

X = P.*sin(TH).*cos(PH);
Y = P.*sin(TH).*sin(PH);
Z = P.*cos(TH);

% (1) 阵列几何 (立体排布, 关于原点对称)
figure(1);
xv = ((0:Nx-1) - (Nx-1)/2) * dx_lambda;            % 相位中心在原点
yv = ((0:Ny-1) - (Ny-1)/2) * dy_lambda;
zv = ((0:Nz-1) - (Nz-1)/2) * dz_lambda;
[xg,yg,zg] = meshgrid(xv, yv, zv);
plot3(xg(:), yg(:), zg(:), 'o', 'MarkerFaceColor',[.15 .45 .85], ...
      'MarkerEdgeColor','k', 'MarkerSize',6, 'LineWidth',0.5); hold on;
plot3(0, 0, 0, 'rp', 'MarkerFaceColor','r', 'MarkerSize',14);   % 原点=相位中心
axis equal; grid on; box on; view(135,25);
xlim([xv(1)-dx_lambda, xv(end)+dx_lambda]);
ylim([yv(1)-dy_lambda, yv(end)+dy_lambda]);
zlim([zv(1)-dz_lambda, zv(end)+dz_lambda]);
xlabel('x / \lambda'); ylabel('y / \lambda'); zlabel('z / \lambda');
title(sprintf('阵列几何 (立体)  UVA %d×%d×%d, d=%.2f\\lambda (对称于原点)', Nx,Ny,Nz,dx_lambda));

% (2) 线性幅度 3D 曲面
figure(2);
surf(X,Y,Z,P,'EdgeColor','none'); hold on;
xa=sin(th0)*cos(ph0); ya=sin(th0)*sin(ph0); za=cos(th0);   % 指向箭头
plot3([0 xa],[0 ya],[0 za],'r-','LineWidth',2);
axis equal; colormap(turbo); colorbar; shading interp;
xlabel('x'); ylabel('y'); zlabel('z'); view(135,25);
title(sprintf('UVA 3D 方向图 (线性)\n%d×%d×%d, d=%.2f\\lambda, 指向(\\theta_0=%d°,\\phi_0=%d°)', ...
      Nx,Ny,Nz,dx_lambda,theta_0,phi_0));

% (3) dB 曲面 (下限 -30dB)
figure(3);
PdB = 20*log10(P+eps);  floorDB = -30;
Pd  = max(PdB, floorDB) - floorDB;             % 平移到 [0, |floor|]
Xd=Pd.*sin(TH).*cos(PH); Yd=Pd.*sin(TH).*sin(PH); Zd=Pd.*cos(TH);
surf(Xd,Yd,Zd,PdB,'EdgeColor','none');
axis equal; colormap(turbo); c=colorbar; c.Label.String='dB';
shading interp; clim([floorDB 0]);
xlabel('x'); ylabel('y'); zlabel('z'); view(135,25);
title(sprintf('UVA 3D 方向图 (dB, 下限-30)\n指向(\\theta_0=%d°,\\phi_0=%d°)', theta_0, phi_0));

% ---- 峰值方向核对 ----
[~, im] = max(P(:));
fprintf('主瓣方向: θ=%.0f°, φ=%.0f°  (期望 %d°, %d°)\n', ...
        rad2deg(TH(im)), rad2deg(PH(im)), theta_0, phi_0);
