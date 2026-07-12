% =====================================================================
% upa_3d_pattern.m  —  均匀平面阵(UPA) 三维方向图
% ---------------------------------------------------------------------
% 标准阵因子 (sin 函数形式, 阵列在 xy 平面, 法向 +z):
%   ψx = kdx·sinθ·cosφ + αx,   αx = -kdx·sinθ_0·cosφ_0
%   ψy = kdy·sinθ·sinφ + αy,   αy = -kdy·sinθ_0·sinφ_0
%   AF = | sin(Nx·ψx/2)/(Nx·sin(ψx/2)) · sin(Ny·ψy/2)/(Ny·sin(ψy/2)) |
%   (θ_0,φ_0) = 波束指向 (仰角, 方位)。αx,αy 为渐进馈电相位。
% 改 Nx,Ny / d / theta_0,phi_0 即可扫描任意指向。
% =====================================================================
clear; clc; close all;

Nx = 8;  Ny = 8;                      % 阵元数 (x×y)
dx_lambda = 0.5;  dy_lambda = 0.5;    % 间距 / λ
theta_0 = 0;   phi_0 = 0;           % 波束指向 (仰角, 方位, 度)

kdx = 2*pi*dx_lambda;  kdy = 2*pi*dy_lambda;
th0 = deg2rad(theta_0);  ph0 = deg2rad(phi_0);

% ---- 渐进馈电相位 (使主瓣指向 θ_0,φ_0) ----
ax_ = -kdx*sin(th0)*cos(ph0);
ay_ = -kdy*sin(th0)*sin(ph0);

% ---- 球面网格 (上半球) ----
th_v = deg2rad(0:1:180);
ph_v = deg2rad(0:1:360);
[TH, PH] = meshgrid(th_v, ph_v);

% ---- 标准阵因子 (sin 形式) ----
psi_x = kdx*sin(TH).*cos(PH) + ax_;
psi_y = kdy*sin(TH).*sin(PH) + ay_;
AF_x = sin(Nx*psi_x/2)./(Nx*sin(psi_x/2));  AF_x(abs(sin(psi_x/2))<1e-9)=1;
AF_y = sin(Ny*psi_y/2)./(Ny*sin(psi_y/2));  AF_y(abs(sin(psi_y/2))<1e-9)=1;
P = abs(AF_x .* AF_y);
P = P / max(P(:));               % 归一化, 峰值=1

X = P.*sin(TH).*cos(PH);
Y = P.*sin(TH).*sin(PH);
Z = P.*cos(TH);

% (1) 阵列几何 (俯视, 关于原点对称)
figure(1);
xv = ((0:Nx-1) - (Nx-1)/2) * dx_lambda;            % 相位中心在原点
yv = ((0:Ny-1) - (Ny-1)/2) * dy_lambda;
[xg, yg] = meshgrid(xv, yv);
plot(xg(:), yg(:), 'o', 'MarkerFaceColor',[.15 .45 .85], ...
     'MarkerEdgeColor','k', 'MarkerSize',9, 'LineWidth',1); hold on;
plot(0, 0, 'rp', 'MarkerFaceColor','r', 'MarkerSize',13);     % 原点=相位中心
axis equal; grid on; box on;
xlim([xv(1)-dx_lambda, xv(end)+dx_lambda]);
ylim([yv(1)-dy_lambda, yv(end)+dy_lambda]);
xlabel('x / \lambda'); ylabel('y / \lambda');
title(sprintf('阵列几何 (俯视)  UPA %d×%d, d=%.2f\\lambda (对称于原点)', Nx, Ny, dx_lambda));

% (2) 线性幅度 3D 曲面
figure(2);
surf(X,Y,Z,P,'EdgeColor','none'); hold on;
xa=sin(th0)*cos(ph0); ya=sin(th0)*sin(ph0); za=cos(th0);   % 指向箭头
plot3([0 xa],[0 ya],[0 za],'r-','LineWidth',2);
axis equal; colormap(turbo); colorbar; shading interp;
xlabel('x'); ylabel('y'); zlabel('z'); view(135,25);
title(sprintf('UPA 3D 方向图 (线性)\n%d×%d, d=%.2f\\lambda, 指向(\\theta_0=%d°,\\phi_0=%d°)', ...
      Nx,Ny,dx_lambda,theta_0,phi_0));

% (3) dB 曲面 (下限 -30dB)
figure(3);
PdB = 20*log10(P+eps);  floorDB = -30;
Pd  = max(PdB, floorDB) - floorDB;             % 平移到 [0, |floor|]
Xd=Pd.*sin(TH).*cos(PH); Yd=Pd.*sin(TH).*sin(PH); Zd=Pd.*cos(TH);
surf(Xd,Yd,Zd,PdB,'EdgeColor','none');
axis equal; colormap(turbo); c=colorbar; c.Label.String='dB';
shading interp; clim([floorDB 0]);
xlabel('x'); ylabel('y'); zlabel('z'); view(135,25);
title(sprintf('UPA 3D 方向图 (dB, 下限-30)\n指向(\\theta_0=%d°,\\phi_0=%d°)', theta_0, phi_0));

% ---- 峰值方向核对 ----
[~, im] = max(P(:));
fprintf('主瓣方向: θ=%.0f°, φ=%.0f°  (期望 %d°, %d°)\n', ...
        rad2deg(TH(im)), rad2deg(PH(im)), theta_0, phi_0);
