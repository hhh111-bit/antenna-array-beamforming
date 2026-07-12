% =====================================================================
% ula_polar_pattern.m  —  均匀线阵(ULA) 方向图
% ---------------------------------------------------------------------
% 标准阵因子 (sin 函数形式, 阵列沿 x 轴, 相位中心在原点):
%   ψ  = kd·sinθ + α,   α = -kd·sinθ_0     (渐进馈电相位)
%   AF = sin(N·ψ/2) / ( N·sin(ψ/2) )        归一化, 峰值=1
%   θ  = 从阵列法向(broadside)量起;  θ_0 = 波束指向角
% 改 N / d_lambda / theta_0 即可扫描任意指向。
% =====================================================================
clear; clc; close all;

N        = 12;          % 阵元数
d_lambda = 0.5;         % 阵元间距 / λ
theta_0  = 0;          % 波束指向角 (度, 相对法向)
kd       = 2*pi*d_lambda;

% ---- 扫描角 (可见区) ----
theta = -90:0.25:90;
th    = deg2rad(theta);
th0   = deg2rad(theta_0);

% ---- 渐进馈电相位 (使主瓣指向 θ_0) ----
alpha = -kd*sin(th0);

% ---- 标准阵因子 (sin 形式) ----
psi = kd*sin(th) + alpha;
AF  = sin(N*psi/2) ./ (N*sin(psi/2));
AF(abs(sin(psi/2)) < 1e-9) = 1;            % 0/0 主瓣极限
P   = abs(AF);                             % 归一化幅度, 峰值=1
PdB = 20*log10(P + eps);

% (1) 阵列几何 (阵元沿 x 轴, 关于原点对称)
figure(1);
xn = ((0:N-1) - (N-1)/2) * d_lambda;               % 相位中心在原点
plot(xn, zeros(1,N), '-', 'Color',[.6 .6 .6], 'LineWidth',1, ...
     'HandleVisibility','off'); hold on;
plot(xn, zeros(1,N), 'o', 'MarkerFaceColor',[.15 .45 .85], ...
     'MarkerEdgeColor','k', 'MarkerSize',10, 'LineWidth',1);
plot(0, 0, 'rp', 'MarkerFaceColor','r', 'MarkerSize',13);     % 原点=相位中心
axis equal; box off;
xlim([xn(1)-d_lambda, xn(end)+d_lambda]);  ylim([-0.8 0.8]);
set(gca, 'YTick', [], 'YColor', 'none');           % 一维阵: 隐藏 y 轴
xlabel('x / \lambda');
title(sprintf('阵列几何 — ULA  N=%d, d=%.2f\\lambda (对称于原点)', N, d_lambda));
legend('阵元','原点/相位中心', 'Location','southoutside', 'Orientation','horizontal');

% (2) 极坐标方向图
figure(2);
polarplot(th, P, 'b', 'LineWidth', 1.5); hold on;
polarplot([th0 th0], [0 1], 'r--', 'LineWidth', 1.5);   % 指向线
ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise';
thetalim([-90 90]); rlim([0 1]);
title(sprintf('ULA 极坐标方向图  指向 \\theta_0=%d°', theta_0));

% (3) 直角坐标 (dB, 下限 -40)
figure(3);
plot(theta, PdB, 'b', 'LineWidth', 1.5); hold on;
xline(theta_0, 'r--', 'LineWidth', 1.5);
xlabel('\theta (°)'); ylabel('归一化增益 (dB)');
xlim([-90 90]); ylim([-40 3]); grid on;
title(sprintf('ULA 直角坐标方向图 (dB)  指向 \\theta_0=%d°', theta_0));

% ---- 峰值方向核对 ----
[~, im] = max(P);
fprintf('主瓣方向: θ=%.1f°  (期望 %d°)\n', theta(im), theta_0);
