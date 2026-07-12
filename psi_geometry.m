% ===================================================
% psi_geometry.m
% 几何示意：远场波程差 = 观察方向 r̂ 在阵元位置 p 上的投影 (r̂·p)
% 用于理解 ψx = k·dx·sinθcosφ+αx , ψy = k·dy·sinθsinφ+αy 的来历
% ===================================================
clear; clc; close all;

% ---- 参数 ----
Nx = 4; Ny = 4;          % 阵元数
d  = 0.5;                % 间距 d/λ
theta = 50;              % 观察方向极角(从 z 轴量)
phi   = 35;              % 观察方向方位角(从 x 轴量)

% ---- 观察方向单位矢量 r̂ ----
rhat = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];

% ---- 阵元坐标(xy 平面) ----
[mx, ny] = meshgrid(0:Nx-1, 0:Ny-1);
Px = mx(:)*d;  Py = ny(:)*d;  Pz = zeros(size(Px));

% ---- 选一个代表阵元 p 演示投影 ----
p = [3*d; 2*d; 0];                 % (m=3, n=2)
projLen = dot(rhat, p);            % r̂·p = 波程差
foot    = projLen * rhat;          % p 在 r̂ 方向上的投影点

figure('Position',[80 80 1150 520],'Color','w');

% ================= 左图：3D 投影关系 =================
subplot(1,2,1); hold on; grid on; axis equal;
% 阵元
plot3(Px,Py,Pz,'ko','MarkerFaceColor',[.7 .7 .7],'MarkerSize',6);
plot3(0,0,0,'ks','MarkerFaceColor','k','MarkerSize',9);           % 参考原点
plot3(p(1),p(2),p(3),'ro','MarkerFaceColor','r','MarkerSize',9);  % 代表阵元 p

L = 3.0;                                                          % r̂ 画长一点
quiver3(0,0,0, rhat(1)*L, rhat(2)*L, rhat(3)*L, 0,'b','LineWidth',2,'MaxHeadSize',0.3);
text(rhat(1)*L, rhat(2)*L, rhat(3)*L,'  $\hat{r}$ (观察方向)','Interpreter','latex','Color','b','FontSize',12);

% 位置矢量 p
quiver3(0,0,0, p(1),p(2),p(3), 0,'r','LineWidth',1.5,'MaxHeadSize',0.5);
text(p(1),p(2),p(3),'  p','Color','r','FontSize',12);

% 投影段(原点→投影点)= 波程差 r̂·p
plot3([0 foot(1)],[0 foot(2)],[0 foot(3)],'g','LineWidth',3);
% 从 p 到投影点的垂线(虚线)
plot3([p(1) foot(1)],[p(2) foot(2)],[p(3) foot(3)],'k--','LineWidth',1);

text(foot(1)/2,foot(2)/2,foot(3)/2+0.15, ...
     sprintf('  r\\cdotp = %.3f\\lambda',projLen),'Color',[0 .5 0],'FontSize',12,'FontWeight','bold');

xlabel('x/\lambda'); ylabel('y/\lambda'); zlabel('z/\lambda');
title('波程差 = r̂ 在 p 上的投影 (r̂·p)');
view(35,20); xlim([-.5 3]); ylim([-.5 3]); zlim([0 2.5]);

% ================= 右图：两条平行波前 =================
subplot(1,2,2); hold on; grid on; axis equal;
plot3(Px,Py,Pz,'ko','MarkerFaceColor',[.7 .7 .7],'MarkerSize',6);
plot3(0,0,0,'ks','MarkerFaceColor','k','MarkerSize',9);
plot3(p(1),p(2),p(3),'ro','MarkerFaceColor','r','MarkerSize',9);

% 构造垂直于 r̂ 的两个基矢量,用于画波前平面
tmp = [0;0;1]; if abs(dot(tmp,rhat))>0.95, tmp=[0;1;0]; end
e1 = cross(rhat,tmp); e1=e1/norm(e1);
e2 = cross(rhat,e1);
s = 1.6;
planeAt = @(c,col) patch(c(1)+s*[-e1(1)-e2(1), e1(1)-e2(1), e1(1)+e2(1), -e1(1)+e2(1)], ...
                          c(2)+s*[-e1(2)-e2(2), e1(2)-e2(2), e1(2)+e2(2), -e1(2)+e2(2)], ...
                          c(3)+s*[-e1(3)-e2(3), e1(3)-e2(3), e1(3)+e2(3), -e1(3)+e2(3)], ...
                          col,'FaceAlpha',0.15,'EdgeColor',col);
planeAt([0;0;0],'b');       % 过原点的波前
planeAt(foot,'r');          % 过阵元 p 的波前(平移了 r̂·p)

% 平行射线
quiver3(0,0,0, rhat(1)*2.5,rhat(2)*2.5,rhat(3)*2.5,0,'b','LineWidth',1.2,'MaxHeadSize',0.3);
quiver3(p(1),p(2),p(3), rhat(1)*2.5,rhat(2)*2.5,rhat(3)*2.5,0,'r','LineWidth',1.2,'MaxHeadSize',0.3);
plot3([0 foot(1)],[0 foot(2)],[0 foot(3)],'g','LineWidth',3);

xlabel('x/\lambda'); ylabel('y/\lambda'); zlabel('z/\lambda');
title(sprintf('两平行波前之间距 = r̂·p = %.3f\\lambda\n\\phi(相位差)=k·r̂·p=%.1f°',projLen,rad2deg(2*pi*projLen)));
view(35,20); xlim([-.5 3]); ylim([-.5 3]); zlim([0 2.5]);

sgtitle(sprintf('远场波程差几何 (\\theta=%d°, \\phi=%d°)：r̂·p = m·d·sin\\theta cos\\phi + n·d·sin\\theta sin\\phi',theta,phi),...
        'FontSize',13,'FontWeight','bold');

% ---- 数值核对：投影 vs 公式 ----
formula = p(1)*sind(theta)*cosd(phi) + p(2)*sind(theta)*sind(phi);
fprintf('r̂·p(点积)      = %.4f λ\n', projLen);
fprintf('公式 m·d·u+n·d·v = %.4f λ  (应相等)\n', formula);
fprintf('=== 几何示意图已生成 ===\n');
