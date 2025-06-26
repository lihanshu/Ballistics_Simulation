%% Documentation
%Program writen by: Brian Wade
%Date 9 Sept 2015

%% Description and setup
%This program solves the 12 simultatious ordinary differential equations
%that describe the ballistic motion of a projectile. These equations are
%based on those found in McCoy, 1998 (see citation below).  The input
%aerodynamic coefficeints are from McCoy's book as well.
% 求解描述弹丸弹道运动的 12 个联立常微分方程，同时说明了方程和输入数据的来源，

%Input data from: 
%1. McCoy, RL, Modern Exterior Ballistics: The Launch and Flight Dynamics 
%of Symmetric Projectiles, Schiffer Military History, Atglen, PA, 1998.


%% Program setup
clear
clc
close all
start_time=tic; %Timer

%% Setup
%Declare global variables needed in EoM.m code
global cdo_wpn
global clo_wpn
global cmo_wpn
global cmqao_wpn
global cd2_wpn
global cl2_wpn
global cm2_wpn
global cmqa2_wpn
global std_atm
global d
global It
global m
global R

%Declare global variables needed in EoM.m code.
global gravity
gravity=9.81; %gravity in metric units (m/s^2)

%Standard Measurements
R=6371220; %Radius of earth in meters

%read standard atmosphere tables
std_atm=readmatrix('std_atm.csv','Range','A2:k43');
    
% Weapon charicteristics - These are the properties of the submunissions.
% These weapons charicteristics are from McCoy, 1998.  See citation above.
% 声明在 EoM.m 中会用到的全局变量，对重力加速度、地球半径等标准参数初始化，
% 读取标准大气数据和武器特性数据（包括空气动力学系数）。
% 弹丸特性参数（来源于McCoy 1998）
d = 119.56 / 1000;  % 弹丸直径（米，转换自毫米）
Ip = 0.02335;       % 轴向转动惯量（kg·m²）
It = 0.23187;       % 横向转动惯量（kg·m²）
m = 13.585;         % 弹丸质量（kg）

% 读取气动系数（从Excel文件中读取，用于计算空气阻力、升力等）
cdo_wpn = xlsread('Aerodynamic_Char_120mm_Mortar.xlsx', 'A5:B11');  % 零升阻力系数
cd2_wpn = xlsread('Aerodynamic_Char_120mm_Mortar.xlsx', 'A15:B22'); % 阻力系数修正项
% （其余系数类似：cl为升力系数，cm为力矩系数，后缀2表示与攻角平方相关的项）
clo_wpn=xlsread('Aerodynamic_Char_120mm_Mortar.xlsx','A26:B30');
cl2_wpn=xlsread('Aerodynamic_Char_120mm_Mortar.xlsx','A34:B41');
cmo_wpn=xlsread('Aerodynamic_Char_120mm_Mortar.xlsx','A45:B51');
cm2_wpn=xlsread('Aerodynamic_Char_120mm_Mortar.xlsx','A55:B63');
cmqao_wpn=xlsread('Aerodynamic_Char_120mm_Mortar.xlsx','A67:B72');
cmqa2_wpn=xlsread('Aerodynamic_Char_120mm_Mortar.xlsx','A76:B83');



%% Initial conditions
% 初始速度与角度
Vo_set = 100;       % 初速度（m/s）
phi_0_set = 45;     % 发射仰角（度，向上为正）
theta_0_set = 15;   % 发射方位角（度，向右为正）

% 初始角速度
w_z0_set = 1;       % 初始俯仰角速度（rad/s，抬头为正）
w_y0_set = 0.5;     % 初始偏航角速度（rad/s，左偏为正）

% 初始姿态角
alpha_0_set = 2;    % 初始俯仰角（度）
beta_0_set = -0.5;  % 初始偏航角（度）

% 初始位置（惯性坐标系下，原点为发射点）
x_0 = 0;  % x轴：射程方向
y_0 = 0;  % y轴：高度方向
z_0 = 0;  % z轴：横程方向



%% Run program
t_max = 300;  % 最大仿真时间（秒）

% 复制初始条件到临时变量（避免直接修改原始参数）
Vo = Vo_set; phi_0 = phi_0_set; theta_0 = theta_0_set;
w_z0 = w_z0_set; w_y0 = w_y0_set; alpha_0 = alpha_0_set; 
beta_0 = beta_0_set;
p = 0;  % 初始自旋角速度（rad/s，此处设为0）

%% Intermediate Calcs
% 计算弹丸初始轴向单位向量（弹体坐标系相对于惯性坐标系的指向）
X1o = cosd(phi_0 + alpha_0) * cosd(theta_0 + beta_0);  % x分量
X2o = sind(phi_0 + alpha_0) * cosd(theta_0 + beta_0);  % y分量
X3o = sind(theta_0 + beta_0);                          % z分量

% 计算轴向单位向量的初始变化率（用于后续角速度转换）
Q = (sind(theta_0 + beta_0))^2 + (cosd(theta_0 + beta_0))^2 * (cosd(phi_0 + alpha_0))^2;
dx_1o = (1/sqrt(Q)) * (-w_z0 * (cosd(theta_0 + beta_0))^2 * sind(phi_0 + alpha_0) * cosd(phi_0 + alpha_0) + w_y0 * sind(theta_0 + beta_0));
dx_2o = (1/sqrt(Q)) * (w_z0 * (cosd(theta_0 + beta_0))^2 * (cosd(phi_0 + alpha_0))^2 + w_z0 * (sind(theta_0 + beta_0))^2);
dx_3o = (1/sqrt(Q)) * (-w_z0 * sind(theta_0 + beta_0) * cosd(theta_0 + beta_0) * sind(phi_0 + alpha_0) - w_y0 * cosd(theta_0 + beta_0) * cosd(phi_0 + alpha_0));


% 状态向量x的定义（共12个元素，描述弹丸的完整运动状态）：
% x(1): x方向速度（m/s）；x(2): y方向速度（m/s）；
% x(3): z方向速度（m/s）
% x(4): 滚转角速度（rad/s）；x(5): 俯仰角速度（rad/s）；
% x(6): 偏航角速度（rad/s）
% x(7)-x(9): 轴向单位向量的x、y、z分量（描述姿态）
% x(10): x方向位置（射程，m）；x(11): y方向位置（高度，m）；
% x(12): z方向位置（横程，m）

% 初始化状态向量
x0(1) = Vo * cosd(phi_0) * cosd(theta_0);  % 初始x方向速度
x0(2) = Vo * sind(phi_0) * cosd(theta_0);  % 初始y方向速度
x0(3) = Vo * sind(theta_0);                % 初始z方向速度
x0(4) = (Ip*p/It)*X1o + X2o*dx_3o - X3o*dx_2o;  % 初始滚转角速度
x0(5) = (Ip*p/It)*X2o + X1o*dx_3o + X3o*dx_1o;  % 初始俯仰角速度
x0(6) = (Ip*p/It)*X3o + X1o*dx_2o + X2o*dx_1o;  % 初始偏航角速度
x0(7:9) = [X1o, X2o, X3o];  % 初始轴向单位向量
x0(10:12) = [x_0, y_0, z_0];  % 初始位置

% 求解ODE：调用ode45求解器，积分运动方程EoM.m
tspan = [0 t_max];  % 仿真时间范围
Opt = odeset('Events', @myEvent);  % 设置事件函数（落地时终止仿真）
[t, x] = ode45(@EoM, tspan, x0, Opt);  % t为时间序列，x为状态向量随时间的变化

%% Intermediate Calcs（续）
% 计算姿态角（俯仰角alpha、偏航角beta）
alpha = acosd(x(:,2) ./ sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2));
beta = acosd(x(:,3) ./ sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2));

% 找到弹道最高点（最大高度及对应索引）
[max_ht, I] = max(x(:,11));

% 插值计算落地时间、射程、横程（从最高点到落地的区间内插值，避免发射点干扰）
impact_time = interp1(x(I:end,11), t(I:end), 0);  % 落地时间
range = interp1(x(I:end,11), x(I:end,10), 0);      % 射程（x方向）
cross_range = interp1(x(I:end,11), x(I:end,12), 0);  % 横程（z方向）

% 计算落地时的姿态角、速度向量及总速度
impact_alpha = interp1(x(I:end,11), alpha(I:end), 0);
impact_beta = interp1(x(I:end,11), beta(I:end), 0);
vel_x_imp = interp1(x(I:end,11), x(I:end,1), 0);  % 落地x方向速度
vel_y_imp = interp1(x(I:end,11), x(I:end,2), 0);  % 落地y方向速度
vel_z_imp = interp1(x(I:end,11), x(I:end,3), 0);  % 落地z方向速度
impact_vel = sqrt(vel_x_imp^2 + vel_y_imp^2 + vel_z_imp^2);  % 落地总速度

% 计算落地角度（与竖直方向的夹角）
impactVect = [cosd(impact_beta)*cosd(impact_alpha), sind(impact_beta)*cosd(impact_alpha), sind(impact_alpha)];
vert = [0,1,0];  % 竖直向下向量
impact_angle = acosd(dot(vert, impactVect)/(norm(vert)*norm(impactVect))) - 90;  % 落地角度

% 总飞行距离（地面平面内的射程与横程合成）
total_dis = sqrt(x(:,10).^2 + x(:,12).^2);
total_dis_impact = interp1(x(I:end,11), total_dis(I:end), 0);  % 落地总距离

%% Outputs
% 打印关键结果到命令行
disp(['总飞行距离 = ', num2str(total_dis_impact), ' 米']);
disp(['落地角度 = ', num2str(impact_angle), ' 度']);
disp(['落地速度 = ', num2str(impact_vel), ' m/s']);
disp(['射程（x方向） = ', num2str(range), ' 米']);
disp(['横程（z方向） = ', num2str(cross_range), ' 米']);

% 绘制弹道可视化图
figure();
subplot(1,3,1);  % 3D弹道图（射程-横程-高度）
plot3(x(:,10), x(:,12), x(:,11)); grid on;
xlabel('射程 (m)'); ylabel('横程 (m)'); zlabel('高度 (m)'); axis tight;

subplot(1,3,2);  % 总距离-高度图
plot(total_dis, x(:,11)); grid on;
xlabel('总距离 (m)'); ylabel('高度 (m)');

subplot(1,3,3);  % 射程-横程图（地面轨迹）
plot(x(:,10), x(:,12)); grid on;
xlabel('射程 (m)'); ylabel('横程 (m)');

% 计算程序运行时间
end_time = toc(start_time);

%% Termination function for ODE
function [value, isterminal, direction] = myEvent(t, y)
    value = y(11);          % 检测变量：高度（y方向位置）
    isterminal = 1;         % 1=满足条件时终止仿真
    direction = -1;         % 只检测高度从正变负（落地）的情况
end

    
    