function dx=EoM(t,x)
%% Documentation
%Program writen by: Brian Wade
%Date 9 Sept 2015

%% Description and setup
%This program runs as a subprogram of the TBM_Flight.m program.  This
%function solves the 12 simultatious ordinary differential equations that
%describe the ballistic motion of a projectile. These equations are based 
%on the equations found in McCoy, 1998 (see citation below).  The input
%aerodynamic coefficeints are from McCoy's book as well.

% EoM.m 是弹道仿真程序的核心函数，用于求解描述弹丸运动的12 个常微分方程。这些方程基于牛顿力学和刚体动力学，
% 考虑了空气阻力、升力、重力和力矩的影响，完整描述了弹丸在三维空间中的运动状态。
% 输入参数
% t：当前时间
% x：状态向量（12 维），包含速度、角速度、姿态和位置信息
% 输出参数
% dx：状态向量的导数（12 维），即每个状态变量的变化率

%Input data from: 
%1. McCoy, RL, Modern Exerioor Balistics: The Launch and Flight Dynamics 
%of Symmetric Projectiles, Schiffer Military History, Atglen, PA, 1998.

%% Program
% 获取主程序中的全局变量
global cdo_wpn clo_wpn cmo_wpn cmqao_wpn cd2_wpn cl2_wpn cm2_wpn cmqa2_wpn ...
       std_atm d It m R gravity

%Find atmospheric variables
rho=interp1(std_atm(:,1),std_atm(:,7),x(11)/1000); %atmpospheric density in
%kg/m^3
% rho：根据弹丸当前高度 x(11) 对标准大气数据 std_atm 进行线性插值，得到当前高度的大气密度。
%
% a：同样根据弹丸当前高度进行线性插值，得到当前高度的声速。

a=interp1(std_atm(:,1),std_atm(:,8),x(11)/1000); %speed of sound in m/s

%Find refrence area of projectile.
% 根据弹丸直径 d 计算弹丸的参考面积 S。
S=(pi/4)*d^2; %refrence area in m^2

%% Body Forces
% 计算弹丸质心总速度大小（惯性坐标系下x、y、z方向速度的合速度）
V = sqrt(x(1)^2 + x(2)^2 + x(3)^2);  % V单位：m/s

% 计算总攻角（alpha_t）：弹丸轴向与速度方向的夹角
% x(7)-x(9)是弹丸轴向单位向量，点乘速度向量再除以V得到夹角余弦值
alpha_t = acos((x(1)*x(7) + x(2)*x(8) + x(3)*x(9))/V);  % 单位：rad

% 计算马赫数（气流速度与音速的比值，反映压缩性影响）
mach = V/a;  % 无量纲

% 计算动压系数（q）：空气动力的关键参数，与密度、速度平方、参考面积相关
q = rho*V*S;  % 单位：kg/(m·s)（后续与气动系数结合得到力的大小）

% 限制马赫数范围：若当前马赫数超过气动系数表的最大值，取表中最大马赫数
if mach > cdo_wpn(size(cdo_wpn,1),1)
    M = cdo_wpn(size(cdo_wpn,1),1);  % 表中最大马赫数
else
    M = mach;  % 当前马赫数
end

% 根据马赫数M插值获取基础气动系数（从Excel表中读取的离散数据）
cdo = interp1(cdo_wpn(:,1), cdo_wpn(:,2), M);  % 零升阻力系数（马赫数相关）
cd2 = interp1(cd2_wpn(:,1), cd2_wpn(:,2), M);  % 阻力系数的攻角修正项
clo = interp1(clo_wpn(:,1), clo_wpn(:,2), M);  % 零攻角升力系数
cl2 = interp1(cl2_wpn(:,1), cl2_wpn(:,2), M);  % 升力系数的攻角修正项
cmo = interp1(cmo_wpn(:,1), cmo_wpn(:,2), M);  % 零攻角俯仰力矩系数
cm2 = interp1(cm2_wpn(:,1), cm2_wpn(:,2), M);  % 力矩系数的攻角修正项
cmqao = interp1(cmqao_wpn(:,1), cmqao_wpn(:,2), M);  % 零攻角翻转力矩系数
cmqa2 = interp1(cmqa2_wpn(:,1), cmqa2_wpn(:,2), M);  % 翻转力矩系数的攻角修正项

% 综合气动系数（加入攻角影响，sin(alpha_t)^2体现非线性关系）
cd = cdo + cd2*(sin(alpha_t))^2;  % 总阻力系数
cl = clo + cl2*(sin(alpha_t))^2;  % 总升力系数
cm = cmo + cm2*(sin(alpha_t))^2;  % 总俯仰力矩系数
cmqa = cmqao + cmqa2*(sin(alpha_t))^2;  % 总翻转力矩系数

% 将气动系数转换为单位质量/转动惯量的力或力矩（便于直接代入运动方程）
Cd = (q*cd)/(2*m);  % 单位质量的阻力加速度（m/s²）
Cl = (q*cl)/(2*m);  % 单位质量的升力加速度（m/s²）
Cm = (q*d*cm)/(2*It);  % 单位转动惯量的俯仰角加速度（rad/s²）
Cmq = (q*d^2*cmqa)/(2*It);  % 单位转动惯量的翻转角加速度（rad/s²）

% 计算惯性坐标系下的重力分量（修正地球曲率影响）
g1 = -gravity * x(10)/R;  % x方向重力分量（因地球曲率，指向地心的分力）
g2 = -gravity * (1 - 2*x(11)/R);  % y方向重力分量（高度越高，重力越小）
g3 = 0;  % z方向无重力分量（假设对称）
% 注：注释掉的部分是简化模型（仅y方向有重力），实际使用修正后的分量
% 中间变量：弹体自旋角速度在轴向的投影（用于修正角加速度计算）
IpP_It = x(4)*x(7) + x(5)*x(8) + x(6)*x(9);  % 无量纲

% 以下为12个状态变量的导数（dx(i,1)即dx/dt）
% ---------------------------
% 1-3：速度分量的导数（加速度）
dx(1,1) = -Cd*x(1) + Cl*(V^2*x(7) - V*x(1)*cos(alpha_t)) + g1;
% 物理意义：x方向加速度 = 阻力加速度（-Cd*x(1)） + 升力加速度（Cl*...） + x方向重力分量（g1）
% 升力项解析：(V²x7 - Vx1 cosα)是升力方向的速度相关项，体现升力垂直于速度的特性

dx(2,1) = -Cd*x(2) + Cl*(V^2*x(8) - V*x(2)*cos(alpha_t)) + g2;  % y方向加速度
dx(3,1) = -Cd*x(3) + Cl*(V^2*x(9) - V*x(3)*cos(alpha_t)) + g3;  % z方向加速度

% ---------------------------
% 4-6：角速度分量的导数（角加速度）
dx(4,1) = Cm*(x(2)*x(9) - x(3)*x(8)) + Cmq*(x(4) - IpP_It*x(7));
% 物理意义：x方向角加速度 = 俯仰力矩产生的角加速度（Cm*...） + 翻转力矩产生的角加速度（Cmq*...）
% 力矩项解析：(x2x9 - x3x8)是速度与姿态向量的叉乘项，体现力矩方向与姿态的关系

dx(5,1) = Cm*(x(3)*x(7) - x(1)*x(9)) + Cmq*(x(5) - IpP_It*x(8));  % y方向角加速度
dx(6,1) = Cm*(x(1)*x(8) - x(2)*x(7)) + Cmq*(x(6) - IpP_It*x(9));  % z方向角加速度

% ---------------------------
% 7-9：轴向单位向量的导数（姿态变化率）
dx(7,1) = x(5)*x(9) - x(6)*x(8);
% 物理意义：姿态向量x分量的变化率 = 俯仰角速度与z分量的叉乘 - 偏航角速度与y分量的叉乘
% 本质：刚体姿态运动学方程（单位向量随角速度的旋转关系）

dx(8,1) = x(6)*x(7) - x(4)*x(9);  % 姿态向量y分量的变化率
dx(9,1) = x(4)*x(8) - x(5)*x(7);  % 姿态向量z分量的变化率

% ---------------------------
% 10-12：位置分量的导数（速度）
dx(10,1) = x(1);  % x方向位置变化率 = x方向速度
dx(11,1) = x(2) + (x(1)^2)/(2*R);  % y方向位置变化率 = y方向速度 + 地球曲率修正项
dx(12,1) = x(3);  % z方向位置变化率 = z方向速度

end
