# 6-DOF Ballistics Simulation

> MATLAB Simulation of a six degree of freedom (6-DOF) ballistic flight of a 120mm mortar round.

![Sample Output](/images/output.png)

![Sample Output](/images/console_output.png)

## Files

The simulation is setup, run, and controlled from the Mortar_Sim.m file. This main file uses the ODE45 solver in MATLAB with the equations of motion (EoM) specified for each time step in the EoM.m file. The file reads input data from a standard atmosphere table and a table of aerodynamic coefficients vs. Mach number (McCoy, 1998). The 6-DOF model was developed based on the methodology in McCoy, chapter 9.

1. [Mortar_Sim.m](Mortar_Sim.m) - Main program
2. [EoM.m](EoM.m) - Equations of Motion
3. [std_atm.csv](std_atm.csv) - Table of the standard atmosphere
4. [Aerodynamic_Char_120mm_Mortar.xlsx](Aerodynamic_Char_120mm_Mortar.xlsx) - Aerodynamic coefficients of the 120mm for different Mach number. Data extracted from McCoy, 1998 on page 220.

## Setup and Run

The simulation is setup from within the main program ([Mortar_Sim.m](Mortar_Sim.m)). The simulation parameters are set in lines 68-81. These lines set the initial conditions of the mortar round as it exits the mortar tube and enters free flight.

![initial_conditions](/images/initial_conditions_set.png)

The first variable (Vo_set) is the initial velocity at the barrel exit in meters per second. The next variable (phi_0_set) is the vertical angle of departure with respect to the inertial frame in degrees (this is the angle of the mortar tub measured from the horizontal ground plane). The third variable (theta_0_set) is the horizontal angle of departure in degrees (this is the angle of the mortar tub measured from a vertical plan perpendicular to the ground and along the x-axis).

The next two variables (w_z_0_set and w_y0_set) are the initial pitch and yaw of the mortar round as it exists the barrel both in radians per second. The sign convention is positive for a nose-up or left-yaw rate.

The third set of variables (alpha_0_set and beta_0_set) are the pitch and yaw angles of the round at the barrel exit. These are measures of the angular difference between the mortar body axis and the inertial (gun tube) axis. Both angles are in degrees. These use the same sign convention as above where a nose-up or left-yaw are positive. Note that these angles establish an initial off-axis flight of the round.

The final set of variables (x_0, y_0, and z_0) are the x, y, and z locations of the end of the mortar barrel (where the round enters free flight). These are all measured in meters.

## MATLAB version and add-ons

This simulation was built using MATLAB 2019b. No other MATLAB add ons should be needed to run the simulation.

## References

1. McCoy, RL, Modern Exterior Ballistics: The Launch and Flight Dynamics of Symmetric Projectiles, Schiffer Military History, Atglen, PA, 1998.

# 120 毫米迫击炮弹丸六自由度（6-DOF）弹道飞行的 MATLAB 仿真

仿真通过 Mortar_Sim.m 文件进行设置、运行和控制。主程序调用 MATLAB 的 ODE45 求解器，结合 EoM.m 文件中定义的逐时间步长运动方程（EoM）进行计算。仿真读取标准大气数据表和不同马赫数下的空气动力系数表（McCoy, 1998），其六自由度模型基于 McCoy 第 9 章的方法构建。

Mortar_Sim.m - 主程序

EoM.m - 运动方程

std_atm.csv - 标准大气数据表

Aerodynamic_Char_120mm_Mortar.xlsx - 120 毫米迫击炮弹在不同马赫数下的空气动力系数表（数据源自 McCoy, 1998 第 220 页）

设置与运行
仿真参数在主程序 Mortar_Sim.m 中设置（第 68-81 行），用于定义迫击炮弹丸出膛进入自由飞行时的初始条件。


初始条件设置截图

Vo_set：炮口初始速度（米 / 秒）
phi_0_set：相对于惯性坐标系的垂直射角（度）（即迫击炮身与水平面的夹角）
theta_0_set：水平射角（度）（即迫击炮身与垂直于地面并沿 x 轴的垂直平面的夹角）

w_z_0_set 和 w_y0_set：弹丸出膛时的初始俯仰角速度和偏航角速度（弧度 / 秒）（抬头或左偏航为正）

alpha_0_set 和 beta_0_set：弹丸出膛时的俯仰角和偏航角（度）（弹体轴线与惯性坐标系（炮管）轴线的夹角，抬头或左偏航为正，用于设定初始离轴飞行状态）

x_0, y_0, z_0：迫击炮管末端（弹丸进入自由飞行点）的三维坐标（米）

MATLAB 版本与附加组件:本仿真基于 MATLAB 2019b 开发，运行无需额外附加组件。

参考文献:McCoy, RL，《现代外弹道学：对称弹丸的发射与飞行动力学》，Schiffer Military History，Atglen, PA，1998。