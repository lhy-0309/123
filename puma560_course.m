function puma560_course
% PUMA 560 机器人运动分析与仿真
% 需要 Peter Corke Robotics Toolbox for MATLAB 及其依赖。
% 将本文件所在目录加入 MATLAB 路径后，执行 puma560_course。
% 计算量使用 SI 单位。显示的角度曲线使用度制。
% 构型图、轨迹图、力矩图以及采样 CSV 写入本文件旁的 puma560_results。

if exist('mdl_puma560', 'file') ~= 2 || exist('jtraj', 'file') ~= 2
    error('请先配置 Peter Corke Robotics Toolbox for MATLAB 及其依赖。');
end
mdl_puma560;
qA = [0 -45 90 0 45 0] * pi/180;
qB = [45 -30 60 30 45 20] * pi/180;
qlim = p560.qlim;
outputDir = fullfile(fileparts(mfilename('fullpath')), 'puma560_results');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

% 1. 正运动学与构型。位姿均以基座坐标系表达。
TA = p560.fkine(qA);
TB = p560.fkine(qB);
TA_matrix = matrix4(TA);
TB_matrix = matrix4(TB);
fprintf('\n1. 起始位姿 TA 和目标位姿 TB\n');
format short g
disp(TA_matrix);
disp(TB_matrix);
fig = figure('Name', 'PUMA 560 - Endpoint Configurations', 'Color', 'w', ...
    'Position', [100 100 1100 475]);
ax = subplot(1,2,1, 'Parent', fig);
drawPose(ax, p560, qA, TA_matrix, 'A');
ax = subplot(1,2,2, 'Parent', fig);
drawPose(ax, p560, qB, TB_matrix, 'B');
print(fig, fullfile(outputDir, 'configurations.png'), '-dpng', '-r200');

% 2. 目标位姿的八个解析分支，逐一正解验证及检查关节限位。
branches = {'lun'; 'luf'; 'ldn'; 'ldf'; 'run'; 'ruf'; 'rdn'; 'rdf'};
angles = nan(numel(branches), 6);
poseError = nan(numel(branches), 2);
withinLimits = false(numel(branches), 1);
tol = 1e-9;
fprintf('\n2. 逆运动学解，角度单位为 deg\n');
fprintf('分支      q1        q2        q3        q4        q5        q6     位置误差(m)  旋转误差(F)  限位内\n');
for k = 1:numel(branches)
    qRaw = p560.ikine6s(TB, branches{k});
    if any(~isfinite(qRaw))
        continue
    end
    angles(k,:) = wrapRadians(qRaw);
    errorMatrix = matrix4(p560.fkine(angles(k,:))) - TB_matrix;
    poseError(k,:) = [norm(errorMatrix(1:3,4)), norm(errorMatrix(1:3,1:3), 'fro')];
    % 检查每关节能否通过添加整数周转角落入其限位区间。
    lowerTurns = ceil((qlim(:,1).' - angles(k,:) - tol)/(2*pi));
    upperTurns = floor((qlim(:,2).' - angles(k,:) + tol)/(2*pi));
    withinLimits(k) = all(lowerTurns <= upperTurns);
    fprintf('%s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f   %11.3e   %11.3e   %d\n', ...
        branches{k}, angles(k,:)*180/pi, poseError(k,:), withinLimits(k));
end
valid = all(isfinite(angles), 2) & all(poseError < tol, 2);
assert(any(valid), '没有通过正运动学复核的逆解。');
cost = vecnorm(wrapRadians(angles-qB), 2, 2);
cost(~valid) = inf;
[~, chosen] = min(cost);
fprintf('与题设 qB 最接近的分支为 %s，其正解为\n', branches{chosen});
disp(matrix4(p560.fkine(angles(chosen,:))));
fprintf('满足关节限位的已验证分支数为 %d。\n', nnz(withinLimits & valid));

% 3. 基座雅可比矩阵与给定关节速度下的末端速度。
J = p560.jacob0(qA);
jointVelocity = [0.1; 0.2; 0; 0; 0; 0];
twist = J * jointVelocity;
fprintf('\n3. 起始构型的雅可比矩阵，秩为 %d\n', rank(J));
disp(J);
fprintf('线速度 (m/s)   [% .8f % .8f % .8f]\n', twist(1:3));
fprintf('角速度 (rad/s) [% .8f % .8f % .8f]\n', twist(4:6));

% 4. 关节空间五次插值。静止起停，与第 3 题的瞬时速度独立。
t = (0:0.05:3).';
[q, qd, qdd] = jtraj(qA, qB, t);
assert(size(q,1) == 61 && norm(q(1,:)-qA) < tol && norm(q(end,:)-qB) < tol);
assert(max(max(abs([qd([1 end],:); qdd([1 end],:)]))) < tol);
% 将边界的数值舍入残差置零，避免误触发库仑摩擦。
q([1 end],:) = [qA; qB];
qd([1 end],:) = 0;
qdd([1 end],:) = 0;
outOfLimits = q < (qlim(:,1).'-tol) | q > (qlim(:,2).'+tol);
fprintf('\n4. 超出模型限位的采样点数（按关节）\n');
disp(sum(outOfLimits,1));
fprintf('各关节的采样峰值速度和加速度（绝对值，SI）\n');
disp([max(abs(qd),[],1); max(abs(qdd),[],1)]);
fig = figure('Name', 'PUMA 560 - Joint Trajectory', 'Color', 'w', ...
    'Position', [100 100 920 720]);
curveData = {q, qd, qdd};
curveUnits = {'Angle (deg)', 'Velocity (deg/s)', 'Acceleration (deg/s^2)'};
curveTitles = {'Joint positions', 'Joint velocities', 'Joint accelerations'};
for k = 1:3
    ax = subplot(3,1,k, 'Parent', fig);
    drawCurves(ax, t, curveData{k}*180/pi, 1:6, curveUnits{k}, curveTitles{k}, k == 1);
end
xlabel(ax, 'Time (s)');
print(fig, fullfile(outputDir, 'joint_trajectory.png'), '-dpng', '-r200');

% 5. 关节端逆动力学。使用默认摩擦及电机惯量，无外部负载。
% rne 中的重力参数为静力补偿所用的基座向量 [0 0 9.81]。
g = [0 0 9.81];
tau = p560.rne(q, qd, qdd, 'gravity', g, 'slow');
tauGravity = p560.rne(q, zeros(size(q)), zeros(size(q)), 'gravity', g, 'slow');
assert(all(isfinite(tau(:))), '力矩存在无效数值。');
assert(max(max(abs(tau([1 end],:)-tauGravity([1 end],:)))) < tol);
[~, peakRow] = max(abs(tau),[],1);
peakIndex = sub2ind(size(tau), peakRow, 1:6);
fprintf('\n5. 驱动力矩 (N*m) 及其峰值采样时刻\n');
fprintf('关节        起始        终点         最小         最大       峰值带符号      时刻(s)\n');
for j = 1:6
    fprintf('%2d   %11.6f %11.6f %11.6f %11.6f %11.6f    %5.2f\n', ...
        j, tau(1,j), tau(end,j), min(tau(:,j)), max(tau(:,j)), tau(peakIndex(j)), t(peakRow(j)));
end
fig = figure('Name', 'PUMA 560 - Inverse Dynamics', 'Color', 'w', ...
    'Position', [130 130 920 580]);
ax = subplot(2,1,1, 'Parent', fig);
drawCurves(ax, t, tau, 1:3, 'Torque (N m)', 'Base, shoulder and elbow', true);
ax = subplot(2,1,2, 'Parent', fig);
drawCurves(ax, t, tau, 4:6, 'Torque (N m)', 'Wrist joints', true);
xlabel(ax, 'Time (s)');
print(fig, fullfile(outputDir, 'joint_torque.png'), '-dpng', '-r200');

headers = {'t_s'};
quantities = {'q_rad', 'qd_rad_s', 'qdd_rad_s2', 'tau_Nm', 'tau_gravity_Nm'};
for k = 1:numel(quantities)
    for j = 1:6
        headers{end+1} = sprintf('%s_%d', quantities{k}, j); %#ok<AGROW>
    end
end
writetable(array2table([t,q,qd,qdd,tau,tauGravity], 'VariableNames', headers), ...
    fullfile(outputDir, 'trajectory_and_torque.csv'));
save(fullfile(outputDir, 'matlab_kinematics.mat'), 'TA_matrix', 'TB_matrix', ...
    'angles', 'poseError', 'withinLimits', 'branches', 'chosen', 'J', 'twist', 'qlim', '-v7');
fprintf('\n输出目录 %s\n', outputDir);
end

function T = matrix4(pose)
% 在矩阵与 SE3 两种位姿返回值之间统一输出。
if isnumeric(pose)
    T = pose;
else
    T = pose.T;
end
end

function q = wrapRadians(q)
q = mod(q + pi, 2*pi) - pi;
end

function drawPose(ax, robot, q, T, label)
% 用模型各连杆坐标系的原点绘制骨架和末端坐标轴。
F = matrix4(robot.base);
points = zeros(3,7);
points(:,1) = F(1:3,4);
for j = 1:6
    F = F * matrix4(robot.links(j).A(q(j)));
    points(:,j+1) = F(1:3,4);
end
hold(ax, 'on');
plot3(ax, points(1,:), points(2,:), points(3,:), '-o', ...
    'Color', [0.08 0.34 0.60], 'LineWidth', 2.5, 'MarkerFaceColor', [0.08 0.34 0.60]);
origin = T(1:3,4);
axisColors = [0.74 0.25 0.15; 0.16 0.50 0.36; 0.45 0.33 0.64];
axisLabels = 'xyz';
for j = 1:3
    u = 0.095*T(1:3,j);
    quiver3(ax, origin(1),origin(2),origin(3), u(1),u(2),u(3), 0, ...
        'Color', axisColors(j,:), 'LineWidth', 1.8, 'MaxHeadSize', 0.5);
    p = origin+1.18*u;
    text(ax, p(1),p(2),p(3), axisLabels(j), 'Color', axisColors(j,:));
end
axis(ax, 'equal');
xlim(ax, [-0.20 0.50]); ylim(ax, [-0.34 0.36]); zlim(ax, [-0.41 0.29]);
view(ax, -54, 24);
grid(ax, 'on');
xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)'); zlabel(ax, 'Z (m)');
title(ax, ['Configuration ', label]);
end

function drawCurves(ax, t, data, joints, units, caption, showLegend)
colors = [0.082 0.345 0.600; 0.765 0.310 0.125; 0.204 0.502 0.380; ...
          0.459 0.271 0.600; 0.635 0.404 0.000; 0.075 0.471 0.541];
styles = {'-', '--', '-.', '-', '--', '-.'};
hold(ax, 'on');
for j = joints
    plot(ax, t, data(:,j), 'Color', colors(j,:), 'LineStyle', styles{j}, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Joint %d', j));
end
xlim(ax, [t(1), t(end)]);
grid(ax, 'on');
ylabel(ax, units);
title(ax, caption);
if showLegend
    legend(ax, 'show', 'Location', 'best');
end
end
