clc;
clear;

%% 1. 干涉图生成
% --------------------------
% 说明：
%   该部分加载雷达图像数据，并根据主图像与从图像计算干涉图，
%   继而提取干涉相位，并以伪彩色图显示原始干涉图相位。
%
% 数据加载
load('C:\Users\栗浩宇\Desktop\Insar滤波\ers_vesuvius.mat');

% 从加载的数据中提取主图像（im_m）与从图像（im_s）
masterImage = im_m;
slaveImage  = im_s;

% 计算复数形式的干涉图：主图像乘以从图像的共轭
interferogram = masterImage .* conj(slaveImage);

% 提取干涉相位（以弧度表示）
interferogram_phase = angle(interferogram);

% 显示原始干涉图相位
figure;
imagesc(interferogram_phase);
colormap('jet');
colorbar;
title('一. 原始干涉图');

%% 2. 去平地效应（消除平地效应）
% --------------------------
% 说明：
%   为消除因卫星几何结构产生的平地效应，
%   本节计算每个像素对应的距离变化，并由此生成相位修正，
%   再将其应用于原始干涉相位，确保结果相位位于 [0, 2π) 区间。

% 计算卫星至地面平均距离（考虑入射角23°）
averageDistance = 780000 / cos(23 * pi / 180);

% 计算单个像素对应的距离变化（依据雷达波在真空中传播速度是光速299792458m/s及采样频率18.93MHz）
distancePerPixel = 299792458 / (2 * 18930000);

% 获取图像的列数
numberOfColumns = size(im_m, 2);

% 初始化距离修正数组（用于存储每列的距离偏移）
distanceCorrection = zeros(1, numberOfColumns);

% 对图像每一列计算距离修正值：
%   - 左半部分为正方向修正
%   - 右半部分为负方向修正
for i = 1:(numberOfColumns / 2)
    distanceCorrection(1, numberOfColumns/2 + 1 - i) = distancePerPixel * i;  % 正向修正
    distanceCorrection(1, numberOfColumns/2 + i)     = -distancePerPixel * i; % 负向修正
end
% 将一维距离修正数组扩展为与图像大小一致的二维矩阵
distanceCorrection = repmat(distanceCorrection, size(im_m, 1), 1);

% 计算平地效应相位修正：
%   使用参数：波长 = 0.0566m，基线 = 251m，入射角 = 23°
flatPhaseCorrection = (4 * pi * 251 * distanceCorrection) ...
                      / (0.0566 * averageDistance * tan(23 * pi / 180));

% 应用相位修正，获得校正后的干涉相位
correctedInterferogramPhase = interferogram_phase - flatPhaseCorrection;
% 确保相位值落在 [0, 2π) 区间内
correctedInterferogramPhase = mod(correctedInterferogramPhase, 2 * pi);

% 显示去平地效应后的干涉图相位
figure;
imagesc(correctedInterferogramPhase);
colormap('jet');
colorbar;
title('二. 去平地效应图');

% 组合幅值与校正后的相位，构造新的复数干涉图
correctedInterferogram = abs(interferogram) .* exp(1i * correctedInterferogramPhase);

%% 3. Goldstein 滤波处理
% --------------------------
% 说明：
%   使用Goldstein滤波方法对校正后的干涉图进行相位滤波，
%   以减小噪声并改善相位连续性。
%
% 设置滤波参数：
alpha       = 0.5;  % Goldstein滤波的幂指数参数
window_size = 5;    % 滑动窗口的尺寸
step_size   = 1;    % 滑动窗口的步长

% 调用 Goldstein 滤波函数，返回滤波后的复数干涉图
filteredInterferogram_Goldstein = goldstein_filter(correctedInterferogram, alpha, window_size, step_size);

% 显示Goldstein滤波后干涉图的相位
figure;
imagesc(angle(filteredInterferogram_Goldstein));
colormap('jet');
colorbar;
title('三. Goldstein 滤波后的干涉图');

%% 4. 均值滤波处理
% --------------------------
% 说明：
%   采用均值滤波对校正后的干涉图进行处理，通过计算滑动窗口内
%   实部和虚部的均值来抑制噪声。
%
% 设置均值滤波的窗口大小（可与Goldstein滤波的窗口大小相同或不同）
mean_window_size = 5;

% 调用均值滤波函数，返回滤波后的复数干涉图
filteredInterferogram_Mean = mean_filter(correctedInterferogram, mean_window_size);

% 显示均值滤波后干涉图的相位
figure;
imagesc(angle(filteredInterferogram_Mean));
colormap('jet');
colorbar;
title('四. 均值滤波后的干涉图');

%% 5. 中值滤波处理
% --------------------------
% 说明：
%   使用中值滤波对校正后的干涉图进行处理，
%   分别对实部和虚部采用中值滤波以有效去除椒盐噪声。
%
% 调用中值滤波函数，返回滤波后的复数干涉图
filteredInterferogram_Median = median_filter(correctedInterferogram, mean_window_size);

% 显示中值滤波后干涉图的相位
figure;
imagesc(angle(filteredInterferogram_Median));
colormap('jet');
colorbar;
title('五. 中值滤波后的干涉图');

%% 6. 自适应中值滤波处理
% --------------------------
% 说明：
%   自适应中值滤波根据局部统计特性自动调整窗口尺寸，
%   分别对干涉图的实部和虚部进行处理，
%   以在保留细节的同时有效抑制噪声。
%
% 设置自适应中值滤波的参数：
wmin = 5;   % 最小窗口尺寸
wmax = 7;   % 最大窗口尺寸

% 调用自适应中值滤波函数，返回滤波后的复数干涉图
filteredInterferogram_AdaptiveMedian = adaptive_median_filter(correctedInterferogram, wmin, wmax);

% 显示自适应中值滤波后干涉图的相位
figure;
imagesc(angle(filteredInterferogram_AdaptiveMedian));
colormap('jet');
colorbar;
title('六. 自适应中值滤波后的干涉图');

%% 7. 定量评价 - 计算MSE和STD指标
% --------------------------
% 说明：
%   为比较各滤波方法的效果，本节采用以下评价指标：
%   1) 均方误差（MSE）——衡量滤波前后相位的整体误差；
%   2) 全局相位残差标准差（STD）——反映全图相位一致性；
%   3) 局部相位标准差——通过滑动窗口计算局部区域相位波动，
%      进而得到改善率。
%
% 评估窗口大小设置：
eval_window_size = 5;

% 使用“去平地后但尚未滤波”的干涉图作为参考图，
% 分别对四种滤波结果进行评估。

% 1. Goldstein滤波的评估
[mse_goldstein, std_global_goldstein, std_local_before_goldstein, std_local_after_goldstein] = ...
    evaluate_filter_performance(correctedInterferogram, filteredInterferogram_Goldstein, eval_window_size);
fprintf('Goldstein 滤波评估结果:\n');
fprintf('  - MSE: %.4f\n', mse_goldstein);
fprintf('  - 全局相位残差 STD: %.4f rad\n', std_global_goldstein);
fprintf('  - 局部相位 STD: 前 %.4f rad → 后 %.4f rad\n', std_local_before_goldstein, std_local_after_goldstein);
fprintf('  - STD 改善率: %.2f%%\n', (1 - std_local_after_goldstein/std_local_before_goldstein) * 100);

% 2. 均值滤波的评估
[mse_mean, std_global_mean, std_local_before_mean, std_local_after_mean] = ...
    evaluate_filter_performance(correctedInterferogram, filteredInterferogram_Mean, eval_window_size);
fprintf('\n均值滤波评估结果:\n');
fprintf('  - MSE: %.4f\n', mse_mean);
fprintf('  - 全局相位残差 STD: %.4f rad\n', std_global_mean);
fprintf('  - 局部相位 STD: 前 %.4f rad → 后 %.4f rad\n', std_local_before_mean, std_local_after_mean);
fprintf('  - STD 改善率: %.2f%%\n', (1 - std_local_after_mean/std_local_before_mean) * 100);

% 3. 中值滤波的评估
[mse_median, std_global_median, std_local_before_median, std_local_after_median] = ...
    evaluate_filter_performance(correctedInterferogram, filteredInterferogram_Median, eval_window_size);
fprintf('\n中值滤波评估结果:\n');
fprintf('  - MSE: %.4f\n', mse_median);
fprintf('  - 全局相位残差 STD: %.4f rad\n', std_global_median);
fprintf('  - 局部相位 STD: 前 %.4f rad → 后 %.4f rad\n', std_local_before_median, std_local_after_median);
fprintf('  - STD 改善率: %.2f%%\n', (1 - std_local_after_median/std_local_before_median) * 100);

% 4. 自适应中值滤波的评估
[mse_adaptive_median, std_global_adaptive, std_local_before_adaptive, std_local_after_adaptive] = ...
    evaluate_filter_performance(correctedInterferogram, filteredInterferogram_AdaptiveMedian, eval_window_size);
fprintf('\n自适应中值滤波评估结果:\n');
fprintf('  - MSE: %.4f\n', mse_adaptive_median);
fprintf('  - 全局相位残差 STD: %.4f rad\n', std_global_adaptive);
fprintf('  - 局部相位 STD: 前 %.4f rad → 后 %.4f rad\n', std_local_before_adaptive, std_local_after_adaptive);
fprintf('  - STD 改善率: %.2f%%\n', (1 - std_local_after_adaptive/std_local_before_adaptive) * 100);

% 5. 绘制各滤波方法评估结果对比图
figure;

% 定义各滤波方法名称与对应评价指标数据
filter_names = {'Goldstein', '均值', '中值', '自适应中值'};
mse_values          = [mse_goldstein, mse_mean, mse_median, mse_adaptive_median];
std_global_values   = [std_global_goldstein, std_global_mean, std_global_median, std_global_adaptive];
std_improvement     = [
    (1 - std_local_after_goldstein/std_local_before_goldstein) * 100, ...
    (1 - std_local_after_mean    /std_local_before_mean)    * 100, ...
    (1 - std_local_after_median  /std_local_before_median)  * 100, ...
    (1 - std_local_after_adaptive/std_local_before_adaptive) * 100
];

% 子图1：MSE对比（越小越好）
subplot(1, 3, 1);
bar(mse_values);
set(gca, 'XTickLabel', filter_names);
title('MSE对比 (越小越好)');
ylabel('均方误差');
grid on;

% 子图2：全局相位残差STD对比（越小越好）
subplot(1, 3, 2);
bar(std_global_values);
set(gca, 'XTickLabel', filter_names);
title('全局相位残差STD对比 (越小越好)');
ylabel('标准差 (rad)');
grid on;

% 子图3：局部STD改善率对比（越大越好）
subplot(1, 3, 3);
bar(std_improvement);
set(gca, 'XTickLabel', filter_names);
title('局部STD改善率对比 (越大越好)');
ylabel('改善率 (%)');
grid on;

% 调整整体图形尺寸与布局
set(gcf, 'Position', [100, 100, 1200, 400]);
sgtitle('各滤波方法性能评估对比');

%% 均值滤波函数
% 说明：
%   对输入的复数干涉图分别对其实部与虚部采用滑动窗口均值滤波，
%   并对边界点进行扩充处理，确保滤波结果的连续性与稳定性。
%
% 输入参数：
%   input_data  - 待滤波的干涉图（复数矩阵）
%   window_size - 滑动窗口尺寸（正整数）
%
% 输出：
%   out_data    - 均值滤波后的干涉图（复数矩阵）
function out_data = mean_filter(input_data, window_size)
    [rows, cols] = size(input_data);
    mid = floor(window_size / 2);
    out_data = input_data;  % 初始化输出矩阵

    % 创建数据副本，将NaN值置为0，防止影响均值计算
    data_copy = input_data;
    data_copy(isnan(data_copy)) = 0;
    
    % 对图像进行边界扩充，采用复制边界法
    expanded_data = padarray(data_copy, [mid mid], 'replicate', 'both');
    
    % 分离扩充数据的实部与虚部
    real_part = real(expanded_data);
    imag_part = imag(expanded_data);

    % 对图像每个像素点采用滑动窗口计算均值
    for i = 1:rows
        for j = 1:cols
            real_window = real_part(i:i+2*mid, j:j+2*mid);
            imag_window = imag_part(i:i+2*mid, j:j+2*mid);
            
            real_mean = mean(real_window(:));
            imag_mean = mean(imag_window(:));
            
            out_data(i, j) = complex(real_mean, imag_mean);
        end
    end
    
    % 处理特殊情况：若原始相位为0或NaN，则对应输出也置为0或NaN
    idx = angle(input_data) == 0;
    out_data(idx) = 0;
    idx = isnan(angle(input_data));
    out_data(idx) = nan;
end

%% 中值滤波函数
% 说明：
%   分别对输入复数干涉图的实部与虚部进行滑动窗口中值滤波，
%   采用边界扩充法处理图像边缘，降低噪声影响。
%
% 输入参数：
%   input_data  - 待滤波的干涉图（复数矩阵）
%   window_size - 滑动窗口尺寸（正整数）
%
% 输出：
%   out_data    - 中值滤波后的干涉图（复数矩阵）
function out_data = median_filter(input_data, window_size)
    [rows, cols] = size(input_data);
    mid = floor(window_size / 2);
    out_data = input_data;  % 初始化输出矩阵

    % 创建数据副本，将NaN值置为0
    data_copy = input_data;
    data_copy(isnan(data_copy)) = 0;
    
    % 对图像进行边界扩充（复制边界）
    expanded_data = padarray(data_copy, [mid mid], 'replicate', 'both');
    
    % 分离实部与虚部
    real_part = real(expanded_data);
    imag_part = imag(expanded_data);

    % 对每个像素使用滑动窗口分别计算实部与虚部的中值
    for i = 1:rows
        for j = 1:cols
            real_window = real_part(i:i+2*mid, j:j+2*mid);
            imag_window = imag_part(i:i+2*mid, j:j+2*mid);
            
            real_median = median(real_window(:));
            imag_median = median(imag_window(:));
            
            out_data(i, j) = complex(real_median, imag_median);
        end
    end
    
    % 同样处理特殊情况：相位为0或NaN的位置
    idx = angle(input_data) == 0;
    out_data(idx) = 0;
    idx = isnan(angle(input_data));
    out_data(idx) = nan;
end

%% 自适应中值滤波函数
% 说明：
%   自适应中值滤波针对每个像素根据局部窗口内的统计特性动态调整窗口尺寸，
%   分别对实部和虚部进行处理，以更好地保留图像细节并抑制噪声。
%
% 输入参数：
%   input_data - 待滤波的干涉图（复数矩阵）
%   wmin       - 最小窗口尺寸（奇数）
%   wmax       - 最大窗口尺寸（奇数）
%
% 输出：
%   out_data   - 自适应中值滤波后的干涉图（复数矩阵）
function out_data = adaptive_median_filter(input_data, wmin, wmax)
    [rows, cols] = size(input_data);
    out_data = input_data;  % 初始化输出矩阵

    % 处理NaN值
    data_copy = input_data;
    data_copy(isnan(data_copy)) = 0;
    
    % 分离实部与虚部
    real_part = real(data_copy);
    imag_part = imag(data_copy);
    
    % 为保证窗口完全覆盖，计算扩充尺寸
    pad_size = (wmax-1)/2;
    expand_real = padarray(real_part, [pad_size pad_size], 0, 'both');
    expand_imag = padarray(imag_part, [pad_size pad_size], 0, 'both');
    
    % 对图像每个像素分别自适应处理中、实部与虚部
    for i = 1:rows
        for j = 1:cols
            real_pixel_replaced = false;
            imag_pixel_replaced = false;
            real_result = real_part(i, j);
            imag_result = imag_part(i, j);
            
            % 对实部进行自适应中值滤波
            for n = wmin:2:wmax
                half_window = (n-1)/2;
                S_real = expand_real(i : i+2*half_window, j : j+2*half_window);
                real_max = max(S_real(:));
                real_min = min(S_real(:));
                real_med = median(S_real(:));
                
                % 若中值位于极值之间，则判断当前像素是否为异常值
                if real_med > real_min && real_med < real_max
                    if real_part(i,j) <= real_min || real_part(i,j) >= real_max
                        real_result = real_med;
                        real_pixel_replaced = true;
                    end
                    break;
                end
            end
            % 若未替换，则采用最大窗口的中值
            if ~real_pixel_replaced
                half_window = (wmax-1)/2;
                S_real = expand_real(i : i+2*half_window, j : j+2*half_window);
                real_result = median(S_real(:));
            end
            
            % 对虚部进行自适应中值滤波（方法同上）
            for n = wmin:2:wmax
                half_window = (n-1)/2;
                S_imag = expand_imag(i : i+2*half_window, j : j+2*half_window);
                imag_max = max(S_imag(:));
                imag_min = min(S_imag(:));
                imag_med = median(S_imag(:));
                
                if imag_med > imag_min && imag_med < imag_max
                    if imag_part(i,j) <= imag_min || imag_part(i,j) >= imag_max
                        imag_result = imag_med;
                        imag_pixel_replaced = true;
                    end
                    break;
                end
            end
            if ~imag_pixel_replaced
                half_window = (wmax-1)/2;
                S_imag = expand_imag(i : i+2*half_window, j : j+2*half_window);
                imag_result = median(S_imag(:));
            end
            
            % 合成处理后的复数结果
            out_data(i,j) = complex(real_result, imag_result);
        end
    end
    
    % 特殊情况处理：相位为0或NaN的位置
    idx = angle(input_data) == 0;
    out_data(idx) = 0;
    idx = isnan(angle(input_data));
    out_data(idx) = nan;
end

%% Goldstein滤波函数
% 说明：
%   Goldstein滤波是一种基于频域加权的方法，本函数展示一种简化的
%   反距离加权Goldstein滤波思路。实际应用中可能需更复杂处理。
%
% 输入参数：
%   input_data  - 待滤波的干涉图复数矩阵
%   alpha       - 滤波参数，用于控制加权因子的幂次
%   window_size - 滑动窗口尺寸
%   step_size   - 滑动窗口步长
%
% 输出：
%   out_data    - Goldstein滤波后的干涉图复数矩阵
function out_data = goldstein_filter(input_data, alpha, window_size, step_size)
    % 构造简易卷积核，并通过傅里叶变换获取其频域形式
    K = ones(3, 3);
    K = K / sum(K(:));
    K = fftshift(fft2(K));

    [rows, cols] = size(input_data);
    out_data = zeros(rows, cols);  % 初始化输出矩阵

    % 若window_size为奇数，则将其调整为偶数（便于窗口分割）
    if mod(window_size,2) ~= 0
        window_size = window_size - 1;  
    end
    
    % 构造简化版的权重矩阵
    x = 1:window_size/2;
    [X, Y] = meshgrid(x, x);
    X = X + Y;
    weight = [X, fliplr(X)];
    weight = [weight; flipud(weight)];

    % 滑动窗口遍历整个图像
    for ii = 1:step_size:rows
        for jj = 1:step_size:cols
            mm = ii + window_size - 1;
            if mm > rows, mm = rows; end
            nn = jj + window_size - 1;
            if nn > cols, nn = cols; end
            
            % 取当前窗口数据
            window = input_data(ii:mm, jj:nn);
            ww = weight(1:(mm-ii+1), 1:(nn-jj+1));

            % 计算窗口数据的傅里叶变换，并移到频域中心
            H = fft2(window);
            H = fftshift(H);

            % 对幅值进行局部卷积平滑，再进行alpha加权
            S = conv2(abs(H), K, 'same');
            S = S / max(S(:) + eps);
            S = S .^ alpha;

            % 加权后反变换回时域
            H = H .* S;
            H = ifftshift(H);
            window_filt = ifft2(H);

            % 累计滤波结果（考虑权重叠加）
            out_data(ii:mm, jj:nn) = out_data(ii:mm, jj:nn) + window_filt .* ww;
        end
    end
    
    % 特殊情况处理：相位为0或NaN的像素
    idx = angle(input_data) == 0;
    out_data(idx) = 0;
    idx = isnan(angle(input_data));
    out_data(idx) = nan;
end

%% 评估函数：计算MSE和STD指标
% 说明：
%   该函数对比原始（去平地但未滤波）干涉图与滤波后干涉图，
%   分别计算：
%     1) 均方误差（MSE），衡量整体相位误差；
%     2) 全局相位残差标准差（global_std）；
%     3) 局部相位标准差（local_std），通过滑动窗口求局部统计平均，
%        用于量化噪声抑制效果。
%
% 输入参数：
%   original_interferogram - 去平地后但未滤波的干涉图（复数矩阵）
%   filtered_interferogram - 滤波后的干涉图（复数矩阵）
%   window_size            - 计算局部标准差的窗口尺寸
%
% 输出：
%   mse             - 均方误差（MSE）
%   global_std      - 全局相位残差的标准差
%   local_std_before- 滤波前局部相位标准差的均值
%   local_std_after - 滤波后局部相位标准差的均值
function [mse, global_std, local_std_before, local_std_after] = evaluate_filter_performance(original_interferogram, filtered_interferogram, window_size)
    % 提取原始与滤波后干涉图的相位信息
    original_phase = angle(original_interferogram);  
    filtered_phase = angle(filtered_interferogram);  

    % 1. 计算均方误差（MSE）
    mse = mean((original_phase(:) - filtered_phase(:)).^2, 'omitnan');

    % 2. 计算全局相位残差标准差（对相位差进行包裹到 [-pi, pi] 范围）
    phase_diff = wrapToPi(filtered_phase - original_phase);
    global_std = std(phase_diff(:), 'omitnan');

    % 3. 计算局部相位标准差（采用stdfilt函数在滑动窗口内计算标准差后取均值）
    local_std_before = stdfilt(original_phase, true(window_size));
    local_std_after  = stdfilt(filtered_phase, true(window_size));

    local_std_before = mean(local_std_before(:), 'omitnan');
    local_std_after  = mean(local_std_after(:), 'omitnan');
end
