function qualitycheck(Data, Params)     % 函数定义：信号质量检查。输入 Data(IQ数据) 与 Params(全局参数)，无输出；只画四张诊断图
    if Params.plotting                   % 只有参数开关 plotting 为真才执行；不画图时本函数直接跳过，不影响解码主流程
        % 读取 IQ 数据：wav 双声道拆成 I/Q 两路，合成复数基带信号
        I = single(Data.raw(:,1));       % 取 raw 的第 1 列(左声道)作同相分量 I，single 单精度省内存、够用
        Q = single(Data.raw(:,2));       % 取 raw 的第 2 列(右声道)作正交分量 Q，同样转单精度
        x = (I + 1j*Q);                  % 合成复数基带信号 x = I + jQ，这就是后续所有处理的输入
        fs = Data.fs;                    % 采样率(Hz)，供功率谱/频谱图横轴刻度使用
        if numel(x) < 4                  % 样本点少于 4 个，做 FFT/画图都没意义
            return                       % 提前退出函数，避免后续报错
        end
        
        % 基带信号时域波形（I / Q 分量）
        figure;                          % 新建图窗
        plot(real(x)); hold on;          % 画实部 I 随时间变化，hold on 保留以便叠加 Q
        plot(imag(x)); hold off;         % 在同一张图叠加虚部 Q，hold off 结束叠加
        ylim([-0.02, 0.02]);             % 固定纵轴到 ±0.02，便于观察幅度是否被限幅/饱和
        xlabel("Sample")                 % 横轴标签：采样点序号
        ylabel("Amplitude")              % 纵轴标签：幅度
        legend("Real Part", "Imaginary Part")  % 图例：实部(I)与虚部(Q)
        grid();                          % 显示网格线
        if Params.export                 % 若开启导出
            exportgraphics(gcf, "data/plots/qc_baseband.pdf")  % 把当前图窗存为 PDF
        end

        x = x(1:min(100000,numel(x)));   % 截断到最多 10 万点：后面三张图累，取前段足够看质量且画得快
        % 基带信号散点图（IQ 平面）
        figure;                          % 新建图窗
        plot(real(x), imag(x), '.', LineStyle='none');  % IQ 平面散点：横轴 I、纵轴 Q，点号小圆点、不连线
        ylim([-0.02, 0.02]);             % 先限定纵轴范围
        xlabel("I-Part")                 % 横轴标签：I 分量
        ylabel("Q-Part")                 % 纵轴标签：Q 分量
        axis equal;                      % 横纵轴等比例，保证星座点不被拉伸变形
        axis([-0.015 0.015 -0.015 0.015]);  % 再精确锁定坐标范围 ±0.015，聚焦信号簇
        grid();                          % 显示网格线
        if Params.export                 % 若开启导出
            exportgraphics(gcf, "data/plots/qc_scatter.pdf")  % 存为 PDF
        end

        % 功率谱密度（Welch 法，中心化频率轴）
        [pxx,f] = pwelch(x,[],[],[],fs,'centered');  % Welch PSD：窗/重叠/FFT点数都用默认，fs 定频率轴，centered 把 0 频移到中间
        psdDb = 10 .* log10(pxx);        % 功率转 dB：10·log10(P)，便于看清动态范围
        figure;                          % 新建图窗
        plot(f/1000, psdDb);             % 画 PSD，频率由 Hz 转 kHz
        xlabel("Frequency [kHz]");       % 横轴标签：频率(kHz)
        ylabel("Magnitude [dB]");        % 纵轴标签：功率(dB)
        grid();                          % 显示网格线
        if Params.export                 % 若开启导出
            exportgraphics(gcf, "data/plots/qc_psd.pdf")  % 存为 PDF
        end

        % 瀑布图（短时傅里叶变换，观察频率随时间变化）
        nFFT    = min(1024,numel(x));    % 每段 FFT 点数，最多 1024，数据不足时取实际长度
        overlap = floor(nFFT/2);         % 相邻段重叠 50%(512 点)，提高时间分辨率、减少缝隙
        window  = hanning(nFFT);         % 加汉宁窗，抑制频谱泄漏
        [S, f, t] = spectrogram(x, window, overlap, nFFT, Data.fs, "centered", "psd");  % 短时傅里叶：返回 S(各时刻频谱)、f(频率)、t(时间)；centered 居中频率轴，psd 输出功率谱密度
        U   = sum(window.^2);                    % 窗的能量和，作归一化分母(相干增益校正)
        Pxx = (abs(S).^2) / (Data.fs * U);       % 由 STFT 幅度平方换算成功率谱密度，使量纲与 PSD 一致
        sdB = 10*log10(Pxx + eps);               % 转 dB；加 eps 防止 Pxx=0 时 log 出 -Inf
        figure;                          % 新建图窗
        imagesc(f/1e3, t, sdB.');        % 画瀑布图：横轴频率(kHz)、纵轴时间(s)、颜色=功率；sdB.' 转置使行列对齐坐标
        axis xy;                         % 原点设在左下角(imagesc 默认左上，纠正为时间向上增长)
        cb = colorbar;                   % 显示颜色条
        cb.Label.String = 'Power [dB]';  % 给颜色条加标签：功率(dB)
        xlabel("Frequency [kHz]");       % 横轴标签：频率(kHz)
        ylabel("Time [s]");              % 纵轴标签：时间(s)
        if Params.export                 % 若开启导出
            exportgraphics(gcf, "data/plots/qc_waterfall.pdf")  % 存为 PDF
        end
   end
end
