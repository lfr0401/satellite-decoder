function symbols = demod(Data, Params, Rcc)
% 物理层核心：复数 IQ 采样 → QPSK 符号，供下一级 constellation.m 软判决
% 输入 Data.raw(:,1)=I, Data.raw(:,2)=Q; Data.fs=原始采样率
%      Params=重采样/分块/绘图参数; Rcc=RRC 滤波器参数
% 输出 symbols：复数 QPSK 符号列向量
% 处理链：限幅归一化 → 4 倍重采样 → RRC 匹配滤波 → 粗频偏 → 定时 → 精相位
% 说明：下面 4 个对象均来自 Communications Toolbox 的 comm 包

    % ① 接收端根升余弦(RRC)匹配滤波器
    % 发射端 RRC × 接收端 RRC = 完整升余弦 = 奈奎斯特滤波器 → 符号时刻零 ISI
    % 同时是白噪声匹配滤波，输出 SNR 最大，是软判决性能的前提
    % RolloffFactor：滚降系数(=0.35)，必须与发射端一致，否则 ISI 不为零
    % FilterSpanInSymbols：冲激响应跨度，单位"符号"，配置值 10；越长越理想但延迟越大
    % InputSamplesPerSymbol：输入 4 sps，必须和 resample 后的倍率一致
    % DecimationFactor=1：不抽取，输出仍是 4 sps，抽取交给 symSync
    % 注意：rrcRx / cfo / symSync 都是固定尺寸对象，输入长度变化前需 release
    rrcRx = comm.RaisedCosineReceiveFilter( ...
                'RolloffFactor', Rcc.rollOff, ...
                'FilterSpanInSymbols', Rcc.spanSym, ...
                'InputSamplesPerSymbol', Params.targetSps, ...
                'DecimationFactor', 1 ...
                );

    % ② 粗频偏补偿：QPSK 四次方环
    % 信号四次方后四相位差消去，谱峰位于 4fc，峰位置/4 即频偏
    % 补偿范围约 ±fs/4，先拉掉卫星多普勒+本振偏差，再交 carrierSync 精调
    cfo = comm.CoarseFrequencyCompensator(Modulation='qpsk', SampleRate=Data.fs);

    % ③ 符号定时同步：判决反馈过零点 TED + 插值器
    % 比较过零点实际位置与理想符号中心，反馈环调整采样相位
    % 输入 4 sps → 输出 1 sps，抽取在这里完成；低 SNR 下抖动偏大
    % 三级同步顺序是硬约束：频偏 → 定时 → 相位，交换任一级都锁不住
    symSync = comm.SymbolSynchronizer( ...
                'TimingErrorDetector', 'Zero-Crossing (decision-directed)', ...
                'SamplesPerSymbol', Params.targetSps ...
                );

    % ④ 载波同步：QPSK 锁相环，精相位恢复，把星座旋到标准象限
    % 四次方鉴相；残留 90° 相位模糊本模块解不掉，留给 constellation.m 试 4 种映射
    % SamplesPerSymbol=1：输入已是 symSync 输出的 1 sps 符号
    carrierSync = comm.CarrierSynchronizer( ...
                    'Modulation','QPSK', ...    
                    'SamplesPerSymbol',1 ...
                    );
    
    % 读取 IQ 数据：wav 第 1 列=I，第 2 列=Q；单精度省内存
    I = single(Data.raw(:,1));
    Q = single(Data.raw(:,2));

    % clip 限幅：掐掉超过门限的脉冲干扰尖峰，避免归一化时把正常信号压小
    % 门限取 Params.minClip/maxClip（约 ±0.01），太小会削平真信号→星座畸变
    I = clip(I, Params.minClip, Params.maxClip);
    Q = clip(Q, Params.minClip, Params.maxClip);

    % 按各自峰值归一化到 ±1；eps('single') 是防全零数据除零的极小值
    % 注意 I/Q 独立归一化：若两路峰值差很多，说明硬件 I/Q 幅度不平衡
    I = I / max(max(abs(I)), eps('single'));
    Q = Q / max(max(abs(Q)), eps('single'));

    % 合成复基带信号，后续全部针对这个复序列处理
    x = (I + 1j*Q);
    
    % FIR 重采样（多相滤波）：fs → 288 kHz = 72 ksym/s × 4
    % p/q 由 rat(288e3/fs) 得到；4 sps 是定时同步器插值判决的前提
    xResampled = resample(x, Params.p, Params.q);
    
    % RRC 匹配滤波：与发射端根升余弦成对，最大化 SNR、消除 ISI
    % 对整段数据一次滤波，不分块，避免块边界的滤波器暂态污染
    yFilteredAll = rrcRx(xResampled);
    
    i = 1;
    symbols = [];

    hPlot = [];
    if Params.plotting
        figure; clf;
        hPlot = plot(nan, nan, '.', 'MarkerSize', 8);
        axis equal;
        axis([-2.5 2.5 -2.5 2.5]);
        grid on;
        xlabel('I');
        ylabel('Q');
        title('QPSK Constellation');
    end
    

    % 分块循环：每块 blockSize=98304 个采样点，逐块跑三级同步
    while i <= numel(yFilteredAll)
        if i+Params.blockSize-1 <= numel(yFilteredAll)
            yFiltered = yFilteredAll(i:i+Params.blockSize-1);
        else
            % 末尾数据块不足一整块时，取到文件结尾
            yFiltered = yFilteredAll(i:end);
        end
    
        % 末尾残块长度不固定，先释放并重置固定尺寸对象
        % 否则 cfo/symSync 会因输入长度与之前不同而报错
        if numel(yFiltered) < Params.blockSize
            if numel(yFiltered) < Params.targetSps
                break   % 不足一个符号，直接结束
            end
            release(cfo);
            release(symSync);
            % carrierSync 吃 1 sps 符号，非固定尺寸，无需 release
        end

        % 第一级：粗载波频偏补偿（QPSK 四次方环），输出仍是 4 sps
        yCfo = cfo(yFiltered);
    
        % 第二级：符号定时恢复（判决反馈过零点检测），4 sps → 1 sps
        ySync = symSync(yCfo);
    
        % 第三级：精载波相位恢复（QPSK 锁相环），输出落在标准星座点
        symBlock = carrierSync(ySync);
    
        % 累加本块符号到总输出
        symbols = [symbols; symBlock];

        i = i + Params.blockSize;
    
        % 实时刷新星座图：drawnow limitrate 限帧率，避免每块重绘拖慢速度
        if Params.plotting
            set(hPlot, 'XData', real(symBlock), 'YData', imag(symBlock));
            drawnow limitrate;
        end
    end
    % 导出 demodulated.pdf（矢量图），仅在 plotting 与 export 都开启时
    if Params.plotting && Params.export
        exportgraphics(gcf, "data/plots/demodulated.pdf")
    end
end