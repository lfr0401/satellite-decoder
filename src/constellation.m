function softBitsAll = constellation(recursive, symbols, Params)
% QPSK 星座解映射 + 编码域定帧
% 输入 symbols：demod.m 输出的复数 QPSK 符号
%      recursive：是否递归调用（处理首尾星座映射异常的帧）
% 输出 softBitsAll：LLR 软比特流（每符号 2 bit），供 decode.m 做 Viterbi 软译码
% 核心：盲试 4 种星座映射，用同步字互相关选对的那种，同时确定帧头位置

    % 嵌套函数：拼接两段软比特后重新找峰、砍掉无帧头的尾部
    function softBitsConfirmed = findPeaksConstellationError(softBitsCombined, syncAsmBits)
        [corrError, lags] = xcorr(softBitsCombined, syncAsmBits);   % 软比特与同步字互相关
        corrNormalized = abs(corrError)./max(abs(corrError));       % 归一化到 0~1
        signalLevelError = quantile(abs(corrNormalized), 0.999);    % 99.9% 分位数当噪声底
        % 找峰：间距≥16384-1（一帧=16384 软比特），阈值 0.1，峰高>噪声底+0.2
        [pks, locs] = findpeaks(corrNormalized, 'MinPeakDistance',16384-1, 'Threshold', 0.1, 'MinPeakHeight', signalLevelError+0.2); % , 'MinPeakHeight', 40+signalLevel
        if isempty(locs)
            softBitsConfirmed = softBitsCombined;                   % 没找到峰，原样返回
        else
            softBitsConfirmed = softBitsCombined(1:max(0,lags(locs(end))));  % 只留到最后一个帧头
        end
    end

    % 编码域同步字（在软比特上做互相关，用于定帧）
    syncAsm = 'FCA2B63DB00D9794';                      % 8 字节=64 bit，编码域同步字
    syncAsmBytes = sscanf(syncAsm, '%2x').';           % 十六进制字符串→8 个字节
    syncAsmBits = de2bi(syncAsmBytes, 8, 'left-msb');  % 每字节展开成 8 bit，高位在前
    syncAsmBits = reshape(syncAsmBits.', 1, []);       % 拉成 1×64 行向量
    syncAsmBits = 2*double(syncAsmBits)-1;             % 0→-1、1→+1，与 LLR 同域才能互相关
    
    hPlot = [];
    hPeaks = [];
    if Params.plotting && ~recursive    % 只在非递归时建图，避免递归重复开图
        figure; clf;
        hPlot = plot(nan, nan); hold on;      % 相关曲线句柄（nan 占位）
        hPeaks = plot(nan, nan, 'rx');        % 红叉标记检出的帧头
        hold off;
        xlabel('Samples');
        ylabel('Cross-correlation Amplitude');
        grid on;
    end
    
    i = 1;                              % 符号块起始指针
    constellationErrorLast = 0;         % 末尾几帧星座映射异常标志
    constellationErrorFirst = 0;        % 开头几帧星座映射异常标志
    lastFrame = 0;                      % 末帧标志
    softBitsAll = [];                   % 总输出累加器
    while i <= numel(symbols)           % 逐块处理到符号末尾
        softBitsConfirmed = [];         % 本块确认的软比特，空表示尚未确认
        if i+Params.blockSize-1 < numel(symbols)
            symbolsBlock = symbols(i:i+Params.blockSize-1);   % 满块
        else
            lastFrame = 1;              % 不足一块，标记为末帧
            symbolsBlock = symbols(i:end);                     % 取到结尾
        end
        % 遍历 4 种星座旋转/对称方式（解决 90°/180° 相位模糊与 I/Q 镜像），自动选出正确的一种
        for j=1:length(Params.constellations)
            % QPSK 软判决解调，输出对数似然比（LLR）
            % symbolsBlock：复数符号；Params.M=4：QPSK；pi/4：星座点在 ±1±j
            % constellations{j}：当前尝试的映射（默认/旋转90°/虚部反/IQ交换）
            % OutputType="llr"：输出 LLR 而非硬 0/1，符号=判决、幅值=置信度
            softBits = pskdemod( ...
                        symbolsBlock, ...
                        Params.M, ...
                        pi/4, ...
                        Params.constellations{j}, ...
                        OutputType="llr" ...
                        );
            % 每符号出 2 个 LLR，故 softBits 长度 = 符号数 × 2
    
            % 跳过含非法值或全零的静音数据块
            if any(~isfinite(softBits)) || ~any(softBits)
                continue
            end
            % 软比特与同步字互相关，找帧头位置
            [corr, lags] = xcorr(softBits, syncAsmBits);
            corrNormalized = abs(corr)./max(abs(corr));         % 归一化到 0~1
            signalLevel = quantile(abs(corrNormalized), 0.999); % 99.9% 分位数当噪声底
    
            % 噪声底过高（峰值与噪声差距太小），该星座方式不成立
            if signalLevel > 0.6
                continue
            end
    
            % 找相关峰并计算峰间距（帧长应为 16384 软比特）
            [pks, locs] = findpeaks(corrNormalized, 'MinPeakDistance',16384-1, 'Threshold', 0.1, 'MinPeakHeight', signalLevel+0.2); % , 'MinPeakHeight', 40+signalLevel
            locs = locs(lags(locs)>=0);    % 只保留非负 lag 的峰（去掉负半轴镜像峰）
            if isempty(locs)
                continue
            end
            widths = diff(locs)/16;        % 峰间距换算：16384 软比特/16 = 1024，用于判断帧长
    
            % 相关峰为负：说明极性反转（相位差 180°），软比特整体取反
            if corr(locs) < 0
                softBits = -softBits;
            end
    
            if all(widths == 1024)         % 帧长规整时才检查首尾异常帧
                firstPeak = locs(1);       % 第一个帧头位置
                expectedNumPeaks = ceil(numel(corrNormalized(firstPeak:end))/16384);  % 按长度算应有峰数
    
                % 末尾几帧星座映射与其他帧不同，需递归单独处理
                if expectedNumPeaks ~= numel(locs)
                    part1 = softBits(1:min(numel(softBits),lags(locs(end))+16384));   % 正常段软比特
                    part2 = symbolsBlock(lags(locs(end))/2+1+8192+1:end);            % 末尾异常段符号
                    softBits = constellation(1, part2, Params);   % 递归解映射末尾段
                    constellationErrorLast = 1; 
    
                % 开头几帧星座映射与其他帧不同，需递归单独处理
                elseif lags(locs(1)) > 16384+1
                    part2 = softBits(lags(locs(1))-1:end);        % 正常段软比特
                    part1 = symbolsBlock(1:lags(locs(1))/2-1);    % 开头异常段符号
                    softBits = constellation(1, part1, Params);   % 递归解映射开头段
                    constellationErrorFirst = 1; 
                end
            end
            
            % 末尾几帧星座映射不同：拼接两段后重新确认帧头
            if constellationErrorLast
                constellationErrorLast = 0;                        % 标志复位
                softBitsCombined = [part1; softBits];              % 正常段 + 递归解出的末尾段
                softBitsConfirmed = findPeaksConstellationError(softBitsCombined, syncAsmBits);
                break

            % 开头几帧星座映射不同：拼接两段后重新确认帧头
            elseif constellationErrorFirst
                constellationErrorFirst = 0;                       % 标志复位
                softBitsCombined = [softBits; part2];              % 递归解出的开头段 + 正常段
                softBitsConfirmed = findPeaksConstellationError(softBitsCombined, syncAsmBits);
                break
            
            % 递归分支：已确认星座方式，直接采纳
            elseif all(~mod(widths, 1024)) && recursive
                softBitsConfirmed = softBits;
                break
            
            % 最后一帧：保留本块全部软比特
            elseif all(~mod(widths, 1024)) && lastFrame
                softBitsConfirmed = softBits;
                break
                
            % 正常情况：只保留到最后一个帧头为止，避免跨块残缺数据进入译码
            elseif all(~mod(widths, 1024)) && ~lastFrame && ~recursive
                softBitsConfirmed = softBits(1:lags(locs(end)));
                break
            end
        end
        if isempty(softBitsConfirmed)
            % 未确认星座时补零，保持后续帧的时间轴对齐
            softBitsConfirmed = zeros(2*numel(symbolsBlock),1,'like',real(symbolsBlock));
        end
        softBitsAll = [softBitsAll; softBitsConfirmed];                 % 累加本块结果
        i = i + numel(softBitsConfirmed)/2;                             % 软比特数/2 = 符号推进量

        % 绘制软比特与同步字的互相关曲线，标注检出的帧头
        if Params.plotting && ~recursive && ~isempty(hPlot) && exist('corrNormalized','var')
            set(hPlot, 'XData', lags, 'YData', corrNormalized);
            set(hPeaks, 'XData', lags(locs), 'YData', corrNormalized(locs));
            xlim([0 max(lags)]);
            title(sprintf('Cross-correlation of FCA2B63DB00D9794 and encoded Softbits using %s', mat2str(Params.constellations{j})));  % 标题显示胜出的映射
            drawnow;
        end
    end
    % 导出 corr_encoded.pdf（矢量图），仅非递归且开启 export 时
    if Params.plotting && Params.export && ~recursive
        exportgraphics(gcf, "data/plots/corr_encoded.pdf")
    end
end
