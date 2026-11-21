function cvcdus = decode(softBits, Viterbi, Descrambler, Params)
    if isempty(softBits)
        cvcdus = zeros(0,1020,'uint8');
        return
    end

    % 维特比译码器：(171,133) 卷积码，1/2 码率，K=7
    trellis = poly2trellis(Viterbi.constLen, Viterbi.codeGenPoly);
    vDec = comm.ViterbiDecoder( ...
            'TrellisStructure', trellis, ...
            'InputFormat', 'Unquantized',...
            'TracebackDepth', Viterbi.tblen...
            );
    % 非量化输入约定：比特 0 对应正 LLR
    decodedBits = vDec(softBits);
    decodedBits = 2*double(decodedBits)-1;
    
    % 解码域同步字（ASM），在译码后的比特流上定帧
    syncAsm = '1ACFFC1D';
    syncAsmBytes = sscanf(syncAsm, '%2x').';
    syncAsmBits = de2bi(syncAsmBytes, 8, 'left-msb');
    syncAsmBits = reshape(syncAsmBits.', [], 1);
    syncAsmBits = 2*double(syncAsmBits)-1;
    
    [corr, lags] = xcorr(decodedBits, syncAsmBits);
    [pks, locs] = findpeaks(abs(corr), 'MinPeakHeight', 29);
    
    decodedBits = double(decodedBits+1)/2;
    
    % 逐帧剥离同步字后解扰，同时补齐缺失帧的占位数据
    payloads = {};
    isPlaceholder = logical([]);
    for i = 1:numel(locs)
        startIdx = lags(locs(i))+length(syncAsmBits)+1;
        expectedStopIdx = startIdx+8160-1;

        if expectedStopIdx > numel(decodedBits)
            break
        end
        if i == numel(locs)
            nextLag = Inf;
        else
            nextLag = lags(locs(i+1));
        end

        % 完整帧：按固定长度 8160 比特截取
        if nextLag >= expectedStopIdx
            stopIdx = expectedStopIdx;

        % 短帧：不足部分补零
        elseif nextLag < expectedStopIdx
            stopIdx = nextLag;
        end

        payload = zeros(8160, 1);
        payload(1:stopIdx-startIdx+1) = decodedBits(startIdx:stopIdx);
        payloads{end+1} = payload;
        isPlaceholder(end+1) = false;

        % 检测不到的帧：插入全零占位行，避免时间轴错位
        nMissing = 0;
        if isfinite(nextLag)
            nMissing = round((nextLag - lags(locs(i))) / 8192) - 1;
        end
        for m = 1:nMissing
            payloads{end+1} = zeros(8160, 1);
            isPlaceholder(end+1) = true;
        end
    end

    % 解扰：与 255 长伪随机序列按位异或（占位帧跳过）
    numFrames = numel(payloads);
    frameLenBytes = 1024;
    cvcdus = zeros(numFrames, frameLenBytes-4, 'uint8');

    for i = 1:numFrames
        if isPlaceholder(i)
            continue  % 占位帧保持全零，不与伪随机序列异或
        end
        payload = payloads{i};
        payload = reshape(payload, 8, []).';
        payload = uint8(bi2de(payload, 'left-msb'));
        idx = mod(0:numel(payload)-1, 255) + 1;
        payloadDescrambled = bitxor(payload, reshape(Descrambler.pn(idx),[],1));
        cvcdus(i,:) = payloadDescrambled;
    end
    
    % 绘制 ASM 互相关与检出的帧头位置
    if Params.plotting
        figure;
        plot(lags, abs(corr)); hold on;
        plot(lags(locs), pks, 'rx'); hold off;
        xlim([0 length(corr)/2]);
        xlabel('Samples');
        ylabel('Cross-correlation');
        title('Cross-correlation of 1ACFFC1D and decoded Softbits');
        grid on;
        if Params.export
            exportgraphics(gcf, "data/plots/corr_decoded.pdf")
        end
    end
end