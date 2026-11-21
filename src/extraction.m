function [mcus, qualityFactors, apids, Meta] = extraction(cvcdus, Params)
    function remaining = readPackets(bytes, multiplePackets)
        remaining = uint8([]);
        pos = 1;
        while pos <= numel(bytes)
            % 不足 6 字节的包头留到下一帧再拼
            if numel(bytes)-pos+1 < 6
                remaining = bytes(pos:end);
                return
            end
            lenBytes = bytes(pos+4:pos+5);
            lenDec = double(lenBytes(1))*256 + double(lenBytes(2));
            totalLen = 6 + lenDec + 1;
            if bitshift(bytes(pos),-5) ~= 0
                return  % 首指针残留，转到下一帧重新开始
            end
            if pos+totalLen-1 > numel(bytes)
                remaining = bytes(pos:end);
                return
            end

            tempPP = bytes(pos:pos+totalLen-1);
            if totalLen >= 21 && tempPP(1) == 8 && ismember(tempPP(2), [64 65 68]) && ...
                    tempPP(18) == 255 && tempPP(19) == 240 && ...
                    tempPP(15) <= 182 && mod(double(tempPP(15)),14) == 0 && ...
                    tempPP(20) >= 1 && tempPP(20) <= 100
                pp{end+1} = tempPP;
            end
            pos = pos + totalLen;
            if ~multiplePackets
                return
            end
        end
    end

    function counter = calcCounter(Header)
        counterPP1 = Header(3).';
        counterPP1 = int2bit(counterPP1.', 8).';
        counterPP2 = Header(4).';
        counterPP2 = int2bit(counterPP2.', 8).';
        counterBit = [counterPP1 counterPP2];
        counterBit = counterBit(:,3:end);
        counter = int16(bi2de(counterBit, 'left-msb'));
    end
    
    % 提取帧头信息：剥离尾部的 RS 校验区，拆出 MPDU 指针与载荷
    vcdus = cvcdus(:,1:end-128);
    mpdus = vcdus(:,9:end);
    mpdusPayload = mpdus(:,3:end);
    mpdusHeader = mpdus(:,1:2);
    mpduPointerDec = double(bitand(uint16(mpdusHeader(:,1))*256 + ...
                                 uint16(mpdusHeader(:,2)), uint16(2047)));
    [nRows, nCols] = size(mpdusPayload);

    pp = {};
    partialPP = uint8([]);
    lastCounter = [];
    lastId = [];

    for row = 1:nRows
        frame = cvcdus(row,:);
        if ~any(frame) || bitshift(frame(1),-6) ~= 1
            partialPP = uint8([]);
            lastCounter = [];
            continue
        end

        % 帧丢失或通道切换时，不跨帧拼接分包
        counter = double(frame(3))*65536 + double(frame(4))*256 + double(frame(5));
        vcduId = double(frame(1))*256 + double(frame(2));
        if ~isempty(lastCounter) && (mod(counter-lastCounter,2^24) ~= 1 || vcduId ~= lastId)
            partialPP = uint8([]);
        end
        lastCounter = counter;
        lastId = vcduId;

        idx = mpduPointerDec(row) + 1;
        if idx == 2048
            % 本帧无新包头指针，继续拼接上一帧的残包
            if ~isempty(partialPP)
                partialPP = readPackets([partialPP, mpdusPayload(row,:)], false);
            end
        elseif idx <= nCols
            if ~isempty(partialPP)
                readPackets([partialPP, mpdusPayload(row,1:idx-1)], false);
            end
            % 指针指向新包起点：即使旧包未收完整也强制切开
            partialPP = readPackets(mpdusPayload(row,idx:end), true);
        else
            partialPP = uint8([]);
        end
    end

    % 清理残包：按计数器与 MCU 编号筛选出完整的缩略图分包
    nPP = numel(pp);
    validApid = [64 65 68 70];
    % 预分配标记数组，缩略图校验通过才置为 true
    keepPP = false(1, nPP);
    
    expectedCounter   = 0:14:182;
    nPerThumb  = numel(expectedCounter);
    
    % 只保留非空分包
    nonEmptyIdx = find(~cellfun(@isempty, pp));   
    
    % 提取所有分包的 APID（应用进程标识）
    apidAll = cellfun(@(p) p(2), pp(nonEmptyIdx));  
    
    for a = 1:numel(validApid)
        apidVal = validApid(a);  
        % 当前 APID 对应的非空分包索引
        idxApid = nonEmptyIdx(apidAll == apidVal);        
    
        mcuCounters = cellfun(@(p) double(p(15)), pp(idxApid));
        % MCU 计数器归零处即为一张缩略图的起点
        startThumbnail = find(mcuCounters == 0);
    
        for k = 1:numel(startThumbnail)
    
            startIdx = startThumbnail(k);
            endIdx = startIdx;
        
            % 计数器严格递增才连续归属同一张图（182 → 0 表示换图） 
            while endIdx + 1 <= numel(mcuCounters) && mcuCounters(endIdx + 1) > mcuCounters(endIdx)
                endIdx = endIdx + 1;
            end
            thumbnailCounters = mcuCounters(startIdx:endIdx);
        
            % 长度与计数器序列完全匹配才保留（14 个包，计数器 0/14/…/182）
            packetCounters = cellfun(@calcCounter, pp(idxApid(startIdx:endIdx)));
            if numel(thumbnailCounters) == nPerThumb && isequal(thumbnailCounters(:).', expectedCounter) && ...
                    all(mod(diff(double(packetCounters)),16384) == 1)
                keepPP(idxApid(startIdx:endIdx)) = true;
            end
        end
    end
    
    % 保留部分扫描行（后续按时间戳与 MCU 号定位，残行也有价值）
    if isfield(Params, 'keepPartialScans') && Params.keepPartialScans
        keepPP(:) = true;
    end

    ppClean = pp(keepPP);
    nPP = numel(ppClean);
    mcus = cell(1, nPP);
    qualityFactors = zeros(1, nPP);
    apids = zeros(1, nPP);
    Meta.mcu = zeros(1, nPP);
    Meta.time = zeros(1, nPP);
    Meta.sequence = zeros(1, nPP);
    
    for i = 1:nPP
        p = ppClean{i};
        apids(i) = p(2);
        qualityFactors(i) = p(20);
        mcusDec = p(21:end);
        mcus{i} = int2bit(mcusDec.', 8).';
        Meta.mcu(i) = p(15);
        Meta.sequence(i) = double(calcCounter(p));
        Meta.time(i) = (double(p(7))*256+double(p(8)))*86400000 + ...
                      double(p(9))*2^24+double(p(10))*65536+double(p(11))*256+double(p(12));
    end
    
    Meta.apid = apids;

    if Params.plotting
        fprintf("Extracted %d MCUs.\n", nPP);
    end
end
