close all; clear; clc; addpath("src")
tic

% 读取全部配置参数
[Data, Params, Rcc, Viterbi, Descrambler, ReedSolomon, Huffman, DCT] = getconfig();

% 预处理：检查 IQ 数据质量（时域 / 散点 / 功率谱 / 瀑布图）
qualitycheck(Data, Params);

% 解调：重采样 → RRC 匹配滤波 → 频偏 / 定时 / 相位恢复
symbols = demod(Data, Params, Rcc);

% 星座映射：QPSK 软判决 + 同步字互相关定帧
softBits = constellation(0, symbols, Params);

% 信道译码：Viterbi 软译码 + CCSDS 解扰 + 帧同步
cvcdus = decode(softBits, Viterbi, Descrambler, Params); 

% RS 纠错：GF(256) 上的 RS(255,223) 译码
[cvcdus, rsErrors] = reedsolomon(cvcdus, ReedSolomon, Params);

% 提取：从 VCDU 中拆出图像数据块（MCU）及其量化因子
[mcus, qualityFactors, apids, Meta] = extraction(cvcdus, Params);

% JPEG 解码：哈夫曼 → 反量化 → 离散余弦逆变换，重建三通道图像
Images = jpegdecoding(mcus, qualityFactors, apids, Huffman, DCT, Params, Meta);

toc
