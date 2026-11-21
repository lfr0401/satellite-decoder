# 卫星信号解码器 —— LRPT（Meteor-M2）

本仓库包含一个 LRPT 信号的软件接收机原型，用于接收俄罗斯 Meteor-M2 系列对地观测卫星下发的信号。
该项目是柏林工业大学（Technische Universität Berlin）"Funksysteme"课程的课程设计成果，完整演示了从复基带采样到还原出对地观测图像的整条数字卫星接收链路。
实现重点是物理层与数据链路层处理步骤的透明化 MATLAB 实现，遵循 CCSDS 标准。

## 项目范围

LRPT（Low Rate Picture Transmission，低速率图像传输）是 Meteor-M2 卫星用于下传多光谱对地图像的数字链路。
本项目实现了从录制的 LRPT 信号中恢复图像数据所需的核心信号处理与解码步骤。

已实现的功能包括：
- 载波与符号同步
- QPSK 解调
- 遵循 CCSDS 标准的信道译码
- 帧同步与解析
- 类 JPEG 图像重建
- 中间结果与最终结果的可视化

---

## 信号处理链路

接收机的处理流程如下：

1. **IQ 输入**
   - 复基带采样数据（例如来自 SDR 的录音）
2. **载波恢复与滤波**
   - 频率与相位校正
   - 根升余弦（Root-Raised-Cosine）滤波
3. **解调**
   - QPSK 解调
   - 符号定时恢复
4. **帧同步**
   - 检测附加同步码（ASM, Attached Sync Marker）
5. **信道译码**
   - Viterbi 译码（卷积码）
   - Reed-Solomon 译码
   - CCSDS 去随机化（derandomization）
6. **包处理与图像处理**
   - LRPT 帧解析
   - 类 JPEG 解压缩

## 标准与参考资料

实现基于以下标准与规范：
- CCSDS 131.0-B —— 遥测同步与信道编码
- CCITT T.81 / ISO/IEC 10918-1 —— JPEG 基线 DCT
- 公开的 Meteor-M2 LRPT 协议文档

所有相关标准与理论背景在配套的项目报告中有详细讨论。

## 仓库结构

```text
.
├── src/
│   ├── constellation.m    星座图与相关同步
│   ├── decode.m           Viterbi 译码与解扰
│   ├── demod.m            QPSK 解调与同步
│   ├── extraction.m       帧与包提取
│   ├── getconfig.m        参数配置（采样率、符号率、滤波器等）
│   ├── huffman.m          Huffman 解码
│   ├── jpegdecoding.m     类 JPEG 图像解码
│   ├── qualitycheck.m     信号质量检查（功率谱、瀑布图等）
│   └── reedsolomon.m      Reed-Solomon 纠错
│
├── data/
│   ├── meteor_m2_72k.wav  示例录音（IQ 采样）
│   ├── report.pdf         项目报告
│   ├── presentation.pdf   答辩演示文稿
│   └── plots/             各处理阶段的输出图
│
├── main.m                 程序入口
└── README.md
```

## 环境要求

- MATLAB（建议 R2022b 或更新版本）
- Image Processing Toolbox
- Communications Toolbox

可选：
- SDR 硬件（如 RTL-SDR），用于自行录制 IQ 数据

> 注：仓库中不包含示例录音 `data/meteor_m2_72k.wav`，需自行准备 IQ 数据并放到 `data/` 目录下。

## 使用方法

1. 将 IQ 录音文件放到 `data/` 目录
2. 运行主程序 `main.m`

所有相关参数（采样率、符号率、滤波器等）都在 `getconfig.m` 中有说明，可在源码中按需修改。

## 作者

Anton Valentin Dilg, Ramin Leon Neymeyer, Ramon Rennert, Steffen August Sigwart

柏林工业大学 电气工程与计算机科学学院

## 许可证

本项目采用 MIT 许可证。
