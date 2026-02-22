function [y,isfftout] = helperOTFSmod(x)
%   Modulation of the input signal in the delay-Doppler domain
%
%   [Y,ISFFT] = HELPEROTFSMOD(X,PADLEN,PADTYPE) 使用矩形脉冲成形窗口对输入
%   Signal X is modulated using OTFS, and the output is Y. X is an array of real or complex numbers of size M x N.
%   M represents the number of subcarriers, and N represents the number of OTFS subsymbols.
%
%   Input：
%   X - OTFS grid (M x N)
%   Output：
%   Y - Time-domain output vector

M = size(x,1);
y = ifft(x.').' / M; % inverse Zak transform
isfftout = fft(y); % ISFFT
end