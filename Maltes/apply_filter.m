function data = apply_filter(data, ftype, soft, hard, xlim)
%GET_FILTER Filter seismic data
%   Filter seismic data with different filter settings. Applys brickwall
%   filter, linear slope filter, or sinusoidal filter between hard and
%   soft filter frequency.

% sets frequencys below soft filter to 0
data(1:soft,1:xlim) = 0;

% brickwall filter
% sets freqs below hard filter to 0
if ftype == 'bri'
    data(1:hard,1:xlim) = 0;

% linear filter
% sets freqs between soft and hard filter to linear slope
elseif ftype == 'lin'
    data(soft:hard,1:xlim) ...
    = data(soft:hard,1:xlim) ...
    .* ((1:1+hard-soft)/soft).';

% sinusoidal filter
% sets freqs between soft and hard filter to sinusoidal slope
elseif ftype == 'sin'
    data(soft:hard,1:xlim) ...
    = data(soft:hard,1:xlim) ...
    .* sin(((1:1+hard-soft)*pi/(soft*2))).';
end
end
