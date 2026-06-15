function [utcMid, MT_pos_all] = GMTI_result_show_CPP(varargin)
% GMTI_result_show_CPP 批量显示并保存 GMTI 结果
% [utcMid, MT_pos_all] = GMTI_result_show_CPP()
% [utcMid, MT_pos_all] = GMTI_result_show_CPP(result_dir)
% [utcMid, MT_pos_all] = GMTI_result_show_CPP(result_dir, idx_range)
%
% 默认:
%   result_dir = '/Volumes/cz/20251108/Mission068/result3';
%   idx_range  = 1:99;  % 对应 GMTI1~GMTI99

    %% 输入参数
    % close all;
    if numel(varargin) < 1
        result_dir = '/home/shy/AIR/小长/GMTI_Data/Mission068/result3';
    else
        result_dir = varargin{1};
    end

    if numel(varargin) < 2
        idx_range = 13:23;   % 需要的序号范围
    else
        idx_range = varargin{2};
    end
    tracklen = 3;

    % 输出：每个序号一个 cell
    utcMid      = cell(1, numel(idx_range));
    MT_pos_all  = cell(1, numel(idx_range));

    %% 循环处理每个 GMTIxx
    % for kk = 1:numel(idx_range)
    for kk = 1:1
        idx = idx_range(kk);

        % 文件基名: GMTIxx
        baseName      = sprintf('GMTI%02d', idx);  % 若你是两位数: sprintf('GMTI%02d',idx)
        img_path      = fullfile(result_dir, [baseName '.png']);
        corners_txt   = fullfile(result_dir, [baseName '.txt']);
        MTresults_bin = fullfile(result_dir, [baseName '.bin']);

        % 检查文件是否存在
        if ~isfile(img_path) || ~isfile(corners_txt) || ~isfile(MTresults_bin)
            warning('GMTI_result_show_CPP:MissingFile', ...
                '序号 %d 文件缺失，跳过。(png/txt/bin 中至少有一个不存在)', idx);
            continue;
        end

        %% 读取并准备图像
        dbs_img = imread(img_path);
        [H, W, ~] = size(dbs_img);

        % 图像角点像素坐标（右上、右下、左下、左上）
        img_corners = [...
            W - 1, 0;      % 右上
            W - 1, H - 1;  % 右下
            0,     H - 1;  % 左下
            0,     0];     % 左上

        %% 读取四角经纬度 [lon, lat]，需与上面像素角点顺序一致
        ll_gps = get_coner_gps(corners_txt);  % 期望 4x2，列为 [lon, lat]

        % 构建 GPS->[像素] 的投影变换（经纬度当作平面坐标使用，局部范围可行）
        tform = fitgeotrans(ll_gps, img_corners, 'projective');

        %% 读取 GMTI 结果
        MT_combined = struct('id',[], 'lon',[], 'lat',[], 'utcMid',[], 'count',0);
        for idx = idx_range
            baseName = sprintf('GMTI%02d', idx);
            binPath = fullfile(result_dir, [baseName '.bin']);
            if isfile(binPath)
                out = read_mt_result_bin(binPath);
                MT_combined.id = [MT_combined.id, out.id];
                MT_combined.lon = [MT_combined.lon, out.lon];
                MT_combined.lat = [MT_combined.lat, out.lat];
                MT_combined.utcMid = [MT_combined.utcMid, repmat(out.utcGlobal, 1, out.count)];
                MT_combined.count = MT_combined.count + out.count;
            end
        end
        [S, ~] = group_mt_by_utcmid(MT_combined);

        % % 保存 utcMid（每个目标一个向量）
        % utcMid{kk} = MT_out.utcMid;

        %% 画图（使用不可见 figure，画完直接保存）
        fig = figure('Visible','off');
        imshow(dbs_img); hold on; axis on;
        title(sprintf('Multiple Period MT Positions on DBS (GMTI%02d)', idx), ...
              'Interpreter','none');
        ax = gca;

        % 颜色/标记循环
        markers = {'o','s','d','^','v','>','<','p','h','x','+'};
        nMk = numel(markers);

        % 注意：这里改名为 MT_pos_this，避免和函数输出变量冲突
        MT_pos_this = cell(1, numel(S));
        % legend_entries = gobjects(numel(S),1);

        for ii = 1:numel(S)
            temp = zeros(S(ii).count, 4);
            for i = 1:S(ii).count
                [temp(i, 1), temp(i, 2)] = Gaussp3(S(ii).lat(i), S(ii).lon(i), 117);
                temp(i, 3:4) = [S(ii).lat(i), S(ii).lon(i)];
            end
            MT_pos_this{ii} = temp;

            % 经纬度 [lat, lon]
            % pts_gps = [S(ii).lat(:),S(ii).lon(:)];

            % 投影到图像坐标
            % [pts_x, pts_y] = transformPointsForward(tform, pts_gps(:,2), pts_gps(:,1));

            % 只画散点
            % mk = markers{mod(ii-1, nMk)+1};
            % legend_entries(ii) = plot(pts_x, pts_y, mk, 'MarkerSize', 6, 'LineWidth', 1.2);
        end

        % 轨迹卡尔曼滤波 + 绘制
        tracks = kalman_track(MT_pos_this, 6.65625, 50, 250, 3, 1);

        plot_tracks_on_dbs(ax, tracks, tform, struct( ...
            'LineWidth', 1.0, ...
            'Marker', 'o', ...
            'MarkerSize', 5, ...
            'ShowEndpoints', true, ...
            'LegendBase', 'Track ' ...
        ));

        hold off;

        % 保存轨迹结果
        MT_pos_all{kk} = MT_pos_this;

        % %% 保存图像为 AGMTIxx.png
        % out_png = fullfile(result_dir, sprintf('AGMTI%02d.png', idx)); % 如要两位数可改为 %02d
        % try
        %     exportgraphics(fig, out_png, 'Resolution', 200);
        % catch
        %     % 旧版 MATLAB 没有 exportgraphics，用 saveas 兜底
        %     saveas(fig, out_png);
        % end
        % 
        % % close(fig);
        % % fprintf('已生成: %s\n', out_png);
    end

end



function ll_gps = get_coner_gps(filename)
% 读取txt文件
lines = readlines(filename);

% 初始化存储
coords = struct();

% 遍历每一行解析
for i = 1:length(lines)
    str = strtrim(lines(i)); % 去掉首尾空格
    if str == ""  % 跳过空行
        continue;
    end
    tokens = regexp(str, '(\w+)\s*=\s*([\d\.\-]+)', 'tokens'); % 匹配键值
    if ~isempty(tokens)
        key = tokens{1}{1};
        val = str2double(tokens{1}{2});
        coords.(key) = val; % 存到结构体
    end
end

% 访问数据
B0 = coords.B0;
B1 = coords.B1;
B2 = coords.B2;
B3 = coords.B3;
L0 = coords.L0;
L1 = coords.L1;
L2 = coords.L2;
L3 = coords.L3;

ll_gps = [L0, B0;L2, B2;L1, B1;L3, B3];
% 打印检查
disp(coords);
end

%% ---------- 辅助函数：读取单个 period 文件 ----------
function out = read_mt_result_bin(binFile)
% read_mt_result_bin  读取 GMTI 目标结果二进制（固定 755 字节）
% 返回:
%   out.count      : 目标数（已截断到<=40）
%   out.id         : [1 x N] 目标编号（uint8 -> double）
%   out.lon        : [1 x N] 经度（double）
%   out.lat        : [1 x N] 纬度（double）
%   out.utcMid     : [1 x N] 每目标 utcMid（double）
%   out.utcGlobal  : 标头 681 处全局 utcMid（double）

    MAX_TGT    = 40;
    REC_BYTES  = 17;
    OFFSET_681 = 1 + MAX_TGT*REC_BYTES;   % = 681
    FILE_BYTES = OFFSET_681 + 4 + 70;     % = 755

    fid = fopen(binFile,'rb','ieee-le');  % 按小端读取（与写入一致）
    assert(fid>0, '无法打开文件: %s', binFile);
    c = onCleanup(@() fclose(fid));

    % 可选：检查文件长度（不强制）
    fseek(fid, 0, 'eof');
    fsz = ftell(fid);
    fseek(fid, 0, 'bof');
    if fsz < FILE_BYTES
        warning('文件长度(%d)小于预期 755 字节，尝试按已有内容解码。', fsz);
    end

    % [0] 目标总数
    N = fread(fid, 1, 'uint8=>double');
    if isempty(N), error('文件内容不足：缺少总数字节'); end
    N = min(N, MAX_TGT);

    % 预读 40 个槽位（逐条取）
    id     = zeros(1, N);
    lon    = zeros(1, N);
    lat    = zeros(1, N);
    utcMid = zeros(1, N);

    for i = 1:N
        % 每条记录固定 17 字节
        id(i)     = fread(fid, 1, 'uint8=>double');
        lon(i)    = fread(fid, 1, 'int32=>double')*8.38191e-8;
        lat(i)    = fread(fid, 1, 'int32=>double')*8.38191e-8;
        utcMid(i) = fread(fid, 1, 'single=>double');
        % 跳过 4 字节保留
        fseek(fid, 4, 'cof');

        % 若文件提前结束，fread 会返回空；做兜底处理
        if any(isempty([id(i),lon(i),lat(i),utcMid(i)]))
            % 用已有长度截断
            N = min(N, i-1);
            break;
        end
    end

    % 偏移 681：全局 utcMid
    fseek(fid, OFFSET_681, 'bof');
    utcGlobal = fread(fid, 1, 'single=>double');
    if isempty(utcGlobal), utcGlobal = NaN; end

    % 输出（按实际 N 截断）
    out.count     = N;
    out.id        = id(1:N);
    out.lon       = lon(1:N);
    out.lat       = lat(1:N);
    out.utcMid    = utcMid(1:N);
    out.utcGlobal = utcGlobal;
end

%%
function hList = plot_tracks_on_dbs(ax, tracks, tform, opts)
% plot_tracks_on_dbs  在DBS图像上绘制轨迹（连线+点）
%
% 输入:
%   ax     : 坐标轴句柄（imshow 后的坐标轴）
%   tracks : 轨迹集合，支持以下任一形式：
%            1) cell 数组，每个单元为 N×2 或 N×3 矩阵，列为 [lat lon (alt)]
%            2) 结构体数组，每个元素含字段：
%               - LL : N×2 或 N×3（[lat lon (alt)]）
%               （如字段名不同，可在下方自行改成你的字段名）
%   tform  : fitgeotrans 得到的投影（经纬度→像素），调用方式 transformPointsForward(tform, lon, lat)
%   opts   : 可选绘图选项结构体，字段（均可省略）：
%               LineWidth     (默认 1.5)
%               Marker        (默认 'o')
%               MarkerSize    (默认 5)
%               ColorOrder    (默认自动 from axes)
%               ShowEndpoints (默认 true；起点圈、终点叉)
%               LegendBase    (默认 'Track ')
%
% 输出:
%   hList  : 每条轨迹主线的图形句柄数组（便于做图例）

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'LineWidth'),      opts.LineWidth = 1.5; end
    if ~isfield(opts,'Marker'),         opts.Marker = 'o';    end
    if ~isfield(opts,'MarkerSize'),     opts.MarkerSize = 5;  end
    if ~isfield(opts,'ShowEndpoints'),  opts.ShowEndpoints = true; end
    if ~isfield(opts,'LegendBase'),     opts.LegendBase = 'Track '; end

    axes(ax); hold(ax, 'on');
    % 统一把 tracks 变成 cell，每个 cell 是 N×2/3 矩阵
    if isstruct(tracks)
        % 假设字段名为 LL（如不同，改这里）
        tracks_cell = arrayfun(@(s) s.pos, tracks, 'UniformOutput', false);
        % tracks_cell = arrayfun(@(s) s.latlon, tracks, 'UniformOutput', false);
    elseif iscell(tracks)
        tracks_cell = tracks;
    else
        error('tracks 需为 cell 或 struct 数组。');
    end

    % 颜色循环
    co = get(ax, 'ColorOrder');
    Cn = size(co,1);

    hList = gobjects(numel(tracks_cell),1);
    shortestLength = 2;
    for i = 1:numel(tracks_cell)
        Ti = tracks_cell{i};
        if isempty(Ti), continue; end
        if size(Ti,2) < shortestLength || size(Ti,1) < shortestLength
            warning('第 %d 条轨迹列数不足 %d ，跳过。', i, shortestLength);
            continue;
        end
        disp("yeah!");

        % 只取 [lat lon]
        lat = Ti(:,1);
        lon = Ti(:,2);
        % disp(lat);disp(lon);
        for ii = 1:size(Ti, 1)
            [lat(ii), lon(ii)] = Gaussp3RV(Ti(ii, 1), Ti(ii, 2), 117);
            disp([lat(ii), lon(ii)]);
        end

        % 经纬度 -> 图像像素坐标
        [x, y] = transformPointsForward(tform, lon, lat);

        % 修改后的主线绘制：连线 + 实心点
        ci = co(mod(i-1, Cn)+1, :);
        plot(ax, x, y, '-o', ...              % '-' 是线，'o' 是圆圈点
            'LineWidth',  opts.LineWidth, ...
            'MarkerSize', opts.MarkerSize, ... % 中间点的大小
            'MarkerFaceColor', ci, ...        % 使点变成实心的
            'Color',      ci, ...
            'DisplayName', sprintf('%s%d', opts.LegendBase, i));

        % 起止点标注（保持不变，用于强化显示起点圆、终点叉）
        if opts.ShowEndpoints && numel(x) >= 1
            % 起点：大圆圈
            plot(ax, x(1),  y(1),  'o', 'MarkerSize', opts.MarkerSize+4, ...
                'LineWidth', 2.5, 'Color', ci, 'HandleVisibility', 'off');
            % 终点：叉号
            plot(ax, x(end), y(end), 'x', 'MarkerSize', opts.MarkerSize+5, ...
                'LineWidth', 2.5, 'Color', ci, 'HandleVisibility', 'off');
        end
    end
end
