function [utcMid, MT_pos] = GMTI_return_EN(varargin)
% 点击图像像素，返回经纬度（四角经纬度做投影配准）
% 光标：红色十字准星，回车结束

    if numel(varargin) < 1 || isempty(varargin{1})
        result_dir = '/Volumes/cz/20251106/Mission080/result';
    else
        result_dir = varargin{1};
    end

    img_path    = fullfile(result_dir, 'GMTI07.png');
    corners_txt = fullfile(result_dir, 'GMTI07.txt');

    dbs_img = imread(img_path);
    [H, W, ~] = size(dbs_img);

    % 像素角点顺序：右上、右下、左下、左上
    img_corners = [W-1, 0;
                   W-1, H-1;
                   0,   H-1;
                   0,   0];

    ll_gps = get_coner_gps(corners_txt); % [L,B]，与上面顺序一致

    % 经纬度 -> 像素 的投影（局部近似平面）
    tform = fitgeotrans(ll_gps, img_corners, 'projective');

    % ==== 显示图像并启用红色十字准星 ====
    hfig = figure('Name','红色十字准星：点击选点，回车结束');
    imshow(dbs_img, 'Border','tight'); axis on; hold on;

    % 红色十字准星（两条线，随鼠标移动）
    hx = line([NaN NaN],[NaN NaN],'Color','r','LineWidth',1.2,'HitTest','off');
    hy = line([NaN NaN],[NaN NaN],'Color','r','LineWidth',1.2,'HitTest','off');

    % 存储状态
    S.W = W; S.H = H; S.hx = hx; S.hy = hy;
    S.xs = []; S.ys = []; S.ax = gca; S.done = false;
    guidata(hfig,S);

    set(hfig, 'WindowButtonMotionFcn', @onMove);
    set(hfig, 'WindowButtonDownFcn',  @onClick);
    set(hfig, 'KeyPressFcn',          @onKey);
    uiwait(hfig);              % 等待回车结束

    if ~ishandle(hfig) || ~isvalid(hfig)
        MT_pos = []; utcMid = [];
        return;
    end
    S = guidata(hfig);
    xs = S.xs; ys = S.ys;

    % ==== 像素 -> 经纬度 ====
    if isempty(xs)
        MT_pos = []; utcMid = [];
        return;
    end
    try
        [L, B] = transformPointsInverse(tform, xs, ys);
    catch
        tformInv = invert(tform);
        [L, B] = transformPointsForward(tformInv, xs, ys);
    end
    ll_click = [L(:), B(:)];

    % 标注（点用红色，文本保持黄色便于对比）
    for i = 1:numel(xs)
        plot(xs(i), ys(i), 'x', 'Color','r', 'MarkerSize',10, 'LineWidth',1.5);
        txt = sprintf('L=%.6f, B=%.6f', ll_click(i,1), ll_click(i,2));
        text(xs(i)+6, ys(i)+6, txt, 'FontSize',10, 'Color','y', 'Interpreter','none');
    end
    drawnow;

    fprintf('--- 点击结果（经度L，纬度B） ---\n');
    for i = 1:size(ll_click,1)
        fprintf('%2d) L=%.8f, B=%.8f\n', i, ll_click(i,1), ll_click(i,2));
    end

    utcMid = [];
    MT_pos = ll_click;

    % ===== 内部回调 =====
    function onMove(src, ~)
        if ~ishandle(src), return; end
        S = guidata(src); if isempty(S) || ~isfield(S,'ax'), return; end
        cp = get(S.ax, 'CurrentPoint'); x = cp(1,1); y = cp(1,2);
        if ~isfinite(x) || ~isfinite(y), return; end
        set(S.hx, 'XData', [1 S.W], 'YData', [y y]);   % 水平线
        set(S.hy, 'XData', [x x],   'YData', [1 S.H]); % 垂直线
    end

    function onClick(src, ~)
        S = guidata(src); cp = get(S.ax,'CurrentPoint');
        x = cp(1,1); y = cp(1,2);
        if x>=1 && x<=S.W && y>=1 && y<=S.H
            S.xs(end+1) = x; S.ys(end+1) = y;
            plot(x, y, 'x', 'Color','r', 'MarkerSize',10, 'LineWidth',1.5);
            guidata(src, S);
        end
    end

    function onKey(src, evt)
        if any(strcmp(evt.Key, {'return','enter','escape'}))
            uiresume(src);
        end
    end
end

function ll_gps = get_coner_gps(filename)
    lines = readlines(filename);
    coords = struct();
    for i = 1:length(lines)
        str = strtrim(lines(i));
        if str == "", continue; end
        tokens = regexp(str, '(\w+)\s*=\s*([\-+]?\d+(\.\d+)?)', 'tokens');
        if ~isempty(tokens)
            key = tokens{1}{1}; val = str2double(tokens{1}{2});
            coords.(key) = val;
        end
    end
    needKeys = {'B0','B1','B2','B3','L0','L1','L2','L3'};
    for k = 1:numel(needKeys)
        if ~isfield(coords, needKeys{k})
            error('缺少字段 %s，请检查 %s', needKeys{k}, filename);
        end
    end
    % 顺序：右上、右下、左下、左上（与 img_corners 一致）
    ll_gps = [coords.L0, coords.B0;
              coords.L2, coords.B2;
              coords.L1, coords.B1;
              coords.L3, coords.B3];
end
