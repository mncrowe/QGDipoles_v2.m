function fig = splot(f, x, y, options)
% Plots f(x, y, n) as a multi-panel plot with n panels
%
% - f: f(x, y, n), array for n sub-figures
% - x: x array for each sub-figure (default: 1:Nx)
% - y: y array for each sub-figure (default: 1:Ny)
% - options:
%   - xlabel: x caption, string
%   - ylabel: y caption, string
%   - title: string array of titles for each figure
%   - type:   0 - surface plot
%             1 - pcolor plot (default)
%             n - (> 1) contour plot with n contours
%   - fig_size: [s_x s_y] for width s_x and height s_y, [0 0] gives default size
%   - close_all: true - close all current figures, false - do nothing (default)
%   - panel_array: [n_x n_y] for n_x by n_y array of subfigures
%   - new_fig: true - open new figure (default), false - use current figure window
%   - font_size: size of font for axis labels (default: 12)
%   - show_contours: number of contours on plot (default: 0)
%   - line_width: width of all plot lines and contours (default: 0.5)
%   - xlim: x axis limits (default: [min(x) max(x)]
%   - xlim: y axis limits (default: [min(y) max(y)]
%   - clim: colormap limits, array of size [n 2], if 0 limits set automatically
%   - circle_radius: if nonzero, add circle to figure of specified radius
%   - contour_colour: colour of contours, character array (default: 'k')
%   - shading: shading style, character array (default: 'interp')
%   - interpreter: interpreter for text, 'none', 'tex' (default) or 'latex'
%   - interior_labels: use xlabels and ylabels inside subplot array (default: true)
%   - label_subplots: labels subplots with 'a).', 'b).', etc (default: false)
%   - bold_labels: subplot labels are shown in bold (default: false)
%
% -------------------------------------------------------------------------
% Note: The array f should have first dimension corresponding to x, second
%       dimension to y. If f is entered with x and y flipped f will be 
%       transposed unless length(x) = length(y).
% -------------------------------------------------------------------------

arguments
    f                       (:,:,:) double
    x                       (:,1)   double  = 1:size(f, 1)
    y                       (:,1)   double  = 1:size(f, 2)
    options.xlabel          (1,1)   string  = ""
    options.ylabel          (1,1)   string  = ""
    options.title           (1,:)   string  = ""
    options.type            (1,:)   double  = 1
    options.fig_size        (1,2)   double  = [0, 0]
    options.close_all       (1,1)   logical = false
    options.panel_array     (1,:)   int64   = 0
    options.new_fig         (1,1)   logical = true
    options.font_size       (1,1)   int64   = 12
    options.show_contours   (1,:)   double  = 0
    options.line_width      (1,1)   double  = 0.5
    options.xlim            (1,2)   double  = [min(x) max(x)]
    options.ylim            (1,2)   double  = [min(y) max(y)]
    options.clim            (:,2)   double  = zeros(size(f, 3), 2);
    options.circle_radius   (1,1)   double  = 0
    options.contour_colour  (1,:)   char    = 'k'
    options.shading         (1,:)   char    = 'interp'
    options.interpreter     (1,:)   char    = 'tex'
    options.interior_labels (1,1)   logical = true
    options.plot_label      (1,:)   string  = ""
    options.label_subplots  (1,1)   logical = false
    options.bold_labels     (1,1)   logical = false
    options.flip_axis       (1,2)   logical = [false false]
    options.x_ticks         (1,:)   double  = []
    options.y_ticks         (1,:)   double  = []
    options.label_start     (1,1)   double  = 1
end

% Close existing figures and open new figure if required:

if options.close_all; close all; end
if options.new_fig; figure; end

fig = gcf;

% Get size of f:

[Nx, Ny, Np] = size(f);

% Check size of x and y match size of first two dimensions of f:

if length(x) ~= Nx; error('Length of x must match the first dimension of f'); end
if length(y) ~= Ny; error('Length of y must match the second dimension of f'); end

% Set clim array to match number of subpanels if only one set of limits given:

if size(options.clim, 1) < Np; options.clim = repmat(options.clim(1, :), Np, 1); end

% Set array entries outside xlim and ylim to NaN to avoid incorrect colorbar:

X_lim = options.xlim + abs(x(2) - x(1)) * [-1 1];
Y_lim = options.ylim + abs(y(2) - y(1)) * [-1 1];

exterior = x < X_lim(1) | x > X_lim(2) | y' < Y_lim(1) | y' > Y_lim(2);

f(repmat(exterior, 1, 1, Np)) = NaN;

% Get size of subfigure array (required only for Np > 1):

if isequal(options.panel_array, 0)
    [nx, ny] = sq_fac(Np);
else
    nx = double(options.panel_array(1)); ny = double(options.panel_array(2));
end

% Set figure window size and swap nx and ny if the figure is longer vertically:

if ~isequal(options.fig_size, [0, 0])
    pos = get(gcf, 'position');
    pos(3) = options.fig_size(1);
    pos(4) = options.fig_size(2);
    set(gcf, 'position', pos)
    if pos(4) > pos(3); [nx, ny] = deal(ny, nx); end
end

% Extend options.type to all subfigures:

if length(options.type) < Np; options.type = repmat(options.type(1), [1 Np]); end

% Create figure:

if Np == 1

    % Creat surface plot, pseudocolour plot or contour plot depending on type:

    hold on

    switch options.type
        case 0
            surf(x, y, f')
        case 1
            pcolor(x, y, f')
        otherwise
            contour(x, y, f', options.type, EdgeColor = options.contour_colour, LineWidth = options.line_width)
    end

    if options.show_contours > 0 | numel(options.show_contours) > 1
        contour(x, y, f', options.show_contours, EdgeColor = options.contour_colour, LineWidth = options.line_width)
    end

    hold off

    % Align plot, set axis limits, add colourbar and interpolated shading:

    view(0, 90)
    shading(options.shading)
    cb = colorbar;
    xlim(options.xlim)
    ylim(options.ylim)
    colormap(cmap())
    if any(options.clim); clim(options.clim); end
    set(gca, TickLabelInterpreter = options.interpreter)
    set(cb, TickLabelInterpreter = options.interpreter)

    if numel(options.x_ticks) > 0; set(gca, 'XTick', options.x_ticks); end
    if numel(options.y_ticks) > 0; set(gca, 'YTick', options.y_ticks); end

    if options.flip_axis(1); set(gca, 'XDir','reverse'); end
    if options.flip_axis(2); set(gca, 'YDir','reverse'); end

    % Add x and y axis labels and title:

    if options.xlabel ~= ""; xlabel(options.xlabel, Interpreter = options.interpreter); end
    if options.ylabel ~= ""; ylabel(options.ylabel, Interpreter = options.interpreter); end
    if options.title ~= ""; title(options.title, Interpreter = options.interpreter); end

    % Add plot label:

    if options.plot_label ~= ""

        if options.bold_labels
            if strcmp(options.interpreter, 'latex')
                options.plot_label = ['\textbf{' options.plot_label '}'];
            end
            bold_or_normal = 'bold';
        else
            bold_or_normal = 'normal';
        end

        text(-0.1, 1.05, options.plot_label, Units = 'normalized', Interpreter = options.interpreter, ...
            FontWeight = bold_or_normal, FontSize = options.font_size, HorizontalAlignment = 'center')

    end

    % Set font and line size, add bounding rectangle:

     box off
     set(gca, FontSize = options.font_size);
     rectangle(Position = [options.xlim(1) options.ylim(1) options.xlim(2)-options.xlim(1) ...
        options.ylim(2)-options.ylim(1)], LineWidth = options.line_width)

     if options.circle_radius > 0
        rectangle(position = options.circle_radius * [-1 -1 2 2], curvature = [1 1], ...
            LineStyle = "--", LineWidth = options.line_width)
     end

else

    hold on;

    % Create Np subplots:

    for n = 1:Np

        subplot(ny, nx, n);

        % Creat surface plot, pseudocolour plot or contour plot depending on type:

        hold on

        switch options.type(n)
            case 0
                surf(x, y, f(:, :, n)')
            case 1
                pcolor(x, y, f(:, :, n)')
            otherwise
                contour(x, y, f(:, :, n)', options.type(n), EdgeColor = options.contour_colour, LineWidth = options.line_width)
        end

        if options.show_contours > 0 | numel(options.show_contours) > 1
            contour(x, y, f(:, :, n)', options.show_contours, EdgeColor = options.contour_colour, LineWidth = options.line_width)
        end

        hold off

        % Align plot, set axis limits, add colourbar and interpolated shading:

        view(0, 90)
        shading(options.shading)
        cb = colorbar;
        xlim(options.xlim)
        ylim(options.ylim)
        colormap(cmap())
        if any(options.clim(n, :)); clim(options.clim(n, :)); end
        set(gca, TickLabelInterpreter = options.interpreter)
        set(cb, TickLabelInterpreter = options.interpreter)

        if numel(options.x_ticks) > 0; set(gca, 'XTick', options.x_ticks); end
        if numel(options.y_ticks) > 0; set(gca, 'YTick', options.y_ticks); end

        if options.flip_axis(1); set(gca, 'XDir','reverse'); end
        if options.flip_axis(2); set(gca, 'YDir','reverse'); end

        % Add x and y axis labels and title:

        if options.xlabel ~= "" && (options.interior_labels || n > (ny - 1) * nx)
            xlabel(options.xlabel, Interpreter = options.interpreter);
        end

        if options.ylabel ~= "" && (options.interior_labels || rem(n, nx) == 1)
            ylabel(options.ylabel, Interpreter = options.interpreter);
        end

        if length(options.title) == Np && options.title(n) ~= ""
            title(options.title(n), Interpreter = options.interpreter);
        end

        % Add subplot labels:

        if options.label_subplots

            lab = [num2str(char(96 + n + options.label_start - 1)) ').'];

            if options.bold_labels
                if strcmp(options.interpreter, 'latex')
                    lab = ['\textbf{' lab '}'];
                end
                bold_or_normal = 'bold';
            else
                bold_or_normal = 'normal';
            end

            if ny == 1; y_pos = 1.05; else; y_pos = 1.1; end

            text(-0.3, y_pos, lab, Units = 'normalized', Interpreter = options.interpreter, ...
                FontWeight = bold_or_normal, FontSize = options.font_size, ...
                HorizontalAlignment = 'center')

        end
        
        % Set font and line size, add bounding rectangle:
        
        box off
        set(gca, FontSize = options.font_size);
        rectangle(Position = [options.xlim(1) options.ylim(1) options.xlim(2)-options.xlim(1) ...
            options.ylim(2)-options.ylim(1)], LineWidth = options.line_width)

        if options.circle_radius > 0
            rectangle(position = options.circle_radius * [-1 -1 2 2], curvature = [1 1], ...
                LineStyle = "--", LineWidth = options.line_width)
        end

    end

    % Add title to whole figure:

    if length(options.title) == 1; sgtitle(options.title, Interpreter = options.interpreter); end

    hold off;

end

if nargout == 0; clear fig; end

end

function [a, b] = sq_fac(n)
% returns the two factors, [a b], closest to sqrt(n), where a*b=n
%
% - n: number to factorise
%
% ----------------------------------------------------------------------------
% Note: cleverer methods should be used when dealing with large numbers
% ----------------------------------------------------------------------------

i = floor(sqrt(n));
while n / i ~= floor(n / i)
    i = i - 1; 
end
a = n / i;
b = i;

end

