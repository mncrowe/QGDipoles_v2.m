function PlotModon(modon, field, grid, params, options)
% Plots the specified field of a modon strucutre
%
% - modon: modon structure created using CreateModon functions
% - field: field to plot, options are: 'psi' (default), 'q', 'b', 'Psi'
% - grid: grid structure, include if not a field of modon
% - params: params structure, include if not a field of modon
% - options:
%   - xlim: x axis limits (default: [min(x) max(x)]
%   - ylim: y axis limits (default: [min(y) max(y)]
%   - contours: nunber/heights of contours to plot (default: 0)
%   - individual: show each layer in a separate window (default: false)
%   - index: layers to plot (default: 1:Nz)
%   - plane: plane to plot through, 'xy' or 'yz'
%
% -------------------------------------------------------------------------
% Note: here Psi is defined as: Psi = psi + (U - U_n)*y. For SQG or 3DQG
%       modons U_n = 0. Hence Psi represents the flow in the vortex frame.
% -------------------------------------------------------------------------

arguments
    modon                    struct
    field              (1,:) char    = 'psi'
    grid                     struct  = modon.grid
    params                   struct  = modon.params
    options.xlim       (1,2) double  = [min(grid.x) max(grid.x)]
    options.ylim       (1,2) double  = [min(grid.y) max(grid.y)]
    options.clim       (1,2) double  = [0 0]
    options.contours   (1,:) double  = 0
    options.individual (1,1) logical = false
    options.index      (1,:) double  = 1:params.Nz
    options.plane      (1,:) char    = 'xy'
    options.show_title (1,1) logical = true
end

% Set plotting field:

switch field
    case 'psi'
        F = modon.psi;
        title = '\psi&';
    case 'Psi'
        if isfield(params, 'U_n')
            U = params.U - params.U_n;
            title = '\Psi& = \psi& + (U - U_n)y';
        else
            U = params.U;
            title = '\Psi& = \psi& + Uy';
        end
        F = modon.psi + grid.y .* reshape(U, 1, 1, []);
    case 'q'
        F = modon.q;
        title = 'q&';
    case 'b'
        F = modon.b;
        title = 'b&';
    otherwise
        error('Unrecognised field.')
end

% Plot using 'splot.m':

switch options.plane
    case 'xy'

        % Create titles for xy plane plotting:

        if options.show_title
            titles = cell(length(options.index), 1);
            if isfield(grid, 'z')   % if z grid exists, use (x, y, z = z_value)
                for i = 1:length(options.index)
                    titles{i} = strrep(title, '&', ['(x, y, z = ' num2str(grid.z(options.index(i))) ')']);
                end
            else                    % if no z grid, use layer number, no number shown for 1 layer models
                for i = 1:length(options.index)
                    if params.Nz == 1
                        titles{i} = [strrep(title, '&', '(x, y)')];
                    else
                        titles{i} = [strrep(title, '&', '(x, y)') ' (Layer ' num2str(options.index(i)) ')'];
                    end
                end
            end
        end

        if options.individual   % plot figures in separate windows
            for i = 1:length(options.index)
                if options.show_title
                    title = titles{i};
                else
                    title = "";
                end
                splot(F(:, :, options.index(i)), grid.x, grid.y, xlabel = 'x', ylabel = 'y', ...
                    circle_radius = params.a, xlim = options.xlim, ylim = options.ylim, ...
                    show_contours = options.contours, title = title, clim = options.clim)
            end
        else                    % plot figure as a panel array in single window
            if options.show_title
                title = titles;
            else
                title = "";
            end
            splot(F(:, :, options.index), grid.x, grid.y, xlabel = 'x', ylabel = 'y', ...
                circle_radius = params.a, xlim = options.xlim, ylim = options.ylim, ...
                show_contours = options.contours, title = title, clim = options.clim)
        end

    case 'yz'

        % Create titles for yz plotting:

        if isfield(grid, 'z')
            z = grid.z;
            title = strrep(title, '&', '(x = 0, y, z)');
            z_label = 'z';
        else
            z = 1:params.Nz;
            title = strrep(title, '&', '(x = 0, y, n)');
            z_label = 'n (Layer)';
        end

        % Plot figure in single window:

        if ~options.show_title
            title = "";
        end

        splot(squeeze(F(params.Nx/2 + 1, :, :)), grid.y, z, xlabel = 'y', ylabel = z_label, ...
            xlim = options.ylim, show_contours = options.contours, title = title, ...
            clim = options.clim)

        % Set z axis to layer number for layered models:

        if ~isfield(grid, 'z')
            set(gca, 'YTick', 1:params.Nz)
            set(gca, 'YDir','reverse')
        end

    otherwise

        error('Unrecognised value of ''plane''.')

end

end