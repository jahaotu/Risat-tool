function Riser_analysis_tool2
    % This Riser tool-RisAT was developed by 
    %Engr.E.A.A and Jonadab A
    % this is a beta version of the program as we hope to improve it in the
    % future and include more configuration 
    %MAIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%%%%%
    %
    % Create the main figure without the default menu bar
    fig = figure('Name', 'Riser Analysis', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 720], 'MenuBar', 'none', 'CloseRequestFcn', @saveDataOnClose, 'Color', [0.53, 0.81, 0.98]);
    %fig.Resize ='off';
    % Global variables for curve data
    global xf yf zf x1 y1 z1;
    
    % Initialize curve data
    % free hanging riser
    xf = [-50.3118, -27.0447, -3.6775, 20.0235, 37.7159, 45.0599];
    yf = [-0.6843, -0.7158, -0.7185, -0.7158, -0.7158, -0.7158];
    zf = [0, 1.2701, 2.0734, 3.1281, 7.0029, 10.9186];
    % lazy wave riser curve
    x1 = [42.3894, 40.3865, 30.3719, 26.0323, 17.3531, 9.6752, -0.6731, -36.7254];
    y1 = [-0.7158, -0.7158, -0.7158, -0.7158, -0.7158, -0.7158, -0.7158, -0.7158];
    z1 = [8.1879, 3.6032, 0.9732, 1.0279, 4.0968, 4.8602, 0.8848, 0.8026];
    
    
    % Initialize data storage struct
    appData = struct('Parameters', [], 'Results', [], 'LastSimulation', []);
    if isfile('riser_analysis_data.mat')
        load('riser_analysis_data.mat', 'appData');
    end
    
    % Create the custom menu
    newSimMenu = uimenu(fig, 'Text', '<html><b>New Simulation</b></html>', 'ForegroundColor', 'blue');
    uimenu(newSimMenu, 'Text', 'Create New Project', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @createNewProject);
    uimenu(newSimMenu, 'Text', 'Enter Simulation Parameters', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @enterParameters);
    
    SimType = uimenu(fig, 'Text', '<html><b>Simulation Type</b></html>', 'ForegroundColor', 'blue');
    static = uimenu(SimType, 'Text', 'Static', 'ForegroundColor', 'blue');
    uimenu(static, 'Text', 'Tension', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @(src, event)visualizeResults('tension'));
    uimenu(static, 'Text', 'Bending Moment', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @(src, event)visualizeResults('bendingMoment'));
    uimenu(static, 'Text', 'Shear Force', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @(src, event)visualizeResults('shearForce'));
    uimenu(static, 'Text', 'Modal Analysis', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @modalAnalysisPrompt);
    
    dynamic = uimenu(SimType, 'Text', 'Dynamic', 'ForegroundColor', 'blue');
    uimenu(dynamic, 'Text', 'Displacement', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @(src, event)dynamicAnalysisPrompt('displacement'));
    uimenu(dynamic, 'Text', 'Velocity', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @(src, event)dynamicAnalysisPrompt('velocity'));
    uimenu(dynamic, 'Text', 'Acceleration', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @(src, event)dynamicAnalysisPrompt('acceleration'));
    
    uimenu(fig, 'Text', '<html><b>Corrosion Analysis</b></html>', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @corrosionAnalysis);
    
    Visualization = uimenu(fig, 'Text', '<html><b>Visualizations</b></html>', 'ForegroundColor', 'blue');
    uimenu(Visualization, 'Text', 'Static Plots', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @visualizeStaticPlots);
    uimenu(Visualization, 'Text', 'Dynamic Plots', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @visualizeDynamicPlots);
    uimenu(Visualization, 'Text', 'Restore', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @restoreView);
    
    uimenu(fig, 'Text', '<html><b>Report</b></html>', 'ForegroundColor', 'blue', 'MenuSelectedFcn', @createPDFReport);

    % Axes for visualizing results
    ax = axes('Parent', fig, 'Position', [0.1, 0.5, 0.8, 0.4]);
    
    % Text area for displaying results
    textArea = uicontrol('Parent', fig, 'Style', 'edit', 'Max', 10, 'Min', 1, 'Enable', 'inactive', 'BackgroundColor', 'white', 'HorizontalAlignment', 'left', 'Position', [0.1 * fig.Position(3), 0.1 * fig.Position(4), 0.8 * fig.Position(3), 0.3 * fig.Position(4)], 'String', 'Text Area for simulation results in tabular form for each analysis');

    % Create Simulate Button
    uicontrol('Parent', fig, 'Style', 'pushbutton', 'Position', [0.9 * fig.Position(3) - 80, 0.05 * fig.Position(4), 80, 30], 'String', 'Simulate', 'BackgroundColor', 'green', 'ForegroundColor', 'white', 'Callback', @computeSimulations);

    % Function to handle creating a new project
    function createNewProject(src, event)
        prompt = {'Enter project name:'};
        title = 'New Project';
        dims = [1 35];
        definput = {'Project1'};
        answer = inputdlg(prompt, title, dims, definput);
        
        if ~isempty(answer)
            projectName = answer{1};
            projectFolder = fullfile(pwd, projectName);
            if ~exist(projectFolder, 'dir')
                mkdir(projectFolder);
            end
            
            appData.CurrentProject = projectName;
            appData.ProjectFolder = projectFolder;
            
            % Initialize Parameters and Results
            appData.Parameters = struct();
            appData.Results = struct();
        end
    end

    % Function to enter simulation parameters
    function enterParameters(src, event)
        if ~isfield(appData, 'CurrentProject')
            errordlg('Please create a new project first.', 'Error');
            return;
        end
        
        prompt = {'Lu', 'Lb', 'ks', 'D', 't', 'rho_s', 'rho_w', 'Cn', 'Ctau', 'g', 'P0', 'v0', 'hw', 'E', 'I', 'N', 'A'};
        title = 'Enter Simulation Parameters';
        dims = [1 35];
        definput = {'1800', '600', '4000', '0.2032', '0.0191', '7850', '1030', '1.2', '0.024', '9.81', '200000', '0.5', '2000', '2.1e11', '', '34', ''};
        answer = inputdlg(prompt, title, dims, definput);
        
        if ~isempty(answer)
            for i = 1:length(prompt)
                appData.Parameters.(prompt{i}) = str2double(answer{i});
            end
            save(fullfile(appData.ProjectFolder, 'parameters.mat'), 'appData');
        end
    end

    % Function to handle different types of simulations
    function computeSimulations(src, event)
        if ~isfield(appData, 'CurrentProject') || isempty(appData.Parameters)
            errordlg('Please create a new project and enter parameters first.', 'Error');
            return;
        end
        
        % Simulation code here
        simulateRiser;
    end

    % Function to visualize results
    function visualizeResults(type)
        if ~isfield(appData, 'Results') || isempty(appData.Results)
            errordlg('Please run the simulation first.', 'Error');
            return;
        end
        
        data = appData.Results.(type);
        plot(ax, data.x, data.y,'LineWidth',1.5);
        switch type
            case 'tension'
                ylabel('Tension (N)');
            case 'bendingMoment'
                ylabel('Bending Moment (N.m)');
            case 'shearForce'
                ylabel('Shear Force (N)');
        end
        xlabel('Arc length of the Riser (m)');
        title([type, ' Results']);
    end


    % Function to visualize static plots
    function visualizeStaticPlots(src, event)
        stl_file = 'shipModel.stl';
        if ~isfield(appData, 'Results') || isempty(appData.Results)
            errordlg('Please run the simulation first.', 'Error');
            return;
        end
        
        figure;
        subplot(2, 2, 1);
        plot(appData.Results.tension.x, appData.Results.tension.y,'LineWidth',1.5);
        xlabel('Arc length of the Riser(m)')
        ylabel('Tension of the riser (N)')
        title('Tension along Riser');
        
        
        subplot(2, 2, 2);
        plot(appData.Results.bendingMoment.x, appData.Results.bendingMoment.y,'LineWidth',1.5);
        xlabel('Arc length of the Riser(m)')
        ylabel('Bending Moment[N.m]')
        title('Bending Moment along Riser');
        
        subplot(2, 2, 3);
        plot(appData.Results.shearForce.x, appData.Results.shearForce.y,'LineWidth',1.5);
        xlabel('Arc Length of the Riser (m)')
        ylabel('Shear Force [N]')
        title('Shear Force along Riser ');

        
        subplot(2, 2, 4);
       
       visualize_vessel_with_riser(appData.Results.modalAnalysis.configuration);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        view(3);
        %axis('off')
        hold off;
    end

    % Function to visualize dynamic plots
    function visualizeDynamicPlots(src, event)
        if ~isfield(appData, 'Results') || isempty(appData.Results)
            errordlg('Please run the simulation first.', 'Error');
            return;
        end
        
        % Prompt user to select riser section
        prompt = {'Select riser section (top, middle, bottom):'};
        title = 'Select Riser Section';
        dims = [1 35];
        definput = {'middle'};
        answer = inputdlg(prompt, title, dims, definput);
        
        if ~isempty(answer)
            section = lower(answer{1});
            n_elements = appData.Parameters.N;
            
            % Determine the index based on the selected section
            switch section
                case 'top'
                    idx = 1;
                case 'middle'
                    idx = ceil(n_elements / 2);
                case 'bottom'
                    idx = (2*n_elements) + 2;
                otherwise
                    errordlg('Invalid section selected.', 'Error');
                    return;
            end
            
            figure;
            subplot(3, 1, 1);
            plot(appData.Results.displacement.x, appData.Results.displacement.y(idx, :),'LineWidth',1.5);
            xlabel('time [s]')
            ylabel('Displacement [m]')
            
            subplot(3, 1, 2);
            plot(appData.Results.velocity.x, appData.Results.velocity.y(idx, :),'LineWidth',1.5);
            xlabel('time [s]')
            ylabel('Velocity [m/s]')
            
            subplot(3, 1, 3);
            plot(appData.Results.acceleration.x, appData.Results.acceleration.y(idx, :),'LineWidth',1.5);
            xlabel('time [s]')
            ylabel('Acceleration [m/s^2]')
        end
    end

    % Function to restore the main view
    function restoreView(src, event)
        clf(fig);
        guidata2;
    end

    % Function to create a PDF report
    function createPDFReport(src, event)
        if ~isfield(appData, 'CurrentProject')
            errordlg('Please start a new simulation first.', 'Error');
            return;
        end

        projectName = appData.CurrentProject;
        projectFolder = appData.ProjectFolder;
        reportFile = fullfile(projectFolder, [projectName, '_report.pdf']);

        import mlreportgen.dom.*;
        import mlreportgen.report.*;

        rpt = Report(reportFile, 'pdf');
        titlepg = TitlePage('Title', 'Riser Analysis Report', 'Author', 'Simulation Tool');
        add(rpt, titlepg);
        toc = TableOfContents;
        add(rpt, toc);

        % Parameters Section
        ch = Chapter('Title', 'Simulation Parameters');
        paramTypes = fieldnames(appData.Parameters);
        for i = 1:length(paramTypes)
            p = Paragraph([paramTypes{i}, ': ', num2str(appData.Parameters.(paramTypes{i}))]);
            add(ch, p);
        end
        add(rpt, ch);

        % Results Section
        ch = Chapter('Title', 'Simulation Results');
        resultTypes = fieldnames(appData.Results);
        for i = 1:length(resultTypes)
            sec = Section('Title', resultTypes{i});
            results = appData.Results.(resultTypes{i});
            
            if isfield(results, 'x') && isfield(results, 'y')
                table = Table();
                table.TableEntriesHAlign = 'center';
                table.TableEntriesVAlign = 'middle';
                table.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
                
                % Add table header
                header = TableRow();
                header.Style = {Bold(true)};
                headerEntries = {'Arc Length(m)', 'Value'};
                for j = 1:length(headerEntries)
                    entry = TableEntry(headerEntries{j});
                    append(header, entry);
                end
                append(table, header);
                
                % Add table rows
                for j = 1:length(results.x)
                    row = TableRow();
                    entry1 = TableEntry(num2str(results.x(j)));
                    entry2 = TableEntry(num2str(results.y(j)));
                    append(row, entry1);
                    append(row, entry2);
                    append(table, row);
                end
                append(sec, table);
            elseif isfield(results, 'V') && isfield(results, 'Cr') % Corrosion results
                table = Table();
                table.TableEntriesHAlign = 'center';
                table.TableEntriesVAlign = 'middle';
                table.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
                
                % Add table header
                header = TableRow();
                header.Style = {Bold(true)};
                headerEntries = {'Velocity (m/s)', 'Corrosion Rate (m/s)'};
                for j = 1:length(headerEntries)
                    entry = TableEntry(headerEntries{j});
                    append(header, entry);
                end
                append(table, header);
                
                % Add table rows
                for j = 1:length(results.V)
                    row = TableRow();
                    entry1 = TableEntry(num2str(results.V(j)));
                    entry2 = TableEntry(num2str(results.Cr(j)));
                    append(row, entry1);
                    append(row, entry2);
                    append(table, row);
                end
                append(sec, table);
            end
            append(ch, sec);
        end
        append(rpt, ch);

        % Add plots
        ch = Chapter('Title', 'Simulation Plots');
        
        % Save plots to images and add them to the report
        plotTypes = {'curvature', 'bendingMoment', 'shearForce', 'tension', 'displacement', 'velocity', 'acceleration'};
        for i = 1:length(plotTypes)
            if isfield(appData.Results, plotTypes{i})
                figPlot = figure('Visible', 'off');
                plot(appData.Results.(plotTypes{i}).x, appData.Results.(plotTypes{i}).y,'LineWidth',1.5);
                title([plotTypes{i}, ' Plot']);
                xlabel('Riser Arc Length');
                ylabel(plotTypes{i});
                imgFile = fullfile(tempdir, [plotTypes{i}, '.png']);
                saveas(figPlot, imgFile);
                close(figPlot);
                
                img = Image(imgFile);
                img.Style = {ScaleToFit};
                append(ch, img);
            end
        end

        % Add vessel visualization and modal shapes
        if isfield(appData.Results, 'modalAnalysis')
            imgFileVessel = fullfile(tempdir, 'vessel_visualization.png');
            imgFileModal = fullfile(tempdir, 'modal_shapes.png');
            
            figVessel = figure('Visible', 'off');
            visualize_vessel_with_riser(appData.Results.modalAnalysis.configuration);
            vesselFrame = getframe(figVessel);
            imwrite(vesselFrame.cdata, imgFileVessel);
            close(figVessel);
            
            figModal = figure('Visible', 'off');
            visualize_modal_shapes(appData.Results.modalAnalysis.curve, appData.Results.modalAnalysis.M, appData.Results.modalAnalysis.K, 6, appData.Results.modalAnalysis.frequencies);
            modalFrame = getframe(figModal);
            imwrite(modalFrame.cdata, imgFileModal);
            close(figModal);

            imgVessel = Image(imgFileVessel);
            imgVessel.Style = {ScaleToFit};
            append(ch, imgVessel);
            
            imgModal = Image(imgFileModal);
            imgModal.Style = {ScaleToFit};
            append(ch, imgModal);
        end

        % Add corrosion plots
        if isfield(appData.Results, 'corrosion')
            imgFileCorrosion1 = fullfile(tempdir, 'corrosion_plot1.png');
            imgFileCorrosion2 = fullfile(tempdir, 'corrosion_plot2.png');
            
            figCorrosion1 = figure('Visible', 'off');
            plot(appData.Results.corrosion.V, appData.Results.corrosion.Cr,'LineWidth',1.5);
            xlabel('Velocity [m/s]');
            ylabel('Corrosion rate');
            title('Corrosion Rate vs. Velocity');
            saveas(figCorrosion1, imgFileCorrosion1);
            close(figCorrosion1);
            
            figCorrosion2 = figure('Visible', 'off');
            plot(appData.Results.corrosion.P, appData.Results.corrosion.Cr,'LineWidth',1.5);
            xlabel('Pressure');
            ylabel('Corrosion rate [m/s]');
            title('Corrosion Rate vs. Pressure');
            saveas(figCorrosion2, imgFileCorrosion2);
            close(figCorrosion2);

            imgCorrosion1 = Image(imgFileCorrosion1);
            imgCorrosion1.Style = {ScaleToFit};
            append(ch, imgCorrosion1);
            
            imgCorrosion2 = Image(imgFileCorrosion2);
            imgCorrosion2.Style = {ScaleToFit};
            append(ch, imgCorrosion2);
        end
        
        add(rpt, ch);
        
        % Save and close report
        close(rpt);
        rptview(reportFile);
    end

    % Function to handle corrosion analysis
    function corrosionAnalysis(src, event)
        if ~isfield(appData, 'CurrentProject')
            errordlg('Please start a new simulation first.', 'Error');
            return;
        end
        
        % Prompt user for corrosion analysis parameters
        prompt = {'Enter current density (I) in A:', 'Enter area (A) in m^2:', 'Enter charge number (z):', ...
                  'Enter pipe length (L) in m:', 'Enter Faraday constant (F) in C/mol:', 'Enter density of metal (q1) in kg/m^3:', ...
                  'Enter pipe diameter (d) in m:', 'Enter roughness height (h) in m:', 'Enter kinematic viscosity (U) in m^2/s:', ...
                  'Enter pressure (P) in Pa:', 'Enter density of fluid (q2) in kg/m^3:'};
        title = 'Corrosion Analysis Parameters';
        dims = [1 50];
        definput = {'110e-6', '1', '2', '2350', '96485.34', '7870', '0.305', '0.01', '0.00905', '18400000', '852.8'};
        answer = inputdlg(prompt, title, dims, definput);
        
        if ~isempty(answer)
            % Parse user inputs
            I = str2double(answer{1});
            A = str2double(answer{2});
            z = str2double(answer{3});
            L = str2double(answer{4});
            F = str2double(answer{5});
            q1 = str2double(answer{6});
            d = str2double(answer{7});
            h = str2double(answer{8});
            U = str2double(answer{9});
            P = str2double(answer{10});
            q2 = str2double(answer{11});
            
            % Compute corrosion rates
            V = 0.5:0.1:7;
            i = I / A;
            Cg = (i * A) / (z * F * q1);
            Re = (q2 * V * d) / U;
            Q = ((0.18 * U) / q2) * (((d / 2 - h) ./ (25 * d * Re .^ (-7/8))) .^ 3);
            SC = U ./ (Q * q2);
            f = 0.00137 * (1 + (20000 * h / d + 10000 ./ Re) .^ 0.33);
            sca = SC .^ (1/3);
            con = (1/8) * L * pi;
            co = con * sca;
            Cd = co .* (Q .* f .* Re);
            Ce = (2 * V .^ 3 * q2 .* f) / P;
            Cr = Cd + Ce + Cg;
            P = linspace(9000000, 20000000, size(Cr, 2));

            % Save corrosion analysis results to appData
            appData.Results.corrosion = struct('V', V, 'Cr', Cr, 'P', P);
            
            % Plot corrosion analysis results
            figure;
            subplot(1, 2, 1);
            plot(V, Cr,'LineWidth',1.5);
            xlabel('Velocity [m/s]');
            ylabel('Corrosion rate[mm/year]');
            
            subplot(1, 2, 2);
            plot(P, Cr,'LineWidth',1.5);
            xlabel('Pressure[N/m^2]');
            ylabel('Corrosion rate [mm/year]');
        end
    end

    % Function to save data when closing the app
    function saveDataOnClose(src, event)
        save('riser_analysis_data.mat', 'appData');
        delete(fig);
    end

    % Revised MATLAB Code for Riser Solutions
    function simulateRiser
        % Load Parameters
        params = appData.Parameters;
        
        % Derived parameters
        params.Ab = params.A;
        params.rho_b = params.rho_s - params.rho_w;
        params.wg = params.rho_b * params.g * params.A;
        params.wb = (params.rho_b - params.rho_w) * params.g * params.Ab;
        params.L0 = params.P0 / params.wg * sinh(acosh(params.wg / params.P0 * params.hw + 1));
        
        % Define constants
        pi = 3.141592653589793;
        E = params.E;
        I = pi / 4 * ((params.D / 2)^4 - (params.D / 2 - params.t)^4);
        ds = params.L0 / params.N;

        % Initial guess for the angles
        theta = 2 * pi - atan(params.wg / params.P0 * (0:params.N) * ds);

        % Function definitions
        y = @(i, theta) -ds * sum(sin(theta(i:params.N+1)));
        ksb = @(i, theta) params.ks * ds * (y(i, theta) >= params.hw);
        Ksb = @(i, theta) sum(arrayfun(@(j) ksb(j, theta), 1:i));

        x = @(i, theta) ds * sum(cos(theta(1:i)));

        vx = @(i, theta) params.v0 * (1 - (y(i, theta) / params.hw)^2);

        fn = @(i, theta) -0.5 * params.rho_w * params.Cn * params.D * ds * abs(vx(i, theta) * sin(theta(i))) * vx(i, theta) * sin(theta(i));
        ftau = @(i, theta) 0.5 * params.rho_w * params.Ctau * params.D * ds * abs(vx(i, theta) * cos(theta(i))) * vx(i, theta) * cos(theta(i));

        fx = @(i, theta) -fn(i, theta) * sin(theta(i)) + ftau(i, theta) * cos(theta(i));
        Fx = @(i, theta) sum(arrayfun(@(j) fx(j, theta), 1:i));

        fy = @(i, theta) fn(i, theta) * cos(theta(i)) + ftau(i, theta) * sin(theta(i));
        Fy = @(i, theta) sum(arrayfun(@(j) fy(j, theta), 1:i));

        w = @(i, theta) params.wg + (params.wb * ((i * ds) > (params.L0 - params.Lb - params.Lu) && (i * ds) <= (params.L0 - params.Lb)));
        W = @(i, theta) sum(arrayfun(@(j) w(j, theta), 1:i));

        % Objective function
        func = @(theta) [
            theta(1) - 2 * pi;
            arrayfun(@(i) E * I / ds * ((theta(i-1) - 2 * theta(i) + theta(i+1)) * (i < params.N)) + ...
                     ds * cos(theta(i)) * sum(arrayfun(@(j) min(Ksb(i, theta), Ksb(j, theta)) * ds * sin(theta(j)), 1:params.N)) + ...
                     params.hw * Ksb(i, theta) * ds * cos(theta(i)) + Fx(i, theta) * ds * sin(theta(i)) + ...
                     (Fy(i, theta) + W(i, theta)) * ds * cos(theta(i)) + params.P0 * sin(theta(i)), 2:params.N)'
        ];

        % Solve for lay angles using fsolve
        options = optimoptions('fsolve', 'Display', 'iter', 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000);
        theta = fsolve(func, theta, options);

        % Compute curvature, bending moment, shear force, and tension
        curvature = diff(theta) / ds; % Curvature k(s)
        curvature = [curvature, curvature(end)]; % Ensure same length as theta
        bending_moment = E * I * curvature; % Bending moment M(s)
        shear_force = diff(bending_moment) / ds; % Shear force S(s)
        shear_force = [shear_force, shear_force(end)]; % Ensure same length as theta

        % Compute tension T(s)
        tension = zeros(1, params.N+1);
        tension(1) =params.P0; 
        for i = 2:params.N+1
            fx_val = fx(i, theta);
            fy_val = fy(i, theta);
            w_val = w(i, theta);
            tension(i) = tension(i-1) - ds * (fx_val * cos(theta(i)) - (fy_val + w_val) * sin(theta(i)));
        end

        % Compute HOP for each strip
        hop_s = ds * (2 * pi - theta);
        hop = cumsum(hop_s); % Cumulative sum to get HOP along the riser length

        % Fit the data with spline functions
        spline_hop = linspace(min(hop), max(hop), 1000);
        shear_force_spline = spline(hop, shear_force, spline_hop);
        bending_moment_spline = spline(hop, bending_moment, spline_hop);
        tension_spline = spline(hop, tension, spline_hop);

        % Save results to appData
        appData.Results.curvature = struct('x', hop, 'y', curvature);
        appData.Results.bendingMoment = struct('x', spline_hop, 'y', bending_moment_spline);
        appData.Results.shearForce = struct('x', spline_hop, 'y', shear_force_spline);
        appData.Results.tension = struct('x', spline_hop, 'y', tension_spline);
        appData.Results.angles = struct('x', hop, 'y', theta);
        
        % Modal analysis using MATLAB built-in functions
        L = params.L0; % Length of the riser
        EI = E * I;
        rhoA = params.rho_s * params.A;
        n_elements = params.N;
        element_length = L / n_elements;

        % Beam element mass and stiffness matrices
        M = zeros(2*n_elements+2);
        K = zeros(2*n_elements+2);

        for i = 1:n_elements
            Me = (rhoA * element_length / 420) * ...
                [156, 22*element_length, 54, -13*element_length;
                 22*element_length, 4*element_length^2, 13*element_length, -3*element_length^2;
                 54, 13*element_length, 156, -22*element_length;
                 -13*element_length, -3*element_length^2, -22*element_length, 4*element_length^2];
            Ke = (EI / element_length^3) * ...
                [12, 6*element_length, -12, 6*element_length;
                 6*element_length, 4*element_length^2, -6*element_length, 2*element_length^2;
                 -12, -6*element_length, 12, -6*element_length;
                 6*element_length, 2*element_length^2, -6*element_length, 4*element_length^2];
            
            idx = [2*i-1, 2*i, 2*i+1, 2*i+2];
            M(idx, idx) = M(idx, idx) + Me;
            K(idx, idx) = K(idx, idx) + Ke;
        end

        % Solve the eigenvalue problem for natural frequencies and mode shapes
        [V, D] = eig(K, M);
        natural_frequencies = sqrt(diag(D)) / (2 * pi);

        % Save modal analysis results
        appData.Results.modalAnalysis = struct('x', 0:element_length:L, 'y', V(:, 1:6)); % First 6 modes
        appData.Results.modalAnalysis.curve = [x1; y1; z1]; % Example curve
        appData.Results.modalAnalysis.M = M;
        appData.Results.modalAnalysis.K = K;
        appData.Results.modalAnalysis.configuration = 'Lazy Wave'; % Example configuration
        appData.Results.modalAnalysis.frequencies = natural_frequencies(1:6); % Store first 6 modal frequencies

        % Newmark-beta explicit scheme for dynamic analysis
        beta = 0.25; % Newmark-beta parameter
        gamma = 0.5; % Newmark-beta parameter
        alpha = 0.02; % Damping coefficient
        time = linspace(0, 10, 1000); % Time vector for 10 seconds
        dt = time(2) - time(1); % Time step

        % Damping matrix
        B = M + alpha * K;

        % Initial conditions
        u = zeros(2*n_elements+2, length(time));
        v = zeros(2*n_elements+2, length(time));
        a = zeros(2*n_elements+2, length(time));

        % Initial force vector (using tension)
        F = zeros(2*n_elements+2, length(time));
        F(2:2:end, :) = repmat(tension', 1, length(time)); % Apply tension force

        % Initial acceleration
        a(:, 1) = M \ (F(:, 1) - B * v(:, 1) - K * u(:, 1));

        % Time-stepping loop
        waitbarHandle = waitbar(0, 'Starting simulations...');
        totalSteps = length(time);
        
        for t = 2:totalSteps
            % Newmark-beta equations
            u(:, t) = u(:, t-1) + dt * v(:, t-1) + (dt^2 / 2) * a(:, t-1);
            v(:, t) = v(:, t-1) + dt * a(:, t-1);
            a(:, t) = M \ (F(:, t) - B * v(:, t) - K * u(:, t));
            
            % Update waitbar
            waitbar(t / totalSteps, waitbarHandle, sprintf('Simulation Progress: %0.2f%%', (t / totalSteps) * 100));
        end
        
        close(waitbarHandle);

        % Save dynamic analysis results
        appData.Results.displacement = struct('x', time, 'y', u);
        appData.Results.velocity = struct('x', time, 'y', v);
        appData.Results.acceleration = struct('x', time, 'y', a);
        
        % Create a table for the results
        tension_val =tension
 % Create a table for the results
        resultTable = table(theta', bending_moment', shear_force', tension', 'VariableNames', {'Lay_angle_rad', 'Bending_Moment_Nm', 'Shear_Force_N', 'Tension_N'});
        
        % Format results into a string
     res=formattedDisplayText(resultTable,'NumericFormat','longEng',...
'SuppressMarkup',true,'LineSpacing','compact');
        %resultStr = evalc('disp(resultTable)');
        
        % Update text area with formatted results
        set(textArea, 'String', res);
        
        
        % Save appData
        save(fullfile(appData.ProjectFolder, 'appData.mat'), 'appData');
    end

    % Function to prompt for modal analysis configuration
    function modalAnalysisPrompt(src, event)
        prompt = {'Enter configuration (e.g., Lazy Wave, Steep Wave, etc.):', 'Enter number of modes:'};
        title = 'Modal Analysis Configuration';
        dims = [1 50];
        definput = {'Lazy Wave', '3'};
        answer = inputdlg(prompt, title, dims, definput);
        
        if ~isempty(answer)
            configuration = answer{1};
            num_modes = str2double(answer{2});
            visualize_riser_configuration(configuration, num_modes);
        end
    end

    % Function to prompt for dynamic analysis configuration
    function dynamicAnalysisPrompt(type)
        prompt = {'Select riser section (top, middle, bottom):'};
        title = 'Select Riser Section';
        dims = [1 50];
        definput = {'middle'};
        answer = inputdlg(prompt, title, dims, definput);
        
        if ~isempty(answer)
            section = lower(answer{1});
            n_elements = appData.Parameters.N;
            
            % Determine the index based on the selected section
            switch section
                case 'top'
                    idx = 1;
                case 'middle'
                    idx = ceil(n_elements / 2);
                case 'bottom'
                    idx = (2*n_elements)+2;
                otherwise
                    errordlg('Invalid section selected.', 'Error');
                    return;
            end
            
            figure;
            switch type
                case 'displacement'
                    plot(appData.Results.displacement.x, appData.Results.displacement.y(idx, :),'LineWidth',1.5);
                    xlabel('time [s]')
                    ylabel('Displacement [m]')
                case 'velocity'
                    plot(appData.Results.velocity.x, appData.Results.velocity.y(idx, :),'LineWidth',1.5);
                    xlabel('time [s]')
                    ylabel('Velocity [m/s]')
                case 'acceleration'
                    plot(appData.Results.acceleration.x, appData.Results.acceleration.y(idx, :),'LineWidth',1.5);
                    xlabel('time [s]')
                    ylabel('Acceleration [m/s^2]')
                otherwise
                    errordlg('Invalid type selected.', 'Error');
                    return;
            end
        end
    end

    % Function to visualize riser configuration and modal shapes
    function visualize_riser_configuration(configuration, num_modes)
        % Main script execution
        stl_file = 'shipModel.stl';
    
        % Select configuration
        switch configuration
            case 'Free Hanging Catenary'
                [x, y, z] = spline_fit_3D(xf, yf, zf);
            case {'Lazy Wave', 'Steep Wave', 'Lazy-S', 'Steep-S', 'Pliant Wave'}
                [x, y, z] = spline_fit_3D(x1, y1, z1);
            otherwise
                error('Invalid configuration');
        end
    
        % Plot 3D riser with vessel
        figure;
        draw_vessel_and_seabed_3d(stl_file, x(1), y(1), z(1));
        hold on;
        tubeplot([x; y; z], 0.1, 20); % Adjusted radius for better visibility
        title(configuration);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        view(3);
        hold off;
        
        % Perform modal analysis and visualize modal shapes
        curve = [x; y; z];
        [M, K] = euler_bernoulli_beam(curve);
        visualize_modal_shapes(curve, M, K, num_modes, appData.Results.modalAnalysis.frequencies);
    end

    % Spline function
    function [x_spline, y_spline, z_spline] = spline_fit_3D(x, y, z)
        % Ensure the input vectors are column vectors
        x = x(:);
        y = y(:);
        z = z(:);
    
        % Parameter t
        t = linspace(0, 1, length(x))';
    
        % Perform spline interpolation
        t_fit = linspace(0, 1, 100);
        x_spline = spline(t, x, t_fit);
        y_spline = spline(t, y, t_fit);
        z_spline = spline(t, z, t_fit);
    end

    % Define the tubeplot function with lighting and material properties
    function h = tubeplot(curve, radius, N)
        if nargin < 3
            N = 8;
        end
        if nargin < 2
            radius = 1;
        end
    
        npoints = size(curve, 2);
        X = zeros(N, npoints);
        Y = zeros(N, npoints);
        Z = zeros(N, npoints);
    
        for i = 1:npoints
            t = curve(:, i);
    
            if i == 1
                prev = curve(:, i);
            else
                prev = curve(:, i-1);
            end
            if i == npoints
                next = curve(:, i);
            else
                next = curve(:, i+1);
            end
    
            tangent = next - prev;
            tangent = tangent / norm(tangent);
    
            if i == 1
                normal = [0 0 1]';
            else
                normal = prev - t;
                normal = normal - dot(normal, tangent) * tangent;
                normal = normal / norm(normal);
            end
    
            binormal = cross(tangent, normal);
    
            theta = linspace(0, 2*pi, N);
            circle = radius * [cos(theta); sin(theta)];
    
            X(:, i) = t(1) + normal(1) * circle(1, :) + binormal(1) * circle(2, :);
            Y(:, i) = t(2) + normal(2) * circle(1, :) + binormal(2) * circle(2, :);
            Z(:, i) = t(3) + normal(3) * circle(1, :) + binormal(3) * circle(2, :);
        end
    
        h = surf(X, Y, Z);
        axis off;
        shading interp;
        lighting gouraud;
        material shiny;
        camlight('headlight');
    end

    % Define the vessel and seabed for 3D plot with STL
    function draw_vessel_and_seabed_3d(stl_file, vessel_x, vessel_y, vessel_z)
        hold on;
        vessel = stlread(stl_file);
        vessel_scale_factor = 0.2; % Scale down the vessel for better visibility
        vessel_points_scaled = vessel.Points * vessel_scale_factor;
        trisurf(vessel.ConnectivityList, vessel_points_scaled(:,1) + vessel_x, vessel_points_scaled(:,2) + vessel_y, vessel_points_scaled(:,3) + vessel_z, 'FaceColor', 'blue', 'EdgeColor', 'none');
        lighting gouraud;
        camlight('headlight');
    end

    % Function to perform modal analysis using Euler-Bernoulli beam theory
    function [M, K] = euler_bernoulli_beam(curve)
        n = size(curve, 2);
        M = zeros(3*n, 3*n);
        K = zeros(3*n, 3*n);
        L = diff(curve, 1, 2); % Segment lengths
    
        for i = 1:(n-1)
            len = norm(L(:, i));
            m = 1; % Mass per unit length (example value)
            E = 1e9; % Young's modulus (example value)
            I = 1; % Moment of inertia (example value)
    
            % Mass matrix for beam element
            Me = m * len / 420 * [156 22*len 54 -13*len; 22*len 4*len^2 13*len -3*len^2; 54 13*len 156 -22*len; -13*len -3*len^2 -22*len 4*len^2];
            % Stiffness matrix for beam element
            Ke = E * I / len^3 * [12 6*len -12 6*len; 6*len 4*len^2 -6*len 2*len^2; -12 -6*len 12 -6*len; 6*len 2*len^2 -6*len 4*len^2];
    
            indices = 3*i-2:3*i+1;
            M(indices, indices) = M(indices, indices) + Me;
            K(indices, indices) = K(indices, indices) + Ke;
        end
    
        % Normalize matrices
        M = M / max(M(:));
        K = K / max(K(:));
    end

    % Function to visualize 3D deformations of the modal shapes
    function visualize_modal_shapes(curve, M, K, num_modes, frequencies)
        [V, D] = eig(K, M);
        [~, idx] = sort(diag(D));
        V = V(:, idx);
    
        npoints = size(curve, 2);
        for i = 1:num_modes
            subplot(ceil(num_modes/3), 3, i);
            deformed_curve = curve + 0.1 * reshape(V(1:3*npoints, i), 3, npoints);
            tubeplot(deformed_curve, 2, 20); % Adjusted radius for better visibility
            title(['Mode ', num2str(i), ' (f = ', num2str(frequencies(i)), ' Hz)']);
            axis off;
            %grid on;
            c = colorbar;
            c.Label.String = 'Deformations(m)';
            xlim([-40 40]);
            ylim([-30 30]);
            zlim([-30 30]);
            view(3);
        end
    end

    function visualize_vessel_with_riser(configuration)
      %  global xf yf zf x1 y1 z1;
        
        % Main script execution
        stl_file = 'shipModel.stl';
    
        % Select configuration
        switch configuration
            case 'Free Hanging Catenary'
                [x, y, z] = spline_fit_3D(xf, yf, zf);
            case {'Lazy Wave', 'Steep Wave', 'Lazy-S', 'Steep-S', 'Pliant Wave'}
                [x, y, z] = spline_fit_3D(x1, y1, z1);
            otherwise
                error('Invalid configuration');
        end
    
        % Plot 3D riser with vessel
        draw_vessel_and_seabed_3d(stl_file, x(1), y(1), z(1));
        hold on;
        tubeplot([x; y; z], 0.1, 20); % Adjusted radius for better visibility
        title(configuration);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        view(3);
        hold off;

    end
end
