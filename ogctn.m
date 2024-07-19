classdef ice < matlab.apps.AppBase

	% Properties that correspond to app components
	properties (Access = public)
		UIFigure                    matlab.ui.Figure
		TabGroup2                   matlab.ui.container.TabGroup
		Tab1                        matlab.ui.container.Tab
		UIAxes5                     matlab.ui.control.UIAxes
		UIAxes4                     matlab.ui.control.UIAxes
		UIAxes3                     matlab.ui.control.UIAxes
		UIAxes2                     matlab.ui.control.UIAxes
		Tab2                        matlab.ui.container.Tab
		UIAxes11_2                  matlab.ui.control.UIAxes
		UIAxes10_2                  matlab.ui.control.UIAxes
		UIAxes11                    matlab.ui.control.UIAxes
		UIAxes10                    matlab.ui.control.UIAxes
		UIAxes9                     matlab.ui.control.UIAxes
		UIAxes8                     matlab.ui.control.UIAxes
		UIAxes7                     matlab.ui.control.UIAxes
		UIAxes6                     matlab.ui.control.UIAxes
		Tab3                        matlab.ui.container.Tab
		UIAxes14                    matlab.ui.control.UIAxes
		UIAxes13                    matlab.ui.control.UIAxes
		UIAxes12                    matlab.ui.control.UIAxes
		StatusLamp                  matlab.ui.control.Lamp
		RUNButton                   matlab.ui.control.Button
		TabGroup                    matlab.ui.container.TabGroup
		InputsTab                   matlab.ui.container.Tab
		Panel_2                     matlab.ui.container.Panel
		EVOEditField                matlab.ui.control.NumericEditField
		EVOEditFieldLabel           matlab.ui.control.Label
		IVCEditField                matlab.ui.control.NumericEditField
		IVCEditFieldLabel           matlab.ui.control.Label
		CombStartEditField          matlab.ui.control.NumericEditField
		CombStartEditFieldLabel     matlab.ui.control.Label
		LoadEditField               matlab.ui.control.NumericEditField
		LoadEditFieldLabel          matlab.ui.control.Label
		RPMEditField                matlab.ui.control.NumericEditField
		RPMEditFieldLabel           matlab.ui.control.Label
		n_combEditField             matlab.ui.control.NumericEditField
		n_combEditFieldLabel        matlab.ui.control.Label
		n_mechEditField             matlab.ui.control.NumericEditField
		n_mechEditFieldLabel        matlab.ui.control.Label
		T_wallEditField             matlab.ui.control.NumericEditField
		T_wallEditFieldLabel        matlab.ui.control.Label
		CombDurationEditField       matlab.ui.control.NumericEditField
		CombDurationEditFieldLabel  matlab.ui.control.Label
		P_aEditField                matlab.ui.control.NumericEditField
		P_aEditFieldLabel           matlab.ui.control.Label
		T_aEditField                matlab.ui.control.NumericEditField
		T_aEditFieldLabel           matlab.ui.control.Label
		Panel_Fuel                  matlab.ui.container.Panel
		ADDFUELButton               matlab.ui.control.Button
		AF_ratio                    matlab.ui.control.NumericEditField
		AFLabel                     matlab.ui.control.Label
		LHV                         matlab.ui.control.NumericEditField
		LHVLabel                    matlab.ui.control.Label
		FuelDropDown                matlab.ui.control.DropDown
		FuelDropDownLabel           matlab.ui.control.Label
		EngineGeometryButtonGroup   matlab.ui.container.ButtonGroup
		Option2Button               matlab.ui.control.RadioButton
		Option1Button               matlab.ui.control.RadioButton
		Panel                       matlab.ui.container.Panel
		Image                       matlab.ui.control.Image
		V_cylEditField              matlab.ui.control.NumericEditField
		V_cylEditFieldLabel         matlab.ui.control.Label
		ConnrodlEditField           matlab.ui.control.NumericEditField
		ConnectingrodlLabel         matlab.ui.control.Label
		CompressionRatioEditField   matlab.ui.control.NumericEditField
		CompressionRatioLabel       matlab.ui.control.Label
		BoreStrokeEditField         matlab.ui.control.NumericEditField
		BoreStrokeEditFieldLabel    matlab.ui.control.Label
		V_dEditField                matlab.ui.control.NumericEditField
		V_dEditFieldLabel           matlab.ui.control.Label
		V_cEditField                matlab.ui.control.NumericEditField
		V_cEditFieldLabel           matlab.ui.control.Label
		V_totalEditField            matlab.ui.control.NumericEditField
		V_totalEditFieldLabel       matlab.ui.control.Label
		NumberofCylindersEditField  matlab.ui.control.NumericEditField
		NumberofCylindersEditFieldLabel  matlab.ui.control.Label
		StrokeSEditField            matlab.ui.control.NumericEditField
		StrokeSLabel                matlab.ui.control.Label
		BoreBEditField              matlab.ui.control.NumericEditField
		BoreBLabel                  matlab.ui.control.Label
	end


	properties (Access = private)
%		Oguzalp Cetin		% 150418031
		B;					% Bore [cm]
		S;					% Stroke [cm]
		BS_ratio;			% Bore/Stroke
		l;					% Connecting rod [cm]
		a;					% Crank radius [cm]
		cr;					% Compression ratio
		ca;					% Crank angle vector [degrees]
		V_total;			% Total volume [cc]
		NumberofCylinders;	% Number of cylinders
		V_cyl;				% Single cylinder total volume [cc]
		V_c;				% Single cylinder clearance volume [cc]
		V_d;				% Single cylinder displacement volume [cc]
		V_i;				% Instantaneous cylidner volume, vector [cc]
		dV_i;				% In cylinder volume rate, vector
		s_i;				% Crank center to end of the con. rod, vector [cm]
		Up_i;				% Instantaneous piston speed, vector [cm/rad]
		U_m;				% Mean piston speed [cm/rad]
		A_i;				% Surface area of a cylinder, vector [cm^2]
		
		Q_LHV;				% Lower heating value of fuel [J/kg]
		gamma;	rho;		% Specific heat ratio (cp/cv), air density [kg/m^3]
		AF;		m_f;		% Air fuel ratio (mass), fuel mass [kg]
		Q_f;	dQ_f;		% Heat input, fuel [kJ/kg], Rate of heat release
		R_air;				% Specific gas constant of air: 287.05 [J/(kg*K)]
		mu;					% The cylinder gas viscosity [kg/(m*s)]
		xb;		d_Xb		% Burned mass fraction, single Wiebe approach
		dQ_w;				% Gas-wall heat transfer rate
		dQ;					% Heat transfer rate
		Q;					% Total heat transter per cycle
		RPM;				% Engine speed [rpm]
		Load;				% Engine load - affects inlet pressure, range: [0 - 1]
		IVC;				% Intake valve close time, ca [deg]
		EVO;				% Exhaust valve open time, ca [deg]
		cd;		soc;		% Combustion duration, start of combustion [degrees]
		T_w;				% Cylinder wall temperature [K]
		T_a;	P_a;		% Ambient temperature [K], pressure [kPa]
		n_mech;	n_comb;		% Mechanical efficiency, Combustion efficiency
		P;		T;			% In cylinder pressure [kPa] and temperature [K] vectors
		dP;		dT;			% 
		W;		Wdot;		% Net Indicated Work [N*m], Power [J]
		Torque;
		panel_opt;			% UI variable
		in_progress;		% UI variable

	end
	
	methods (Access = private)

		function	[] = auto_plot(~, ax, x, y, x_lim, tick_dx)
			if x_lim(1) > x_lim(2)
				tmp = x_lim(2);	x_lim(2) = x_lim(1);	x_lim(1) = tmp;
			end
			ax.XLim = x_lim;
			ax.XTick = x_lim(1): tick_dx : x_lim(2);
			r = (max(y) - min(y)) / 2;
			if r ~= 0
				ax.YLim = [min(y) - 0.15 * r, max(y) + 0.15*r];
			end
			clear r; 
			plot(ax, x, y);
			grid(ax, "on");
		end

		function	[] = reset_figures(app)
			ax = [app.UIAxes2, app.UIAxes3, app.UIAxes4, app.UIAxes5, ...
				app.UIAxes6, app.UIAxes7, app.UIAxes8, app.UIAxes9, ...
				app.UIAxes10, app.UIAxes10_2, app.UIAxes11_2, ...
				app.UIAxes11, app.UIAxes12, app.UIAxes13, app.UIAxes14];
			for i = 1:length(ax)
				cla(ax(i));
			end
		end
		
		function	[] = plot_all(app)
			auto_plot(app, app.UIAxes2, app.ca, 1e6*app.V_i, [0 720], 45);
			auto_plot(app, app.UIAxes3, app.ca, 1e4*app.A_i, [0 720], 45);
			auto_plot(app, app.UIAxes4, app.ca, 1e6*app.dV_i, [0 720], 45);
			auto_plot(app, app.UIAxes5, app.ca, 1e2*app.Up_i, [0 720], 45);
			comb = [app.soc, app.soc+app.cd];
			c_t = app.soc:app.soc+app.cd;
			auto_plot(app, app.UIAxes6, c_t, app.gamma(c_t), comb, 20);
			auto_plot(app, app.UIAxes7, c_t, app.mu(c_t), comb, 20);
			auto_plot(app, app.UIAxes8, app.ca, app.xb, comb, 10);
			auto_plot(app, app.UIAxes9, app.ca, app.d_Xb, comb, 10);
			auto_plot(app, app.UIAxes10, app.ca, 1e-3*app.dQ_f, [app.IVC app.EVO], 45);
			auto_plot(app, app.UIAxes10_2, app.ca, 1e-3*app.dQ_w, [app.IVC app.EVO], 45);
			auto_plot(app, app.UIAxes11, app.ca, 1e-3*app.dQ, [0 720], 45);
			plot(app.UIAxes11_2,c_t, 1e-3*app.Q(c_t));hold(app.UIAxes11_2, "on");plot(app.UIAxes11_2, c_t, zeros(1, length(c_t))+1e-3*app.Q_f);
			app.UIAxes11_2.XLim = comb; app.UIAxes11_2.XTick = app.soc:15:app.soc+app.cd;
			auto_plot(app, app.UIAxes12, app.ca, app.T, [0 720], 45);
			auto_plot(app, app.UIAxes13, app.ca, 1e-5*app.P, [0 720], 45);
			semilogy(app.UIAxes14, app.V_i/app.V_c, 1e-5*app.P);
			app.UIAxes14.XLim = [0 app.cr+1]; app.UIAxes14.XTick = 0:app.cr+1;
			app.UIAxes14.YLim = 1e-5*[0.85*min(app.P) 1.15*max(app.P)];
		end

		function	[] = update_edit_fields(app)
			if (app.panel_opt == 1)
				app.BoreStrokeEditField.Value = app.BS_ratio;
				app.V_totalEditField.Value = app.V_total*1e6;
			elseif (app.panel_opt == 2)
				app.StrokeSEditField.Value = app.S * 10^2;
				app.BoreBEditField.Value = app.B * 10^2;
			end
			app.V_cEditField.Value = app.V_c*1e6;
			app.V_dEditField.Value = app.V_d*1e6;
			app.V_cylEditField.Value = app.V_cyl*1e6;
		end

		function	[] = init_app_data(app)
			app.ca = 1 : 720;	% 4 stroke engine
			app.cr = app.CompressionRatioEditField.Value;
			app.NumberofCylinders = app.NumberofCylindersEditField.Value;
			if (app.panel_opt == 1)
				app.B = app.BoreBEditField.Value * 1e-3;
				app.S = app.StrokeSEditField.Value * 1e-3;	% 2*crank radius
				app.BS_ratio = app.B / app.S;
				app.V_d = (pi * app.B^2 / 4) * app.S; % Displacement volume [m^3]
				app.V_c = app.V_d / (app.cr - 1);
				app.V_total = app.NumberofCylinders * (app.V_c + app.V_d);
			elseif (app.panel_opt == 2)
				app.BS_ratio = app.BoreStrokeEditField.Value;
				app.V_total = app.V_totalEditField.Value * 1e-6; % [m^3]
				app.V_c = (app.V_total / app.NumberofCylinders) / app.cr; % [m^3]
				app.V_d = (app.V_total / app.NumberofCylinders) - app.V_c; % [m^3]
				app.B = (4 * app.BS_ratio * app.V_d / pi)^(1/3);
				app.S = app.B / app.BS_ratio;
			end
			app.V_cyl = app.V_total / app.NumberofCylinders;
			update_edit_fields(app);
			app.l = app.ConnrodlEditField.Value * 1e-3;
			app.a = app.S / 2;
			app.RPM = app.RPMEditField.Value;
			app.Load = app.LoadEditField.Value / 100;
			app.IVC = app.IVCEditField.Value;		% IVO set to 1 by default
			app.EVO = app.EVOEditField.Value;		% EVC set to 720 by default
			app.cd = app.CombDurationEditField.Value;
			app.soc = app.CombStartEditField.Value;
			app.R_air = 287.05;
			app.T_a = app.T_aEditField.Value;
			app.P_a = app.P_aEditField.Value * 1e3;
			app.T_w = app.T_wallEditField.Value;
			app.n_comb = app.n_combEditField.Value / 100;
			app.xb = zeros(1, length(app.ca));		app.d_Xb = zeros(1, length(app.ca));
			app.dQ_f = zeros(1, length(app.ca));	app.dQ_w = zeros(1, length(app.ca));
			app.Q = zeros(1, length(app.ca));		app.dQ = zeros(1, length(app.ca));
			app.P = zeros(1, length(app.ca));		app.dP = zeros(1, length(app.ca));
			app.T = zeros(1, length(app.ca));		app.dT = zeros(1, length(app.ca));
			app.rho = zeros(1, length(app.ca));		app.gamma = zeros(1, length(app.ca));
			app.mu = zeros(1, length(app.ca));
			app.W = zeros(1, length(app.ca));	app.Wdot = zeros(1, length(app.ca));
			app.Torque = zeros(1, length(app.ca));
		end

		function	[] = calculate_kinematics(app)
			app.s_i = app.a * cosd(app.ca) + sqrt(app.l^2 - (app.a^2 * sind(app.ca).^2));
			app.Up_i = -app.a * sind(app.ca) + (-app.a^2 * sind(app.ca) .* cosd(app.ca) ./ ...
						(sqrt(app.l^2 - app.a^2 * sind(app.ca).^2)));
			app.V_i = app.V_c + (pi * app.B^2 / 4) * (app.l + app.a - app.s_i);
			app.dV_i = diff([app.V_i app.V_i(1)]);
			app.A_i = 2 * (pi * app.B^2 / 4) + (pi * app.B * (app.l + app.a - app.s_i));
		end

		function	[imep, bmep, fmep, torque, power] = performance_calculations(app, rpm)
			% https://www.philadelphia.edu.jo/academics/nbadarneh/uploads/ICE%20CH4%202018-2019.pdf
			n = 2;	% number of revolutions per event: 2 for four stroke engine
			solver(app, rpm);
			for i = 2 : length(app.ca)
				app.W(i) = app.W(i - 1) + (app.P(i) * app.dV_i(i));	% single cylinder work on expansion [N*m]
				app.Wdot(i) = (app.NumberofCylinders * app.W(i) * (rpm / 60) / n) / 1000;	% Power, [kW]
				app.Torque(i) = (app.Wdot(i) * 1000)/(2 * pi * (rpm / 60));	% Torque, [N*m]
			end
			power = app.Wdot(end);
			imep = app.Wdot(end) * n * 1000 / (app.V_d * 1000 * (rpm / 60));
			% According to Blair, fmep [Pa]
			if app.V_d*1e4 > 500 % in cc
				fmep = 100000 + 350 * app.S * rpm;
			else
				fmep = 100000 + 100 * (500 - app.V_d*1e4) + 350 * app.S * rpm;
			end
			bmep = imep - fmep;
			torque = app.Torque(end);
			app.n_mechEditField.Value = (bmep / imep) * 100;
		end

		function [] = anands_heat_release(app, i, rpm)
			% https://doi.org/10.1016/j.enconman.2007.01.031
			% According to Heywood [2], the cylinder gas viscosity
			app.mu(i) = 7.457e-6 + (4.1547e-8) * app.T(i) - (7.4793e-12) * app.T(i)^2;		% [kg / (m*s)]
			% Below viscosity formula works better for T>2000
				% https://www.omnicalculator.com/physics/kinematic-viscosity-of-air
			% app.mu(i) = (1.458e-6 * app.T(i) ^ (3/2)) / (app.T(i) + 110.4);		% [kg / (m*s)]
			% https://doi.org/10.1016/j.ijhydene.2017.07.016
			app.rho(i) = app.P(i) / (app.R_air * app.T(i));
			Reynolds = app.rho(i) * app.U_m * app.B / app.mu(i);
			Nusselt = 0.49 * Reynolds^0.7;
			k_gas = (6.1944*1e-3) + (7.3814*1e-5)*app.T(i) - (1.2491*1e-8)*app.T(i)^2; % [W / (m*K)]
			h_c = k_gas * Nusselt / app.B; % [W/(K*m^2)]
			h_r = (4.25*1e-9) * ((app.T(i)^4 - app.T_w^4) / (app.T(i) - app.T_w));	% [W/(K*m^2)]
			app.dQ_w(i) = (h_c + h_r) * app.A_i(i) * (app.T(i) - app.T_w) * (60 / (360 * rpm));
			% Overall heat transfer calculations
			app.Q_f = app.n_comb * app.m_f * app.Q_LHV;
			app.dQ_f(i) = app.Q_f * app.d_Xb(i);
			app.dQ(i) = app.dQ_f(i) - app.dQ_w(i);
			app.Q(i+1) = app.Q(i) + app.dQ(i);
		end

		function	[g] = get_gamma(~, T)
		% Specific heat ratio depending on temperature, https://doi.org/10.1016/j.enconman.2004.12.009
			if T <= 1233.5
				% Function of Brunt et al
				g = 1.393 - (1.6 * 10^-4.43)*T - (2.94 * 10^-8) * T^2 + (1.3814 * 10^-11) * T;
			else
				% Function of Egnell
				g = 1.393 - 0.288 * exp(-1100 / T);
			end
		end

		function	[] = set_fractions(app, a, m)
			% Single wiebe function to estimate mass burned fraction rate
			for i = app.soc : (app.soc + app.cd + 1)
				tmp = (i - app.soc) / app.cd;
				app.xb(i) = 1 - exp(-a * (tmp ^ (m+1)));
				app.d_Xb(i) = (a * (m+1) / app.cd) * (tmp ^ m) * exp(-a * tmp^(m+1));	
			end
		end

		function	[] = solver(app, rpm)
			app.U_m = 2 * (rpm / 60) * app.S;
			set_fractions(app, 5, 3);
			app.P(1:app.IVC) = app.P_a * app.Load;
			app.T(1:app.IVC) = app.T_a;
			app.gamma(1 : app.IVC) = get_gamma(app, app.T(app.IVC));
			app.m_f = (app.V_i(app.IVC) * app.P(app.IVC) / (app.R_air * app.T(app.IVC))) / app.AF;
			for i = app.IVC : app.EVO
				app.gamma(i) = get_gamma(app, app.T(i));
				anands_heat_release(app, i, rpm);	% d_Qw, d_Qf and d_Q estimation
				app.dT(i) = app.T(i) * (app.gamma(i) - 1) * ((app.dQ(i) / (app.P(i) * app.V_i(i))) - (app.dV_i(i) / app.V_i(i))); 
				app.T(i+1) = app.T(i) + app.dT(i);
				app.dP(i) = (-app.dV_i(i) * app.P(i) / app.V_i(i)) + (app.P(i) * app.dT(i) / app.T(i));
				app.P(i+1) = app.P(i) + app.dP(i);
			end
			app.P(app.EVO : end) = app.P_a;
			app.T(app.EVO : end) = app.T_a;
			app.gamma(app.EVO : end) = get_gamma(app, app.T(app.EVO));
		end

	end

	% Callbacks that handle component events
	methods (Access = private)

		% Code that executes after component creation
		function startupFcn(app)
			app.in_progress = 0;
			app.panel_opt = 1;
			app.V_totalEditField.Enable = "off";
			app.V_cEditField.Enable = "off";
			app.V_dEditField.Enable = "off";
			app.V_cylEditField.Enable = "off";
			app.BoreStrokeEditField.Enable = "off";
			app.Panel_Fuel.Visible = "off";
			app.Q_LHV = 43.4 * 1e6;			% Default, gasoline
			app.AF = 14.7;					% Default, gasoline
		end

		% Button pushed function: RUNButton
		function RUNButtonPushed(app, event)
			if (~app.in_progress)
				app.in_progress = 1; app.StatusLamp.Color = "r"; drawnow;
				reset_figures(app);
				init_app_data(app);
				calculate_kinematics(app);
				rpm = app.RPM;
				solver(app, rpm);
				plot_all(app);
				app.StatusLamp.Color = "g"; drawnow; app.in_progress = 0;
			end
		end

		% Selection changed function: EngineGeometryButtonGroup
		function EngineGeometryButtonGroupSelectionChanged(app, event)
			if (app.EngineGeometryButtonGroup.SelectedObject == app.Option1Button)
				app.panel_opt = 1;
				app.BoreBEditField.Enable = "on";
				app.StrokeSEditField.Enable = "on";
				app.V_totalEditField.Enable = "off";
				app.BoreStrokeEditField.Enable = "off";
			elseif (app.EngineGeometryButtonGroup.SelectedObject == app.Option2Button)
				app.panel_opt = 2;
				app.BoreBEditField.Enable = "off";
				app.StrokeSEditField.Enable = "off";
				app.V_totalEditField.Enable = "on";
				app.BoreStrokeEditField.Enable = "on";
			end
		end

		% Value changed function: FuelDropDown
		function FuelDropDownValueChanged(app, event)
			if (app.FuelDropDown.Value == "Other")
				app.Panel_2.Visible = "off";
				app.Panel_Fuel.Visible = "on";
			else
				app.Panel_2.Visible = "on";
				if (app.FuelDropDown.Value == "Gasoline")
					app.Q_LHV = 43.4e6;
					app.AF = 14.7;
				elseif (app.FuelDropDown.Value == "Methane")
					app.Q_LHV = 50.0e6;
					app.AF = 17.19;
				elseif (app.FuelDropDown.Value == "Propane")
					app.Q_LHV = 46.4e6;
					app.AF = 15.67;
				end
			end
		end

		% Button pushed function: ADDFUELButton
		function ADDFUELButtonPushed(app, event)
			app.Q_LHV = app.LHV.Value * 1e6;
			app.AF = app.AF_ratio.Value;
			app.Panel_Fuel.Visible = "off";
			app.Panel_2.Visible = "on";
		end
	end

	% Component initialization
	methods (Access = private)

		% Create UIFigure and components
		function createComponents(app)

			% Get the file path for locating images
			pathToMLAPP = fileparts(mfilename('fullpath'));

			% Create UIFigure and hide until all components are created
			app.UIFigure = uifigure('Visible', 'off');
			app.UIFigure.Position = [100 100 1852 951];
			app.UIFigure.Name = 'MATLAB App';

			% Create TabGroup
			app.TabGroup = uitabgroup(app.UIFigure);
			app.TabGroup.Position = [20 151 475 780];

			% Create InputsTab
			app.InputsTab = uitab(app.TabGroup);
			app.InputsTab.Title = 'Inputs';

			% Create Panel
			app.Panel = uipanel(app.InputsTab);
			app.Panel.Position = [23 338 443 382];

			% Create BoreBLabel
			app.BoreBLabel = uilabel(app.Panel);
			app.BoreBLabel.HorizontalAlignment = 'center';
			app.BoreBLabel.Position = [62 323 48 22];
			app.BoreBLabel.Text = 'Bore [B]';

			% Create BoreBEditField
			app.BoreBEditField = uieditfield(app.Panel, 'numeric');
			app.BoreBEditField.ValueDisplayFormat = '%11.4g (mm)';
			app.BoreBEditField.HorizontalAlignment = 'center';
			app.BoreBEditField.Position = [119 323 80 22];
			app.BoreBEditField.Value = 101.6;

			% Create StrokeSLabel
			app.StrokeSLabel = uilabel(app.Panel);
			app.StrokeSLabel.HorizontalAlignment = 'center';
			app.StrokeSLabel.Position = [55 289 58 22];
			app.StrokeSLabel.Text = 'Stroke [S]';

			% Create StrokeSEditField
			app.StrokeSEditField = uieditfield(app.Panel, 'numeric');
			app.StrokeSEditField.ValueDisplayFormat = '%11.4g (mm)';
			app.StrokeSEditField.HorizontalAlignment = 'center';
			app.StrokeSEditField.Position = [119 289 80 22];
			app.StrokeSEditField.Value = 88.4;

			% Create NumberofCylindersEditFieldLabel
			app.NumberofCylindersEditFieldLabel = uilabel(app.Panel);
			app.NumberofCylindersEditFieldLabel.HorizontalAlignment = 'center';
			app.NumberofCylindersEditFieldLabel.Position = [45 116 68 30];
			app.NumberofCylindersEditFieldLabel.Text = {'Number'; 'of Cylinders'};

			% Create NumberofCylindersEditField
			app.NumberofCylindersEditField = uieditfield(app.Panel, 'numeric');
			app.NumberofCylindersEditField.HorizontalAlignment = 'center';
			app.NumberofCylindersEditField.Position = [119 120 80 22];
			app.NumberofCylindersEditField.Value = 4;

			% Create V_totalEditFieldLabel
			app.V_totalEditFieldLabel = uilabel(app.Panel);
			app.V_totalEditFieldLabel.HorizontalAlignment = 'center';
			app.V_totalEditFieldLabel.Position = [68 157 42 22];
			app.V_totalEditFieldLabel.Text = 'V_total';

			% Create V_totalEditField
			app.V_totalEditField = uieditfield(app.Panel, 'numeric');
			app.V_totalEditField.ValueDisplayFormat = '%11.4g (cc)';
			app.V_totalEditField.HorizontalAlignment = 'center';
			app.V_totalEditField.Position = [119 157 80 22];
			app.V_totalEditField.Value = 3200;

			% Create V_cEditFieldLabel
			app.V_cEditFieldLabel = uilabel(app.Panel);
			app.V_cEditFieldLabel.HandleVisibility = 'off';
			app.V_cEditFieldLabel.HorizontalAlignment = 'center';
			app.V_cEditFieldLabel.Position = [83 83 26 22];
			app.V_cEditFieldLabel.Text = 'V_c';

			% Create V_cEditField
			app.V_cEditField = uieditfield(app.Panel, 'numeric');
			app.V_cEditField.ValueDisplayFormat = '%11.4f (cc)';
			app.V_cEditField.HorizontalAlignment = 'center';
			app.V_cEditField.Position = [119 83 80 22];

			% Create V_dEditFieldLabel
			app.V_dEditFieldLabel = uilabel(app.Panel);
			app.V_dEditFieldLabel.HorizontalAlignment = 'center';
			app.V_dEditFieldLabel.Position = [84 50 26 22];
			app.V_dEditFieldLabel.Text = 'V_d';

			% Create V_dEditField
			app.V_dEditField = uieditfield(app.Panel, 'numeric');
			app.V_dEditField.ValueDisplayFormat = '%11.4f (cc)';
			app.V_dEditField.HorizontalAlignment = 'center';
			app.V_dEditField.Position = [119 50 80 22];

			% Create BoreStrokeEditFieldLabel
			app.BoreStrokeEditFieldLabel = uilabel(app.Panel);
			app.BoreStrokeEditFieldLabel.HorizontalAlignment = 'center';
			app.BoreStrokeEditFieldLabel.Position = [42 256 68 22];
			app.BoreStrokeEditFieldLabel.Text = 'Bore/Stroke';

			% Create BoreStrokeEditField
			app.BoreStrokeEditField = uieditfield(app.Panel, 'numeric');
			app.BoreStrokeEditField.HorizontalAlignment = 'center';
			app.BoreStrokeEditField.Position = [119 256 80 22];
			app.BoreStrokeEditField.Value = 1.149;

			% Create CompressionRatioLabel
			app.CompressionRatioLabel = uilabel(app.Panel);
			app.CompressionRatioLabel.HorizontalAlignment = 'center';
			app.CompressionRatioLabel.Position = [8 223 105 22];
			app.CompressionRatioLabel.Text = 'Compression Ratio';

			% Create CompressionRatioEditField
			app.CompressionRatioEditField = uieditfield(app.Panel, 'numeric');
			app.CompressionRatioEditField.HorizontalAlignment = 'center';
			app.CompressionRatioEditField.Position = [119 223 80 22];
			app.CompressionRatioEditField.Value = 8;

			% Create ConnectingrodlLabel
			app.ConnectingrodlLabel = uilabel(app.Panel);
			app.ConnectingrodlLabel.HorizontalAlignment = 'center';
			app.ConnectingrodlLabel.Position = [29 190 81 22];
			app.ConnectingrodlLabel.Text = 'Conn. rod [ l ]';

			% Create ConnrodlEditField
			app.ConnrodlEditField = uieditfield(app.Panel, 'numeric');
			app.ConnrodlEditField.ValueDisplayFormat = '%11.4g (mm)';
			app.ConnrodlEditField.HorizontalAlignment = 'center';
			app.ConnrodlEditField.Position = [119 190 80 22];
			app.ConnrodlEditField.Value = 144.9;

			% Create V_cylEditFieldLabel
			app.V_cylEditFieldLabel = uilabel(app.Panel);
			app.V_cylEditFieldLabel.HorizontalAlignment = 'center';
			app.V_cylEditFieldLabel.Position = [77 17 34 22];
			app.V_cylEditFieldLabel.Text = 'V_cyl';

			% Create V_cylEditField
			app.V_cylEditField = uieditfield(app.Panel, 'numeric');
			app.V_cylEditField.ValueDisplayFormat = '%11.4f (cc)';
			app.V_cylEditField.HorizontalAlignment = 'center';
			app.V_cylEditField.Position = [119 17 80 22];

			% Create Image
			app.Image = uiimage(app.Panel);
			app.Image.Position = [220 17 214 328];
			app.Image.ImageSource = fullfile(pathToMLAPP, 'fig_engine.PNG');

			% Create EngineGeometryButtonGroup
			app.EngineGeometryButtonGroup = uibuttongroup(app.InputsTab);
			app.EngineGeometryButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @EngineGeometryButtonGroupSelectionChanged, true);
			app.EngineGeometryButtonGroup.TitlePosition = 'centertop';
			app.EngineGeometryButtonGroup.Title = 'Engine Geometry';
			app.EngineGeometryButtonGroup.Position = [162 694 165 49];

			% Create Option1Button
			app.Option1Button = uiradiobutton(app.EngineGeometryButtonGroup);
			app.Option1Button.Text = 'Option 1';
			app.Option1Button.Position = [11 3 67 22];
			app.Option1Button.Value = true;

			% Create Option2Button
			app.Option2Button = uiradiobutton(app.EngineGeometryButtonGroup);
			app.Option2Button.Text = 'Option 2';
			app.Option2Button.Position = [90 3 67 22];

			% Create FuelDropDownLabel
			app.FuelDropDownLabel = uilabel(app.InputsTab);
			app.FuelDropDownLabel.HorizontalAlignment = 'right';
			app.FuelDropDownLabel.Position = [142 296 28 22];
			app.FuelDropDownLabel.Text = 'Fuel';

			% Create FuelDropDown
			app.FuelDropDown = uidropdown(app.InputsTab);
			app.FuelDropDown.Items = {'Gasoline', 'Methane', 'Propane', 'Other'};
			app.FuelDropDown.ValueChangedFcn = createCallbackFcn(app, @FuelDropDownValueChanged, true);
			app.FuelDropDown.Position = [185 296 100 22];
			app.FuelDropDown.Value = 'Gasoline';

			% Create Panel_Fuel
			app.Panel_Fuel = uipanel(app.InputsTab);
			app.Panel_Fuel.Visible = 'off';
			app.Panel_Fuel.Position = [31 184 435 97];

			% Create LHVLabel
			app.LHVLabel = uilabel(app.Panel_Fuel);
			app.LHVLabel.HorizontalAlignment = 'right';
			app.LHVLabel.Position = [11 55 28 22];
			app.LHVLabel.Text = 'LHV';

			% Create LHV
			app.LHV = uieditfield(app.Panel_Fuel, 'numeric');
			app.LHV.ValueDisplayFormat = '%11.4g [MJ/kg]';
			app.LHV.Position = [54 55 73 22];

			% Create AFLabel
			app.AFLabel = uilabel(app.Panel_Fuel);
			app.AFLabel.HorizontalAlignment = 'right';
			app.AFLabel.Position = [173 55 25 22];
			app.AFLabel.Text = 'AF';

			% Create AF_ratio
			app.AF_ratio = uieditfield(app.Panel_Fuel, 'numeric');
			app.AF_ratio.Position = [213 55 70 22];

			% Create ADDFUELButton
			app.ADDFUELButton = uibutton(app.Panel_Fuel, 'push');
			app.ADDFUELButton.ButtonPushedFcn = createCallbackFcn(app, @ADDFUELButtonPushed, true);
			app.ADDFUELButton.Position = [183 14 100 23];
			app.ADDFUELButton.Text = 'ADD FUEL';

			% Create Panel_2
			app.Panel_2 = uipanel(app.InputsTab);
			app.Panel_2.Position = [32 24 433 257];

			% Create T_aEditFieldLabel
			app.T_aEditFieldLabel = uilabel(app.Panel_2);
			app.T_aEditFieldLabel.HorizontalAlignment = 'right';
			app.T_aEditFieldLabel.Position = [19 215 26 22];
			app.T_aEditFieldLabel.Text = 'T_a';

			% Create T_aEditField
			app.T_aEditField = uieditfield(app.Panel_2, 'numeric');
			app.T_aEditField.ValueDisplayFormat = '%11.4g [K]';
			app.T_aEditField.HorizontalAlignment = 'center';
			app.T_aEditField.Position = [60 215 57 22];
			app.T_aEditField.Value = 298;

			% Create P_aEditFieldLabel
			app.P_aEditFieldLabel = uilabel(app.Panel_2);
			app.P_aEditFieldLabel.HorizontalAlignment = 'right';
			app.P_aEditFieldLabel.Position = [157 215 26 22];
			app.P_aEditFieldLabel.Text = 'P_a';

			% Create P_aEditField
			app.P_aEditField = uieditfield(app.Panel_2, 'numeric');
			app.P_aEditField.ValueDisplayFormat = '%11.4g [kPa]';
			app.P_aEditField.HorizontalAlignment = 'center';
			app.P_aEditField.Position = [198 215 70 22];
			app.P_aEditField.Value = 101.325;

			% Create CombDurationEditFieldLabel
			app.CombDurationEditFieldLabel = uilabel(app.Panel_2);
			app.CombDurationEditFieldLabel.HorizontalAlignment = 'right';
			app.CombDurationEditFieldLabel.Position = [211 30 89 22];
			app.CombDurationEditFieldLabel.Text = 'Comb. Duration';

			% Create CombDurationEditField
			app.CombDurationEditField = uieditfield(app.Panel_2, 'numeric');
			app.CombDurationEditField.ValueDisplayFormat = '%11.4g [deg]';
			app.CombDurationEditField.HorizontalAlignment = 'center';
			app.CombDurationEditField.Position = [305 30 69 22];
			app.CombDurationEditField.Value = 65;

			% Create T_wallEditFieldLabel
			app.T_wallEditFieldLabel = uilabel(app.Panel_2);
			app.T_wallEditFieldLabel.HorizontalAlignment = 'right';
			app.T_wallEditFieldLabel.Position = [294 215 40 22];
			app.T_wallEditFieldLabel.Text = 'T_wall';

			% Create T_wallEditField
			app.T_wallEditField = uieditfield(app.Panel_2, 'numeric');
			app.T_wallEditField.ValueDisplayFormat = '%11.4g [K]';
			app.T_wallEditField.HorizontalAlignment = 'center';
			app.T_wallEditField.Position = [349 215 71 22];
			app.T_wallEditField.Value = 400;

			% Create n_mechEditFieldLabel
			app.n_mechEditFieldLabel = uilabel(app.Panel_2);
			app.n_mechEditFieldLabel.HorizontalAlignment = 'right';
			app.n_mechEditFieldLabel.Position = [25 133 48 22];
			app.n_mechEditFieldLabel.Text = 'n_mech';

			% Create n_mechEditField
			app.n_mechEditField = uieditfield(app.Panel_2, 'numeric');
			app.n_mechEditField.ValueDisplayFormat = '%11.4g [%%]';
			app.n_mechEditField.HorizontalAlignment = 'center';
			app.n_mechEditField.Position = [88 133 70 22];
			app.n_mechEditField.Value = 95;

			% Create n_combEditFieldLabel
			app.n_combEditFieldLabel = uilabel(app.Panel_2);
			app.n_combEditFieldLabel.HorizontalAlignment = 'right';
			app.n_combEditFieldLabel.Position = [211 133 48 22];
			app.n_combEditFieldLabel.Text = 'n_comb';

			% Create n_combEditField
			app.n_combEditField = uieditfield(app.Panel_2, 'numeric');
			app.n_combEditField.ValueDisplayFormat = '%11.4g [%%]';
			app.n_combEditField.HorizontalAlignment = 'center';
			app.n_combEditField.Position = [274 133 100 22];
			app.n_combEditField.Value = 98;

			% Create RPMEditFieldLabel
			app.RPMEditFieldLabel = uilabel(app.Panel_2);
			app.RPMEditFieldLabel.HorizontalAlignment = 'right';
			app.RPMEditFieldLabel.Position = [32 173 32 22];
			app.RPMEditFieldLabel.Text = 'RPM';

			% Create RPMEditField
			app.RPMEditField = uieditfield(app.Panel_2, 'numeric');
			app.RPMEditField.HorizontalAlignment = 'center';
			app.RPMEditField.Position = [79 173 74 22];
			app.RPMEditField.Value = 1500;

			% Create LoadEditFieldLabel
			app.LoadEditFieldLabel = uilabel(app.Panel_2);
			app.LoadEditFieldLabel.HorizontalAlignment = 'right';
			app.LoadEditFieldLabel.Position = [220 173 32 22];
			app.LoadEditFieldLabel.Text = 'Load';

			% Create LoadEditField
			app.LoadEditField = uieditfield(app.Panel_2, 'numeric');
			app.LoadEditField.ValueDisplayFormat = '%11.4g [%%]';
			app.LoadEditField.HorizontalAlignment = 'center';
			app.LoadEditField.Position = [267 173 100 22];
			app.LoadEditField.Value = 90;

			% Create CombStartEditFieldLabel
			app.CombStartEditFieldLabel = uilabel(app.Panel_2);
			app.CombStartEditFieldLabel.HorizontalAlignment = 'right';
			app.CombStartEditFieldLabel.Position = [230 74 69 22];
			app.CombStartEditFieldLabel.Text = 'Comb. Start';

			% Create CombStartEditField
			app.CombStartEditField = uieditfield(app.Panel_2, 'numeric');
			app.CombStartEditField.ValueDisplayFormat = '%11.4g [deg]';
			app.CombStartEditField.Position = [303 74 71 22];
			app.CombStartEditField.Value = 330;

			% Create IVCEditFieldLabel
			app.IVCEditFieldLabel = uilabel(app.Panel_2);
			app.IVCEditFieldLabel.HorizontalAlignment = 'right';
			app.IVCEditFieldLabel.Position = [85 74 25 22];
			app.IVCEditFieldLabel.Text = 'IVC';

			% Create IVCEditField
			app.IVCEditField = uieditfield(app.Panel_2, 'numeric');
			app.IVCEditField.ValueDisplayFormat = '%11.4g [deg]';
			app.IVCEditField.Position = [124 74 70 22];
			app.IVCEditField.Value = 180;

			% Create EVOEditFieldLabel
			app.EVOEditFieldLabel = uilabel(app.Panel_2);
			app.EVOEditFieldLabel.HorizontalAlignment = 'right';
			app.EVOEditFieldLabel.Position = [82 30 30 22];
			app.EVOEditFieldLabel.Text = 'EVO';

			% Create EVOEditField
			app.EVOEditField = uieditfield(app.Panel_2, 'numeric');
			app.EVOEditField.ValueDisplayFormat = '%11.4g [deg]';
			app.EVOEditField.Position = [126 30 68 22];
			app.EVOEditField.Value = 540;

			% Create RUNButton
			app.RUNButton = uibutton(app.UIFigure, 'push');
			app.RUNButton.ButtonPushedFcn = createCallbackFcn(app, @RUNButtonPushed, true);
			app.RUNButton.Position = [202 54 154 66];
			app.RUNButton.Text = 'RUN';

			% Create StatusLamp
			app.StatusLamp = uilamp(app.UIFigure);
			app.StatusLamp.Position = [313 69 36 36];

			% Create TabGroup2
			app.TabGroup2 = uitabgroup(app.UIFigure);
			app.TabGroup2.Position = [528 108 1304 823];

			% Create Tab1
			app.Tab1 = uitab(app.TabGroup2);
			app.Tab1.Title = 'Tab1';

			% Create UIAxes2
			app.UIAxes2 = uiaxes(app.Tab1);
			title(app.UIAxes2, 'Instantaneous Cylinder Volume')
			xlabel(app.UIAxes2, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes2, 'Volume [cc]')
			zlabel(app.UIAxes2, 'Z')
			app.UIAxes2.Position = [119 424 453 315];

			% Create UIAxes3
			app.UIAxes3 = uiaxes(app.Tab1);
			title(app.UIAxes3, 'Total Heat Transfer Area')
			xlabel(app.UIAxes3, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes3, 'Area [cm^2]')
			zlabel(app.UIAxes3, 'Z')
			app.UIAxes3.Position = [684 424 453 315];

			% Create UIAxes4
			app.UIAxes4 = uiaxes(app.Tab1);
			title(app.UIAxes4, 'Cylinder Volume Rate')
			xlabel(app.UIAxes4, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes4, 'Cylinder Volume Rate [cc/degree]')
			zlabel(app.UIAxes4, 'Z')
			app.UIAxes4.Position = [119 67 453 315];

			% Create UIAxes5
			app.UIAxes5 = uiaxes(app.Tab1);
			title(app.UIAxes5, 'Instantaneous Piston Velocity')
			xlabel(app.UIAxes5, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes5, 'Velocity [cm/rad]')
			zlabel(app.UIAxes5, 'Z')
			app.UIAxes5.Position = [684 67 453 315];

			% Create Tab2
			app.Tab2 = uitab(app.TabGroup2);
			app.Tab2.Title = 'Tab2';

			% Create UIAxes6
			app.UIAxes6 = uiaxes(app.Tab2);
			title(app.UIAxes6, 'Specific heat ratio, γ(T)')
			xlabel(app.UIAxes6, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes6, 'Specific heat ratio, k(T)')
			zlabel(app.UIAxes6, 'Z')
			app.UIAxes6.Position = [44 567 356 197];

			% Create UIAxes7
			app.UIAxes7 = uiaxes(app.Tab2);
			title(app.UIAxes7, 'Gas viscosity, µ(T)')
			xlabel(app.UIAxes7, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes7, 'Viscosity [kg/(m*s)]')
			zlabel(app.UIAxes7, 'Z')
			app.UIAxes7.Position = [44 332 356 197];

			% Create UIAxes8
			app.UIAxes8 = uiaxes(app.Tab2);
			title(app.UIAxes8, 'MFB (Mass Fraction Burned), x_b')
			xlabel(app.UIAxes8, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes8, 'MFB')
			zlabel(app.UIAxes8, 'Z')
			app.UIAxes8.Position = [481 567 356 197];

			% Create UIAxes9
			app.UIAxes9 = uiaxes(app.Tab2);
			title(app.UIAxes9, 'Rate of MFB, dx_b/dθ')
			xlabel(app.UIAxes9, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes9, 'Rate of MFB [1/degree]')
			zlabel(app.UIAxes9, 'Z')
			app.UIAxes9.Position = [481 332 356 197];

			% Create UIAxes10
			app.UIAxes10 = uiaxes(app.Tab2);
			title(app.UIAxes10, 'dQ_f / dθ')
			xlabel(app.UIAxes10, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes10, 'dQ_f / dθ')
			zlabel(app.UIAxes10, 'Z')
			app.UIAxes10.Position = [46 85 356 197];

			% Create UIAxes11
			app.UIAxes11 = uiaxes(app.Tab2);
			title(app.UIAxes11, 'dQ / dθ : (dQ_w / dθ - dQ_f / dθ)')
			xlabel(app.UIAxes11, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes11, 'dQ / dθ')
			zlabel(app.UIAxes11, 'Z')
			app.UIAxes11.Position = [897 86 356 197];

			% Create UIAxes10_2
			app.UIAxes10_2 = uiaxes(app.Tab2);
			title(app.UIAxes10_2, 'dQ_w / dθ')
			xlabel(app.UIAxes10_2, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes10_2, 'dQ_w/ dθ')
			zlabel(app.UIAxes10_2, 'Z')
			app.UIAxes10_2.Position = [481 85 356 197];

			% Create UIAxes11_2
			app.UIAxes11_2 = uiaxes(app.Tab2);
			title(app.UIAxes11_2, 'Q')
			xlabel(app.UIAxes11_2, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes11_2, 'Q [kJ]')
			zlabel(app.UIAxes11_2, 'Z')
			app.UIAxes11_2.Position = [897 441 356 197];

			% Create Tab3
			app.Tab3 = uitab(app.TabGroup2);
			app.Tab3.Title = 'Tab3';

			% Create UIAxes12
			app.UIAxes12 = uiaxes(app.Tab3);
			title(app.UIAxes12, 'In Cylinder Temperature')
			xlabel(app.UIAxes12, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes12, 'Temperature [K]')
			zlabel(app.UIAxes12, 'Z')
			app.UIAxes12.Position = [102 485 470 278];

			% Create UIAxes13
			app.UIAxes13 = uiaxes(app.Tab3);
			title(app.UIAxes13, 'In Cylinder Pressure')
			xlabel(app.UIAxes13, 'Crank Angle, θ [degree]')
			ylabel(app.UIAxes13, 'Pressure [bar]')
			zlabel(app.UIAxes13, 'Z')
			app.UIAxes13.Position = [775 485 470 278];

			% Create UIAxes14
			app.UIAxes14 = uiaxes(app.Tab3);
			title(app.UIAxes14, 'Normalized Pressure-Volume')
			xlabel(app.UIAxes14, 'V/V_c')
			ylabel(app.UIAxes14, 'Pressure [bar]')
			zlabel(app.UIAxes14, 'Z')
			app.UIAxes14.Position = [264 13 670 412];

			% Show the figure after all components are created
			app.UIFigure.Visible = 'on';
		end
	end

	% App creation and deletion
	methods (Access = public)

		% Construct app
		function app = ice

			% Create UIFigure and components
			createComponents(app)

			% Register the app with App Designer
			registerApp(app, app.UIFigure)

			% Execute the startup function
			runStartupFcn(app, @startupFcn)

			if nargout == 0
				clear app
			end
		end

		% Code that executes before app deletion
		function delete(app)

			% Delete UIFigure when app is deleted
			delete(app.UIFigure)
		end
	end
end
