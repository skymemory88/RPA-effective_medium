function  [cVar, freq_total, chi0, chiq, qvec, dscrt_var] = MF_RPA_Yikai(mion, scanMode, dscrt_var, freq_total, theta, phi, gama, hyp, RPA_mode)
% Current version assumes complete symmetrical equivalence among the four spin moments per unit cell
% mion: Magnetic ion type: 'Er', 'Ho'
% scanMode: 'field' or 'temp' scan
% dscrt_var: discrete variable (temperature or field)
% freq_total frequency range [GHz]
% theta: misalignment angle from c-axis [degree]
% phi: inplane angle in ab-plane between the magnetic field and a/b-axis [degree]
% gama: linewidth of hyperfine levels [meV].
% hyp: Hyperfine isotope proportion (0~1)
% RPA_Mode: RPA option (true/false)

Options.plotting = true; % Decide whether or not to plot the data at thes end
Options.unit = 'GHz'; % Energy unit choice: J, GHz, meV (default)
Options.saving = false; % Options to save the susceptibility tensors
Options.scanMode = scanMode; % 1. Field plot with RPA; 2. Temp plot with RPA; 3. wavevector plot with RPA
Options.nZee = true;
Options.RPA = RPA_mode; % Apply random phase approximation (RPA) correction
Options.Kplot = false; % k-dependent plot for RPA susceptibilities

if Options.RPA == false
    Options.Kplot = false;
end
if Options.Kplot == true
    qz = linspace(0,-1, 101)';
    qx = linspace(0, 1, 101)';
    % qx = ones(length(qz),1);
    qy = zeros(length(qx),1);
    % qz = zeros(length(qx),1);
    qvec = [qx qy qz];
    cVar0 = [0.44]; % selected points of the continuous variable for k-plot
else
    qx = 0.0;
    % qx = [0.01 0.1 0.3 0.6 1]';
    qy = zeros(size(qx,1),1);
    qz = zeros(size(qx,1),1);
    qvec = [qx qy qz];
end

% Declare physical constants as global for consistency
const.hbar = 1.05457E-34; % Reduced Planck constant [J.s]
const.J2meV = 6.24151e+21; % convert Joule to meV
const.Gh2mV = const.hbar * 2*pi * 1e9 * const.J2meV; % convert GHz to meV
const.muB = 9.27401e-24; % Bohr magneton [J/T]
const.muN = 5.05078e-27; % Nuclear magneton [J/T]
const.mu0 = 4*pi*1e-7; % [H/m]
const.dpRng = 100; % dipole summation range (number of unit cell)

% Select the final output unit (default: meV)
switch Options.unit
    case 'GHz'
    ConvUnit = 1/const.Gh2mV;
    case 'J'
    ConvUnit = 1/const.J2meV;
    otherwise
    ConvUnit = 1;
end

for ii = 1:length(dscrt_var) 
    ipt = false;
    counter = 0;
    if hyp > 0
        nZee_path = 'Hz_I=1';
    else
        nZee_path = 'Hz_I=0';
    end
    if ispc
        Options.location = ['C:\Users\skyme\OneDrive - Nexus365\Postdoc\Research projects\Li',mion,...
            'F4 project\Data\Simulations\mean field\eigen-states\', nZee_path, '\'];
    else
        Options.location = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Postdoc/Research projects/',...
            'Li', mion,'F4 project/Data/Simulations/mean field/eigen-states/', nZee_path, '/'];
    end
    while ~ipt
        switch Options.scanMode % 1. Field plot with RPA. 2. wavevector plot with RPA
            case 'field'
                filename = strcat(['Hscan_Li',mion,'F4_'],...
                    sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var(ii), theta, phi, hyp),'.mat');
                file = fullfile(Options.location, filename);
                load(file,'-mat','eee','fff','ttt','vvv','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
                file_part = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_%5$.3f-%6$.3fGHz'...
                    , dscrt_var(ii), theta, phi, gama, min(freq_total), max(freq_total));
                cVar = vecnorm(fff,2,1); % choose the continuous variable to be field
                ttt = repelem(ttt, length(cVar));
                ipt = true;
                fprintf('Calculating for T = %.3f K.\n', dscrt_var(ii));
            case 'temp'
                filename = strcat(['Tscan_Li',mion,'F4_'],...
                    sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var(ii), theta, phi, hyp),'.mat');
                file = fullfile(Options.location,filename);
                load(file,'-mat','eee','ttt','vvv','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
                file_part = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_%4$.2e_%5$.3f-%6$.3fGHz'...
                    , dscrt_var(ii), theta, phi, gama, min(freq_total), max(freq_total));
                cVar = ttt; % choose the continuous variable to be temperature
                ipt = true;
                fprintf('Calculating for B = %.3f T.\n', dscrt_var(ii));
            otherwise
                prompt = sprintf('Unknown Scan Mode! Please select a mode between tempscan and fieldscan.\n');
                answer = input(prompt,'s');
                Options.scanMode = lower(answer);
                counter = counter + 1; % input request counter
                if counter >= 3
                    return % terminate the program after three failed requests
                end
        end
    end
    const.elem = find(ion.prop);    
    const.ELEf = ion.gLande(const.elem) * const.muB * const.J2meV; % Lande factor * Bohr magneton (meV/T)
    const.NUCf = ion.nLande(const.elem) * const.muN * const.J2meV; % [meV/T]
    const.gfac = ion.gLande(const.elem)^2 * (const.muB)^2 * (const.mu0/4/pi) * const.J2meV * 10^30; % (gL.muB)^2.mu0/4pi [meV.Ang^3]

    eee = eee - min(eee,[],2); % Normalize the eigen-energies to the ground state
    save_name = strcat('chi_Li',mion,'F4_',file_part, '.mat');
    
    eigenE = eee; % load eigenenergy [meV]
%     eigenE(:,9:end) = 0; % truncate the Hilbert space
    eigenW = vvv; % load eigenstates
%     eigenW(:,9:end,9:end) = 0; % truncate the Hilbert space
    if Options.RPA == true
        if Options.Kplot == true
            bidx = int16.empty(0,length(cVar0));
            for jj = 1:length(cVar0)
                [~,bidx(jj)] = min(abs( cVar - cVar0(jj) ));
            end
            eigenW = vvv(bidx,:,:); % eigen-functions
            eigenE = eee(bidx,:); % eigen-energies [meV]
            cVar = cVar(bidx); % magnetic field or temperature
        end
        [cVar, freq_total, chi0, ~, ~] = linear_response(ion, eigenE, cVar, freq_total, ttt, eigenW, gama, const, Options);
        [~, ~, ~, chiq, ~] = RPA(qvec, cVar, freq_total, ion, chi0, const); % Electronic susceptibilitie
%         chiq = const.ELEf^2 * chiq .* ConvUnit; % [J/T^2 or GHz/T^2 or meV/T^2]
        chiq = chiq .* ConvUnit; % [J/T^2 or GHz/T^2 or meV/T^2]
    else
        [cVar, freq_total, chi0, ~, ~] = linear_response(ion, eigenE, cVar, freq_total, ttt, eigenW, gama, const, Options);
    end
%     chi0 = const.ELEf^2 * chi0 .* ConvUnit; % [J/T^2 or GHz/T^2 or meV/T^2]
    chi0 = chi0 .* ConvUnit; % [J/T^2 or GHz/T^2 or meV/T^2]

    if Options.saving == true % Save the susceptibilities
        if Options.RPA == true
            savefile1 = fullfile(Options.location,save_name);
            save_vars = {dscrt_var, cVar, freq_total, ion, chi0, gama, qvec, chiq};
            save_file(save_vars, savefile1, Options); % save the data w. RPA corrections
        else
            savefile2 = fullfile(Options.location,save_name);
            save_vars = {dscrt_var, cVar, freq_total, ion, chi0, gama};
            save_file(save_vars, savefile2, Options); % Save the data w/o RPA corrections
        end
    end
    
    if Options.plotting == true % Plot the susceptibilities
        plot_var = {cVar, freq_total, chi0, gama, [0,0,0], dscrt_var(ii)};
        if Options.Kplot == true
            chi0_p = repmat(chi0,1,1,1,1,length(qvec));
%             chi0_p = permute(chi0_p,[1 2 3 5 4]); % Add the additional dimension for the plotting purpose
            plot_var{5} = qvec;
            if Options.RPA == true
                chi_bundle.chi0 = chi0_p;
                chi_bundle.chiq = chiq;
                chi_bundle.diff = chiq - chi0_p;
                plot_var{3} = chi_bundle;
                figs(plot_var, Options, const, "\chi", Options.Kplot);
            else
                plot_var{3} = chi0_p;
                figs(plot_var, Options, const, "\chi_0", Options.Kplot);
            end
        else
            if Options.RPA == true
                chi_bundle.chi0 = chi0;
                chi_bundle.chiq = chiq;
                chi_bundle.diff = chiq - chi0;
                plot_var{3} = chi_bundle;
                figs(plot_var, Options, const, "\chi", Options.Kplot);
            else
                figs(plot_var, Options, const, "\chi_0", Options.Kplot);
            end
        end
    end
end
end

function figs(input_var, Options, const, fig_tit, Qplot)
continu_var = input_var{1};
freq_total = input_var{2};
chi = input_var{3};
gama = input_var{4};
qvec = input_var{5};
dscrt_var = input_var{6};

% custom colormap scale
cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
[linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
[ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
cmap = flip(cmap,1);

% convert frequency axis to appropriate choice
if strcmp(Options.unit, 'J')
    freq_total = freq_total*const.Gh2mV/const.J2meV;
    ylab = "Energy (J)"; % yaxis label
elseif strcmp(Options.unit, 'meV')
    freq_total = freq_total*const.Gh2mV;
    ylab = "Energy (meV)";
else
    ylab = "Frequency (GHz)";
end
 
pos0 = [100 300 600 400]; % initial figure position
pos_inc = [150 0 0 0];
if Qplot == true
    qv = vecnorm(qvec,2,2);
    % qv = qvec(:,1); % use qx
    if isstruct(chi) && all(isfield(chi, {'chi0','chiq','diff'}))
        component_data = {chi.chi0, chi.chiq, chi.diff};
        component_labels = ["\chi_0","\chi_{RPA}","\chi_{RPA}-\chi_0"];
    else
        component_data = {chi};
        component_labels = string(fig_tit);
    end

    axis_labels = ["xx","yy","zz"];
    ncols = numel(component_data);
    base_array = component_data{1};
    n_states = max(1, size(base_array,4));

    for nb = 1:n_states
        if numel(continu_var) >= nb
            cont_val = continu_var(nb);
        elseif ~isempty(continu_var)
            cont_val = continu_var(end);
        else
            cont_val = NaN;
        end

        switch Options.scanMode
            case 'field'
                file_part1 = sprintf('T = %.3f K and ', dscrt_var);
                file_part2 = sprintf('B = %1$.2f T.', cont_val);
            case 'temp'
                file_part1 = sprintf('B = %.3f T and ', dscrt_var);
                file_part2 = sprintf('T = %1$.2f K.', cont_val);
            otherwise
                disp("Unknow scan mode!")
                return
        end

        if isvector(freq_total)
            freq_axis = freq_total(:);
        else
            freq_axis = squeeze(freq_total(min(nb,size(freq_total,1)),:));
        end
        freq_axis = freq_axis(:);
        if isempty(freq_axis)
            freq_axis = freq_total(:);
        end
        len_freq = numel(freq_axis);

        for jj = 1:3
            fig = figure;
            fig_pos = pos0 + (jj-1)*pos_inc;
            fig_pos(3) = max(fig_pos(3), 1050);
            fig_pos(4) = max(fig_pos(4), 600);
            set(fig,'position',fig_pos);

            for cc = 1:ncols
                state_idx = min(nb, size(component_data{cc},4));
                data_slice = squeeze(component_data{cc}(jj,jj,:,state_idx,:));
                if isvector(data_slice)
                    data_slice = reshape(data_slice, len_freq, []);
                end
                if size(data_slice,1) ~= len_freq && size(data_slice,2) == len_freq
                    data_slice = data_slice.';
                end
                if size(data_slice,1) ~= len_freq
                    data_slice = reshape(data_slice, len_freq, []);
                end
                len_q = size(data_slice,2);
                q_axis = qv(1:len_q);
                % Show start/mid/end q-points with explicit vector labels
                num_ticks = min(3, len_q);
                tick_idx = unique(round(linspace(1, len_q, num_ticks)));
                tick_idx = tick_idx(:)';
                q_tick_vals = q_axis(tick_idx);
                q_labels = arrayfun(@(idx) sprintf('[%0.1f %0.1f %0.1f]',...
                    qvec(idx,1), qvec(idx,2), qvec(idx,3)), tick_idx, 'UniformOutput', false);

                real_map = mag2db(abs(real(data_slice)));
                imag_map = mag2db(abs(imag(data_slice)));

                ax_real = subplot(2,ncols,cc);
                hp_real = pcolor(q_axis,freq_axis,real_map);
                set(hp_real,'edgeColor','none');
                set(ax_real,'xdir','reverse');
                set(ax_real,'XTick',q_tick_vals,'XTickLabel',q_labels,'XTickLabelRotation',45);
                colormap(ax_real,cmap);
                caxis(ax_real,'auto');
                colorbar(ax_real);
                title(ax_real,sprintf('Re[%s^{%s}]', char(component_labels(cc)), char(axis_labels(jj))));
                if cc == 1
                    ylabel(ax_real,ylab);
                end

                ax_imag = subplot(2,ncols,cc + ncols);
                hp_imag = pcolor(q_axis,freq_axis,imag_map);
                set(hp_imag,'edgeColor','none');
                set(ax_imag,'xdir','reverse');
                set(ax_imag,'XTick',q_tick_vals,'XTickLabel',q_labels,'XTickLabelRotation',45);
                colormap(ax_imag,cmap);
                caxis(ax_imag,'auto');
                colorbar(ax_imag);
                title(ax_imag,sprintf('Im[%s^{%s}]', char(component_labels(cc)), char(axis_labels(jj))));
                if cc == 1
                    ylabel(ax_imag,ylab);
                end
                xlabel(ax_imag,'Q = [h, k, l]');
            end

            sgtitle(sprintf('%s%s (\gamma = %s meV)', file_part1, file_part2, num2str(gama,'%.2e')), 'Interpreter','tex');
        end
    end


else
    for nq = 1:size(qvec,1)
        switch Options.scanMode
            case 'field'
                xlab = 'Magnetic Field (T)';
                file_part1 = sprintf('T = %.3f K and ', dscrt_var);
                file_part2 = sprintf('Q = [%1$.2f %2$.2f %3$.2f]', qvec(nq,1), qvec(nq,2), qvec(nq,3));
            case 'temp'
                xlab = 'Temperature (K)';
                file_part1 = sprintf('T = %.3f K and ', dscrt_var);
                file_part2 = sprintf('Q = [%1$.2f %2$.2f %3$.2f]', qvec(nq,1), qvec(nq,2), qvec(nq,3));
            otherwise
                disp("Unknow scan mode!")
                return
        end

        if isstruct(chi) && all(isfield(chi, {'chi0','chiq','diff'}))
            component_data = {chi.chi0, chi.chiq, chi.diff};
            component_labels = ["\chi_0","\chi_{RPA}","\chi_{RPA}-\chi_0"];
        else
            component_data = {chi};
            component_labels = string(fig_tit);
        end

        ncols = numel(component_data);
        axis_labels = ["xx","yy","zz"];

        for jj = 1:3
            fig = figure;
            fig_pos = pos0 + (jj-1)*pos_inc;
            fig_pos(3) = max(fig_pos(3), 1050);
            fig_pos(4) = max(fig_pos(4), 600);
            set(fig,'position',fig_pos);

            for cc = 1:ncols
                data_slice = squeeze(component_data{cc}(jj,jj,:,:,nq));
                real_map = mag2db(abs(real(data_slice)));
                imag_map = mag2db(abs(imag(data_slice)));

                ax_real = subplot(2,ncols,cc);
                hp_real = pcolor(continu_var(1,:),freq_total,real_map);
                set(hp_real,'edgeColor','none');
                colormap(ax_real,cmap);
                caxis(ax_real,'auto');
                colorbar(ax_real);
                title(ax_real,sprintf('Re[%s^{%s}]', char(component_labels(cc)), char(axis_labels(jj))));
                if cc == 1
                    ylabel(ax_real,ylab);
                end

                ax_imag = subplot(2,ncols,cc + ncols);
                hp_imag = pcolor(continu_var(1,:),freq_total,imag_map);
                set(hp_imag,'edgeColor','none');
                colormap(ax_imag,cmap);
                caxis(ax_imag,'auto');
                colorbar(ax_imag);
                title(ax_imag,sprintf('Im[%s^{%s}]', char(component_labels(cc)), char(axis_labels(jj))));
                if cc == 1
                    ylabel(ax_imag,ylab);
                end
                xlabel(ax_imag,xlab);
            end

            sgtitle(sprintf('%s%s (\gamma = %s meV)', file_part1, file_part2, num2str(gama,'%.2e')), 'Interpreter','tex');
        end
    end
end
end

function save_file(save_vars, savefile, opt)
switch opt.scanMode
    case 'field'
        temp = save_vars{1};
        fields = save_vars{2};
    case 'temp'
        fields = save_vars{1};
        temp = save_vars{2};
end
freq_total = save_vars{3};
ion = save_vars{4};
chi = save_vars{5};
gama = save_vars{6};
unit = opt.unit;
if length(save_vars) == 8
    qvec = save_vars{7};
    chiq = save_vars{8};
    save(savefile,'temp','fields','freq_total','ion','chi','gama','qvec','chiq','unit','-v7.3');
else
    save(savefile,'temp','fields','freq_total','ion','chi','gama','unit','-v7.3');
end
end

function [cVar, freq_total, chi, JIz_exp, gamma] = linear_response(ion, eigenE, cVar, freq_total,...
    temperatures, eigenW, gma, const, Options)
kB = 8.61733e-2; % Boltzmann constant [meV/K]
gamma = ones(size(eigenE,2))*gma; % spin linewidth [meV]

% Declare susceptibility tensor
chi0 = zeros(3,3,length(freq_total(1,:)),size(cVar,2));

%Initiate ionJ operators
ionJ = ion.J(const.elem);
Jz = diag(ionJ:-1:-ionJ); % Jz = -J, -J+1, ... ,J-1,J
Jp = diag(sqrt((ionJ-((ionJ-1):-1:-ionJ) ).*(ionJ+1+( (ionJ-1):-1:-ionJ) )),1); % electronic spin ladder operator
% Jp = diag(sqrt(ionJ*(ionJ+1) - ((ionJ-1):-1:-ionJ).*(ionJ:-1:(-ionJ+1))),1); % alternative format of spin ladder operator
Jm = Jp'; % electronic spin ladder operator

if ion.hyp(const.elem) > 0 % hyperfine interaction option
    %Initiate I operators
    ionI = ion.I(const.elem);
    Iz = diag(ionI:-1:-ionI); % Iz = -I, -I+1,...,I-1,I
%     Iz = diag([ones(1,4)*0.5 ones(1,4)*(-0.5)]); % projection to Ising space
    IhT.z = kron(eye(2*ionJ+1), Iz); % Expand Hilbert space
    Ip = diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1); % Nuclear spin ladder operator
    % Ip = diag(sqrt(ionI*(ionI+1) - ((ionI-1):-1:-ionI).*(ionI:-1:(-ionI+1))),1); % alternative format of spin ladder operator
    Im = Ip'; % Nuclear spin ladder operator
    Iph = kron(eye(2*ionJ+1), Ip); % Expand to match the dimension of Hilbert space
    Imh = kron(eye(2*ionJ+1), Im);
    IhT.x = (Iph+Imh)/2;
    IhT.y = (Iph-Imh)/2i;
    
    Jph = kron(Jp, eye(2*ionI+1)); % Expand to match the dimension of Hilbert space
    Jmh = kron(Jm, eye(2*ionI+1));
    JhT.x = (Jph+Jmh)/2; % Expand J space to include nuclear degree of freedom
    JhT.y = (Jph-Jmh)/2i;
    JhT.z = kron(Jz, eye(2*ionI+1));
else
    JhT.x = (Jp+Jm)/2;
    JhT.y = (Jp-Jm)/2i;
    JhT.z = Jz;

    IhT.x = 0;
    IhT.y = 0;
    IhT.z = 0;
end
JhT.x = JhT.x + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.x; % Hybridized electronuclear spin operator
JhT.y = JhT.y + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.y;
JhT.z = JhT.z + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.z;

% Single out <Jz+Iz> calculations
JIz_exp = double.empty(size(cVar,2), size(JhT.z,1),0); % Expectation value of J-I pseudo-spin
for ii = 1:size(cVar,2) % calculate susceptibility for all fields
    v = squeeze(eigenW(ii,:,:)); % Obtain the corresponding eigen vectors
    ttz  = v'  * JhT.z * v;
    JIz_exp(ii,:,1) = real(diag(ttz));
end

% For Ho, If exists, load spin-linewidths extracted from experiments for field scan simulation
if ion.prop(2)
    if temperatures(1) == temperatures(end) % for fixed temperature calculations
        if Options.nZee && ismember(temperatures(1), [0.15 0.18 0.22 0.30])
            gma_load = load([Options.location, sprintf('Exp_fit_gamma_%umK.mat', 1000*temperatures(1))],'gma');
            gma_load = flip(gma_load.gma)*const.Gh2mV; % [meV]
            for kk = 1:length(gma_load)
                gamma(kk,kk+1) = gma_load(kk);
                gamma(kk+1,kk) = gma_load(kk);
            end
        elseif ~Options.nZee && ismember(temperatures(1), [0.1 0.13 0.15 0.18 0.24 0.25])
            gma_load = load([Options.location, sprintf('Exp_fit_gamma_%umK.mat', 1000*temperatures(1))],'gma');
            gma_load = flip(gma_load.gma)*const.Gh2mV; % [meV]
            for kk = 1:length(gma_load)
                gamma(kk,kk+1) = gma_load(kk);
                gamma(kk+1,kk) = gma_load(kk);
            end
        end
    end
end


for nf = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
    freq = freq_total (nf);
    omega = freq*const.Gh2mV;   % GHz to meV
%     for jj = 1:size(continu_var,2) % for debugging: calculate susceptibility for all fields
    parfor jj = 1:size(cVar,2) % calculate susceptibility for all fields
        v = squeeze(eigenW(jj,:,:)); % Obtain the corresponding eigen vectors
        en = squeeze(eigenE(jj,:)); % Obtain the corresponding eigen energies [meV]
        if temperatures(jj) ~= 0
            beta = 1/(kB*temperatures(jj)); % [meV^-1]
            Z = sum(exp(-beta*en));
            zn = exp(-beta*en)/Z;
            [n,np] = meshgrid(zn,zn);
            NN = n-np;
        else
            Z = zeros(size(en));
            Z(1) = 1; % Occupy the ground state with unity probability
            [n,np] = meshgrid(Z,Z);
            NN = n-np;
        end
        % % pseudo-statistics of occupation number for debugging/testing
        % NN = ones(size(NN));
        
        [ee,eep] = meshgrid(en,en);
        EE = (eep-ee-omega); % [meV]

        % computer separatly the electronic and nuclear spin susceptibilities 
        deno0 = 1 ./ (EE - 1i*gamma); % [meV^-1]
        chi0(:,:,nf,jj) = chi_Mx(JhT, JhT, v, NN, deno0); % Electornic spin operators
    end
end
chi = chi0; % Electronic susceptibilities
end

function chi0 = chi_Mx(varargin)
% Calculation of matrix elements for the real and imaginary susceptibilities
% opr1: operator 1
% opr2: operator 2
% wav: eigenfunction of the hamiltonian
% NN: population difference between two eigenstates
% deno1: denominator for real part of the susceptibility
% deno2: denominator for imaginary part of the susceptibility
opr1 = varargin{1};
opr2 = varargin{2};
wav = varargin{3};
NN = varargin{4};

ttx1  = wav'  * opr1.x * wav;
tty1  = wav'  * opr1.y * wav;
ttz1  = wav'  * opr1.z * wav;

ttx2  = wav'  * opr2.x * wav;
tty2  = wav'  * opr2.y * wav;
ttz2  = wav'  * opr2.z * wav;

if nargin == 6 % Compute separately the real and imaginary parts
    deno1 = varargin{5};
    deno2 = varargin{6};
    % Calculate susceptibilities along a-axis
    chi_tx  = (ttx1) .* (ttx2.') .* NN .* deno1; % Imaginary part
    chi_t1x = (ttx1) .* (ttx2.') .* NN .* deno2; % Real part
    xx = sum(sum(chi_tx));
    xx1 = sum(sum(chi_t1x));
    
    % Calculate susceptibilities along b-axis
    chi_ty  = (tty1) .* (tty2.') .* NN .* deno1;
    chi_t1y = (tty1) .* (tty2.') .* NN .* deno2;
    yy = sum(sum(chi_ty));
    yy1 = sum(sum(chi_t1y));
    
    % Calculate susceptibilities along c-axis
    chi_tz  = (ttz1) .* (ttz2.') .* NN .* deno1;
    chi_t1z = (ttz1) .* (ttz2.') .* NN .* deno2;
    zz = sum(sum(chi_tz));
    zz1 = sum(sum(chi_t1z));
    
    % Calculate susceptibilities along ab-axis
    chi_txy  = (ttx1) .* (tty2.') .* NN .* deno1;
    chi_t1xy = (ttx1) .* (tty2.') .* NN .* deno2;
    xy = sum(sum(chi_txy));
    xy1 = sum(sum(chi_t1xy));
    
    % Calculate susceptibilities along ac-axis
    chi_txz  = (ttx1) .* (ttz2.') .* NN .* deno1;
    chi_t1xz = (ttx1) .* (ttz2.') .* NN .* deno2;
    xz = sum(sum(chi_txz));
    xz1 = sum(sum(chi_t1xz));
    
    % Calculate susceptibilities along ba-axis
    chi_tyx  = (tty1) .* (ttx2.') .* NN .* deno1;
    chi_t1yx = (tty1) .* (ttx2.') .* NN .* deno2;
    yx = sum(sum(chi_tyx));
    yx1 = sum(sum(chi_t1yx));
    
    % Calculate susceptibilities along bc-axis
    chi_tyz  = (tty1) .* (ttz2.') .* NN .* deno1;
    chi_t1yz = (tty1) .* (ttz2.') .* NN .* deno2;
    yz = sum(sum(chi_tyz));
    yz1 = sum(sum(chi_t1yz));
    
    % Calculate susceptibilities along ca-axis
    chi_tzx  = (ttz1) .* (ttx2.') .* NN .* deno1;
    chi_t1zx = (ttz1) .* (ttx2.') .* NN .* deno2;
    zx = sum(sum(chi_tzx));
    zx1 = sum(sum(chi_t1zx));
    
    % Calculate susceptibilities along cb-axis
    chi_tzy  = (ttz1) .* (tty2.') .* NN .* deno1;
    chi_t1zy = (ttz1) .* (tty2.') .* NN .* deno2;
    zy = sum(sum(chi_tzy));
    zy1 = sum(sum(chi_t1zy));
    
    xr = [xx xy xz
          yx yy yz
          zx zy zz]; % Real part of susceptibility
    xi = [xx1 xy1 xz1
          yx1 yy1 yz1
          zx1 zy1 zz1]; % Imaginary part of susceptibility
    chi0 = xr + 1i.*xi;
      
elseif nargin == 5 % Compute the complex matrix element
    deno0 = varargin{5};
    % Calculate susceptibilities along a-axis
    chi_tx  = (ttx1) .* (ttx2.') .* NN .* deno0;
    xx = sum(sum(chi_tx));
    
    % Calculate susceptibilities along b-axis
    chi_ty  = (tty1) .* (tty2.') .* NN .* deno0;
    yy = sum(sum(chi_ty));
    
    % Calculate susceptibilities along c-axis
    chi_tz  = (ttz1) .* (ttz2.') .* NN .* deno0;
    zz = sum(sum(chi_tz));
    
    % Calculate susceptibilities along ab-axis
    chi_txy  = (ttx1) .* (tty2.') .* NN .* deno0;
    xy = sum(sum(chi_txy));
    
    % Calculate susceptibilities along ac-axis
    chi_txz  = (ttx1) .* (ttz2.') .* NN .* deno0;
    xz = sum(sum(chi_txz));
    
    % Calculate susceptibilities along ba-axis
    chi_tyx  = (tty1) .* (ttx2.') .* NN .* deno0;
    yx = sum(sum(chi_tyx));
    
    % Calculate susceptibilities along bc-axis
    chi_tyz  = (tty1) .* (ttz2.') .* NN .* deno0;
    yz = sum(sum(chi_tyz));
    
    % Calculate susceptibilities along ca-axis
    chi_tzx  = (ttz1) .* (ttx2.') .* NN .* deno0;
    zx = sum(sum(chi_tzx));
    
    % Calculate susceptibilities along cb-axis
    chi_tzy  = (ttz1) .* (tty2.') .* NN .* deno0;
    zy = sum(sum(chi_tzy));
       
    chi0 = [xx xy xz
            yx yy yz
            zx zy zz];
end
end

function [qvec, cVar, freq_total, chiq, RPA_deno] = RPA(qvec, cVar, freq_total, ion, chi0, const)
unitN = 4; % Number of magnetic atoms in unit cell
lattice = ion.abc{const.elem};
Vc = sum( lattice(1,:) .* cross(lattice(2,:), lattice(3,:)) ); % Volume of unit cell [Ang^-3]
eins = diag([1 1 1]);
eins = repmat(eins,1,1,4,4);
demagn_t = ellipsoid_demagn(ion.alpha);
demagn = repmat(demagn_t,1,1,4,4);

chiq = zeros(3, 3, length(freq_total(1,:)), size(cVar,2), size(qvec,1));
D = zeros(3, 3, unitN, unitN, size(qvec,1));

% k-space calculation
parfor jj = 1:size(qvec,1)
% for jj = 1:size(qvec,1) % for debugging
    lorz_on = 1; % keep on the Lorentz term
    % Optional: turn off Lorentz term at finite q = [h k l]
    if abs(real(sum(exp(1i*2*pi*qvec(jj,:)*ion.tau'))/size(ion.tau,1))-1) > 1e-10
        lorz_on = 0;
    end
    D(:,:,:,:,jj) = const.gfac * (MF_dipole(qvec(jj,:), const.dpRng, lattice, ion.tau) - lorz_on*4*pi/Vc*(eins/3 - ion.demag*demagn))...
         + exchange(qvec(jj,:), abs(ion.ex(const.elem)), lattice, ion.tau); % [meV]
end

RPA_deno = zeros(3, 3, size(freq_total,1), size(cVar,2)); % RPA correction factor (denominator)
for nb = 1:size(cVar,2) % field/temperature iterator
    for nq = 1:size(qvec,1) % q vector iterator
        Jav = squeeze( sum(sum(D(:,:,:,:,nq),4),3)/unitN ); % [meV] average over the unit cell
        Jq = -diag(ion.renorm(const.elem,:)) .* Jav; % [meV]
        parfor nf = 1:length(freq_total(1,:))
%         for nf = 1:length(freq_total(1,:)) % for debugging
            chi_mf = squeeze(chi0(:,:,nf,nb));
            MM = chi_mf * Jq; % [meV^-1 * meV], non-commuting operation for matrices
            deno = squeeze(eye(size(MM))- MM);
            chiq(:,:,nf,nb,nq) = deno\chi_mf;
            RPA_deno(:,:,nf,nb,nq) = det(deno); % RPA denominator, save for pole analysis
        end
    end
end
% subs = ["xx", "yy", "zz"];
% for ii = 1:3
%     figure
%     hp0 = pcolor(continu_var,freq_total,squeeze((Mx(ii,ii,:,:))));
%     set(hp0, 'edgeColor','none')
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     
%     figure
%     hp0 = pcolor(continu_var,freq_total,squeeze(abs(deno(ii,ii,:,:))));
%     set(hp0, 'edgeColor','none')
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     for nq = 1:size(qvec,1)
%         figure
%         hp0 = pcolor(continu_var,freq_total,real(squeeze(chiq(ii,ii,:,:,nq))));
%         set(hp0, 'edgeColor','none');
%         colorbar
%         xlabel('Magnetic field (T)')
%         ylabel('Frequency (GHz)')
%         title(['Real part of \chi_',char(subs(ii)),'^{RPA}'])
%         figure
%         hp0 = pcolor(continu_var,freq_total,imag(squeeze(chiq(ii,ii,:,:,nq))));
%         set(hp0, 'edgeColor','none');
%         colorbar
%         xlabel('Magnetic field (T)')
%         ylabel('Frequency (GHz)')
%         title(['Imaginary part of \chi_',char(subs(ii)),'^{RPA}'])
%     end
% end
end
