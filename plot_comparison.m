function plot_comparison(cVar, freq_total, chi_ini, chi_emt, scanMode, dscrt_var)
    fprintf('\n=== Creating comparison plots ===\n');

    figure('Position', [100, 100, 1800, 1000]);
    components = [1, 2, 3];
    comp_labels = {'xx', 'yy', 'zz'};

    cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
        [linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
        [ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
    cmap = flip(cmap,1);

    switch scanMode
        case 'field'
            xlab = 'Magnetic Field (T)';
            title_prefix = sprintf('T = %.3f K', dscrt_var);
        case 'temp'
            xlab = 'Temperature (K)';
            title_prefix = sprintf('B = %.3f T', dscrt_var);
        otherwise
            xlab = 'cVar';
            title_prefix = '';
    end

    for jj = 1:3
        comp_idx = components(jj);

        subplot(2, 3, jj);
        chi_rpa_comp = squeeze(chi_ini(comp_idx, comp_idx, :, :, :));
        chi_rpa_avg = mean(chi_rpa_comp, 3);
        data_rpa = mag2db(abs(imag(chi_rpa_avg)));
        hp = pcolor(cVar, freq_total, data_rpa);
        set(hp, 'EdgeColor', 'none');
        colormap(gca, cmap);
        colorbar;
        ylabel('Frequency (GHz)');
        title(sprintf('RPA: Im[\\chi^{%s}] (%s)', comp_labels{jj}, title_prefix), 'Interpreter', 'tex');
        set(gca, 'FontSize', 10);

        subplot(2, 3, jj + 3);
        chi_emt_comp = squeeze(chi_emt(comp_idx, comp_idx, :, :));
        data_emt = mag2db(abs(imag(chi_emt_comp)));
        hp = pcolor(cVar, freq_total, data_emt);
        set(hp, 'EdgeColor', 'none');
        colormap(gca, cmap);
        colorbar;
        xlabel(xlab);
        ylabel('Frequency (GHz)');
        title(sprintf('Eff. Medium: Im[\\chi^{%s}]', comp_labels{jj}), 'Interpreter', 'tex');
        set(gca, 'FontSize', 10);
    end

    sgtitle('Comparison: RPA vs Effective Medium Theory', 'FontSize', 14, 'FontWeight', 'bold');
end