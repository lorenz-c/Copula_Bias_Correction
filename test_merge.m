for i = 1:230
    if ~isnan(C_smap.scan.corr(i)) && ~isnan(C_smos.scan.corr(i)) && ...
            valid_dta_smap.scan(i) > 100 && valid_dta_smos.scan(i) > 100
        
        XY     = [TS1.scan.Data.soil_moisture(i, :)' TS3.scan.Data.soil_moisture(i, :)'];
        XY(XY(:, 1) < 1e-4) = 0;
        XY(isnan(XY(:, 1)), :) = [];
        XY(isnan(XY(:, 2)), :) = [];
       
        
        
        if size(XY, 1) > 50
            X_pred = TS2.scan.Data.soil_moisture(i, :)';
            
            val_dta = find(~isnan(X_pred));
            X_pred(isnan(X_pred)) = [];

        
            [s_dC, theta_rl, sel_family] = copula_merge(XY(:, 1), XY(:, 2), X_pred);
        
            T_out_mn = NaN(1, length(TS2.scan.TimeStamp));
            T_out_md = NaN(1, length(TS2.scan.TimeStamp));

            T_out_mn(val_dta) = nanmean(s_dC{1});
            T_out_md(val_dta) = nanmedian(s_dC{1});

            C_out(i, 1) = nancorr(TS2.scan.Data.soil_moisture(i, :)', ...
                           ismn_at_smos.scan.Data.soil_moisture(i, :)');
            C_out(i, 2) = nancorr(T_out_mn', ...
                           ismn_at_smos.scan.Data.soil_moisture(i, :)');
            C_out(i, 3) = nancorr(T_out_md', ...
                           ismn_at_smos.scan.Data.soil_moisture(i, :)');
            
                       
            D1 = TS2.scan.Data.soil_moisture(i, :)' - ismn_at_smos.scan.Data.soil_moisture(i, :)';
            D2 = T_out_mn' - ismn_at_smos.scan.Data.soil_moisture(i, :)';
            D3 = T_out_md' - ismn_at_smos.scan.Data.soil_moisture(i, :)';
            
            R_out(i, 1) = nanrms(D1, 1)/nanmean(ismn_at_smos.scan.Data.soil_moisture(i, :)');
            R_out(i, 2) = nanrms(D2, 1)/nanmean(ismn_at_smos.scan.Data.soil_moisture(i, :)');
            R_out(i, 3) = nanrms(D3, 1)/nanmean(ismn_at_smos.scan.Data.soil_moisture(i, :)');
            
            C_out(i, :)
            R_out(i, :)
            
            TS_out_mn(i, 1:length(T_out_mn)) = T_out_mn';
            TS_out_md(i, 1:length(T_out_mn)) = T_out_md';
        else
            TS_out_mn(i, 1:length(TS1.scan.TimeStamp)) = NaN;
            TS_out_md(i, 1:length(TS1.scan.TimeStamp)) = NaN;
        
            C_out(i, 1) = NaN;
            C_out(i, 2) = NaN;
            C_out(i, 3) = NaN;
                       
            R_out(i, 1) = NaN;
            R_out(i, 2) = NaN;
            R_out(i, 3) = NaN;
        end
        
    else
        TS_out_mn(i, 1:length(TS1.scan.TimeStamp)) = NaN;
        TS_out_md(i, 1:length(TS1.scan.TimeStamp)) = NaN;
        
        C_out(i, 1) = NaN;
        C_out(i, 2) = NaN;
        C_out(i, 3) = NaN;
                       
        R_out(i, 1) = NaN;
        R_out(i, 2) = NaN;
        R_out(i, 3) = NaN;
    end
    i
end
