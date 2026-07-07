function CubeSatInertia_Revised_Discrete()
% CUBESATINERTIA  Block placement optimizer for 3U CubeSat inertia matching.
%
% PURPOSE
%   Determines the placement and mass of discrete tuning blocks within
%   the 3D-printed 3U CubeSat dynamic model so that the assembled hardware
%   recreates a target inertia tensor within a specified tolerance (1%).
%
% AXIS CONVENTION (matches SolidWorks)
%   X = cross-section width axis
%   Y = long (stacking) axis of the 3U structure
%   Z = cross-section depth axis

    clc; clear; close all;
    %% =====================================================================
    %  SECTION 1: USER INPUTS — Edit these values to match the hardware
    %  =====================================================================
    % --- Target inertia tensor (Cerberus Target) ---
    % AXES FIXED: Y and Z have been swapped so the long axis matches SolidWorks.
    Target_COM_cm      = [0.097, -0.711, -0.283];   % Z origin shifted to geometric center 
    Target_I_diag_gcm2 = [294229.69, 94578.57, 284732.34];
    Target_I_off_gcm2  = [2840.56, 87.86, 3663.86];

    % --- Static hardware properties (from SolidWorks at full density) ---
    Printed_1U_Mass_g  = 277.46;                        
    SW_Mass_g          = 2188.08;                       
    SW_COM_cm          = [0.00, -0.04, -0.01];          
    SW_I_diag_gcm2     = [157212.03, 48215.13, 157208.24];  
    SW_I_off_gcm2      = [0.00, 0.00, -0.47];          

    % --- HARDWARE CONSTRAINTS: Exact Available Blocks ---
    % Enter an array of the EXACT physical block weights (in grams) available. 
    % The solver will only use these specified masses (0g is always included automatically).
    Available_Block_Masses_g = [15, 27]; 

    % --- Slot geometry (interior cavities of the printed shell) ---
    Box_Size_Cross_cm = 1.52;                           
    Box_Size_Stack_cm = 1.56;                           
    Pitch_Cross_cm    = 1.72;                           
    nSlots_X = 5;                                       
    nSlots_Z = 5;                                       

    % --- Accuracy tolerances ---
    MOI_Tolerance_Pct = 1;                            
    Off_Diag_Tol_Pct  = 1;                            
    COM_Tolerance_cm  = 100;                          

    % --- Solver settings ---
    Num_Trials = 24;

    %% =====================================================================
    %  SECTION 1.5: DERIVED CONSTANTS — DO NOT EDIT
    %  =====================================================================
    Printed_3U_Mass_g  = Printed_1U_Mass_g * 3;
    Infill_Ratio       = Printed_3U_Mass_g / SW_Mass_g;
    Static_Mass_g      = Printed_3U_Mass_g;
    Static_COM_cm      = SW_COM_cm;
    Static_I_diag_gcm2 = SW_I_diag_gcm2 * Infill_Ratio;
    Static_I_off_gcm2  = SW_I_off_gcm2 * Infill_Ratio;
    
    % Sort and include 0g in the available mass array
    Allowed_M_g = unique([0, Available_Block_Masses_g]);
    Max_Mass_g  = max(Allowed_M_g);
    Min_Active_Mass_g = min(Allowed_M_g(Allowed_M_g > 0));

    fprintf('=== CONFIG ===\n');
    fprintf('  Axis convention: X = width, Y = stacking (long), Z = depth\n');
    fprintf('  SW Mass: %.2f g | Infill Ratio: %.4f | Static Mass: %.2f g\n', ...
            SW_Mass_g, Infill_Ratio, Static_Mass_g);
    fprintf('  Allowed Masses (g): [%s]\n', join(string(Allowed_M_g), ', '));

    I_target_matrix = [Target_I_diag_gcm2(1), Target_I_off_gcm2(1), Target_I_off_gcm2(2);
                       Target_I_off_gcm2(1), Target_I_diag_gcm2(2), Target_I_off_gcm2(3);
                       Target_I_off_gcm2(2), Target_I_off_gcm2(3), Target_I_diag_gcm2(3)];
    target_eigenvalues = eig(I_target_matrix);
    if any(target_eigenvalues <= 0)
        error('Target inertia tensor is not positive definite.');
    end

    %% =====================================================================
    %  SECTION 2: UNIT CONVERSION (cm, g) -> SI (m, kg) — DO NOT EDIT
    %  =====================================================================
    Target_COM    = Target_COM_cm    / 100;
    Target_I_diag = Target_I_diag_gcm2 / 1e7;
    Target_I_off  = Target_I_off_gcm2  / 1e7;
    S_M      = Static_Mass_g    / 1000;
    S_C      = Static_COM_cm    / 100;
    S_I_diag = Static_I_diag_gcm2 / 1e7;
    S_I_off  = Static_I_off_gcm2  / 1e7;
    
    Allowed_M     = Allowed_M_g / 1000;
    Max_Mass      = Max_Mass_g  / 1000;
    Min_Mass      = Min_Active_Mass_g / 1000;
    
    Box_XZ = Box_Size_Cross_cm / 100;                   
    Box_Y  = Box_Size_Stack_cm / 100;                   

    %% =====================================================================
    %  SECTION 3: GRID GENERATION — DO NOT EDIT
    %  =====================================================================
    W = 0.60; F = 0.20; B = 1.56;
    y_1U = zeros(1, 5);
    for L = 1:5
        y_1U(L) = W + F + (B/2) + (L-1)*(B + F);
    end
    height_1U = 10.0;
    y_3U      = [y_1U, y_1U + height_1U, y_1U + 2*height_1U];
    y_vec     = (y_3U - (height_1U * 3) / 2) / 100;    
    x_vec = ((1:nSlots_X) - (nSlots_X+1)/2) * (Pitch_Cross_cm / 100);
    z_vec = ((1:nSlots_Z) - (nSlots_Z+1)/2) * (Pitch_Cross_cm / 100);
    [X, Y, Z] = ndgrid(x_vec, y_vec, z_vec);
    Fixed_Pos = [X(:), Y(:), Z(:)];

    %% =====================================================================
    %  SECTION 4: MULTI-START OPTIMIZATION — DO NOT EDIT
    %  =====================================================================
    fprintf('\nStarting Optimization (%d Trials)...\n', Num_Trials);
    current_pool = gcp('nocreate');
    if isempty(current_pool)
        fprintf('Initializing parallel pool...\n');
        parpool('local');
    end
    options_strict = optimoptions('lsqnonlin', ...
        'Display', 'none', ...
        'Algorithm', 'trust-region-reflective', ...
        'StepTolerance', 1e-7, ...
        'FunctionTolerance', 1e-7, ...
        'UseParallel', false);
    w = struct('Reg', 0.01, 'C', 0, 'I_diag', 100000, 'I_off', 100000);
    num_slots = size(Fixed_Pos, 1);
    trial_costs = inf(1, Num_Trials);
    trial_vars  = cell(1, Num_Trials);
    D = parallel.pool.DataQueue;
    trials_completed = 0;
    fprintf('Optimization Progress: %5.1f%%', 0);
    afterEach(D, @nUpdateProgress);

    parfor i = 1:Num_Trials
        stream = RandStream('mlfg6331_64', 'Seed', i*7919 + 13);
        initial_masses = rand(stream, 1, num_slots) * Max_Mass;
        lb = zeros(1, num_slots);
        ub = ones(1, num_slots) * Max_Mass;

        curr_vars  = initial_masses;
        rw_weights = ones(1, num_slots);
        rw_epsilon = Max_Mass * 0.05;
        for rw_iter = 1:5
            ResidualFcn_rw = @(v) compute_residuals_rw(v, Fixed_Pos, Target_COM, ...
                Target_I_diag, Target_I_off, S_M, S_C, S_I_diag, S_I_off, ...
                w, Box_XZ, Box_Y, rw_weights);
            [curr_vars, ~] = lsqnonlin(ResidualFcn_rw, curr_vars, lb, ub, options_strict);
            rw_weights = 1 ./ (curr_vars + rw_epsilon);
            rw_weights = rw_weights / max(rw_weights);
        end

        ResidualFcn_clean = @(v) compute_residuals(v, Fixed_Pos, Target_COM, ...
            Target_I_diag, Target_I_off, S_M, S_C, S_I_diag, S_I_off, ...
            w, Box_XZ, Box_Y);
        [curr_vars, ~] = lsqnonlin(ResidualFcn_clean, curr_vars, lb, ub, options_strict);

        max_consolidation_iters = 50;
        iter_count = 0;
        while iter_count < max_consolidation_iters
            iter_count = iter_count + 1;
            dust_idx = find(curr_vars > 1e-4 & curr_vars < Min_Mass/2);
            if isempty(dust_idx), break; end
            [~, sort_order] = sort(curr_vars(dust_idx));
            num_to_kill = max(1, round(0.25 * length(dust_idx)));
            victim_idx  = dust_idx(sort_order(1:num_to_kill));
            lb(victim_idx) = 0;
            ub(victim_idx) = 0;
            curr_vars(victim_idx) = 0;
            [curr_vars, ~] = lsqnonlin(ResidualFcn_clean, curr_vars, lb, ub, options_strict);
        end

        survivors = find(curr_vars > 1e-4);
        [~, removal_order] = sort(curr_vars(survivors));
        survivors = survivors(removal_order);
        for s_idx = 1:length(survivors)
            slot = survivors(s_idx);
            if curr_vars(slot) < 1e-4, continue; end
            test_vars = curr_vars;
            test_vars(slot) = 0;
            lb_test = lb; ub_test = ub;
            lb_test(slot) = 0; ub_test(slot) = 0;
            [test_vars, ~] = lsqnonlin(ResidualFcn_clean, test_vars, lb_test, ub_test, options_strict);
            if check_accuracy(test_vars, Fixed_Pos, Target_COM, Target_I_diag, ...
                              Target_I_off, S_M, S_C, S_I_diag, S_I_off, ...
                              Box_XZ, Box_Y, MOI_Tolerance_Pct, ...
                              COM_Tolerance_cm, Off_Diag_Tol_Pct)
                curr_vars = test_vars;
                lb = lb_test; ub = ub_test;
            end
        end
        final_resnorm  = sum(ResidualFcn_clean(curr_vars).^2);
        trial_costs(i) = final_resnorm;
        trial_vars{i}  = curr_vars;
        send(D, i);
    end
    fprintf('\n');

    %% =====================================================================
    %  SECTION 5: BEST-TRIAL SELECTION — DO NOT EDIT
    %  =====================================================================
    acceptable   = false(1, Num_Trials);
    block_counts = inf(1, Num_Trials);
    for i = 1:Num_Trials
        if isempty(trial_vars{i}), continue; end
        if check_accuracy(trial_vars{i}, Fixed_Pos, Target_COM, Target_I_diag, ...
                          Target_I_off, S_M, S_C, S_I_diag, S_I_off, ...
                          Box_XZ, Box_Y, MOI_Tolerance_Pct, ...
                          COM_Tolerance_cm, Off_Diag_Tol_Pct)
            acceptable(i)   = true;
            block_counts(i) = sum(trial_vars{i} > 1e-4);
        end
    end
    if any(acceptable)
        candidates = find(acceptable);
        min_blocks = min(block_counts(candidates));
        finalists  = candidates(block_counts(candidates) == min_blocks);
        [~, fin_idx] = min(trial_costs(finalists));
        best_idx = finalists(fin_idx);
        fprintf('Selected trial %d: %d blocks, cost %.4e (within tolerances)\n', ...
                best_idx, min_blocks, trial_costs(best_idx));
    else
        [~, best_idx] = min(trial_costs);
        fprintf('WARNING: No trial met tolerances. Using best-cost trial %d (%d blocks).\n', ...
                best_idx, sum(trial_vars{best_idx} > 1e-4));
    end
    Global_Best_Vars = trial_vars{best_idx};

    %% =====================================================================
    %  SECTION 6: EXACT DISCRETIZATION ENGINE — DO NOT EDIT
    %  =====================================================================
    fprintf('\nPhase 2: Snapping to exact allowed inventory blocks...\n');
    Discrete_Masses = zeros(size(Global_Best_Vars));
    for k = 1:length(Global_Best_Vars)
        % Find the absolute closest allowed block mass for each continuous mass
        [~, nearest_idx] = min(abs(Allowed_M - Global_Best_Vars(k)));
        Discrete_Masses(k) = Allowed_M(nearest_idx);
    end

    fprintf('Phase 3: Pruning unnecessary blocks...\n');
    survivors = find(Discrete_Masses > 1e-4);
    [~, removal_order] = sort(Discrete_Masses(survivors));
    survivors = survivors(removal_order);
    removed_count = 0;
    for s_idx = 1:length(survivors)
        slot = survivors(s_idx);
        if Discrete_Masses(slot) < 1e-4, continue; end
        test_masses = Discrete_Masses;
        test_masses(slot) = 0;
        if check_accuracy(test_masses, Fixed_Pos, Target_COM, Target_I_diag, ...
                          Target_I_off, S_M, S_C, S_I_diag, S_I_off, ...
                          Box_XZ, Box_Y, MOI_Tolerance_Pct, ...
                          COM_Tolerance_cm, Off_Diag_Tol_Pct)
            Discrete_Masses = test_masses;
            removed_count = removed_count + 1;
        end
    end
    fprintf('Pruning removed %d additional blocks.\n', removed_count);
    
    fprintf('Phase 4: Reducing oversized blocks...\n');
    downsize_count = 0;
    survivors = find(Discrete_Masses > 1e-4);
    for s_idx = 1:length(survivors)
        slot = survivors(s_idx);
        curr_val = Discrete_Masses(slot);
        
        % Locate current block's position in the allowed inventory array
        idx_match = find(abs(Allowed_M - curr_val) < 1e-6, 1); 
        
        % Try stepping down to the next smallest available block size
        while idx_match > 1 
            test_masses = Discrete_Masses;
            idx_match = idx_match - 1;
            test_masses(slot) = Allowed_M(idx_match);
            
            if check_accuracy(test_masses, Fixed_Pos, Target_COM, Target_I_diag, ...
                              Target_I_off, S_M, S_C, S_I_diag, S_I_off, ...
                              Box_XZ, Box_Y, MOI_Tolerance_Pct, ...
                              COM_Tolerance_cm, Off_Diag_Tol_Pct)
                Discrete_Masses = test_masses;
                downsize_count = downsize_count + 1;
                if idx_match == 1, break; end 
            else
                break; 
            end
        end
    end
    fprintf('Downsizing reduced mass on %d blocks.\n', downsize_count);

    continuous_cost   = trial_costs(best_idx);
    discrete_cost     = sum(compute_residuals(Discrete_Masses, Fixed_Pos, ...
        Target_COM, Target_I_diag, Target_I_off, S_M, S_C, S_I_diag, S_I_off, ...
        w, Box_XZ, Box_Y).^2);
    final_block_count = sum(Discrete_Masses > 1e-4);
    
    fprintf('Final result: %d blocks, total tuning mass %.1f g\n', ...
            final_block_count, sum(Discrete_Masses)*1000);
    
    final_passes = check_accuracy(Discrete_Masses, Fixed_Pos, Target_COM, ...
                                   Target_I_diag, Target_I_off, ...
                                   S_M, S_C, S_I_diag, S_I_off, ...
                                   Box_XZ, Box_Y, MOI_Tolerance_Pct, ...
                                   COM_Tolerance_cm, Off_Diag_Tol_Pct);
    if ~final_passes
        fprintf(['\nNOTE: The final discrete solution does not meet the\n' ...
                 '      %.2f%% MOI tolerance on at least one component.\n'], MOI_Tolerance_Pct);
    end

    %% =====================================================================
    %  SECTION 7: REPORT GENERATION — DO NOT EDIT
    %  =====================================================================
    generate_report(Discrete_Masses, Fixed_Pos, Pitch_Cross_cm/100, y_vec, ...
                    nSlots_X, nSlots_Z, Target_COM, Target_I_diag, Target_I_off, ...
                    S_M, S_C, S_I_diag, S_I_off, Box_XZ, Box_Y, ...
                    MOI_Tolerance_Pct, Off_Diag_Tol_Pct);
    
    function nUpdateProgress(~)
        trials_completed = trials_completed + 1;
        pct = (trials_completed / Num_Trials) * 100;
        fprintf(repmat('\b', 1, 6));
        fprintf('%5.1f%%', pct);
    end
end

%% =========================================================================
%  HELPER FUNCTIONS (Resdiuals, Accuracy, Print, Plot) 
%  =========================================================================
function F = compute_residuals(masses, pos, T_C, T_I_diag, T_I_off, ...
                                S_M, S_C, S_I_diag, S_I_off, w, bXZ, bY)
    masses = masses(:)';
    curr_M = sum(masses) + S_M;
    if curr_M < 1e-6
        F = 1e6 * ones(length(masses) + 9, 1);
        return;
    end
    tuning_moment = sum(pos .* masses', 1);
    static_moment = S_M * S_C;
    curr_C = (tuning_moment + static_moment) / curr_M;
    dx = pos(:,1) - curr_C(1);
    dy = pos(:,2) - curr_C(2);
    dz = pos(:,3) - curr_C(3);
    Ixx_box = (1/12) * masses' * (bY^2  + bXZ^2);      
    Iyy_box = (1/12) * masses' * (bXZ^2 + bXZ^2);      
    Izz_box = (1/12) * masses' * (bXZ^2 + bY^2);       
    Ixx_w = sum(Ixx_box + masses' .* (dy.^2 + dz.^2));
    Iyy_w = sum(Iyy_box + masses' .* (dx.^2 + dz.^2));
    Izz_w = sum(Izz_box + masses' .* (dx.^2 + dy.^2));
    Ixy_w = -sum(masses' .* (dx .* dy));
    Ixz_w = -sum(masses' .* (dx .* dz));
    Iyz_w = -sum(masses' .* (dy .* dz));
    dx_s = S_C(1) - curr_C(1);
    dy_s = S_C(2) - curr_C(2);
    dz_s = S_C(3) - curr_C(3);
    Ixx_s = S_I_diag(1) + S_M * (dy_s^2 + dz_s^2);
    Iyy_s = S_I_diag(2) + S_M * (dx_s^2 + dz_s^2);
    Izz_s = S_I_diag(3) + S_M * (dx_s^2 + dy_s^2);
    Ixy_s = S_I_off(1) - S_M * (dx_s * dy_s);
    Ixz_s = S_I_off(2) - S_M * (dx_s * dz_s);
    Iyz_s = S_I_off(3) - S_M * (dy_s * dz_s);
    Ixx = Ixx_w + Ixx_s; Iyy = Iyy_w + Iyy_s; Izz = Izz_w + Izz_s;
    Ixy = Ixy_w + Ixy_s; Ixz = Ixz_w + Ixz_s; Iyz = Iyz_w + Iyz_s;
    ref_I = max(T_I_diag); if ref_I < 1e-9, ref_I = 1; end
    safe_diag_div = max(T_I_diag, 1e-9);
    F_C      = sqrt(w.C)      * (curr_C - T_C);
    F_I_diag = sqrt(w.I_diag) * ([Ixx, Iyy, Izz] - T_I_diag) ./ safe_diag_div;
    F_I_off  = sqrt(w.I_off)  * ([Ixy, Ixz, Iyz] - T_I_off) / ref_I;
    F_Reg    = sqrt(w.Reg)    * masses(:);
    F = [F_C(:); F_I_diag(:); F_I_off(:); F_Reg(:)];
end

function F = compute_residuals_rw(masses, pos, T_C, T_I_diag, T_I_off, ...
                                   S_M, S_C, S_I_diag, S_I_off, w, bXZ, bY, rw_weights)
    F = compute_residuals(masses, pos, T_C, T_I_diag, T_I_off, ...
                          S_M, S_C, S_I_diag, S_I_off, w, bXZ, bY);
    n = length(masses);
    F(end-n+1:end) = sqrt(w.Reg) * rw_weights(:) .* masses(:);
end

function ok = check_accuracy(masses, pos, T_C, T_I_diag, T_I_off, ...
                              S_M, S_C, S_I_diag, S_I_off, bXZ, bY, ...
                              moi_tol_pct, com_tol_cm, off_tol_pct)
    masses = masses(:)';
    curr_M = sum(masses) + S_M;
    if curr_M < 1e-6, ok = false; return; end
    curr_C = (sum(pos .* masses', 1) + S_M * S_C) / curr_M;
    if any(abs(curr_C - T_C) > com_tol_cm/100)
        ok = false; return;
    end
    dx = pos(:,1) - curr_C(1); dy = pos(:,2) - curr_C(2); dz = pos(:,3) - curr_C(3);
    Ixx_box = (1/12) * masses' * (bY^2  + bXZ^2);
    Iyy_box = (1/12) * masses' * (bXZ^2 + bXZ^2);
    Izz_box = (1/12) * masses' * (bXZ^2 + bY^2);
    Ixx = sum(Ixx_box + masses' .* (dy.^2 + dz.^2));
    Iyy = sum(Iyy_box + masses' .* (dx.^2 + dz.^2));
    Izz = sum(Izz_box + masses' .* (dx.^2 + dy.^2));
    Ixy = -sum(masses' .* (dx .* dy)); Ixz = -sum(masses' .* (dx .* dz)); Iyz = -sum(masses' .* (dy .* dz));
    dx_s = S_C(1) - curr_C(1); dy_s = S_C(2) - curr_C(2); dz_s = S_C(3) - curr_C(3);
    Ixx = Ixx + S_I_diag(1) + S_M * (dy_s^2 + dz_s^2);
    Iyy = Iyy + S_I_diag(2) + S_M * (dx_s^2 + dz_s^2);
    Izz = Izz + S_I_diag(3) + S_M * (dx_s^2 + dy_s^2);
    Ixy = Ixy + S_I_off(1) - S_M * (dx_s * dy_s);
    Ixz = Ixz + S_I_off(2) - S_M * (dx_s * dz_s);
    Iyz = Iyz + S_I_off(3) - S_M * (dy_s * dz_s);
    diag_actual = [Ixx, Iyy, Izz];
    safe_target = max(abs(T_I_diag), 1e-9);
    if any(abs(diag_actual - T_I_diag) ./ safe_target > moi_tol_pct/100)
        ok = false; return;
    end
    ref = max(T_I_diag); off_actual = [Ixy, Ixz, Iyz];
    if any(abs(off_actual - T_I_off) > (off_tol_pct/100) * ref)
        ok = false; return;
    end
    ok = true;
end

function generate_report(masses, pos, pitch_cross, y_vec, nX, nZ, T_C, ...
                         T_I_diag, T_I_off, S_M, S_C, S_I_diag, S_I_off, bXZ, bY, ...
                         moi_tol_pct, off_tol_pct)
    active_mask = masses > 0.001; active_pos  = pos(active_mask, :); active_mass = masses(active_mask);
    curr_M = sum(active_mass) + S_M;
    curr_C = (sum(active_pos .* active_mass', 1) + S_M * S_C) / curr_M;
    dx = active_pos(:,1) - curr_C(1); dy = active_pos(:,2) - curr_C(2); dz = active_pos(:,3) - curr_C(3);
    Ixx_box = (1/12) * active_mass' * (bY^2  + bXZ^2); Iyy_box = (1/12) * active_mass' * (bXZ^2 + bXZ^2); Izz_box = (1/12) * active_mass' * (bXZ^2 + bY^2);
    Ixx_w = sum(Ixx_box + active_mass' .* (dy.^2 + dz.^2)); Iyy_w = sum(Iyy_box + active_mass' .* (dx.^2 + dz.^2)); Izz_w = sum(Izz_box + active_mass' .* (dx.^2 + dy.^2));
    Ixy_w = -sum(active_mass' .* (dx .* dy)); Ixz_w = -sum(active_mass' .* (dx .* dz)); Iyz_w = -sum(active_mass' .* (dy .* dz));
    dx_s = S_C(1) - curr_C(1); dy_s = S_C(2) - curr_C(2); dz_s = S_C(3) - curr_C(3);
    Ixx_s = S_I_diag(1) + S_M * (dy_s^2 + dz_s^2); Iyy_s = S_I_diag(2) + S_M * (dx_s^2 + dz_s^2); Izz_s = S_I_diag(3) + S_M * (dx_s^2 + dy_s^2);
    Ixy_s = S_I_off(1) - S_M * (dx_s * dy_s); Ixz_s = S_I_off(2) - S_M * (dx_s * dz_s); Iyz_s = S_I_off(3) - S_M * (dy_s * dz_s);
    Tensor_gcm2 = [Ixx_w+Ixx_s, Ixy_w+Ixy_s, Ixz_w+Ixz_s; Ixy_w+Ixy_s, Iyy_w+Iyy_s, Iyz_w+Iyz_s; Ixz_w+Ixz_s, Iyz_w+Iyz_s, Izz_w+Izz_s] * 1e7;
    curr_M_g = curr_M * 1000; T_C_cm = T_C * 100; curr_C_cm = curr_C * 100; T_I_diag_gcm2 = T_I_diag * 1e7; T_I_off_gcm2 = T_I_off * 1e7;
    fprintf('\n=================================================================\n');
    fprintf('             3U CUBESAT INERTIA OPTIMIZATION REPORT              \n');
    fprintf('=================================================================\n');
    fprintf('Static Hardware Mass:      %.0f g\n', S_M * 1000);
    fprintf('Active Tuning Blocks:      %d\n\n', length(active_mass));
    fprintf('| %-16s | %-12s | %-12s | %-12s | %-10s | %-6s |\n', 'Metric', 'Target', 'Actual', 'Diff', 'Error (%)', 'Status');
    fprintf('|%s|\n', repmat('-',1,85));
    fprintf('| %-16s | %-12s | %-12.4f | %-12s | %-10s | %-6s |\n', 'Mass (g)', 'INFO', curr_M_g, '-', '-', 'N/A');
    fprintf('| %-16s | %-12.4f | %-12.4f | %-12.4e | %-10s | %-6s |\n', 'COM X (cm)', T_C_cm(1), curr_C_cm(1), abs(T_C_cm(1)-curr_C_cm(1)), '-', 'INFO');
    fprintf('| %-16s | %-12.4f | %-12.4f | %-12.4e | %-10s | %-6s |\n', 'COM Y (cm)', T_C_cm(2), curr_C_cm(2), abs(T_C_cm(2)-curr_C_cm(2)), '-', 'INFO');
    fprintf('| %-16s | %-12.4f | %-12.4f | %-12.4e | %-10s | %-6s |\n', 'COM Z (cm)', T_C_cm(3), curr_C_cm(3), abs(T_C_cm(3)-curr_C_cm(3)), '-', 'INFO');
    moi_tol_frac = moi_tol_pct / 100;
    print_row('I_xx (g*cm^2)', T_I_diag_gcm2(1), Tensor_gcm2(1,1), moi_tol_frac * max(T_I_diag_gcm2(1), 1));
    print_row('I_yy (g*cm^2)', T_I_diag_gcm2(2), Tensor_gcm2(2,2), moi_tol_frac * max(T_I_diag_gcm2(2), 1));
    print_row('I_zz (g*cm^2)', T_I_diag_gcm2(3), Tensor_gcm2(3,3), moi_tol_frac * max(T_I_diag_gcm2(3), 1));
    max_I_ref = max(T_I_diag_gcm2); cross_axis_tol = (off_tol_pct / 100) * max_I_ref; if cross_axis_tol == 0, cross_axis_tol = 1; end
    print_row('I_xy (g*cm^2)', T_I_off_gcm2(1), Tensor_gcm2(1,2), cross_axis_tol, max_I_ref);
    print_row('I_xz (g*cm^2)', T_I_off_gcm2(2), Tensor_gcm2(1,3), cross_axis_tol, max_I_ref);
    print_row('I_yz (g*cm^2)', T_I_off_gcm2(3), Tensor_gcm2(2,3), cross_axis_tol, max_I_ref);
    fprintf('\n--- ASSEMBLED INERTIA TENSOR (g*cm^2) ---\n');
    fprintf('X [%9.2f  %9.2f  %9.2f]\n', Tensor_gcm2(1,1), Tensor_gcm2(1,2), Tensor_gcm2(1,3));
    fprintf('Y [%9.2f  %9.2f  %9.2f]\n', Tensor_gcm2(2,1), Tensor_gcm2(2,2), Tensor_gcm2(2,3));
    fprintf('Z [%9.2f  %9.2f  %9.2f]\n\n', Tensor_gcm2(3,1), Tensor_gcm2(3,2), Tensor_gcm2(3,3));
    
    fprintf('--- BLOCK PLACEMENT (Active Slots Only) ---\n');
    fprintf('| %-10s | %-10s | %-10s | %-10s |\n', 'Mass (g)', 'X (cm)', 'Y (cm)', 'Z (cm)');
    fprintf('|%s|\n', repmat('-', 51, 1));
    for b = 1:length(active_mass)
        fprintf('| %-10.0f | %-10.2f | %-10.2f | %-10.2f |\n', active_mass(b)*1000, active_pos(b,1)*100, active_pos(b,2)*100, active_pos(b,3)*100);
    end
    fprintf('\n');

    figure('Color','w', 'Name', '3U Inertia Optimization Result', 'NumberTitle', 'off');
    hold on; axis equal; grid on; hLx = (bXZ * 100) / 2; hLy = (bY  * 100) / 2; hLz = (bXZ * 100) / 2;
    unique_masses_g = unique(round(active_mass * 1000)); num_unique = length(unique_masses_g);
    color_palette = lines(max(num_unique, 1)); legend_handles = gobjects(num_unique, 1); legend_labels  = cell(num_unique, 1);
    for i = 1:num_unique
        legend_handles(i) = patch(NaN, NaN, color_palette(i,:), 'EdgeColor', 'k', 'FaceAlpha', 0.9);
        legend_labels{i} = sprintf('%d g', unique_masses_g(i));
    end
    for b = 1:length(active_mass)
        xc = active_pos(b,1)*100; yc = active_pos(b,2)*100; zc = active_pos(b,3)*100;
        v = [xc-hLx, zc-hLz, yc-hLy; xc+hLx, zc-hLz, yc-hLy; xc+hLx, zc+hLz, yc-hLy; xc-hLx, zc+hLz, yc-hLy;
             xc-hLx, zc-hLz, yc+hLy; xc+hLx, zc-hLz, yc+hLy; xc+hLx, zc+hLz, yc+hLy; xc-hLx, zc+hLz, yc+hLy];
        f = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
        mass_g = round(active_mass(b) * 1000); color_idx = find(unique_masses_g == mass_g, 1);
        patch('Vertices', v, 'Faces', f, 'FaceColor', color_palette(color_idx,:), 'EdgeColor', 'k', 'FaceAlpha', 0.9);
    end
    if num_unique > 0, lgd = legend(legend_handles, legend_labels, 'Location', 'bestoutside'); title(lgd, 'Block Mass'); end
    pitch_cross_cm = pitch_cross * 100; cross_ticks = ((1:nX) - (nX+1)/2) * pitch_cross_cm;
    set(gca, 'XTick', cross_ticks, 'YTick', cross_ticks, 'ZTick', y_vec*100);
    xl = xlabel('X - Width (cm)'); yl = ylabel('Z - Depth (cm)'); zl = zlabel('Y - Stacking (cm)');
    set(xl, 'FontSize', 10, 'FontWeight', 'bold'); set(yl, 'FontSize', 10, 'FontWeight', 'bold'); set(zl, 'FontSize', 10, 'FontWeight', 'bold');    
    xlim([-5, 5]); ylim([-5, 5]); zlim([-16, 16]);
    plot3(T_C(1)*100, T_C(3)*100, T_C(2)*100, 'rx', 'MarkerSize', 15, 'LineWidth', 3, 'DisplayName', 'Target COM');
    plot3(curr_C_cm(1), curr_C_cm(3), curr_C_cm(2), 'go', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Actual COM');
    title(sprintf('3U Optimized Configuration (%d Blocks)', length(active_mass))); view(135, 25); set(gca, 'FontSize', 8);
end

function print_row(name, target, actual, tol, ref_val)
    if nargin < 5, ref_val = abs(target); end
    diff_val = abs(target - actual); status = 'PASS';
    if diff_val > tol, status = 'FAIL'; end
    if abs(ref_val) > 1e-9
        err_pct = (diff_val / abs(ref_val)) * 100; pct_str = sprintf('%8.3f%%', err_pct);
    else
        pct_str = '    -     ';
    end
    fprintf('| %-16s | %-12.4f | %-12.4f | %-12.4e | %-10s | %-6s |\n', name, target, actual, diff_val, pct_str, status);
end