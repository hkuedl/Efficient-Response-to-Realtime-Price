clc; close all; clear all;

% Set random seed for reproducibility
rng(20);

Total_1 = zeros(4,10); %4: Ground and three comfort; 10: days
Total_2 = zeros(4,10);
Total_3 = zeros(4,10);

Price_samss = [10,50,100];

for opt_day = 1:10 %5 = Feb 5
    for Price_sam_i = 1:3
        for Utili = 3:3  %1-3: "Abs", "Quad","Piece"
            Price_sam = Price_samss(Price_sam_i);

            c_time = 288; % time length
            
            % Parameters
            N_days = 10; %1-6, Feb, 2023
            if Utili == 1
                U_pen = 4.5;
            elseif Utili == 2
                U_pen = 3.5;
            elseif Utili == 3
                U_k = [-10, -8.2, -6.5, -4.2, -1.4, 0.8, 2.8, 4.9, 6.6, 20];
                U_b = [-10, -8, -6, -3, 0, 1, 3, 5, 9];
                U_start = 3*1.4 + 3*4.2 + 2*6.5 + 2*8.2;
                U_s = U_start;
                for i = 2:length(U_b)
                    c_sum = U_start;
                    for ik = 1:i-1
                        c_sum = c_sum + (U_b(ik+1) - U_b(ik)) * U_k(ik+1);
                    end
                    U_s = [U_s, c_sum];
                end
            end
            P_min = 0.0; P_max = 2.0; % MW
            eta_ = 0.93; beta_ = 2.2;
            S0 = 20;
            Sref = [20, repmat(20, 1, 12 * 8), repmat(22, 1, 12 * 4), repmat(21, 1, 12 * 4), repmat(19, 1, 12 * 8)];
            
            % Load data
            data_in = load('Data.mat');
            delta = data_in.delta;
            Pri_rt_Feb = data_in.Pri_rt_Feb;
            Pri_da_Feb = data_in.Pri_da_Feb;
            
            % Initialize result matrices
            Cr_obj_da = zeros(N_days, 3);
            Cr_obj_act = zeros(N_days, 6);
            
            % Optimization for each day
            for d = opt_day:opt_day
                % Define variables
                c_v_tem = sdpvar(c_time+1, 1);
                if Utili == 3
                    c_v_tem_piece = sdpvar(c_time, 1);
                    c_v_tem_piece_i = sdpvar(c_time, length(U_b)-1);
                    c_v_tem_z = binvar(c_time,length(U_b)-1);
                end
                c_v_p = sdpvar(c_time, 1);
                
                % Objective functions
                c_obj1 = 0;
                c_obj1_rt = 0;
                for i = 1:c_time
                    c_obj1 = c_obj1 + (1/12)*Pri_da_Feb((d-1)*288+i,1)*c_v_p(i,1);
                    c_obj1_rt = c_obj1_rt + (1/12)*Pri_rt_Feb((d-1)*288+i,1)*c_v_p(i,1);
                end
                
                % Utility function (Absolute style)
                VT = 0;
                if Utili == 1
                    c_obj2 = sum(U_pen*abs(c_v_tem(2:end,1))) + VT;
                elseif Utili == 2
                    c_obj2 = sum(U_pen*(c_v_tem(2:end,1)).^2) + VT;
                elseif Utili == 3
                    c_obj2 = sum(sum(c_v_tem_piece_i)) + VT;
                end
                c_obj = c_obj1 + c_obj2;
                c_obj_rt = c_obj1_rt + c_obj2;
                
                % Constraints
                cons = [c_v_p <= P_max, c_v_p >= P_min, c_v_tem(1,1) == 0];
                for t = 1:c_time
                    cons = [cons, c_v_tem(t+1,1) == eta_*c_v_tem(t,1) + beta_*c_v_p(t,1) + ...
                           (delta(t,d) - Sref(t+1) + eta_*Sref(t))];
                end
                if Utili == 3
                    MM = 30;
                    for t = 1:c_time
                        for zz = 1:length(U_b)-1
                            cons = [cons, c_v_tem(1+t,1) >= U_b(zz)*c_v_tem_z(t,zz) - MM*(1-c_v_tem_z(t,zz))];
                            cons = [cons, c_v_tem(1+t,1) <= U_b(zz+1)*c_v_tem_z(t,zz) + MM*(1-c_v_tem_z(t,zz))];
                            cons = [cons, c_v_tem_piece_i(t,zz) == (U_s(zz)+U_k(1+zz)*(c_v_tem(1+t,1)-U_b(zz)))*c_v_tem_z(t,zz)];
                        end
                        cons = [cons, sum(c_v_tem_z(t,1:end)) == 1];
                    end
                end
                
                % Solve with day-ahead price
                options = sdpsettings('solver', 'gurobi', 'verbose', 0);
                optimize(cons, c_obj, options);
                %disp(value(c_v_p));
                % Store results
                Cr_obj_da(d,:) = [value(c_obj1), value(c_obj2), value(c_obj)];
                obj_act = 0;
                for i = 1:c_time
                    obj_act = obj_act + (1/12)*Pri_rt_Feb((d-1)*288+i,1)*value(c_v_p(i,1));
                end
                Cr_obj_act(d,1:3) = [obj_act, value(c_obj2)-value(VT), obj_act+value(c_obj2)-value(VT)];
                
                % Solve with real-time price
                optimize(cons, c_obj_rt, options);
                Cr_obj_act(d,4:6) = [value(c_obj1_rt), value(c_obj2)-value(VT), value(c_obj_rt)-value(VT)];
                
                % Store power and temperature
                c_p = value(c_v_p);
                c_tem = value(c_v_tem);
                
                % % Display results
                % disp('solve with real-time price');
                % disp(Cr_obj_act(d,4:6));
                % disp('solve with day-ahead price');
                % disp(Cr_obj_act(d,1:3));
                % 
                % disp(Pri_rt_Feb((d-1)*288+1:(d-1)*288+12,1));
                %disp(c_p);
                %disp(min(c_tem));
                %disp(max(c_tem));
                
            end
            disp(['RT-DA: ' num2str(Cr_obj_act(d,6)) ' ' num2str(Cr_obj_act(d,3))]);
            if Utili == 1
                Total_1(1,opt_day) = Cr_obj_act(d,6);
            elseif Utili == 2
                Total_2(1,opt_day) = Cr_obj_act(d,6);
            elseif Utili == 3
                Total_3(1,opt_day) = Cr_obj_act(d,6);
            end


            % Proposed method: using percentiles to represent the price uncertainty
            quantiles_all = load('Data_quan.mat');
            s_sam_number = 2001;
            quantiles_x = quantiles_all.(['x_' num2str(Price_sam)])(1,2:end-1);
            s_sam = -10 + 0.01*(0:s_sam_number-1)';
            u = zeros(s_sam_number, 2);
            u(:,1) = s_sam;
            if Utili == 1
                u(:,2) = U_pen*sign(u(:,1));
            elseif Utili == 2
                u(:,2) = U_pen*2*u(:,1);
            elseif Utili == 3
                for ii = 1:s_sam_number
                    indices = find(U_b <= u(ii,1));
                    left_index = indices(end);
                    u(ii,2) = U_k(left_index+1);
                end
            end

            v = zeros(N_days, c_time, s_sam_number);
            v_3 = zeros(N_days, c_time, s_sam_number);
            w_huizong = zeros(N_days, c_time, s_sam_number, 2);

            for d = opt_day:opt_day
                for t = c_time:-1:1
                    if t == c_time
                        v(d,t,:) = zeros(1, 1, s_sam_number);
                    else
                        quantile_y = quantiles_all.(['y_' num2str(d-1) '_' num2str(Price_sam)]);
                        quantiles_rt = quantile_y(t+1,:);

                        for ii = 1:s_sam_number
                            [v(d,t,ii), w_x1, w_x2, v_3(d,t,ii)] = ...
                                F_value(s_sam, quantiles_x, quantiles_rt, s_sam(ii), ...
                                       eta_, beta_, u, squeeze(v(d,t+1,:)), P_max, ...
                                       delta(t,d) - Sref(t+1) + eta_*Sref(t));
                            w_huizong(d,t,ii,:) = [w_x1, w_x2];
                        end
                    end
                end
            end

            aa = squeeze(v(d,:,:));
            % figure;
            % plot(aa');
            % xlabel('Time step');
            % ylabel('Value function');
            % figure;
            % plot(1:size(aa,1), max(aa,[],2), 'b-', 1:size(aa,1), min(aa,[],2), 'r-');
            % xlabel('Time step');
            % legend('Max', 'Min');

            aa1 = diff(aa, 1, 2);
            aa2 = aa1 > 0.01;
            % disp(['Number of irnormal value: ', num2str(sum(aa2(:)))]);
            % disp(['Max. irnormal value: ', num2str(max(aa1(:)))])
            % Dynamic programming approach
            DP_obj_act = zeros(N_days, 3);
            DP_v = zeros(N_days, c_time);
            DP_tem = zeros(N_days, c_time);
            s_cand = -10 + 0.01*(0:s_sam_number-1);

            for d = opt_day:opt_day
                tic;
                for t = 1:c_time
                    if t == 1
                        tem_ini = 0;
                    else
                        tem_ini = DP_tem(d,t-1);
                    end

                    v_dt = squeeze(v(d,t,:));
                    v_dt_0 = diff(v_dt);
                    v_dt1 = [v_dt(1); v_dt_0];

                    price_rt = Pri_rt_Feb((d-1)*288+t,1);
                    p_2 = P_max;
                    s_2 = eta_*tem_ini + beta_*P_max + (delta(t,d) - Sref(t+1) + eta_*Sref(t));
                    p_3 = P_min;
                    s_3 = eta_*tem_ini + beta_*P_min + (delta(t,d) - Sref(t+1) + eta_*Sref(t));

                    [~, min_ind] = min(abs(s_cand - s_3));
                    [~, max_ind] = min(abs(s_cand - s_2));

                    s = s_cand(min_ind:max_ind);
                    p = (s - eta_*tem_ini - (delta(t,d) - Sref(t+1) + eta_*Sref(t))) / beta_;

                    if p(1) < P_min
                        s = s(2:end);
                        p = p(2:end);
                    end
                    if p(end) > P_max
                        s = s(1:end-1);
                        p = p(1:end-1);
                    end

                    c = zeros(length(s), 4);
                    for i = 1:length(s)
                        c(i,1) = (1/12)*price_rt*p(i);
                        if Utili == 1
                            c(i,2) = U_pen*abs(s(i));
                        elseif Utili == 2
                            c(i,2) = U_pen*(s(i)^2);
                        elseif Utili == 3
                            indices = find(U_b <= s(i));
                            left_index = indices(end);
                            c(i,2) = U_s(left_index) + (s(i)-U_b(left_index))*U_k(left_index+1);
                        end
                        c(i,3) = -V_cal(v_dt1, s_sam, s(i));
                        c(i,4) = sum(c(i,1:3));
                    end

                    [~, idx] = min(c(:,4));
                    DP_v(d,t) = p(idx);
                    DP_tem(d,t) = s(idx);
                end
                time_elapsed = toc;

                if DP_v(d,t) > P_max || DP_v(d,t) < P_min
                    disp('Wrong!');
                end

                DP_obj_act(d,1) = sum((1/12)*Pri_rt_Feb((d-1)*288+(1:c_time),1).*DP_v(d,1:c_time)');
                if Utili == 1
                    DP_obj_act(d,2) = sum(U_pen*abs(DP_tem(d,1:c_time)));
                elseif Utili == 2
                    DP_obj_act(d,2) = sum(U_pen*(DP_tem(d,1:c_time)).^2);
                elseif Utili == 3
                    DP_obj_act(d,2) = 0;
                    for ikk = 1:c_time
                        DP_i = DP_tem(d,ikk);
                        indices = find(U_b <= DP_i);
                        left_index = indices(end);
                        DP_obj_i = U_s(left_index) + (DP_i-U_b(left_index))*U_k(left_index+1);
                        DP_obj_act(d,2) = DP_obj_act(d,2) + DP_obj_i;
                    end
                end
                DP_obj_act(d,3) = DP_obj_act(d,1) + DP_obj_act(d,2);

                % disp(['time: ' num2str(time_elapsed)]);
                % disp(['solution: ' num2str(DP_obj_act(d,3))]);
                %disp('TCL: '); disp(DP_v(d,1:c_time)');
                % disp('Temperature: '); disp(DP_tem(d,1:c_time)');
                disp(min(DP_tem(d,1:c_time)));
                disp(max(DP_tem(d,1:c_time)));
                disp(min(DP_v(d,1:c_time)));
                disp(max(DP_v(d,1:c_time)));
                
            end

            disp(['Pro: ' num2str(DP_obj_act(d,3))]);
            if Utili == 1
                Total_1(1+Price_sam_i,opt_day) = DP_obj_act(d,3);
            elseif Utili == 2
                Total_2(1+Price_sam_i,opt_day) = DP_obj_act(d,3);
            elseif Utili == 3
                Total_3(1+Price_sam_i,opt_day) = DP_obj_act(d,3);
            end
        end
    end
end



% Functions
function [F_val, w_x1, w_x2, num_filtered] = F_value(s_sam, quantiles_x, quantiles_rt, s, eta_, beta_, u, v, Pmax, delta)
    x1 = eta_*s+ beta_*Pmax + delta;
    x2 = eta_*s + delta;

    [w_x1, w_x2] = deal(v_u(s_sam, v, u, x1), v_u(s_sam, v, u, x2));

    P1 = eta_*w_x1*F_price(quantiles_x, quantiles_rt, 12*beta_*w_x1);
    P2 = eta_*w_x2*(1 - F_price(quantiles_x, quantiles_rt, 12*beta_*w_x2));

    weights0 = diff(quantiles_x);
    weights = [weights0(1), weights0];

    mask = (quantiles_rt >= 12*beta_*w_x1) & (quantiles_rt <= 12*beta_*w_x2);
    filtered_values = quantiles_rt(mask);
    filtered_weights = weights(mask);

    P3 = (eta_/(12*beta_))*sum(filtered_values .* filtered_weights);
    F_val = P1 + P2 + P3;
    num_filtered = length(filtered_values);
end

function p = F_price(quantiles_x, quantiles_rt, x)
    if x < quantiles_rt(1)
        p = 0;
    elseif x > quantiles_rt(end)
        p = 1;
    elseif x == quantiles_rt(end)
        p = quantiles_x(end);
    else
        for i = 1:length(quantiles_x)-1
            if x >= quantiles_rt(i) && x < quantiles_rt(i+1)
                p = quantiles_x(i);
                return;
            end
        end
    end
end

function diff_vu = v_u(s_sam, v, u, x)
    if x < u(1,1)
        index_u = 1;
    elseif x > u(end,1)
        index_u = size(u,1);
    else
        [~, index_u] = min(abs(u(:,1) - x));
    end
    ux = u(index_u,2);

    if x < s_sam(1)
        index_v = 1;
    elseif x > s_sam(end)
        index_v = length(s_sam);
    else
        [~, index_v] = min(abs(s_sam - x));
    end
    vx = v(index_v);

    diff_vu = vx - ux;
end

function Vx = V_cal(v_dt1, s_sam, x)
    Vx = 0;
    for i = 1:length(s_sam)
        Vx = Vx + v_dt1(i)*max(x - s_sam(i), 0);
    end
end
