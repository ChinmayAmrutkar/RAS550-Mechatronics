fname = 'heater_data_plot3-2.xlsx';
T = readtable(fname);
t = T.(1);              % Time
y = T.(2);              % Temperature
u = T.(3);              % duty cycle
plot(t,y,t,u);

dt = t(2) - t(1);
fs = 1/dt;

seg1_index = (1:500)*fs -1;
seg2_index = (500:1000)*fs -1;
seg3_index = (1000:1500)*fs -1;

% plot(t(seg2_index),y(seg2_index),t(seg2_index),u(seg2_index))

%
% Model
% y = y0 + G * step(t - 10 - theta_p) * (1 - exp(-(t - (10 + theta_p))/tau))
% parameters p = [y0, G, tau, theta_p]
model = @(p,tl) p(1) + p(2) .* (tl -10 >= p(4)) .*(1 - exp(-(tl - (10 + p(4))) ./ max(p(3),eps)));

wins = [0 500; 500 1000; 1000 1500];   % [start end) per segment

names = {'y_init','G','tau','theta_p'};
P = zeros(3,4);

figure('Color','w'); tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

for k = 1:3
    w   = wins(k,:);
    idx = (t >= w(1)) & (t < w(2));  
    ts  = t(idx);                       
    ys  = y(idx); 
    us = u(idx);

    % local time axis to stabilize fit
    tl  = ts - ts(1);

    % initial guesse
    y0   = ys(1);
    G0   = ys(end) - y0;
    tau0 = max((tl(end)-tl(1))/5, 5);
    th0  = 5;
    p0   = [y0, G0, tau0, th0];

    % bounds
    lb = [min(ys)-5, -200, 1e-3, -50];
    ub = [max(ys)+5,  200, 5e4,   tl(end)+50];

    % fit using least square method
    opts = optimoptions('lsqcurvefit','Display','off','MaxIterations',2000);
    p = lsqcurvefit(@(pp,xx) model(pp,xx), p0, tl, ys, lb, ub, opts);
    P(k,:) = p;

    %plots
    nexttile;
    plot(tl, ys, 'k-', 'DisplayName','data'); hold on;
    plot(tl, us, 'r--', 'DisplayName','DutyCycle');
    plot(tl, model(p, tl), 'LineWidth',1.6,'lineStyle','-', 'DisplayName','fit');
    grid on; xlabel('t (s, local)'); ylabel('Temp (°C)');
    title(sprintf('Segment %d: [%g, %g) s', k, w(1), w(2)));
    legend('Location','best');
end

ParamTable = array2table(P, 'VariableNames', names, 'RowNames', {'Seg1','Seg2','Seg3'});
disp(ParamTable);
