function PlanetarySimulation(NumPlanets)
    close all    
    %% Make Sun
    figure('color', 'w', 'Units','normalized', 'Position',[0.0,0.05,0.73,0.8]);
    [X,Y,Z] = sphere(20); Xsun = X*2; Ysun = Y*2; Zsun = Z*2;
    surf(Xsun,Ysun,Zsun, 'FaceColor', 'y', 'EdgeAlpha', 0.1);
    daspect([1,1,1]); hold on; ax1 = gca;
    %% Make Planets, their States, Paths, and Period labels
    [Bodies, states, paths, texts, semiMajor, titleHandle] = ...
        CreatePlanets(ax1, NumPlanets);
    %% Extract Start Position and Initialize Period
    startPoints = states(:,[1,2]); Revolutions = zeros(1, NumPlanets); 
    figure('color', 'w',  'Units','normalized', 'Position',[0.73,0.05,0.28,0.8]);
    subplot(2,1,1); Periods = zeros(1, NumPlanets); 
    TR = plot(semiMajor, Periods, "ko"); hold on;  
    xlabel('a',  'Interpreter','latex'); ylabel('T',  'Interpreter','latex');
    sm = linspace(semiMajor(1),semiMajor(end)); 
    TRfit = plot(sm, zeros(1,100), "k", 'LineWidth',2); 
    TRfitText = text(sm(50), 0.5, '$T = 0\times a^{0}$', ....
        'Interpreter','latex', 'FontSize',15, 'HorizontalAlignment','center');
    subplot(2,1,2); axis([0,1,0,1]);
    % Equations
    Planet_s = {'$Planets$'};
    semiMajor_s = {'$a(m)$'};
    Period_s = {'$T(sec)$'};
    Revolution_s = {'$ No~Revs $'};
    
    for n = 1:NumPlanets
        Planet_s = [Planet_s;{['$~~~~~',num2str(n),'$']}];
        semiMajor_s = [semiMajor_s;{['$', num2str(semiMajor(n),'%.4f'),'$']}];
        Period_s = [Period_s;{['$', num2str(Periods(n),'%.4f'),'$']}];
        Revolution_s = [Revolution_s;{['$~~~~~',num2str(Revolutions(n)),'$']}]; 
    end 
    text(0.05, 0.90, Planet_s, 'VerticalAlignment', 'cap', 'FontSize', 12,...
        'interpreter', 'latex')
    text(0.30, 0.90, semiMajor_s, 'VerticalAlignment', 'cap', ....
         'FontSize', 12, 'interpreter', 'latex'); 
    PerEst = text(0.55, 0.90, Period_s, 'VerticalAlignment', 'cap', ....
         'FontSize', 12, 'interpreter', 'latex'); 
    RevCount = text(0.80, 0.90, Revolution_s, 'VerticalAlignment', 'cap', ....
         'FontSize', 12, 'interpreter', 'latex'); axis off
    %% Set Dynamics Handle
    dydt = @(y) Dynamics(y); t = 0; dt = 0.02;
    while(1)
        for i = 1:NumPlanets
            % Get Body, Path, State and Text
            Body = Bodies(i); path = paths(i); txt = texts(i); state = states(i,:); 
            % Compute Update
            dstate = rk4(dydt, state, dt); 
            % Update State, Path amd Text
	        Body.XData = Body.XData + dstate(1); 
	        Body.YData = Body.YData + dstate(2);
            pn = state(1:2); state = state + dstate; 
            states(i,:) = state; pnp1 = state(1:2);
            txt.Position = [pnp1, 1.2]; 
            txt.String = ['Rev = ', num2str(Revolutions(i))];
            if(Revolutions(i) < 1)
                path.XData = [path.XData, state(1)]; 
                path.YData = [path.YData, state(2)];
            end
            % Update Period if necesary
            if(testIntersection(startPoints(i,:), pn, pnp1))
                Revolutions(i) = Revolutions(i) + 1;
                Periods(i) = t/Revolutions(i);
                TR.YData = Periods;
                Period_s(i+1) = {['$', num2str(Periods(i),'%.4f'),'$']};
                PerEst.String = Period_s;
                Revolution_s(i+1) = {['$~~~~~',num2str(Revolutions(i)),'$']};
                RevCount.String = Revolution_s;
                %perform regression
                if(NumPlanets > 1)
                    Model = @(b, X) b(1)*X.^(b(2)); b0 = [1,1.5];
                    b = nlinfit(semiMajor(Revolutions>0), Periods(Revolutions>0),...
                         Model, b0);
                    TRfit.YData = Model(b, sm);
                    TRfitText.Position = [sm(50), TRfit.YData(50)];
                    TRfitText.String = ['$T = ', num2str(b(1)), ...
                                  '\times a^{', num2str(b(2)),'}$'];
                end
            end
        end
        t = t + dt; titleHandle.String = ['time = ', num2str(t,'%.1f')]; drawnow; 
    end
end

function [Bodies, states, paths, texts, semiMajor, titleHandle] = CreatePlanets(ax, NumPlanets)
    Rp = cumsum([2.5,rand(1,NumPlanets-1)]); Ra = cumsum([10,10*rand(1,NumPlanets-1)]);
    a  = (Rp+Ra)/2; e = (a - Rp)./a; b  = a.*sqrt(1 - e.^2);
    c  = [e.*a; zeros(1,NumPlanets)];
    states = zeros(NumPlanets,4); paths = []; texts = [];
    for n = 1:NumPlanets
        i = randi(1000); t = 2*pi*i/1000;
        x = c(1,n) + a(n)*cos(t); y = c(2,n) + b(n)*sin(t);
        l = plot(x,y,'LineWidth',2); paths = [paths, l];
        r = sqrt(x.^2 + y.^2); v = sqrt(20*pi^2*(2/r - 1/a(n)));
        vx = v*y/r; vy = -v*x/r; states(n,:) = [x,y,vx,vy]; 
        txt = text(x,y,1.2,'Rev = 0'); texts = [texts, txt];
    end
    [Bodies] = MakeSystem(ax, paths); semiMajor = a;
    titleHandle = title('time = 0', 'Interpreter','latex');
    RaMax = max(Ra);
    axis([-RaMax,RaMax,-RaMax,RaMax,-3,3]); 
end

function bl = testIntersection(startPoint, pn, pnp1)
    bl = Intersect([0,0], 1.1*startPoint, pn, pnp1);
end

function ds = Dynamics(y)
    ds = [y(3), y(4), -20*pi^2*y(1)/norm(y([1,2]))^3, -20*pi^2*y(2)/norm(y([1,2]))^3];
end

function dy = rk4(dydt, y, dt)
    k1 = dydt(y); k2 = dydt(y + dt*k1/2);
    k3 = dydt(y + dt*k2/2); k4 = dydt(y + dt*k3);
    dy = dt*(k1+2*k2+2*k3+k4)/6;
end
  
function [Bodies] = MakeSystem(ax, paths)
    Bodies = [];  axes(ax); hold on
    [Xs,Ys,Zs] = sphere(20);
    for i = 1:numel(paths)
        factor = 0.5+0.5*rand;
        Z2 = Zs*factor;
        X2 = Xs*factor + paths(i).XData; 
        Y2 = Ys*factor + paths(i).YData; 
        Body = surf(X2,Y2,Z2, 'EdgeAlpha', 0.1, ...
                 'FaceColor', paths(i).Color);
        Bodies = [Bodies, Body];
    end
end

function bl = Intersect(a, b, x, y)
    den0 = (b(1) - a(1)) * (x(2) - y(2)) - (b(2) - a(2)) * (x(1) - y(1));
    num1 = (x(1) - y(1)) * (a(2) - x(2)) - (x(2) - y(2)) * (a(1) - x(1));
    num2 = (a(1) - x(1)) * (b(2) - a(2)) - (a(2) - x(2)) * (b(1) - a(1));
    num1 = num1*(abs(num1) > 1e-10); num2 = num2*(abs(num2) > 1e-10);
    s = num1 / den0; t = num2 / den0; bl = (0 < s && s < 1) && (0 < t && t < 1);
end