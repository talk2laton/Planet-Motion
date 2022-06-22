function PlanetarySimulationAreaLaw(f1, f2, f3, f4)
    close all    
    %% Make Sun
    figure('color', 'w', 'Units','normalized', 'Position',[0.0,0.05,0.73,0.8]);
    [X,Y,Z] = sphere(20); Xsun = X*2; Ysun = Y*2; Zsun = Z*2;
    surf(Xsun,Ysun,Zsun, 'FaceColor', 'y', 'EdgeAlpha', 0.1);
    daspect([1,1,1]); hold on; ax1 = gca;
    %% Make Planets, their States, Paths, and Period labels
    [Body, state, path, txt, Period, titleHandle] = ...
        CreatePlanets(ax1);
    %% Extract Start Position and Initialize Period
    startPoints = state(:,[1,2]); Revolutions = 0; 
    poly1 = fill(0, 0, 'g', 'FaceAlpha',0.5); 
    poly2 = fill(0, 0, 'r', 'FaceAlpha',0.5);
    t1 = f1*Period; t2 = f2*Period; t3 = f3*Period; t4 = f4*Period;
    %% Set Dynamics Handle
    dydt = @(y) Dynamics(y); t = 0; dt = 0.1;
    while(1)
        % Compute Update
        dstate = rk4(dydt, state, dt); 
        % Update State, Path amd Text
        Body.XData = Body.XData + dstate(1); 
        Body.YData = Body.YData + dstate(2);
        pn = state(1:2); 
        state = state + dstate; 
        pnp1 = state(1:2);
        txt.Position = [pnp1, 1.2]; 
        txt.String = ['Rev = ', num2str(Revolutions(1))];
        if(Revolutions(1) < 1)
            path.XData = [path.XData, state(1)]; 
            path.YData = [path.YData, state(2)];
        end
        % Update Period if necesary
        if(testIntersection(startPoints(1,:), pn, pnp1))
            Revolutions(1) = Revolutions(1) + 1;
        end
        t = t + dt; titleHandle.String = ['time = ', num2str(t)]; drawnow;
        if(t1 <= t && t <= t2)
            poly1.XData = [poly1.XData; state(1)];
            poly1.YData = [poly1.YData; state(2)];
        end
        if(t3 <= t && t <= t4)
            poly2.XData = [poly2.XData; state(1)];
            poly2.YData = [poly2.YData; state(2)];
        end
    end
end

function [Body, state, path, txt, Period, titleHandle] =....
    CreatePlanets(ax)
    Rp = 2+2*rand; Ra = 10+20*rand; a  = (Rp+Ra)/2; e = (a - Rp)./a; 
    b  = a.*sqrt(1 - e.^2); c  = [e.*a; 0]; Period = sqrt(a^3/5);
    i = randi(1000); t = 2*pi*i/1000;
    x = c(1) + a*cos(t); y = c(2) + b*sin(t);
    path = plot(x,y,'LineWidth',2); 
    r = sqrt(x.^2 + y.^2); v = sqrt(20*pi^2*(2/r - 1/a));
    vx = v*y/r; vy = -v*x/r; state = [x,y,vx,vy]; 
    txt = text(x,y,1.2,'Rev = 0'); 
    Body = MakeSystem(ax, state(:,[1,2]));
    titleHandle = title('time = 0', 'Interpreter','latex');
    RaMax = max(Ra);
    axis([-RaMax,RaMax,-RaMax,RaMax,-3,3]); 
end

function bl = testIntersection(startPoint, pn, pnp1)
    x0 = startPoint(1); y0 = startPoint(2);
    bl = (pn(1)-x0)*(pnp1(1)-x0) < 0 && (pn(2)-y0)*(pnp1(2)-y0) < 0;
end

function ds = Dynamics(y)
    ds = [y(3), y(4), -20*pi^2*y(1)/norm(y([1,2]))^3, -20*pi^2*y(2)/norm(y([1,2]))^3];
end

function dy = rk4(dydt, y, dt)
    k1 = dydt(y); k2 = dydt(y + dt*k1/2);
    k3 = dydt(y + dt*k2/2); k4 = dydt(y + dt*k3);
    dy = dt*(k1+2*k2+2*k3+k4)/6;
end
  
function [Bodies] = MakeSystem(ax, loc)
    Bodies = [];  axes(ax); hold on
    [Xs,Ys,Zs] = sphere(20);
    for i = 1:size(loc,1)
        factor = 0.5+0.5*rand;
        Z2 = Zs*factor;
        X2 = Xs*factor + loc(i,1); 
        Y2 = Ys*factor + loc(i,2); 
        Body = surf(X2,Y2,Z2, 'EdgeAlpha', 0.1, ...
                 'FaceColor', rand(1,3));
        Bodies = [Bodies, Body];
    end
end