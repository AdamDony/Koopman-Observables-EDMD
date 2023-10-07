function c = max_c(x1_data, x2_data, V, V_dot)
    
    figure
    cont = contour(x1_data, x2_data, V_dot, [0 0], 'r');
    
    numLines = 0;
    index = 1;
    shapes = [];
    while index < length(cont)
        data = [];
        numPoints = cont(2, index);
        numLines = numLines + 1;
        data = cont(:, index+1:index+numPoints);
        pgon = polyshape(data.');
        if pgon.NumRegions ~= 0
            shapes = [shapes, pgon];
        end
%         plot(pgon)
%         hold on

        index = index + numPoints + 1;
    end
    
    % Initial try to get max_c
    norm_vec = zeros(1, length(cont));
    for i = 1:length(cont)
        norm_vec(i) = norm(cont(:, i));
    end
    
    normv = norm_vec(norm_vec~=0);
    r = min(norm_vec(norm_vec~=0));
    nx = norm_vec;
    ind = find(norm_vec==r);
    p = cont(:, ind(1));
    
    syms x1 x2
    c = subs(V, [x1 x2], [p(1) p(2)]);

    % Check if it works and that V is enclosed within where V_dot is
    % negative
    fp = fimplicit(V == c);
    V_shape = polyshape(fp.XData, fp.YData);
    plot(V_shape)

    tf = zeros(1, length(shapes))
    for s = 1:length(shapes)
        % Should be zero if V = c doesnt overlap with V_dot > 0
        tf(s) = overlaps(V_shape, shapes(s));
    end
    
    % If V = c cuts into where V_dot is positive, we need to decrease c
    if any(tf)
        c = find_c(V, c, shapes, c/2, 1);
    end
    
    close

    % Generate Figure
    figure('Position', [100, 100, 1100, 400]);
    
    % Plot 1 - V_dot(x)
    subplot(1, 2, 1)
    surf(x1_data, x2_data, V_dot) %Plot the surface
    hold on
    contour(x1_data, x2_data, V_dot, [0 0], 'r');
    xlabel("x_1")
    ylabel("x_2")
    zlabel('$\dot{V}(x)$', 'Interpreter','latex')
    
    % Plot 2 - V_dot(x) = 0, V(x) < c
    subplot(1, 2, 2)
    cont = contour(x1_data, x2_data, V_dot, [0 0], 'r');
    grid
    hold on
    xline(0)
    yline(0)
    xlabel("x_1")
    ylabel("x_2")
    fimplicit(V == c)

end

function c_opt = find_c(V, c, shapes, diff, iterations)
    fp = fimplicit(V == c);
    V_shape = polyshape(fp.XData, fp.YData);

    tf = zeros(1, length(shapes));
    for s = 1:length(shapes)
        tf(s) = overlaps(V_shape, shapes(s));
    end
    if any(tf)
        c_opt = find_c(V, c - diff, shapes, diff/2, iterations + 1);
    else
        if iterations > 10
            c_opt = c;
        else
            c_opt = find_c(V, c + diff, shapes, diff/2, iterations + 1);
        end
    end
end


function point = getClosestBoundary(p, rad)
    x = p(1)
    y = p(2)
    % if x coordinate is closer to boundary
    if abs(abs(x) - rad) < abs(abs(y) - rad)
        % if x is closer to +rad
        if abs(rad - x) < abs(-rad - x)
            point = [rad; y];
        % if x is close to -rad
        else
            point = [-rad; y];
        end
    % if y coordinate is close to boundary
    else
        % if y is closer to +rad
        if abs(rad - y) < abs(-rad - y)
            point = [x; rad];
        % if y is close to -rad
        else
            point = [x; -rad];
        end
    end
end



