%newton raphson done
%secant done
%bisection done
%fixed point done
%false position done
%muller done
%bairstow done
prompt = 'Is equation a polynomial?: ';
t = input(prompt,'s');
if strcmp(t,'Y') 
    func = input('Enter a polynomial: ','s');
    f = str2sym(func);
    prompt = 'What is your choice of method?: ';
    reply = input(prompt, 's');
    if strcmp(reply,'Muller')
        y = matlabFunction(str2sym(func));
        x0 = input('Enter x0: ');
        x1= input('Enter x1: ');
        x2 = input('Enter x2: ');
        relApproxErr = input('What would you like as the maximum relative approximate error(number you enter will be taken as percentage): ');
        maxIter = input('What should be the maximum iteration number the program should run: ');
        Z = zeros;
        flag = 0;
        for i = 1: maxIter
            c = y(x2);
            b = (y(x2) - y(x1))/(x2 - x1) + (x2-x1)*(((y(x2) - y(x1))/(x2 - x1)) - ((y(x1) - y(x0))/(x1 - x0)))/(x2-x0);
            a = (((y(x2) - y(x1))/(x2 - x1)) - ((y(x1) - y(x0))/(x1 - x0)))/(x2-x0);
            D = (b^2 - 4*a*c)^(1/2);
            if abs(b-D) < abs(b+D)
                E = b + D;
            else
                E = b - D;
            end
            h = -2*c/E;
            p = x2 + h;            
            errPercent = abs(((p-x2)/p)*100);
            if errPercent > relApproxErr
                x0 = x1;
                x1 = x2;
                x2 = p;
                Z(i) = errPercent;
            end
            if errPercent < relApproxErr
                flag = 1;
            break;
            end
        end
        if flag == 1
            fprintf('Convergence criterion for maximum relative error is met\n');
        end
        if flag == 0
            fprintf('Convergence criterion for maximum iteration is met\n');
        end
        fprintf('The root of the function: %f\n', p);
        figure(1);
        fplot(f,[(p - 5) (p + 5)]);
        figure(2);
        plot(Z);
    end
    
    %Bairstow
    if strcmp(reply,'Bairstow')
        r = input('Enter starting point r: ');
        s = input('Enter starting point s: ');
        relApproxErr = input('What would you like as the maximum relative approximate error(number you enter will be taken as percentage): ');
        maxIter = input('What should be the maximum iteration number the program should run: ');
        y = matlabFunction(f);
        n = polynomialDegree(f) + 1;
        m = n - 1;
        a = sym2poly(f);
        a = fliplr(a);
        Z = zeros;
        Y = zeros;
        while n > 0
            for j = 1: maxIter
                b = zeros;
                c = zeros;
                b(n) = a(n);
                b(n-1) = a(n-1) + r*b(n);
                c(n) = b(n);
                c(n-1) = b(n-1) + r*c(n);
                for i = n-2: -1: 1 
                    b(i) = a(i) + r*b(i+1) + s*b(i+2);
                    c(i) = b(i) + r*c(i+1) + s*c(i+2);
                end              
                syms x y 
                eqn1 = c(2)*x + c(3)*y == -b(1);
                eqn2 = c(3)*x + c(4)*y == -b(2);              
                sol = solve([eqn1, eqn2], [x, y]);
                delr = sol.x;
                dels = sol.y;       
                errPercentr = abs((delr/(r+delr))*100);       
                errPercents = abs((dels/(r+dels))*100);
                Z(j) = errPercentr;
                Y(j) = errPercents;
                if errPercents > relApproxErr || errPercentr > relApproxErr 
                    r = r + delr;
                    s = s + dels;
                end 
                if errPercents < relApproxErr && errPercentr < relApproxErr
                    r = r + delr;
                    s = s + dels;
                break;
                end
            end
            m = m - 2;
            if (r^2 + 4*s) < 0
                x1 = r/2;
                x2 = (-(r^2 + 4*s))/2;
                fprintf('The two roots of the functions are: %f + %fi, %f - %fi\n', x1, x2, x1, x2);
            else    
                x1 = (r + (r^2 + 4*s)^(1/2))/2;
                x2 = (r - (r^2 + 4*s)^(1/2))/2;    
                fprintf('The two roots of the functions are: %f, %f\n', x1, x2);
            end 
            if m > 2
                a = ones;
                for k = n : -1: 3
                   a(k-2) = b(k);
                end  
            end
            n = n - 2;
            if m == 2
                if (b(4)^2 - 4*b(5)*b(3)) < 0
                    x3 = -b(4)/2*b(5);
                    x4 = (-(b(4)^2 - 4*b(5)*b(3)))^(1/2)/2*b(5);
                    fprintf('The other two roots of the functions are: %f + %fi, %f - %fi\n', x3, x4, x3, x4);
                else
                    x3 = (-b(4) + (b(4)^2 - 4*b(5)*b(3))^(1/2))/2*b(5);
                    x4 = (-b(4) - (b(4)^2 - 4*b(5)*b(3))^(1/2))/2*b(5);
                    fprintf('The other two roots of the functions are: %f, %f\n', x3, x4);
                end
            break;
            end
            if m == 1
                x1 = -b(3)/b(4);
                fprintf('The other root of the functions is: %f\n', x1);
            break;
            end
        end
        figure(1);
        ezplot(f);
        grid on
        figure(2);
        plot(Z);
        figure(3);
        plot(Y);
    end
end
if strcmp(t,'N')
    prompt = 'What is your choice of method?: ';
    reply = input(prompt, 's');
    
    %BISECTION
    if strcmp(reply,'Bisection')
        func = input('Enter a function: ','s');
        x1 = input('Enter x1: ');
        x2 = input('Enter x2: ');
        relApproxErr = input('What would you like as the maximum relative approximate error(number you enter will be taken as percentage): ');
        funcConverg = input('What would you like as the convergence criterion for function value(number you enter will be taken as percentage): ');
        maxIter = input('What should be the maximum iteration number the program should run: ');
        y = matlabFunction(str2sym(func));
        Z = zeros;
        flag = 0;
        if y(x1)*y(x2) > 0
            fprintf('No roots exist within the given interval\n');
        end
        if y(x1) == 0
            fprintf('x1 is one of the roots\n')
        elseif y(x2) == 0
            fprintf('x2 is one of the roots\n')
        end
        for i = 1: maxIter
            xh = (x1+x2)/2;
            if abs(y(xh)) < funcConverg/100
                flag = 2;
            break ;
            end
            if y(x1)*y(xh) < 0
                errPercent = abs(((x1-xh)/xh)*100);
                x2 = xh;
                Z(i) = errPercent;
            else
                errPercent = abs(((x2-xh)/xh)*100);
                x1 = xh;
                Z(i) = errPercent;
            end
            if errPercent < relApproxErr
                flag = 1;
            break ;
            end    
        end  
        if flag == 2
            fprintf('Convergence criterion for function value is met\n');
        end
        if flag == 1
            fprintf('Convergence criterion for maximum relative error is met\n');
        end
        if flag == 0
            fprintf('Convergence criterion for maximum iteration is met\n');
        end
            
        fprintf('The root of the function: %f\n', xh);
        figure(1);
        ezplot(func);
        grid on
        figure(2);
        plot(Z);
    end
    
    %False-position
    if strcmp(reply,'False-position')
        func = input('Enter a function: ','s');
        x1 = input('Enter x1: ');
        x2 = input('Enter x2: ');
        relApproxErr = input('Enter relative error limit: ');
        funcConverg = input('Enter function value limit: ');
        maxIter = input('Enter iteration number limit: ');
        flag = 0;
        y = matlabFunction(str2sym(func));
        Z = zeros;
        if y(x1)*y(x2) > 0
            fprintf('No roots exist within the given interval\n');
        end
        if y(x1) == 0
            fprintf('x1 is one of the roots\n')
        elseif y(x2) == 0
            fprintf('x2 is one of the roots\n')
        end
        for i = 1: maxIter
            xh = x1 - y(x1) * (x2 - x1)/(y(x2) - y(x1));
            if abs(xh-x1) > abs(xh-x2)
                errPercent = abs(((x2-xh)/xh)*100);
            else
                errPercent = abs(((x1-xh)/xh)*100);
            end
            if abs(y(xh)) < funcConverg/100
                flag = 2;
            break;
            end
            if y(x1)*y(xh) < 0
                x2 = xh;
                Z(i) = errPercent;
            else
                x1 = xh;
                Z(i) = errPercent;
            end
            if errPercent < relApproxErr
                flag = 1;
            break;
            end    
        end  
        if flag == 2
            fprintf('Convergence criterion for function value is met\n');
        end
        if flag == 1
            fprintf('Convergence criterion for maximum relative error is met\n');
        end
        if flag == 0
            fprintf('Convergence criterion for maximum iteration is met\n');
        end
            
        fprintf('The root of the function: %f\n', xh);
        figure(1);
        ezplot(func);
        grid on
        figure(2);
        plot(Z);
    end 
    
    %Fixed-point
    if strcmp(reply,'Fixed-point')
        func = input('Enter a function g(x): ','s');
        f = input('Enter a function f(x): ','s');
        x1 = input('Enter x1: ');
        relApproxErr = input('Enter relative error limit: ');
        funcConverg = input('Enter function value limit: ');
        maxIter = input('Enter iteration number limit: ');
        flag = 0;
        y = matlabFunction(str2sym(func));
        Z = zeros;
        if y(x1) == 0
            fprintf('x1 is one of the roots\n')
        end
        for i = 1 : maxIter
            xh = y(x1);
            errPercent = abs(((x1-xh)/xh)*100);
            if errPercent > relApproxErr
                x1 = xh;
                Z(i) = errPercent;
            end
            if abs(y(xh)) < funcConverg/100
                flag = 2;
            break;
            end
            if errPercent < relApproxErr
                flag = 1;
            break;
            end
        end
        if flag == 2
            fprintf('Convergence criterion for function value is met\n');
        end
        if flag == 1
            fprintf('Convergence criterion for maximum relative error is met\n');
        end
        if flag == 0
            fprintf('Convergence criterion for maximum iteration is met\n');
        end
        fprintf('The root of the function: %f\n', xh);
        figure(1);
        ezplot(f);
        grid on
        figure(2);
        plot(Z);
    end 
    
    %Newton-Raphson
    if strcmp(reply,'Newton-Raphson')
        func = input('Enter a function: ','s');
        x1 = input('Enter x1: ');
        relApproxErr = input('Enter relative error limit: ');
        funcConverg = input('Enter function value limit: ');
        maxIter = input('Enter iteration number limit: ');
        flag = 0;
        y = matlabFunction(str2sym(func));
        f = str2sym(func);
        f = diff(f);
        y1 = matlabFunction(f);
        Z = zeros;
        if y(x1) == 0
            fprintf('x1 is one of the roots\n')
        end
        for i = 1 : maxIter
            xh = x1 - y(x1)/y1(x1);
            errPercent = abs(((x1-xh)/xh)*100);
            if errPercent > relApproxErr
                x1 = xh;
                Z(i) = errPercent;
            end
            if abs(y(xh)) < funcConverg/100
                flag = 2;
            break;
            end
            if errPercent < relApproxErr
                flag = 1;
            break;
            end
        end
        if flag == 2
            fprintf('Convergence criterion for function value is met\n');
        end
        if flag == 1
            fprintf('Convergence criterion for maximum relative error is met\n');
        end
        if flag == 0
            fprintf('Convergence criterion for maximum iteration is met\n');
        end
        fprintf('The root of the function: %f\n', xh);
        figure(1);
        ezplot(func);
        grid on
        figure(2);
        plot(Z);
    end 
    
    %Secant
    if strcmp(reply,'Secant')
        func = input('Enter a function: ','s');
        x1 = input('Enter x1: ');
        x2 = input('Enter x2: ');
        relApproxErr = input('Enter relative error limit: ');
        funcConverg = input('Enter function value limit: ');
        maxIter = input('Enter iteration number limit: ');
        flag = 0;
        f = matlabFunction(str2sym(func));
        Z = zeros;
        if y(x1) == 0
            fprintf('x1 is one of the roots\n')
        end
        if y(x2) == 0
            fprintf('x2 is one of the roots\n')
        end
        for i = 1 : maxIter
            xh = x2 - y(x2)*((x2 - x1)/(y(x2) - y(x1)));
            errPercent = abs(((x2-xh)/xh)*100);
            if errPercent > relApproxErr
                x1 = x2;
                x2 = xh;
                Z(i) = errPercent;
            end
            if abs(y(xh)) < funcConverg/100
                flag = 2;
            break;
            end
            if errPercent < relApproxErr
                flag = 1;
            break;
            end
        end
        if flag == 2
            fprintf('Convergence criterion for function value is met\n');
        end
        if flag == 1
            fprintf('Convergence criterion for maximum relative error is met\n');
        end
        if flag == 0
            fprintf('Convergence criterion for maximum iteration is met\n');
        end
        fprintf('The root of the function: %f\n', xh);
        figure(1);
        ezplot(func);
        grid on
        figure(2);
        plot(Z);
    end 
end
