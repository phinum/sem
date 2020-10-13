a = 8.779; %(L^2bar/mol^2)
b = 0.08445; %(L/mol)
R = 0.08314; %((L*bar)/(K*mol)) 
P = 0;
A = 4.53678;
B = 1149.36;
C = 24.906;
x = zeros; %different values of pressure
y = zeros; %different values of temperature
for T = 200 : 1 : 400
    x(T - 199) = 10^(A - (B / (T + C))); % Antoine's Equation
    y(T -199) = T;
end
r = zeros(201, 3);
for i = 1 : 201
    coeff = [x(i), -(x(i) * b + R * y(i)), a, -a*b]; % Van der Waals equation in coefficient form
    r(i, :) = roots(coeff);
end
s = zeros(201,3);
for i = 1 : 201 % check for real data
    c = isreal(r(i,2));
    if c == 1
        s(i,:) = r(i,:);
    end
end
t = 201;
while t > 0 % eliminate imaginary roots and clean data
    if s(t, 1) == 0
        s(t, :) = [];
        x(t) = [];
        y(t) = [];
    end
    t = t - 1;
end
m = zeros(140,2);
for i = 1 : 140 % initiale min and max root
    m(i,2) = max(r(i,:));
    m(i,1) = min(r(i,:));
end
Tc = y(140) + 1;  % critical temperature form the graph % error calculation
TcT = (8 * a) / (27 * R * b); % Theoretical critical temperature
error = (TcT - Tc)/TcT;
disp(error);
for i = 1 : 15 : 140 % plot for (T < Tc)
    fplot(@(V) ((R * y(i)) / (V-b)) - (a/(V^2)), 'blue');
    title('PV plot');
    xlabel('Volume');
    ylabel('Pressure');
    line([m(i,1) m(i,2)], [x(i) x(i)], 'color', 'black');
    hold on;
    plot(m,x,'k--');
    hold on;
    grid on;
end
for i = 370 % plot for (T = Tc)
    fplot(@(V) ((R * i) / (V-b)) - (a/(V^2)), 'green');
    title('PV plot');
    xlabel('Volume');
    ylabel('Pressure');
    hold on;
    plot(m,x,'k--');
    hold on;
    grid on;
end
for i = 385 : 15 : 490 % plot for (T > Tc)
    fplot(@(V) ((R * i) / (V-b)) - (a/(V^2)), 'red');
    hold on;
    title('PV plot');
    xlabel('Volume');
    ylabel('Pressure');
    hold on;
    plot(m,x,'k--');
    hold on;
    grid on;
end
