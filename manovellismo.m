clc
close all
clear all

%% Inserimento dati iniziali

corretto = 1;
while corretto == 1
    a = input('Raggio manovella in [m]: ');
    b = input('Lunghezza biella in [m]: ');

    if a < b
        corretto = 0;
    else
        disp('La manovella non può essere più lunga della biella.')
    end
end

% Definizione della discretizzazione dell'analisi

alfa_d = 0:0.1:36000;                                                       
alfa = alfa_d * (pi / 180);                                               

alfap = input('Assegna velocità angolare manovella [rad/s]: ');  
alfapp = input('Assegna accelerazione angolare manovella [rad/s^2]: ');                                                                

%% Rappresentazione della configurazione iniziale

beta_a = asin(-(a * sin(alfa(451))) / b);
c_a = a * cos(alfa(451)) + b * cos(beta_a);
gamma = 0;
OA = a * exp(1i * alfa(451));
AB = b * exp(1i * beta_a);
OB = c_a * exp(1i * gamma);
figure()

% Rappresentazione delle cerniere
plot(0, 0, 'ko')
xlabel('x [m]')
ylabel('y [m]')
title('Configurazione a \alpha = 45°')
hold on

% Posizioni di O, A e B
O = 0;
A = OA;
B = OB;

% Rappresentazione delle posizioni O, A e B
plot(real(O), imag(O), 'ko')
plot(real(A), imag(A), 'ko')
plot(real(B), imag(B), 'ko')

% Etichette
offset = 0.05;
text(real(O), imag(O) + offset, 'O', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
text(real(A), imag(A) + offset, 'A', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
text(real(B) + offset, imag(B) + offset, 'B', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')

% Rappresentazione delle aste
plot([0 real(A)], [0 imag(A)], 'r-', 'LineWidth', 3)
plot([real(A) real(A + AB)], [imag(A) imag(A + AB)], 'g-', 'LineWidth', 3)
plot([0 real(B)], [0 imag(B)], 'b-', 'LineWidth', 3)

% Estensione degli assi
offset_extend = 0.1;
xlim([min([0 real(O) real(A) real(A + AB) real(B)]) - offset_extend, max([0 real(O) real(A) real(A + AB) real(B)]) + offset_extend])
ylim([min([0 imag(O) imag(A) imag(A + AB) imag(B)]) - offset_extend, max([0 imag(O) imag(A) imag(A + AB) imag(B)]) + offset_extend])

axis equal
grid on

%% Analisi della posizione

beta = asin(-(a * sin(alfa)) / b);
c = a * cos(alfa) + b * cos(beta);

%% Analisi della velocità

for i = 1:length(alfa)
    A = [1, b * sin(beta(i));
         0, -b * cos(beta(i))];

    B = [-alfap * a * sin(alfa(i));
          alfap * a * cos(alfa(i))];

    xp = inv(A) * B;

    cp(i) = xp(1);
    betap(i) = xp(2);
end

%% Analisi dell'accelerazione

for i = 1:length(alfa)
    A = [1, b * sin(beta(i));
         0, -b * cos(beta(i))];

    B = [-alfapp * a * sin(alfa(i)) - alfap ^ 2 * a * cos(alfa(i)) - betap(i) ^ 2 * b * cos(beta(i));
          alfapp * a * cos(alfa(i)) - alfap ^ 2 * a * sin(alfa(i)) - betap(i) ^ 2 * b * sin(beta(i))];

    xpp = inv(A) * B;

    cpp(i) = xpp(1);
    betapp(i) = xpp(2);
end

%% Visualizzazione risultati

visualizza = 1;
while visualizza == 1
    k = menu('Visualizzazione risultati', 'Grafici posizione', ...
             'Grafici velocità', 'Grafici accelerazione', ...
             'Animazione del moto', 'Fine');

    % Grafico posizione
    if k == 1
        figure(1)
        subplot(211)
        plot(alfa * (180 / pi), c)
        xlabel('Rotazione manovella [°]')
        ylabel('Posizione piede di biella [m]')
        ax = axis;
        axis([0, 360, ax(3), ax(4)])
        grid on

        subplot(212)
        plot(alfa * (180 / pi), beta * (180 / pi))
        xlabel('Rotazione manovella [°]')
        ylabel('Rotazione biella [rad]')
        ax = axis;
        axis([0, 360, ax(3), ax(4)])
        grid on

    % Grafico velocità
    elseif k == 2
        figure(2)
        subplot(211)
        plot(alfa * (180 / pi), cp)
        xlabel('Rotazione manovella [°]')
        ylabel('Velocità piede di biella [m/s]')
        ax = axis;
        axis([0, 360, ax(3), ax(4)])
        grid on 

        subplot(212)
        plot(alfa * (180 / pi), betap)
        xlabel('Rotazione manovella [°]')
        ylabel('Velocità angolare biella [rad/s]')
        ax = axis;
        axis([0, 360, ax(3), ax(4)])
        grid on

    % Grafico accelerazione
    elseif k == 3
        figure(3)
        subplot(211)
        plot(alfa * (180 / pi), cpp)
        xlabel('Rotazione manovella [°]')
        ylabel('Accelerazione piede di biella [m/s^2]')
        ax = axis;
        axis([0, 360, ax(3), ax(4)])
        grid on

        subplot(212)
        plot(alfa * (180 / pi), betapp)
        xlabel('Rotazione manovella [°]')
        ylabel('Accelerazione angolare biella [rad/s^2]')
        ax = axis;
        axis([0, 360, ax(3), ax(4)])
        grid on

    % Animazione del moto
    elseif k == 4
        figure(4)
        set(gcf, 'DoubleBuffer', 'on')

        ngiri = 3;
        nxframe = 2;
        nframe = ngiri * 360 / nxframe;
        
        for i = 1:nframe
            igrado = nxframe * i;
            n = floor(igrado / 360);
            igrado = igrado - 360 * n;

            ic = find(alfa * 180 / pi >= igrado);

            x_manov = a * cos(alfa(ic(1)));
            y_manov = a * sin(alfa(ic(1)));

            x_biell = c(ic(1));
            y_biell = 0;

            plot([0, x_manov], [0, y_manov], 'r', ...
                 [x_manov, x_biell], [y_manov, y_biell], 'b',...
                 [0, x_manov, x_biell], [0, y_manov, y_biell], 'ok')
                 
            axis([-a - b, a + b, -a - b, a + b])
            axis equal;
            title(sprintf('Ciclo: %d', n + 1))
            xlabel('[m]')
            ylabel('[m]') 
            xlim([-a - 0.1, a + b + 1])
            ylim([-a - 0.1, a + 0.1])
            grid on
            hold on;
            viscircles([x_biell, y_biell], 0.5, 'Color', 'g', 'LineWidth', 0.5);
            hold off;
            drawnow
        end
    else
        visualizza = 0;
        close all
    end
end
