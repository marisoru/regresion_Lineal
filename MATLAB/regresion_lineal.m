function varargout = regresion_lineal(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @regresion_lineal_OpeningFcn, ...
                   'gui_OutputFcn',  @regresion_lineal_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before regresion_lineal is made visible.
function regresion_lineal_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes regresion_lineal wait for user response (see UIRESUME)
% uiwait(handles.figure1);
axes(handles.axes7);
img=imread('logocitedi.png');
axis off;
imshow(img);


% --- Outputs from this function are returned to the command line.
function varargout = regresion_lineal_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on selection change in display.
function display_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function display_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
%%Se definen las variables de los selectores.

popup_sel_index = get(handles.display, 'Value');
popup_sel_index2=get(handles.verbose,'Value');
clc
%% Cargar los datos
% x-values corresponden a edades de niños
% y-values corresponden a la altura de los niños en metros
x = load('ex2x.dat');
y = load('ex2y.dat');

%% Graficar los x vs y
axes(handles.axes2);
cla;
plot(x, y, 'o');
ylabel('Altura en metros')
xlabel('Edad en años')

%% Total de número de muestras
% Equivale a la cantidad de datos sumistrados (50)
m = length(y); 
% Agregar una columna de '1's a x
x = [ones(m, 1), x]; 

%% Inicializar parámetros
theta = zeros(size(x,2),1);
alpha = 0.07;
grad = ones(size(theta));
theta0=0;
theta1=0;

%% Iniciar iteraciones 
% 0.00001 Como criterio de toleracia para detener las iteraciones
while abs(max(grad(:))) > 0.01
    h = sum(x * theta,2);
    error = h - y;
    grad = x' * error / m;
    if popup_sel_index2==1
    theta0=theta(1,1)
    theta1=theta(2,1)
    theta = theta - alpha * grad;
    end
    if popup_sel_index2==2
    theta0=theta(1,1);
    theta1=theta(2,1);
    theta = theta - alpha * grad;
    end
  if popup_sel_index==1
  
    hold on 
plot(x(:,2), x*theta, '-')
legend('Datos de entrada', 'Regrensión Lineal')
drawnow;


    end

end
theta0
theta1

    
 %% Graficar los nuevos datos
if popup_sel_index==2
hold on 
plot(x(:,2), x*theta, '-')
legend('Datos de entrada', 'Regrensión Lineal')
end

%% Valores de J
% Matriz inicial de 100x100 de 0's
J_vals = zeros(100, 100);
% Valores alrededor de los cuales será calculado J
% Valores de Theta_0 Desde -3 a 3 
theta0_vals = linspace(-3, 3, 100);
% Valores de Theta_1 Desde -1 a 1
theta1_vals = linspace(-1, 1, 100);

%% Calcular valores de J - Fución de costo
if popup_sel_index2==1
for i = 1:length(theta0_vals)
  for j = 1:length(theta1_vals)
    t = [theta0_vals(i); theta1_vals(j)];
    h = sum(x * t);
    % J = (1/2M)*(Theta*X-Y)'*(Theta*X-Y);
    
    J_vals(i,j) = (0.5/m) .* (x * t - y)' * (x * t - y);
    
   
  end
end
J_vals
end
if popup_sel_index2==2
for i = 1:length(theta0_vals)
  for j = 1:length(theta1_vals)
    t = [theta0_vals(i); theta1_vals(j)];
    h = sum(x * t);
    % J = (1/2M)*(Theta*X-Y)'*(Theta*X-Y);
    
    J_vals(i,j) = (0.5/m) .* (x * t - y)' * (x * t - y);
    
   
  end
end
J_vals(100,100)
end

%% Graficar J 
J_vals = J_vals';
axes(handles.axes6);
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1')
% --- Executes on selection change in verbose.
function verbose_Callback(hObject, eventdata, handles)

function verbose_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
