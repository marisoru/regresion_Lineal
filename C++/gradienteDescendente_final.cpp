#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <iostream>
#ifdef WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif // win32

using namespace std;

double max_element(vector<double> toCompare);
vector <double> linspace(double a, double b, int N) ;
void sleepcp(int milliseconds);
double h_theta(double theta0, double theta1, double x);

int main() {
	//Declaración de varibles
	double xNum, yNum, alpha, hNum, temp0, temp1;
	string line;
	vector<double> xValues;
	vector<double> yValues;
	vector<double> theta_0_ToPrint;
	vector<double> theta_1_ToPrint;
	vector<double> J_valsToPrint;
	vector<double> f_theta;
	int m;
	int max_iter = 1000;
	int verbose = 0;
	int display = 0;
        
	//Archivos que contedrá los puntos x e y
 	FILE *epoints;
    epoints = fopen( "epoints.txt" , "w" );
	//Archivo que contendrá la información para pasarla a GNUplot
    FILE *toPlot;
    toPlot = fopen( "toPlot.dat" , "w" );
    //Archivo que contendrá la información para pasarla a GNUplot
    FILE *lineToPlot;
    lineToPlot = fopen( "lineToPlot.txt" , "w" );
    //Archivo que contendrá la información para pasarla a GNUplot
    FILE *jpoints;
    jpoints = fopen( "jpoints.txt" , "w" );
    //Archivo que contendrá la información para pasarla a GNUplot
    FILE *toPlotJ;
    toPlotJ = fopen( "toPlotJ.dat" , "w" );
	//Archivo para crear un gif
	FILE *archivo;                       /// inicialización de variable para guardar en
	archivo = fopen("creategif.dat","w");                    /// crea o abre el archivo correspondiente
	//Archivo para crear un gif
	FILE *archivo1;                       /// inicialización de variable para guardar en
	archivo1 = fopen("creategif1.dat","w");                    /// crea o abre el archivo correspondiente	

	//Obtener datos en ex2x en un vector para usar 
	ifstream infileX("ex2x.dat");
	while (getline(infileX, line))
       {
       istringstream issx(line);
       if (issx >>xNum)
           xValues.push_back(xNum);
    	}
	
	//Obtener datos en ex2x en un vector para usar 
	ifstream infileY("ex2y.dat");
	while (getline(infileY, line))
    	{
       istringstream issy(line);
       if (issy >>yNum)
           yValues.push_back(yNum);
    	}
		
	//Total de número de muestras. Equivale a la cantidad de datos suministrados
	m = yValues.size();
	
	//Vector de 1's
	vector<double> onesVector (yValues.size(),1);
	
	//El vector theta contiene los valores de Theta_0 y Theta_1
	//Se inician en 0's
	vector<double> theta (2,0);
	
	//Valor de alpha asignado
	//alpha = 0.07;
	//Valor de alpha ingresado por el usuario
	cout << "Ingrese el valor de alpha"<<endl;
	cin >> alpha;
	cout << endl;

	//El vector grad contiene el valor del gradiente que será calculado.
	vector<double> grad (theta.size(),1);
	
	//vetor h de 50 ceros
	vector<double> h (xValues.size(),0);
	//vetor error de 50 ceros
	vector<double> error (xValues.size(),0);
	
	//Mensaje en panatalla para usuario
	cout << "Desea imprimir uno a uno los valores Theta?" << endl;
	cout << "   "<<"0 - NO " << endl;
	cout << "   "<<"1 - SI " << endl;
	cin >> verbose;
	
	//while ( fabs(max_element(grad)) > 0.00001){
	while (max_iter != 0){
		for (int i = 0; i < xValues.size(); i++){
			h.at(i) = onesVector.at(i) * theta.at(0) + xValues.at(i) * theta.at(1);
			error.at(i) = h.at(i) - yValues.at(i);
			//Primera componente del gradiente
			temp0 = onesVector.at(i) * error.at(i);	
			grad.at(0) = grad.at(0) + temp0;
			//Segunda componente del gradiente
			temp1 = xValues.at(i) * error.at(i);	
			grad.at(1) = grad.at(1) + temp1;
		}
		
		grad.at(0) = grad.at(0)/m;
		grad.at(1) = grad.at(1)/m;
		
		//Valores de theta
		theta.at(0) = theta.at(0) - alpha * grad.at(0);
		theta.at(1) = theta.at(1) - alpha * grad.at(1);
		
		theta_0_ToPrint.push_back(theta.at(0));
		theta_1_ToPrint.push_back(theta.at(1));
		
		max_iter--;
	}
	
	
///////////////////////////////////////////////////////////////////////////////

   //Matriz inicial de 0's
   double J_vals[100][100];
   for (int i = 0 ; i < 100; i++){
   	for(int j = 0; j < 100; j++){
   		J_vals[i][j] = 0;
	   }
   }
   
   //Valores alrededor de los cuales será calculado J
   //Valores de Theta_0 desde -3 a 3
   vector<double> Theta_0;
   Theta_0 = linspace (-3, 3, 100);
   //Valores de Theta_1 Desde -1 a 1
   vector<double> Theta_1;
   Theta_1 = linspace (-1, 1, 100);
    
    //Calcular los valores de la función de costo
    vector <double> t (2,0); 
    
    vector<double> h_y;    
    double sum_hy = 0;
    vector<double> hfinal;    
	
    for ( int i = 0; i < Theta_0.size(); i++ ){
    	for ( int j = 0; j < Theta_1.size(); j++){
    		t.at(0) = Theta_0.at(i);
    		t.at(1) = Theta_1.at(j);
    		
    		hfinal = h;
			h.clear();
            for (int s = 0; s < xValues.size(); s++){
    			hNum = onesVector.at(s) * t.at(0) + xValues.at(s) * t.at(1); 
				h.push_back(hNum);
			}
			
			sum_hy = 0;
			h_y.clear();
			for (int k = 0 ; k < yValues.size(); k++){
  				h_y.push_back( h.at(k) - yValues.at(k)) ;
  				h_y.at(k) = pow( h_y.at(k), 2 );
   				sum_hy = sum_hy + h_y.at(k);
			}
			
			//J = (1/2M)*(Theta*X-Y)'*(Theta*X-Y);
			J_vals[i][j] = sum_hy / ( 2 * m );
			J_valsToPrint.push_back(J_vals[i][j]);
            
		}
	}
	
	if ( verbose == 0 ){
		cout << "__________________________________" << endl;
		cout << "______Resultados finales__________" << endl;
		cout << endl;
		cout << "Theta_0 final = " << theta_0_ToPrint.back() << endl;
		cout << "Theta_1 final = " << theta_1_ToPrint.back() << endl;
		cout << "Valor final funcion de costo J = " << J_valsToPrint.back() << endl;
                cout << "__________________________________" << endl;
	}else{
		cout << "__________________________________" << endl;
		cout << "____Resultados por iteracio'n_____" << endl;
		sleepcp(250);
		for (unsigned ii=0; ii<theta_0_ToPrint.size(); ii++){
			cout << "__________________________________" << endl;
			cout << endl;
			cout << "Theta_0 = " << theta_0_ToPrint.at(ii) << endl;
			cout << "Theta_1 = " << theta_1_ToPrint.at(ii) << endl;
			cout << "Valor funcion de costo J = " << J_valsToPrint.at(ii) << endl;
			sleepcp(250);
		}
	}
	
	
	//llenar el archivo epoints con los puntos de ex y de ey
	for (int i = 0 ; i < yValues.size(); i ++){
          fprintf(epoints, "%f %f\n", xValues.at(i), yValues.at(i)) ;
	}
	fclose(epoints);
	
	//Para graficar J	
 	//llenar el archivo JToPlot con los puntos de ex y de ey
	for (int i = 0 ; i < yValues.size(); i ++){
       fprintf(jpoints, "%f %f %f\n", theta_0_ToPrint.at(i), theta_1_ToPrint.at(i),J_valsToPrint.at(i)) ;
	}
	fclose(jpoints);
	
	//Mensaje en pantalla para el usuario
	cout << "Desea ver una a una las rectas obtenidas?" << endl;
	cout << "   "<<"0 - NO " << endl;
	cout << "   "<<"1 - SI " << endl;
	cin >> display;

	//Si display es 0, solo se mostrarán los puntos y la recta final
	if(display == 0){
	   
	   fprintf(toPlot, "set title 'GRADIENTE DESCENDENTE'\n");
	   fprintf(toPlot, "set xlabel 'Edad en anos'\n");
	   fprintf(toPlot, "set ylabel 'Altura en metros'\n");
   	   fprintf(toPlot, "plot \"epoints.txt\"\n");
       fprintf(toPlot, "pause 1\n");
	   fprintf(toPlot, "replot %f + %f*x\n",theta_0_ToPrint.back(), theta_1_ToPrint.back());    
       fprintf(toPlot, "pause 5\n");
	   fclose(toPlot);
	   system("gnuplot \"toPlot.dat\"");

       fprintf(archivo,"reset\n");                              /// escribe en el archivo el texto correspondiente
   	   fprintf(archivo,"set term gif animate\n");               /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"set terminal gif animate delay 50\n");  /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"set output \"final.gif\"\n");         /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"n = 20\ni = 0\n");                      /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"load \"toPlot.dat\"\n");               /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"set output\n");                         /// escribe en el archivo el texto correspondiente
       fclose(archivo);                                         /// cierre del archivo de salida
 	   system("gnuplot \"creategif.dat\"");  /// guarda la grafica en un archivo .gif
   	  
       //Escribir el archivo para ver J
	   fprintf(toPlotJ, "set title 'Grafica de J'\n");
	   fprintf(toPlotJ, "set xlabel 'Theta_0'\n");
	   fprintf(toPlotJ, "set ylabel 'Theta_1'\n");
	   fprintf(toPlotJ, "set dgrid 30,30\n");
	   fprintf(toPlotJ, "set hidden3d\n");
   	   fprintf(toPlotJ, "splot \"jpoints.txt\" u 1:2:3 with lines\n");
       fprintf(toPlotJ, "pause 5\n");
	   fclose(toPlotJ);
	   system("gnuplot \"toPlotJ.dat\"");

       fprintf(archivo,"reset\n");                              /// escribe en el archivo el texto correspondiente
	   fprintf(archivo,"set term gif animate\n");               /// escribe en el archivo el texto correspondiente
       fprintf(archivo1,"set terminal gif animate delay 50\n");  /// escribe en el archivo el texto correspondiente
       fprintf(archivo1,"set output \"jota.gif\"\n");         /// escribe en el archivo el texto correspondiente
       fprintf(archivo1,"n = 20\ni = 0\n");                      /// escribe en el archivo el texto correspondiente
       fprintf(archivo1,"load \"toPlotJ.dat\"\n");               /// escribe en el archivo el texto correspondiente
       fprintf(archivo1,"set output\n");                         /// escribe en el archivo el texto correspondiente
       fclose(archivo1);                                         /// cierre del archivo de salida
 	   system("gnuplot \"creategif1.dat\"");  /// guarda la grafica en un archivo .gif

	   
	}else{
	   fprintf(toPlot, "set yrange [0:3.5]\n\n");
	   fprintf(toPlot, "set title 'GRADIENTE DESCENDENTE'\n");
	   fprintf(toPlot, "set xlabel 'Edad en anos'\n");
	   fprintf(toPlot, "set ylabel 'Altura en metros'\n");
   	   for (unsigned ii=0; ii<theta_0_ToPrint.size(); ii++){
              fprintf(toPlot, "plot \"epoints.txt\"\n");
       	      //fprintf(toPlot, "pause 1\n");
	      fprintf(toPlot, "replot %f + %f*x\n",theta_0_ToPrint.at(ii), theta_1_ToPrint.at(ii));    
              fprintf(toPlot, "pause 0.2\n");
	      fprintf(toPlot, "refresh\n\n");
 	   }
	   fclose(toPlot);
	   system("gnuplot \"toPlot.dat\"");

	   //GIF para ver todas las iteraciones
	   fprintf(archivo,"reset\n");                              /// escribe en el archivo el texto correspondiente
	   fprintf(archivo,"set term gif animate\n");               /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"set terminal gif animate delay 50\n");  /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"set output \"itearciones.gif\"\n");         /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"n = 20\ni = 0\n");                      /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"load \"toPlot.dat\"\n");               /// escribe en el archivo el texto correspondiente
       fprintf(archivo,"set output\n");                         /// escribe en el archivo el texto correspondiente
       fclose(archivo);                                         /// cierre del archivo de salida
 	   system("gnuplot \"creategif.dat\"");  /// guarda la grafica en un archivo .gif
	}	
 	
}


double max_element(vector<double> toCompare){
	double max;
	if ( toCompare.at(0) == toCompare.at(1) )
		max = toCompare.at(0);
	else{
		if ( toCompare.at(0) > toCompare.at(1) )
		   max = toCompare.at(0);
		else
		   max = toCompare.at(1);
	}
	return max;
}

vector <double> linspace(double a, double b, int N) {
    double h = 0;
	h = (b - a) / static_cast<double>(N-1);
    vector<double> xs(N);
    typename std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

void sleepcp(int milliseconds) // Cross-platform sleep function
{
    #ifdef WIN32
        Sleep(milliseconds);
    #else
        usleep(milliseconds * 1000);
    #endif // win32
}

double h_theta(double theta0, double theta1, double x){
   return theta0 + theta1*x;
}
