#include<iostream> 
#include<cmath>
#include<fstream>
using namespace std;

static double a = 0.2, k = 0.8*3.1415;
double f(double x) { return x - a*sin(k*x); }

double derivative(  double x, 
                    double function(double),
                    double dx = pow(10,-12)){
    return ( f(x+dx)-f(x) )/dx;
}

double newton(  double x, 
                double function(double), 
                double c,
                double epsilon = pow(10,-12)) {

    int i = 0;
    double xn = x - (function(x) - c)/derivative(x,function);
    while ( abs(xn-x) > epsilon ) {
        x = xn;
        xn = x - (function(x) - c)/derivative(x,function);
        i++;
    }
    return xn;
}

int main() {
    ofstream out;
    out.open("./newton.out");
    for (double i = 0; i<4; i+=0.01){
        out << i
            << "\t"
            << newton(i,f,i)    
            << endl;
    }

    return 0;
}
