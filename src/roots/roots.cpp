#include <iostream>
#include <cmath>
#include "roots.hpp"

const double epsilon = 1e-6;
const int max_iteration = 1e6;

bool bisection(std::function<double(double)> f, double a, double b, double *root) {

    if (std::abs(f(a)) < epsilon){ // If a is already a root, returns true
        *root = a;
        return true;
    }
    if (std::abs(f(b)) < epsilon){ // If b is already a root, returns true
        *root = b;
        return true;
    }
    if (f(a) * f(b) > 0){
        return false;
    }
    for(int i = 0; i < max_iteration; i++){
        *root = (a + b) / 2.0;
        if (std::abs(f(*root)) < epsilon || (b-a) / 2.0 < epsilon){ // checks if interval is close enough to zero
            return true;
        }

        if (f(a) * f(*root) < 0){
            b = *root;
        }
        else a = *root;
    }
    return false;
}

bool regula_falsi(std::function<double(double)> f, double a, double b, double *root){

    if (std::abs(f(a)) < epsilon){ // If a is already a root, returns true
        *root = a;
        return true;
    }
    if (std::abs(f(b)) < epsilon){ // If b is already a root, returns true
        *root = b;
        return true;
    }

    if (f(a) * f(b) > 0){
        return false;
    }

    for (int i = 0; i < max_iteration; i++){
        *root = a - ((f(a) * (b-a)) / (f(b) - (f(a))));

        if (std::abs(f(*root)) < epsilon){
            return true;
        }
        if (f(a) * f(*root) < 0){
            b = *root;
        }
        else a = *root;
    }
    return false;
}

bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root){
    double x = c;
    for (int i = 0; i < max_iteration; i++){
        double df = g(x);
        if (std::abs(df) < 1e-12){ // returns false if derivative is zero
            return false;
        }
        x = x - f(x) / df;
        if (x < a || x > b){ // returns false if not in interval
            return false;
        }
        if (std::abs(f(x)) < epsilon){
            *root = x;
            return true;
        }
    }
    return false;
}

bool secant(std::function<double(double)> f, double a, double b, double c, double *root){
    
    double x0 = a; // first point
    double x1 = c; // second point which is guess

    for (int i = 0; i < max_iteration; i++){
        double fx1 = f(x1);
        double fx0 = f(x0);

        if (std::abs(fx1 - fx0) < 1e-12){ // prevents divison by zero
            return false;
        }

        double x_next = x1 - (fx1 * (x1 - x0) / (fx1 - fx0));
        
        if (std::abs(f(x_next)) < epsilon){
            *root = x_next;
            return true;
        }
        x0 = x1;
        x1 = x_next;
    }
    return false;
}

