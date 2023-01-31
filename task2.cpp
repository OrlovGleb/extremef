#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/array.hpp>
#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double , 4> state_type;

state_type newton = {{-1000.0, -1000.0, -1000.0, -1000.0}};
vector <double> u;
vector <double> x2_0;
vector <double> t_int;

static double alpha = 0.1; 

void func(const state_type &x , state_type &dxdt , const double t) { // система ODE
/*x1'*/    dxdt[0] = x[1]; //x[1]=x2
/*x2'*/    dxdt[1] = x[3]  - x[0] * exp(-alpha * t);//x[2]=p1 x[0]=x1
/*p1'*/    dxdt[2] = x[3] * exp(-alpha * t);//x[3]=p2
/*p2'*/    dxdt[3] = -x[2];//x[2]=p1
}

void write_out(const state_type &x , const double t){
    const double temp = M_PI / 2.0;
    if (t == temp) {
        newton[0] = x[0];
        newton[1] = x[1];
        newton[2] = x[2];
        newton[3] = x[3];
    }
}

void final_out(const state_type &x , const double t){
    u.push_back((x[2] / 2.0) * (x[2] / 2.0));
    t_int.push_back(t);
    x2_0.push_back(x[2] / 2.0);
    if (t < 0.00001) {
        cout << "t = " << t << '\n' << "x1 = " << x[0] << '\n' << "p2 = " << x[3]  << endl;
    }
}


state_type calculate(double a1, double a2) {
    double b = M_PI / 2.0; //конец отрезка
    double t = 0.0; // начало отрезка
    double dt = 0.01; // шаг
    double delta = 1e-6;
    double F1, F2, F1x, F2x, F1y, F2y, det, f_norm;
    
    state_type F;
    state_type x3; // краевые условия для задачи
    
    typedef runge_kutta_dopri5 <state_type> dopri5_type;
    typedef controlled_runge_kutta <dopri5_type> controlled_dopri5_type;
    typedef dense_output_runge_kutta <controlled_dopri5_type> dense_output_dopri5_type;
    
    dense_output_dopri5_type dopri5 = make_dense_output( 1E-9 , 1E-9 , dopri5_type() );
         //x1(0),x2(0),p1(0),p2(0)
    x3 = {{ a1, 0, 0, a2}}; // краевые условия для задачи
    integrate_adaptive(dopri5 , func , x3 , t , b , dt, write_out);
    F1 = newton[0] ;
    F2 = newton[1] + M_PI/2;
          //x1(0),x2(0),p1(0),p2(0)
    x3 = {{a1 + delta, 0, 0, a2}}; // краевые условия для задачи
    integrate_adaptive(dopri5 , func , x3 , t , b , dt, write_out);
    F1x = (newton[0] - F1) / delta;
    F2x = (newton[1] + M_PI/2 - F2) / delta;
           //x1(0),x2(0),p1(0),p2(0)
    x3 = {{ a1, 0, 0, a2 + delta}}; // краевые условия для задачи
    integrate_adaptive(dopri5 , func , x3 , t , b , dt, write_out);
    F1y = (newton[0] - F1) / delta;
    F2y = (newton[1] + M_PI/2 - F2) / delta;
    
    det = F1x * F2y - F1y * F2x;
    f_norm = sqrt(F1 * F1 / (F1x * F1x + F1y * F1y) + F2 * F2 / (F2x * F2x + F2y * F2y));
    F = {{(F2y * F1 - F1y * F2) / det, (F1x * F2 - F2x * F1) / det, f_norm}};
    return F;
}

int controller(double a1, double a2, double eps, double al) {
    double b = M_PI / 2.0; //конец отрезка
    double t = 0.0; // начало отрезка
    double dt = 0.01; // шаг
    double integral = 0;
    state_type params, x; // краевые условия для задачи
    typedef runge_kutta_dopri5 <state_type> dopri5_type;
    typedef controlled_runge_kutta <dopri5_type> controlled_dopri5_type;
    typedef dense_output_runge_kutta <controlled_dopri5_type> dense_output_dopri5_type;
    dense_output_dopri5_type dopri5 = make_dense_output( 1E-9 , 1E-9 , dopri5_type() );

    alpha = al;
    for (int i = 1; i < 200; i++) {
        params = calculate(a1, a2);
        
        if(params[2] < eps) {
            break;
        }
        else {
            a1 -= params[0];
            a2 -= params[1];
        }
    }
    
    cout << "alpha = " << alpha << endl;
    cout << "error = " << params[2] << endl;
         //x1(0),x2(0),p1(0),p2(0)
    x = {{a1, 0, 0, a2}}; // краевые условия для задачи
    x2_0.clear();
    u.clear();
    t_int.clear();
    integrate_adaptive(dopri5 , func , x , t , b , dt, final_out);
    for (int i = 0; i < u.size() - 1; i++) {
        integral += (u[i + 1] + u[i]) * (t_int[i + 1] - t_int[i]) / 2;
    }
    cout << "Value = " << integral << endl;
    cout << endl;
    return alpha;
}

int main(int argc, const char * argv[]) {
    cout.precision(8);
    
 /*x1(0)*/   double a1 = M_PI/2; //начальное условие
 /*p2(0)*/   double a2 = 0.0; //начальное условие
    double eps = 1e-9;
    boost::array <double, 4> alphas = {{0.0, 0.1,  1.5, 10.0}};
    for (auto al: alphas) {
        controller(a1, a2, eps, al);
    }
}