#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <functional>
#include <stdexcept>
#include <utility>
#include <cstdlib>
#include <ctime>  


using namespace std;
int func_eval = 0;
double R = 0.1;
int C = 10;
double pi = M_PI;

// //1st function
// // actual objective function
// double obj_fun(const vector<double>x){
//     double sum=0;
//     sum+=pow(x[0]-10,3)+pow(x[1]-20,3);
//     return sum;
// }
// // penalized function needed to use for the conjugate method
// double objective_function(const vector<double>& x) {
//     func_eval++;
//     double sum=0;
//     double sum1=0,sum2=0;
//     double g1x=pow(x[0]-5,2)+pow(x[1]-5,2)-100;
//     double g2x=82.81-pow(x[0]-6,2)-pow(x[1]-5,2);
//     if(g1x<0){
//         sum1=R*pow(g1x,2);
//     }else{sum1=0;}
//     if(g2x<0){
//         sum2=R*pow(g2x,2);
//     }else{sum2=0;}
//     sum+=obj_fun(x)+sum1+sum2;
//     return sum;
// }


//2nd function
// actual objective function
double obj_fun(const vector<double>x){
    double sum=0;
    sum+=-(pow(sin(2*pi*x[0]),3)*sin(2*pi*x[1]))/(pow(x[0],3)*(x[0]+x[1]));
    return sum;
}
// penalized function needed to use for the conjugate method
double objective_function(const vector<double>& x) {
    func_eval++;
    double sum=0;
    double sum1=0,sum2=0;
    double g1x=x[1]-pow(x[0],2)-1;
    double g2x=x[0]-pow(x[1]-4,2)-1;
    if(g1x<0){
        sum1=R*pow(g1x,2);
    }else{sum1=0;}
    if(g2x<0){
        sum2=R*pow(g2x,2);
    }else{sum2=0;}
    sum+=obj_fun(x)+sum1+sum2;
    return sum;
}


// //3rd function
// // actual objective function
// double obj_fun(const vector<double>x){
//     double sum=0;
//     sum+=x[0]+x[1]+x[2];
//     return sum;
// }
// // penalized function needed to use for the conjugate method
// double objective_function(const vector<double>& x) {
//     func_eval++;
//     double sum=0;
//     double sum1=0,sum2=0,sum3=0,sum4=0,sum5=0,sum6=0;
//     double g1x=1-0.0025*(x[3]+x[5]);
//     double g2x=1-0.0025*(-x[3]+x[4]+x[6]);
//     double g3x=1-0.01*(-x[5]+x[7]);
//     double g4x=x[0]*x[5]-100*x[0]-833.33252*x[3]+83333.333;
//     double g5x=x[1]*x[6]-x[1]*x[3]+1250*x[3]+1250*x[4];
//     double g6x=x[2]*x[7]-x[2]*x[4]+2500*x[4]-1250000;
//     if(g1x<0){
//         sum1=R*pow(g1x,2);
//     }else{sum1=0;}
//     if(g2x<0){
//         sum2=R*pow(g2x,2);
//     }else{sum2=0;}
//     if(g3x<0){
//         sum3=R*pow(g3x,2);
//     }else{sum3=0;}
//     if(g4x<0){
//         sum4=R*pow(g4x,2);
//     }else{sum4=0;}
//     if(g5x<0){
//         sum5=R*pow(g5x,2);
//     }else{sum5=0;}
//     if(g6x<0){
//         sum6=R*pow(g6x,2);
//     }else{sum6=0;}
//     sum+=obj_fun(x)+sum1+sum2+sum3+sum4+sum5+sum6;
//     return sum;
// }


vector<double> gradient(const vector<double>& x) {
    int n = x.size();
    vector<double> grad(n,1.0);
    double h = 0.00001;
    for (int i = 0; i < n; ++i) {
        vector<double> a_up = x;
        a_up[i] += h;
        vector<double> a_down = x;
        a_down[i] -= h;
        grad[i] = (objective_function(a_up) - objective_function(a_down)) / (2 * h);
    }
    return grad;
}

vector<double> normalize(const vector<double>& v) {
    int n = v.size();
    double norm = sqrt(inner_product(v.begin(), v.end(), v.begin(), 0.0));
    if (norm < 1e-10) {  
        throw runtime_error("Normalization error: Vector norm is too small.");
    }
    vector<double> v_normalized(n, 0.0);
    for (int i = 0; i < n; ++i) {
        v_normalized[i] = v[i] / norm;
    }
    return v_normalized;
}

bool check_linear_independence(const vector<double>& s_k, const vector<double>& s_k_1) {
    double y = inner_product(s_k.begin(), s_k.end(), s_k_1.begin(), 0.0);
    return abs(y) < 1e-3;
}

pair<double, double> bounding_phase_method(function<double(double)> f, double a, double b, double delta, int max_iter = 1000000) {

    int input = 1;
    double x = 1.0;
    while(input == 1){
        if(max_iter <= 0){
            throw runtime_error("Bounding Phase Method failed: max_iter reached.");
        }
        
        x = a + (b - a) * ((double) rand() / RAND_MAX);

        double f_x_minus = f(x - fabs(delta));
        double f_x = f(x);
        double f_x_plus = f(x + fabs(delta));

        if (f_x_minus >= f_x && f_x >= f_x_plus) {
            delta = fabs(delta);
            break;
        } else if (f_x_minus <= f_x && f_x <= f_x_plus) {
            delta = -fabs(delta);
            break;
        } else {
            max_iter--;
        }
    }

    int k =0;
    while (k < max_iter) {
        double step = pow(2, k) * delta;
        double x_next = x + step;
        
        
        if(x_next < a || x_next > b){
            double x_high = min(x + pow(2, k-1) * delta, b);
            double x_low = max(x - pow(2, k-1) * delta, a);
            return {x_low, x_high};
        }

        double f_x_next = f(x_next);
        double f_x = f(x);
        if (f_x_next < f_x) {
            x = x_next;
            f_x = f_x_next;
            k++;
        } else {
            double x_high = x_next;
            double x_low = x - pow(2, k-1) * delta;
           
            x_low = max(x_low, a);
            x_high = min(x_high, b);
            return {x_low, x_high};
        }
    }

    throw runtime_error("Bounding Phase Method failed to find a suitable interval within max iterations.");
}

double golden_section_search(function<double(double)> f, double a, double b, double epsilon, int max_iter = 1000000) {
    const double golden_ratio = (1 + sqrt(5)) / 2;
    const double inv_golden_ratio = 1 / golden_ratio;

    double a_omega = 0;
    double b_omega = 1;
    double L_omega = 1;
    int k = 1;

    double omega1 = a_omega + (0.618) * L_omega;
    double omega2 = b_omega - (0.618) * L_omega;

    double x1 = a + omega1 * (b - a);
    double x2 = a + omega2 * (b - a);

    double f_x1 = f(x1);
    double f_x2 = f(x2);

    while (L_omega > epsilon && k < max_iter) {
        if (f_x1 > f_x2) {
            b_omega = omega1;
            L_omega = b_omega - a_omega;
            omega1 = a_omega + (0.618) * L_omega;
            omega2 = b_omega - (0.618) * L_omega;
            x1 = a + omega1 * (b - a);
            x2 = a + omega2 * (b - a);
            f_x1 = f(x1);
            f_x2 = f(x2);
        } else {
            a_omega = omega2;
            L_omega = b_omega - a_omega;
            omega1 = a_omega + (0.618) * L_omega;
            omega2 = b_omega - (0.618) * L_omega;
            x1 = a + omega1 * (b - a);
            x2 = a + omega2 * (b - a);
            f_x1 = f(x1);
            f_x2 = f(x2);
        }
        k++;
    }

    if (f_x1 < f_x2)
        return x1;
    else
        return x2;
}

double uni_search(function<double(double)> f, double delta, double epsilon, double a, double b) {
    pair<double, double> interval = bounding_phase_method(f, a, b, delta);
    return golden_section_search(f, interval.first, interval.second, epsilon);
}


void conjugate_gradient_optimization_variable(vector<double>& x, double epsilon1, double epsilon2, double epsilon3, double delta, vector<double> lo, vector<double> up) {
    int n = x.size();
    vector<double> grad = gradient(x);
    vector<double> s = normalize(grad);
    vector<double> s_y = grad;
    for(int i = 0; i < n; ++i) {
        s[i] = -s[i];
        s_y[i] = -s_y[i];
    }

    while (true) {
        
        auto f_y = [&x, &s](double y) -> double {
            int n = x.size();
            vector<double> x_new(n,0.0);
            for (int i = 0; i < n; ++i) {
                x_new[i] = x[i] + y * s[i];
            }
            return objective_function(x_new);
        };

        int n = x.size();
        for(int i=0;i<n;i++){
            x[i]=clamp(x[i],lo[i],up[i]);
        }

        // Finding optimal y range for unidirectional search
        double lower=0.0;
        double upper=0.0;
        if(abs(s_y[0])<1e-14){
            throw runtime_error("y Limiting error: Normalized direction is too small.");
        }
        double minab = (lo[0]-x[0])/s_y[0],maxab=(up[0]-x[0])/s_y[0];
        for(int i=1;i<n;++i){
            if(abs(s_y[i])<1e-14){
                throw runtime_error("y Limiting error: Normalized direction is too small.");
            }
            lower = (lo[i]-x[i])/s_y[i];
            upper = (up[i]-x[i])/s_y[i];
            if(lower>minab){
                minab = lower;
            }
            if(upper<maxab){
                maxab = upper;
            }
        }
        double c = minab;
        double d = maxab;
        


        double y;
        y = uni_search(f_y, delta, epsilon1, c, d);
        

        vector<double> x_new(n, 0.0);
        for (int i = 0; i < n; ++i) {
            x_new[i] = x[i] + y * s[i];
        }

        
        double x_diff_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = x_new[i] - x[i];
            x_diff_norm += diff * diff;
        }
        x_diff_norm = sqrt(x_diff_norm);

        double x_norm = sqrt(inner_product(x_new.begin(), x_new.end(), x_new.begin(), 0.0));

        if ((x_norm > 1e-10 && (x_diff_norm / x_norm <= epsilon2)) || (sqrt(inner_product(grad.begin(), grad.end(), grad.begin(), 0.0)) <= epsilon3)) {
            x = x_new;
            break;
        }

        x = x_new;

        
        vector<double> grad_new;
        grad_new = gradient(x);
        

        double grad_new_dot = inner_product(grad_new.begin(), grad_new.end(), grad_new.begin(), 0.0);
        double grad_dot = inner_product(grad.begin(), grad.end(), grad.begin(), 0.0);

        double beta = grad_new_dot / grad_dot;

        vector<double> s_new(n, 0.0);
        for (int i = 0; i < n; ++i) {
            s_new[i] = -grad_new[i] + beta * s[i];
        }

        // Restart conjugate gradient method for faster convergence
        if (check_linear_independence(s_new, s)) {
            grad_new = gradient(x);
            s_new = grad_new;
            for (int i = 0; i < n; ++i) {
                s_new[i] = -s_new[i];
            }
        }

        s = normalize(s_new);
        grad = grad_new;
    }
}

void bracket_operator_method(vector<double>& xprev,vector<double>& xnext,double epsilon1,double epsilon2,double epsilon3,double delta,vector<double> lo,vector<double> up){
    while(1){
    xprev = xnext;
    conjugate_gradient_optimization_variable(xnext, epsilon1, epsilon2, epsilon3, delta, lo, up);
    if(abs(objective_function(xnext)-objective_function(xprev))<epsilon3){
        break;
    }
    R = C*R;
    }  
}

int main() {
    srand(static_cast<unsigned int>(time(NULL))); 

    int n;
    cout << "Comment or Uncomment the objective function as neccesary" << endl;
    cout << "Enter number of variables: ";
    cin >> n;

    vector<double> x(n);

    double a , b;
    vector<double> lo(n,0.0);
    vector<double> up(n,0.0);
    char yes;
    cout << "Are all Xi bounds same (Y/N): ";
    cin >> yes;
    
    cout << "Input all bounds corresponding to Xi: (Note: Press enter after each 'a b' input): " << endl;
    for(int i=0;i<n;i++){
        cin >> lo[i] >> up[i];
    }
    

    

    char YN;
    int input = 1;
    if(yes == 'Y' || yes =='y'){
        while(input==1){
            cout << "Initial guess for x--> Y: Manual guess  N: Random Guess --> ";
            cin >> YN;
            if(YN == 'Y' || YN == 'y'){
                cout << "Enter initial guess for x: ";
                for (int i = 0; i < n; ++i) {
                    cin >> x[i];
                }
                input = 0;
            }
            else if(YN == 'N' || YN == 'n'){
                for (int i = 0; i < n; ++i) {
                    x[i] = a + (b - a) * ((double) rand() / RAND_MAX);
                }
                input = 0;
            }
            else{
                cout << "Invalid input, try again." << endl;
            }
        }
    }else{
        while(input==1){
            cout << "Initial guess for Xi--> Y: Manual guess  N: Random Guess --> ";
            cin >> YN;
            if(YN == 'Y' || YN == 'y'){
                cout << "Enter initial guess for x: ";
                for (int i = 0; i < n; ++i) {
                    cin >> x[i];
                }
                input = 0;
            }
            else if(YN == 'N' || YN == 'n'){
                for (int i = 0; i < n; ++i) {
                    x[i] = lo[i] + (up[i] - lo[i]) * ((double) rand() / RAND_MAX);
                }
                input = 0;
            }
            else{
                cout << "Invalid input, try again." << endl;
            }
        }
    }

    cout << "Initial guess for x: ";
    for (double xi : x) {
        cout << xi << " ";
    }
    cout << endl;

    double epsilon1 = 0.0001, epsilon2 = 0.0001, epsilon3 = 0.0001, delta = 0.0001;
    cout << "Manually enter epsilon1 (for line search), epsilon2, epsilon3, and delta:(Y/N): ";
    cin >> YN;
    if(YN == 'y' || YN == 'Y'){
        cout << "Enter values respectively : ";
        cin >> epsilon1 >> epsilon2 >> epsilon3 >> delta;
    }else{cout << "The values respectively are :" << epsilon1 <<' '<< epsilon2 <<' '<< epsilon3 <<' '<< delta << endl;}
    
    conjugate_gradient_optimization_variable(x, epsilon1, epsilon2, epsilon3, delta, lo, up);
    R=C*R;
    vector<double> xprev = x;
    vector<double> xnext = x;
    bracket_operator_method(xprev,xnext,epsilon1, epsilon2, epsilon3, delta, lo, up);
    cout << "Optimized variables: ";
    for (double xi : xnext) {
        cout << xi << " ";
    }
    cout << endl;
    cout << "Optimized objective function value: " << obj_fun(xnext) << endl;
    cout << "Total function evaluations: " << func_eval << endl;

    return 0;
}
