#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>

// Objective Function
double objective_function(double x, int choice, int& func_eval) {
    func_eval++;
    switch (choice) {
        case 1:
            return pow(x*x - 1, 3) - pow(2*x - 5, 4);
        case 2:
            return 2*x + 2*exp(x) - pow(x, 3) - 8;
        case 3:
            return -4*x * sin(x);
        case 4:
            return 2*pow(x - 3, 2) + exp(0.5*pow(x, 2));
        case 5:
            return pow(x, 2) - 10*exp(0.1*x);
        case 6:
            return 15*pow(x, 2) - 20*sin(x);
        default:
            throw std::invalid_argument("Invalid function choice");
    }
}

// Bounding Phase Method
std::pair<double, double> bounding_phase_method(std::function<double(double)> f, double a, double b, double delta, int& func_eval, int max_iter = 1000) {
    double x = a + (b - a) * ((double) rand() / RAND_MAX);  // Random initial guess in the interval
    int k = 0;

    double f_x_minus = f(x - fabs(delta));
    double f_x = f(x);
    double f_x_plus = f(x + fabs(delta));

    if (f_x_minus >= f_x && f_x <= f_x_plus) {
        delta = fabs(delta);
    } else if (f_x_minus <= f_x && f_x >= f_x_plus) {
        delta = -fabs(delta);
    } else {
        return bounding_phase_method(f, a, b, delta, func_eval);  // Retry with a new random guess
    }

    while (k < max_iter) {
        double x_next = x + pow(2, k) * delta;
        if (f(x_next) < f_x) {
            x = x_next;
            f_x = f(x);
            k++;
        } else {
            double x_high = x_next;
            double x_low = x - pow(2, k-1) * delta;
            return {x_low, x_high};
        }
    }

    throw std::runtime_error("Bounding Phase Method failed to find a suitable interval within max iterations.");
}

// Golden Section Search Method
std::pair<double, double> golden_section_search(std::function<double(double)> f, double a, double b, double epsilon, int& func_eval, int max_iter = 1000) {
    const double golden_ratio = (1 + sqrt(5)) / 2;
    const double inv_golden_ratio = 1 / golden_ratio;

    // Normalize the interval [a, b] to [0, 1]
    double a_omega = 0;
    double b_omega = 1;
    double L_omega = b_omega - a_omega;
    int k = 1;

    // Initial points
    double omega1 = a_omega + (1 - inv_golden_ratio) * L_omega;
    double omega2 = a_omega + inv_golden_ratio * L_omega;

    // Rescale points
    double x1 = a + omega1 * (b - a);
    double x2 = a + omega2 * (b - a);

    double f_x1 = f(x1);
    double f_x2 = f(x2);

    while (L_omega > epsilon && k < max_iter) {
        if (f_x1 < f_x2) {
            b_omega = omega2;
            omega2 = omega1;
            f_x2 = f_x1;
            L_omega = b_omega - a_omega;
            omega1 = a_omega + (1 - inv_golden_ratio) * L_omega;
            x1 = a + omega1 * (b - a);
            f_x1 = f(x1);
        } else {
            a_omega = omega1;
            omega1 = omega2;
            f_x1 = f_x2;
            L_omega = b_omega - a_omega;
            omega2 = a_omega + inv_golden_ratio * L_omega;
            x2 = a + omega2 * (b - a);
            f_x2 = f(x2);
        }

        k++;
    }

    // Return the best point
    if (f_x1 < f_x2) {
        return {x1, f_x1};
    } else {
        return {x2, f_x2};
    }
}

int main() {
    std::srand(std::time(0));  // Seed the random number generator

    std::cout << "Choose an objective function:" << std::endl;
    std::cout << "1. (x^2 - 1)^3 - (2x - 5)^4" << std::endl;
    std::cout << "2. 2x + 2e^x - x^3 - 8" << std::endl;
    std::cout << "3. -4x * sin(x)" << std::endl;
    std::cout << "4. 2(x - 3)^2 + e^(0.5x^2)" << std::endl;
    std::cout << "5. x^2 - 10e^(0.1x)" << std::endl;
    std::cout << "6. 15x^2 - 20sin(x)" << std::endl;

    int choice;
    std::cout << "Enter the number of the function you want to use: ";
    std::cin >> choice;

    if (choice < 1 || choice > 6) {
        std::cerr << "Invalid choice." << std::endl;
        return 1;
    }

    double a, b, delta, epsilon;
    std::cout << "Enter the lower limit of the interval a: ";
    std::cin >> a;
    std::cout << "Enter the upper limit of the interval b: ";
    std::cin >> b;
    std::cout << "Enter the increment delta for Bounding Phase Method: ";
    std::cin >> delta;
    std::cout << "Enter the small number epsilon for Golden Section Search: ";
    std::cin >> epsilon;

    int func_eval = 0;

    // Bounding Phase Method
    auto [x_low, x_high] = bounding_phase_method(
        [choice, &func_eval](double x) { return objective_function(x, choice, func_eval); }, a, b, delta, func_eval
    );
    std::cout << "Bounding Phase Method -> Minimum point lies in the interval: [" << x_low << ", " << x_high << "]" << std::endl;

    // Golden Section Search Method with the interval found by Bounding Phase Method
    auto [x_min, f_min] = golden_section_search(
        [choice, &func_eval](double x) { return objective_function(x, choice, func_eval); }, x_low, x_high, epsilon, func_eval
    );
    std::cout << "Golden Section Search Method Result:" << std::endl;
    std::cout << "The minimum value found is at x = " << x_min << " with function value f(x) = " << f_min << std::endl;
    std::cout << "Total number of function evaluations: " << func_eval << std::endl;

    return 0;
}