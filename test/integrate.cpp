#include <integrate.hpp>
#include "constants.hpp"
#include <iostream>
#include <cmath>

std::size_t g_n_total_calls = 0;
double sine_wave(double t, double frequency, double amplitude, double phase) {
    ++g_n_total_calls;
    return amplitude * std::sin(2 * constants::pi_d * frequency * t + phase);
}

double box(double x, double height) {
    ++g_n_total_calls;
    if (x > -0.5 && x < 0.5) {
        return height;
    }
    return 0.0;
}

double constant(double x, double height) {
    (void)(x);
    ++g_n_total_calls;
    return height;
}

double difficult_function(double x) {
    ++g_n_total_calls;
    return std::sin(1000.0 * x) / (1.0 + 10000.0 * (x - 0.37) * (x - 0.37));
}

double sharp_peak(double x) {
    ++g_n_total_calls;
    return 1.0 / (1.0 + 1'000'000.0 * (x - 0.5) * (x - 0.5));
}

double endpoint_singularity(double x) {
    ++g_n_total_calls;
    return 1.0 / std::sqrt(x);
}

int main() {
    {
        g_n_total_calls = 0;
        Integrate::IntegrationResult<double> ret = Integrate::integrate_gauss_kronrod(
            sine_wave, 0., 0.5, 1e-7, 1e-7, 1.0, 2.0, 0.0
        );
        std::cout << ret << std::endl;
        std::cout << "Number of function calls: " << g_n_total_calls << std::endl;
    }

    {
        g_n_total_calls = 0;
        Integrate::IntegrationResult<double> ret = Integrate::integrate_gauss_kronrod(
            box, -constants::pi_d, constants::e_d, 1e-7, 1e-7, 1.0
        );
        std::cout << ret << std::endl;
        std::cout << "Number of function calls: " << g_n_total_calls << std::endl;
    }
    
    {
        g_n_total_calls = 0;
        Integrate::IntegrationResult<double> ret = Integrate::integrate_gauss_kronrod(
            constant, -constants::pi_d, constants::e_d, 1e-7, 1e-7, 1.0
        );
        std::cout << ret << std::endl;
        std::cout << "Number of function calls: " << g_n_total_calls << std::endl;
    }

    {
        g_n_total_calls = 0;
        Integrate::IntegrationResult<double> ret = Integrate::integrate_gauss_kronrod(
            difficult_function, 0., 1., 1e-7, 1e-7
        );
        std::cout << ret << std::endl;
        std::cout << "Number of function calls: " << g_n_total_calls << std::endl;
    }

    {
        g_n_total_calls = 0;
        Integrate::IntegrationResult<double> ret = Integrate::integrate_gauss_kronrod(
            sharp_peak, 0., 1., 1e-7, 1e-7
        );
        std::cout << ret << std::endl;
        std::cout << "Number of function calls: " << g_n_total_calls << std::endl;
    }

    {
        g_n_total_calls = 0;
        Integrate::IntegrationResult<double> ret = Integrate::integrate_gauss_kronrod(
            endpoint_singularity, 0., 1., 1e-7, 1e-7
        );
        std::cout << ret << std::endl;
        std::cout << "Number of function calls: " << g_n_total_calls << std::endl;
    }
}