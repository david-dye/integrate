#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <ostream>
#include <queue>
#include <string>

class Integrate {
public:
    Integrate() = delete;
    ~Integrate() = delete;

    static constexpr std::size_t MAX_N_INTERVALS = 100'000;

    template <typename T>
    requires std::floating_point<T>
    struct IntegrationResult {
        
        IntegrationResult(
            T value,
            T absolute_error_estimate
        ) : 
            value(value),
            absolute_error_estimate(absolute_error_estimate)
        {}
        
        IntegrationResult(
            T value,
            T absolute_error_estimate,
            bool success,
            const std::string& details
        ) :
            value(value),
            absolute_error_estimate(absolute_error_estimate),
            success(success),
            details(details)
        {}
        
        T value;
        T absolute_error_estimate;
        bool success = true;
        std::string details = "success";
        
        friend std::ostream& operator<<(std::ostream& os, const IntegrationResult& result) {
            os << "Success: " << (result.success ? "true" : "false") 
               << "\nMessage: " << result.details 
               << "\nValue:   " << result.value << " +/- " << result.absolute_error_estimate;
            return os;
        }
    };


    /*
    Applies the Gauss-Kronrod quadrature formula to integrate
    an input function over some finite interval. This uses the
    (G7, K15) rule.
    */
    template <typename T, typename F>
    requires std::floating_point<T> && std::invocable<F&, T>
    static IntegrationResult<T> integrate_gauss_kronrod(
        F&& fn,
        T x_min,
        T x_max
    ) {
        return integrate_gauss_kronrod(
            fn, x_min, x_max, 1e-7, 1e-7
        );
    }


    template <typename T, typename F, typename ...Args>
    requires std::floating_point<T> && std::invocable<F&, T, Args&...>
    static IntegrationResult<T> integrate_gauss_kronrod(
        F&& fn,
        T x_min,
        T x_max,
        T abs_convergence_threshold,
        T rel_convergence_threshold,
        Args&&... params
    ) {
        const std::array<T, 2> init = compute_integral(fn, x_min, x_max, params...);
        T integral = init[0];
        T error = init[1];
        std::priority_queue<Interval<T>> pq;
        pq.push(Interval<T>(x_min, x_max, integral, error));
        
        while (!has_converged(
            integral,
            error,
            abs_convergence_threshold,
            rel_convergence_threshold
        )) {
            if (pq.size() >= MAX_N_INTERVALS) {
                return IntegrationResult<T>(
                    integral,
                    error,
                    false,
                    "Convergence Error: integration maximum interval count reached."
                );
            }
            const Interval<T>& worst_interval = pq.top();
            const T x_min = worst_interval.x_min;
            const T x_max = worst_interval.x_max;
            const T x_mid = std::midpoint(x_min, x_max);
            integral -= worst_interval.value;
            error -= worst_interval.error;
            const std::array<T, 2> left = compute_integral(fn, x_min, x_mid, params...);
            const std::array<T, 2> right = compute_integral(fn, x_mid, x_max, params...);
            pq.pop();
            pq.push(Interval<T>(x_min, x_mid, left[0], left[1]));
            pq.push(Interval<T>(x_mid, x_max, right[0], right[1]));
            integral += left[0] + right[0];
            error += left[1] + right[1];
        }

        return IntegrationResult<T>(integral, error);
    }

private:

    template <typename T>
    requires std::floating_point<T>
    struct Interval {
        
        Interval(
            T x_min,
            T x_max,
            T value,
            T error
        ) : 
            x_min(x_min),
            x_max(x_max),
            value(value),
            error(error) 
        {}
        
        T x_min;
        T x_max;
        T value;
        T error;

        // Define an ordering on intervals based on error
        bool operator<(const Interval& other) const {
            return error < other.error;
        }

    };

    template <typename T>
    static bool has_converged(T integral, T error, T abs_tol, T rel_tol) {
        return error < std::max(abs_tol, rel_tol * std::abs(integral));
    }

    template <typename T, typename F, typename ...Args>
    requires std::floating_point<T> && std::invocable<F&, T, Args&...>
    static std::array<T, 2> compute_integral(
        F& fn,
        T x_min,
        T x_max,
        Args&... params
    ) {
        constexpr std::array<std::array<T, 2>, 7> g7 = {{
            {-0.949107912342759, 0.129484966168870},
            {-0.741531185599394, 0.279705391489277},
            {-0.405845151377397, 0.381830050505119},
            {0.0,                0.417959183673469},
            {0.405845151377397,  0.381830050505119},
            {0.741531185599394,  0.279705391489277},
            {0.949107912342759,  0.129484966168870}
        }};

        constexpr std::array<std::array<T, 2>, 15> k15 = {{
            {-0.991455371120813, 0.022935322010529},
            {-0.949107912342759, 0.063092092629979},
            {-0.864864423359769, 0.104790010322250},
            {-0.741531185599394, 0.140653259715525},
            {-0.586087235467691, 0.169004726639267},
            {-0.405845151377397, 0.190350578064785},
            {-0.207784955007898, 0.204432940075298},
            {0.0,                0.209482141084728},
            {0.207784955007898,  0.204432940075298},
            {0.405845151377397,  0.190350578064785},
            {0.586087235467691,  0.169004726639267},
            {0.741531185599394,  0.140653259715525},
            {0.864864423359769,  0.104790010322250},
            {0.949107912342759,  0.063092092629979},
            {0.991455371120813,  0.022935322010529}
        }};

        T kronrad_val = 0;
        T gauss_val = 0;
        const T half_width = (x_max - x_min) / 2;
        for (std::size_t i = 0; i < k15.size(); ++i) {
            const T x = x_min + half_width * (k15[i][0] + 1);
            const T w_k = k15[i][1] * half_width;
            const T fn_out = fn(x, params...);
            kronrad_val += w_k * fn_out;
            if (i % 2 == 1) {
                // Kronrad point is also a Gauss point
                const T w_g = g7[(i - 1) / 2][1] * half_width;
                gauss_val += w_g * fn_out; 
            }
        }

        T abs_error = std::abs(kronrad_val - gauss_val);

        return std::array<T, 2>({kronrad_val, abs_error});
    }
};

