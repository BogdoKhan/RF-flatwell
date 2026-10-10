#define NOMINMAX

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <matplot/matplot.h>
#include <rapidcsv.h>


constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
constexpr double ROOT_TOLERANCE = 1.0e-10;
constexpr double CURVATURE_TOLERANCE = 1.0e-8;

double rad2deg(double angle) { return angle * 180.0 / PI; }
double deg2rad(double angle) { return angle * PI / 180.0; }

double wrap_angle(double angle) {
    angle = fmod(angle, TWO_PI);    //remainder of division angle/2pi, if angle > 2pi, return angle - 2pi, otherwise if angle < 0, angle + 2pi
    return angle < 0.0 ? angle + TWO_PI : angle;
}

double minimum_allowed_f(double r) { return r > 0.5 ? -PI / 2.0 : 0.0; } //if r>0.5, we can assume sync phase < 0

double U1p(double f, double f2s, double f1s, double r) {    //U' = V(phi) - Vs
    return sin(f) + r * sin(2.0 * f + f2s) - sin(f1s) - r * sin(2.0 * f1s + f2s);
}

double U(double f, double f2s, double f1s, double r) {              //potential
    const double Vs =  sin(f1s) + r * sin(2.0 * f1s + f2s);         //sync voltage
    return cos(f1s) - cos(f) + 0.5 * r * (cos(2.0 * f1s + f2s) - cos(2.0 * f + f2s)) - (f - f1s) * Vs;  //integral of (V(phi)-Vs) from f1s to f.
}

double U1p_q(double f, double f2s, double f1s, double r) {      //U' BUT trivial zero U'(f1s)=0 is excluded, used for W-shaped well
    return cos((f + f1s) / 2.0) + 2.0 * r * cos(f + f1s + f2s) * cos((f - f1s) / 2.0);
}

double U2p(double f, double f2s, double r) {        //U'', if <0 gives W-shaped well, if > 0 gives single well, if == 0 gives flat well
    return cos(f) + 2.0 * r * cos(2.0 * f + f2s);
}

double U3p(double f, double f2s, double r) {        //U''' if 0 gives symmetric well, important!!!
    return -sin(f) - 4.0 * r * sin(2.0 * f + f2s);
}

double U4p(double f, double f2s, double r) {        //U'''' > 0 gives well instead of hill
    return -cos(f) - 8.0 * r * cos(2.0 * f + f2s);
}

double quartic_r(double f1s) { return 0.25 * sqrt(1.0 + 3.0 * cos(f1s) * cos(f1s)); }   //gives exactly flat well for given f1s

double quartic_f2s(double f1s) {
    const double theta = atan2(-sin(f1s), -2.0 * cos(f1s));                         //f2s for flat well, works with quartic_r
    return wrap_angle(theta - 2.0 * f1s);
}

struct Bucket {
    double f1s = 0.0;
    double f2s = PI;
    double r = 0.0;
    double fc = 0.0;
    bool accepted = false;
};

struct Separatrix {
    std::vector<double> phase_deg;
    std::vector<double> upper;
    std::vector<double> lower;
    double left = 0.0;
    double right = 0.0;
    double energy = 0.0;
    double area = 0.0;
    bool outer_w = false;
    bool valid = false;
};

struct TablePoint {
    double f2s = std::numeric_limits<double>::quiet_NaN();
    double area = std::numeric_limits<double>::quiet_NaN();
    double fc = std::numeric_limits<double>::quiet_NaN();
};

double bisect_root(
    double left, double right, double target,
    double (*function)(double, double, double, double),
    double f2s, double f1s, double r) {
    double value_left = function(left, f2s, f1s, r) - target;

    for (int iteration = 0; iteration < 70; ++iteration) {
        const double middle = 0.5 * (left + right);
        const double value_middle = function(middle, f2s, f1s, r) - target;

        if (value_left * value_middle <= 0.0) right = middle;
        else {left = middle; value_left = value_middle;}
    }
    return 0.5 * (left + right);
}

bool find_adjacent_stationary(
    double center, int direction, bool find_minimum,
    double search_range, double f2s, double f1s,
    double r, double& stationary_point) {
    constexpr double step = PI / 1440.0; // 0.125 degree
    constexpr double offset = 1.0e-7;
    const int maximum_steps =
        static_cast<int>(ceil(search_range / step));

    double previous = center + direction * offset;
    double previous_value = U1p(previous, f2s, f1s, r);

    for (int i = 1; i <= maximum_steps; ++i) {
        const double current = center + direction * i * step;
        const double current_value = U1p(current, f2s, f1s, r);

        if (fabs(current_value) < 1.0e-12 || previous_value * current_value < 0.0) {
            double root = current;

            if (fabs(current_value) >= 1.0e-12) {
                root = bisect_root(std::min(previous, current), std::max(previous, current), 0.0, U1p, f2s, f1s, r);
            }

            const double u2 = U2p(root, f2s, r);
            const bool correct_type =  find_minimum ? u2 > CURVATURE_TOLERANCE : u2 < -CURVATURE_TOLERANCE;

            if (correct_type) {
                stationary_point = root; return true;
            }
        }
        previous = current;
        previous_value = current_value;
    }

    return false;
}

Separatrix build_separatrix(const Bucket& b, bool plotting = false) {

    Separatrix result;
    if (!b.accepted) return result;
    result.outer_w = U2p(b.fc, b.f2s, b.r) < -CURVATURE_TOLERANCE;
    double left_maximum = 0.0;
    double right_maximum = 0.0;

    if (!find_adjacent_stationary(b.fc, -1, false, TWO_PI,
        b.f2s, b.f1s, b.r, left_maximum)) return result;
    if (!find_adjacent_stationary(b.fc, +1, false, TWO_PI,
        b.f2s, b.f1s, b.r, right_maximum)) return result;

    result.energy = std::min(U(left_maximum, b.f2s, b.f1s, b.r), U(right_maximum, b.f2s, b.f1s, b.r));
    if (result.energy <= U(b.fc, b.f2s, b.f1s, b.r) + ROOT_TOLERANCE) return result;
    result.left = bisect_root(left_maximum, b.fc, result.energy, U, b.f2s, b.f1s, b.r);
    result.right = bisect_root(b.fc, right_maximum, result.energy, U, b.f2s, b.f1s, b.r);

    constexpr int samples = 1600;

    double df = (result.right - result.left) / samples;
    double energy_tolerance = 1.0e-9 * std::max(1.0, fabs(result.energy));

    if (plotting) {
        result.phase_deg.resize(samples + 1);
        result.upper.resize(samples + 1);
        result.lower.resize(samples + 1);
    }

    double previous_p = 0.0;
    for (int i = 0; i <= samples; ++i) {
        double f = result.left + i * df;
        double kinetic = result.energy - U(f, b.f2s, b.f1s, b.r);
        if (kinetic < -energy_tolerance) return Separatrix{};
        kinetic = std::max(0.0, kinetic);
        double p = sqrt(2.0 * kinetic);

        if (plotting) {
            result.phase_deg[i] = rad2deg(f);
            result.upper[i] = p;
            result.lower[i] = -p;
        }
        if (i > 0) result.area += (previous_p + p) * df;
        previous_p = p;
    }
    result.valid = true; return result;
}

Eigen::Vector2d equations(
    const Eigen::Vector2d& x, double f1s, double r) {
    const double f = x(0);
    const double f2s = x(1);
    Eigen::Vector2d y;
    y << U1p_q(f, f2s, f1s, r), U3p(f, f2s, r);
    return y;
}

bool solve(Eigen::Vector2d& x, double f1s, double r) {
    constexpr double h = 1.0e-6;

    for (int iteration = 0; iteration < 50; ++iteration) {
        const Eigen::Vector2d y = equations(x, f1s, r);
        if (y.norm() < ROOT_TOLERANCE) return true;

        Eigen::Matrix2d J;
        for (int column = 0; column < 2; ++column) {
            Eigen::Vector2d xp = x; Eigen::Vector2d xm = x;
            xp(column) += h;        xm(column) -= h;
            J.col(column) = (equations(xp, f1s, r) - equations(xm, f1s, r)) / (2.0 * h);
        }

        Eigen::Vector2d step = J.colPivHouseholderQr().solve(-y);
        if (!step.allFinite()) return false;
        if (step.norm() > 0.5) step *= 0.5 / step.norm();
        x += step;
    }

    return equations(x, f1s, r).norm() < ROOT_TOLERANCE;
}

Bucket find_symmetric_bucket(double f1s, double r) {
    Bucket result;
    result.f1s = f1s;
    result.r = r;

    if (r <= 0.0) return result;
    double start_r = quartic_r(f1s);
    Eigen::Vector2d root(f1s, quartic_f2s(f1s));

    if (fabs(f1s) > ROOT_TOLERANCE) {                                   // Continue from the flat solution to the requested r.
        int steps = std::max(1, (int)(ceil(fabs(r - start_r)/0.01)));
        for (int i = 1; i <= steps; ++i) {
            double current_r = start_r + (r - start_r) * i / steps;
            if (!solve(root, f1s, current_r)) return result;
        }
    }

    result.fc = root(0);
    result.f2s = wrap_angle(root(1));
    if (result.fc < minimum_allowed_f(r) - ROOT_TOLERANCE || result.fc > PI / 2.0 + ROOT_TOLERANCE) return result;

    double u2 = U2p(result.fc, result.f2s, r);
    double u4 = U4p(result.fc, result.f2s, r);

    if (u2 > CURVATURE_TOLERANCE || u4 <= 0.0) return result;
    result.accepted = true;
    return result;
}

Bucket find_best_single_well(double f1s, double r) {
    Bucket result;
    result.f1s = f1s;
    result.r = r;
    result.fc = f1s;

    if (r < 0.0 || f1s < minimum_allowed_f(r) ||  f1s > PI / 2.0) return result;

    if (r < ROOT_TOLERANCE) result.f2s = 0.0;   // The second-harmonic phase is irrelevant when r = 0.
    else {
        double wanted_sin = std::clamp(-sin(f1s) / (4.0 * r), -1.0, 1.0);
        double angle = PI - asin(wanted_sin); // lower-curvature stable candidate.
        result.f2s = wrap_angle(angle - 2.0 * f1s);
    }

    if (U2p(f1s, result.f2s, r) <= CURVATURE_TOLERANCE) return result;

    result.accepted = true;
    return result;
}

bool select_bucket(
    double f1s, double r,
    Bucket& selected, Separatrix& separatrix) {
    const double flat_threshold = quartic_r(f1s);

    if (r >= flat_threshold - ROOT_TOLERANCE) selected = find_symmetric_bucket(f1s, r);  //flat well @threshold, above W-shaped
    else selected = find_best_single_well(f1s, r);                                      //otherwise well with U'' > 0

    separatrix = build_separatrix(selected);
    return selected.accepted && separatrix.valid;
}

TablePoint calculate_table_point(double f1s, double r) {
    Bucket selected;
    Separatrix separatrix;

    if (!select_bucket(f1s, r, selected, separatrix)) return {};
    return {selected.f2s, separatrix.area, selected.fc};
}

std::vector<std::vector<TablePoint>> make_table() {
    constexpr int maximum_f1s_deg = 85;
    constexpr int maximum_r_index = 50;

    std::vector<std::vector<TablePoint>> table(maximum_f1s_deg + 1, std::vector<TablePoint>(maximum_r_index + 1));

    for (int f1s_deg = 0;
        f1s_deg <= maximum_f1s_deg; ++f1s_deg) {
        const double f1s = deg2rad(static_cast<double>(f1s_deg));

        for (int r_index = 0;
            r_index <= maximum_r_index; ++r_index) {
            const double r = static_cast<double>(r_index) / 100.0;
            table[f1s_deg][r_index] = calculate_table_point(f1s, r);
        }

        std::cout << "\rCalculated f1s = " << f1s_deg << " deg" << std::flush;
    }

    std::cout << '\n';
    return table;
}

bool save_table(
    const std::vector<std::vector<TablePoint>>& table,
    const char* filename) {
    std::ofstream file(filename);
    if (!file) return false;

    file << "f1s_deg,r,f2s_deg,area_normalized,fc_deg\n";
    file << std::fixed << std::setprecision(12);

    for (std::size_t f1s_index = 0; f1s_index < table.size(); ++f1s_index) {
        for (std::size_t r_index = 0; r_index < table[f1s_index].size(); ++r_index) {
            const TablePoint& point = table[f1s_index][r_index];

            file << f1s_index << ',' << static_cast<double>(r_index) / 100.0
                << ',' << rad2deg(point.f2s) << ',' << point.area << ',' << rad2deg(point.fc) << '\n';
        }
    }

    return static_cast<bool>(file);
}

void print_bucket(const Bucket& b, const char* title) {
    const double u2 = U2p(b.fc, b.f2s, b.r);

    std::cout << "\n" << title << '\n'
        << std::fixed << std::setprecision(9)
        << "f = " << rad2deg(b.fc)
        << " deg, f2s = " << rad2deg(b.f2s)
        << " deg, r = " << b.r << '\n'
        << std::scientific
        << "U = " << U(b.fc, b.f2s, b.f1s, b.r)
        << ", U1p = " << U1p(b.fc, b.f2s, b.f1s, b.r)
        << ", q = " << U1p_q(b.fc, b.f2s, b.f1s, b.r)
        << ", U2p = " << u2
        << ", U3p = " << U3p(b.fc, b.f2s, b.r)
        << ", U4p = " << U4p(b.fc, b.f2s, b.r)
        << '\n'
        << (b.accepted ? "Accepted\n" : "Rejected\n");
}

void plot_bucket(const Bucket& b, const Separatrix& separatrix) {
    if (!separatrix.valid) return;

    std::vector<double> v_f1s;
    std::vector<double> v_u;
    std::vector<double> v_hsep;

    for (int i = -180; i < 180; ++i) {
        const double f = deg2rad(static_cast<double>(i));

        v_f1s.push_back(static_cast<double>(i));
        v_u.push_back(U(f, b.f2s, b.f1s, b.r));
        v_hsep.push_back(separatrix.energy);
    }

    matplot::plot(v_f1s, v_u);
    matplot::hold(matplot::on);
    matplot::plot(v_f1s, v_hsep);
    matplot::xlabel("phase f [deg]");
    matplot::ylabel("normalized potential U");
    matplot::title("Potential well and separatrix energy");
    matplot::hold(matplot::off);
    matplot::show();
}

void print_separatrix(const Separatrix& s) {
    if (!s.valid) {
        std::cout << "\nNo connected separatrix was found.\n";
        return;
    }

    std::cout << "\n"
        << (s.outer_w ? "Outer W-bucket separatrix"
            : "Single-well separatrix")
        << '\n'
        << std::fixed << std::setprecision(9)
        << "left = " << rad2deg(s.left)
        << " deg, right = " << rad2deg(s.right)
        << " deg, width = "
        << rad2deg(s.right - s.left) << " deg\n"
        << std::scientific
        << "Hsep = " << s.energy
        << ", normalized phase-space area = "
        << s.area << '\n';
}

void plot_separatrix(const Separatrix& s) {
    if (!s.valid) return;

    matplot::plot(s.phase_deg, s.upper);
    matplot::hold(matplot::on);
    matplot::plot(s.phase_deg, s.lower);
    matplot::xlabel("phase f [deg]");
    matplot::ylabel("normalized momentum P");
    matplot::title(
        s.outer_w ? "Outer W-bucket separatrix"
        : "Single-well separatrix");
    matplot::hold(matplot::off);
    matplot::show();
}

double interpolate(double a, double b, double fraction) {
    return a + fraction * (b - a);
}

double interpolate_angle(
    double a, double b, double fraction) {
    double difference = b - a;

    if (difference > PI) difference -= TWO_PI;
    if (difference < -PI) difference += TWO_PI;

    return wrap_angle(a + fraction * difference);
}

TablePoint table_point_at_r(const std::vector<std::vector<TablePoint>>& table, int f1s_deg, double r) {
    double r_position = 100.0 * r;
    int r_left = static_cast<int>(r_position);

    // Use columns 0.49 and 0.50 when r = 0.50.
    if (r_left == 50) r_left = 49;

    int r_right = r_left + 1;
    double fraction = r_position - r_left;

    TablePoint point;
    point.area = interpolate(table[f1s_deg][r_left].area, table[f1s_deg][r_right].area, fraction);
    point.fc = interpolate(table[f1s_deg][r_left].fc, table[f1s_deg][r_right].fc, fraction);

    if (r <= 0.01) point.f2s = table[f1s_deg][1].f2s;
    else point.f2s = interpolate(table[f1s_deg][r_left].f2s, table[f1s_deg][r_right].f2s, fraction);

    return point;
}

bool find_acceleration_phases(const std::vector<std::vector<TablePoint>>& table, double AcptStat, double C, double r,
    double& f1s, double& f2s, double& fc) {

    if (AcptStat <= 0.0 || C <= 0.0) return false;
    if (r < 0.0 || r > 0.5) return false;

    TablePoint left = table_point_at_r(table, 0, r);
    double target_area = left.area * C / AcptStat;

    // Scan upward and return the lowest solution.
    for (int f = 0; f < 85; ++f) {
        TablePoint right = table_point_at_r(table, f + 1, r);
        double minimum_area = std::min(left.area, right.area);
        double maximum_area = std::max(left.area, right.area);

        if (target_area >= minimum_area && target_area <= maximum_area) {
            double area_difference = right.area - left.area;
            double fraction = 0.0;

            // If the interval is horizontal, its lowest phase is f.
            if (area_difference != 0.0) fraction = (target_area - left.area) / area_difference;

            f1s = deg2rad(f + fraction);
            f2s = interpolate(left.f2s, right.f2s, fraction);
            fc = interpolate(left.fc, right.fc, fraction);

            return true;
        }

        left = right;
    }

    return false;
}


std::vector<std::vector<TablePoint>> load_table(const std::string& filename) {
    rapidcsv::Document csv(filename, rapidcsv::LabelParams(0, -1));
    std::vector<std::vector<TablePoint>> table(85+1, std::vector<TablePoint>(51));

    for (std::size_t row = 0; row < csv.GetRowCount(); ++row) {
        int f1s_deg = csv.GetCell<int>("f1s_deg", row);
        int r_index = static_cast<int>(std::lround(100.0 * csv.GetCell<double>("r", row)));
        // Ignore values outside the required table.
        if (f1s_deg > 85 || r_index > 50) continue;
        table[f1s_deg][r_index] = {
            deg2rad(csv.GetCell<double>("f2s_deg", row)),
            csv.GetCell<double>("area_normalized", row), deg2rad(csv.GetCell<double>("fc_deg", row)) };
    }

    return table;
}

int main() {
    const char* filename = "rf_well_table.csv";

    const std::vector<std::vector<TablePoint>> table = make_table();

    if (!save_table(table, filename)) {
        std::cerr << "Could not write " << filename << ".\n";
        return 1;
    }

    //std::vector<std::vector<TablePoint>> table =
    //    load_table("rf_well_table.csv");

    double C = 0.16664;
    double A = 0.39037;
    double r = 0.32;
    double f1s, f2s, fc;

    find_acceleration_phases(table, A, C, r, f1s, f2s, fc);

    std::cout << "f1s: " << f1s * 180 / PI << " f2s: " << f2s * 180 / PI << ".\n";
    return 0;
}

