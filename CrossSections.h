//
// Created by Trevor on 8/28/2026.
//

#pragma once
#import <iostream>
#import <vector>
#import <cmath>
#import <numbers>
#import <numeric>
#import <fstream>

#ifndef NE684A2_CROSSSECTIONS_H
#define NE684A2_CROSSSECTIONS_H

#define ENERGY_MIN 1e-5
#define ENERGY_MAX 2e7
#include <filesystem>
#include <ranges>

/**
 * @brief Holds isotope data used to generate cross-section data.
 */
struct Fissionable_Isotope {
    float N{};
    float g_J{};
    double E_lambda{};
    double gamma_n{};
    double gamma_gamma{};
    double gamma_f{};
    double a_c{};
    double mass{};
    bool do_scattering_override = false;
    double rho_0 = 0.002196807122623 / 2.0;
};


class CrossSections0K {
public:
    /**
     * @brief Calculates the cross section data for a single isotope at a given energy.
     * @param energy The energy the collision even happened at.
     * @param isotope An isotope the collision may have happened with.
     * @return A vector of cross section values ordered by scatter, capture, fission.
     */
    static std::vector<double> getFissionableCrossSections(const double energy, const Fissionable_Isotope *isotope) {
        std::vector<double> XSec(3); // ordered scatter-capture-fission

        const double k = isotope->rho_0 * std::sqrt(energy);
        const double gamma_n = isotope->gamma_n * std::sqrt(energy / isotope->E_lambda);
        const double gamma_t = gamma_n + isotope->gamma_gamma + isotope->gamma_f;
        const double d = std::pow(energy - isotope->E_lambda, 2.0) + std::pow(gamma_t / 2.0, 2);

        XSec[0] = 4.0 * std::numbers::pi * std::pow(isotope->a_c, 2.0) * isotope->N +
                isotope->N * (isotope->g_J * std::numbers::pi / d) * (
                std::pow(gamma_n / k, 2) +
                4.0 * isotope->a_c * (energy - isotope->E_lambda) * gamma_n / k  +
                -2.0 * std::pow(isotope->a_c, 2.0) * gamma_n * gamma_t
            );

        XSec[1] = isotope->N * (isotope->g_J * std::numbers::pi / d) * gamma_n * isotope->gamma_gamma / std::pow(k, 2);
        XSec[2] = isotope->N * (isotope->g_J * std::numbers::pi / d) * gamma_n * isotope->gamma_f     / std::pow(k, 2);

        return XSec;
    }

    static double getFissionableCrossSection(const unsigned int interaction_index, const double energy, const Fissionable_Isotope *isotope) {
        const double k = isotope->rho_0 * std::sqrt(energy);
        const double gamma_n = isotope->gamma_n * std::sqrt(energy / isotope->E_lambda);
        const double gamma_t = gamma_n + isotope->gamma_gamma + isotope->gamma_f;
        const double d = std::pow(energy - isotope->E_lambda, 2.0) + std::pow(gamma_t / 2.0, 2);

        switch (interaction_index) {
            case 2:
                return isotope->N * (isotope->g_J * std::numbers::pi / d) * gamma_n * isotope->gamma_f / std::pow(k, 2);
            case 1:
                return isotope->N * (isotope->g_J * std::numbers::pi / d) * gamma_n * isotope->gamma_gamma / std::pow(k, 2);
            case 0:
                return 4.0 * std::numbers::pi * std::pow(isotope->a_c, 2.0) * isotope->N +
                        isotope->N * (isotope->g_J * std::numbers::pi / d) * (
                        std::pow(gamma_n / k, 2) +
                        4.0 * isotope->a_c * (energy - isotope->E_lambda) * gamma_n / k  +
                        -2.0 * std::pow(isotope->a_c, 2.0) * gamma_n * gamma_t
                    );
            default:
                throw std::invalid_argument("Invalid interaction index");
        }
    }
};

class DopplerBroadener {
public:
    static double broaden(const double energy, const double temperature, const unsigned int i, const Fissionable_Isotope *isotope) {
        constexpr double deviations = 6.0;
        const double beta = isotope->mass / (k_B * temperature);
        const double gauss_width = deviations / std::sqrt(beta);

        const double v_low = std::max(0.0, std::sqrt(energy) - gauss_width);
        const double v_high = std::sqrt(energy) + gauss_width;

        return std::sqrt(beta / std::numbers::pi) / (2.0 * energy) * composite_simpson(
            v_low, v_high, i, std::sqrt(energy), beta, isotope);
    }

    static double composite_simpson(const double a, const double b, const unsigned int interaction, const double velocity, const double beta, const Fissionable_Isotope *isotope) {
        constexpr unsigned int N = 100;
        const double h = (b - a) / (2.0 * N);
        double sum = internal_function(interaction, velocity, a, beta, isotope);

        for (unsigned int i = 1; i < N; ++i) {
            sum += 4 * internal_function(interaction, velocity, a + h * (2.0 * i - 1.0), beta, isotope);
            sum += 2 * internal_function(interaction, velocity, a + h * (2.0 * i),       beta, isotope);
        }

        //correction factor so as f(x_n) is not doubled in composite simpson
        sum -= internal_function(interaction, velocity, b, beta, isotope);

        sum *= h/3.0;
        return sum;
    }

    static double internal_function(const unsigned int interaction, const double velocity, const double v_prime, const double beta, const Fissionable_Isotope *isotope) {
        if (v_prime <= 0.0) return 0.0;
        double xSec_0K = CrossSections0K::getFissionableCrossSection(interaction, v_prime, isotope);
        if (isotope->do_scattering_override) {
            xSec_0K = 1000;
        }
        //std::sqrt(e_prime) => 2.0 * std::pow(e_prime, 2.0)
        return xSec_0K * 2.0 * std::pow(v_prime, 2.0) * (
                std::exp(-beta * std::pow(v_prime - velocity, 2.0)) -
                std::exp(-beta * std::pow(v_prime + velocity, 2.0))
            );
    }

private:
    static constexpr double k_B = 8.61733326e-5; //ev-K^-1
};

/**
 * @brief Holds continuous-energy cross section data.
 */
class DopplerCrossSections {
public:
    explicit DopplerCrossSections() = default;

    void LoadCrossSections(double temperature) {
        LoadedCrossSections.clear();
        bool exists = std::filesystem::exists(std::format("cross_sections_{}K.csv", static_cast<int>(temperature)));
        if (!exists) {
            write_cross_sections(temperature);
        }
        const std::string file_name = std::format("cross_sections_{}K.csv", static_cast<int>(temperature));
        std::ifstream file(file_name);
        if (!file.is_open()) {
            throw std::runtime_error(std::format("Could not open file: {}", file_name));
        }
        std::string line;
        while (std::getline(file, line)) {
            std::vector<double> row;
            std::stringstream ss(line);
            std::string cell;

            while (std::getline(ss, cell, ',')) {
                row.push_back(std::stod(cell));
            }
            LoadedCrossSections.push_back(row);
        }
        file.close();
    }

    /**
     * @brief Calculates the cross sections for a collision at a given energy.
     * @param energy The energy the collision even happened at.
     * @return A vector of cross section values ordered by Pu:(scatter, capture, fission), U:(scatter, capture, fission), Unknown:(scatter), Total.
     */
    [[nodiscard]] std::vector<double> getCrossSections(const double energy) const {
        std::vector<double> XSec(8);
        const double logEmin = std::log(ENERGY_MIN);
        const double logEmax = std::log(ENERGY_MAX);
        const double delta_logE = (logEmax - logEmin) / static_cast<double>(LoadedCrossSections.size() - 1);

        const double position = (std::log(energy) - logEmin) / delta_logE;
        const auto index = static_cast<size_t>(std::floor(position));
        const double fraction = position - static_cast<double>(index);

        for (unsigned int i = 0; i < XSec.size(); ++i) {
            XSec[i] = LoadedCrossSections[index][i + 1] * (1.0 - fraction) +
                LoadedCrossSections[index + 1][i + 1] * fraction;
        }

        return XSec;
    }

    const std::vector<float> neutronsFromFission = {2.88f, 0.0f, 0.0f};
    const std::vector<float> targetMasses = {239.0f, 238.0f, 12.0f};

private:
    void write_cross_sections(const double temperature) const {
        ///starting at e_min, walk through vector of cross sections at this energy and append to file

        // Step 1: open and clear file.
        const std::string file_name = std::format("cross_sections_{}K.csv", static_cast<int>(temperature));
        std::ofstream outfile(file_name);
        if (!outfile.is_open()) {
            throw std::runtime_error(std::format("Could not open file: {}", file_name));
        }
        outfile.clear();

        // Step 2: determine the point spacing in lethargy
        constexpr unsigned int N = 10000;
        const double delta_l = std::log(ENERGY_MAX / ENERGY_MIN) / N;

        // Step 3a: cycle through each energy point
        for (unsigned int i = 0; i< N; ++i) {
            double energy = std::exp(std::log(ENERGY_MIN) + i * delta_l);
            outfile << energy << ",";

                // Step 3b and 3c should be compressible by making a vector of my isotopes
            // Step 3b: write each Xsec value to a position in a vector
            std::vector<double> XSec(8);
            XSec[0] = DopplerBroadener::broaden(energy, temperature, 0, &Pu239);
            XSec[1] = DopplerBroadener::broaden(energy, temperature, 1, &Pu239);
            XSec[2] = DopplerBroadener::broaden(energy, temperature, 2, &Pu239);
            XSec[3] = DopplerBroadener::broaden(energy, temperature, 0, &U238);
            XSec[4] = DopplerBroadener::broaden(energy, temperature, 1, &U238);
            XSec[5] = DopplerBroadener::broaden(energy, temperature, 2, &U238);
            XSec[6] = DopplerBroadener::broaden(energy, temperature, 0, &C12);
            XSec[7] = std::accumulate(XSec.begin(), XSec.end(), 0.0);

            // Step 3c: write line to file
            for (size_t j = 0; j < XSec.size(); ++j) {
                outfile << XSec[j];
                if (j < XSec.size() - 1) {
                    outfile << ",";
                }
            }
            outfile << "\n";
        }
        outfile.close();
    }

    const Fissionable_Isotope Pu239 = {
        .N = 1.0f, .g_J = 3.0f/4.0f, .E_lambda = 2.956243e-1, .gamma_n = 7.947046e-5,
        .gamma_gamma = 3.982423e-2, .gamma_f = 5.619673e-2, .a_c = 9.41e-4, .mass = 239
    };

    const Fissionable_Isotope U238 = {
        .N = 0.124954f, .g_J = 1.0f, .E_lambda = 6.674280e0, .gamma_n = 1.4923e-2,
        .gamma_gamma = 2.2711e-2, .gamma_f = 9.88e-9, .a_c = 9.48e-4, .mass = 238
    };

    const Fissionable_Isotope C12 = {
        0.0f, 0.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, true
    };

    std::vector<std::vector<double>> LoadedCrossSections;
};

#endif //NE684A2_CROSSSECTIONS_H
