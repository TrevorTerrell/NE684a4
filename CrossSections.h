//
// Created by Trevor on 8/28/2026.
//

#pragma once
#import <vector>
#import <cmath>
#import <numbers>
#import <numeric>
#import <fstream>

#ifndef NE684A2_CROSSSECTIONS_H
#define NE684A2_CROSSSECTIONS_H

#define ENERGY_MIN 1e-5
#define ENERGY_MAX 2e7

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
    double rho_0 = 0.002196807122623 / 2.0;
};

/**
 * @brief Holds continuous-energy cross section data.
 */
class CrossSections {
public:
    explicit CrossSections() = default;

    /**
     * @brief Calculates the cross sections for a collision at a given energy.
     * @param energy The energy the collision even happened at.
     * @return A vector of cross section values ordered by Pu:(scatter, capture, fission), U:(scatter, capture, fission), Unknown:(scatter), Total.
     */
    [[nodiscard]] std::vector<double> getCrossSections(const double energy) const {
        std::vector<double> XSec;
        XSec.reserve(8);

        std::vector<double> XSec_Pu = getFissionableCrossSections(energy, &Pu239);
        XSec.insert(XSec.end(), XSec_Pu.begin(), XSec_Pu.end());

        std::vector<double> XSec_U = getFissionableCrossSections(energy, &U238);
        XSec.insert(XSec.end(), XSec_U.begin(), XSec_U.end());

        XSec.push_back(0.1668101);
        XSec.push_back(std::accumulate(XSec.begin(), XSec.end(), 0.0));

        return XSec;
    }

    const std::vector<float> neutronsFromFission = {2.88f, 0.0f, 0.0f};

private:
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

    const Fissionable_Isotope Pu239 = {
        .N = 1.0f, .g_J = 3.0f/4.0f, .E_lambda = 2.956243e-1, .gamma_n = 7.947046e-5,
        .gamma_gamma = 3.982423e-2, .gamma_f = 5.619673e-2, .a_c = 9.41e-4
    };

    const Fissionable_Isotope U238 = {
        .N = 0.124954f, .g_J = 1.0f, .E_lambda = 6.674280e0, .gamma_n = 1.4923e-3,
        .gamma_gamma = 2.2711e-2, .gamma_f = 9.88e-9, .a_c = 9.48e-4
    };
};

/**
 * Exports the cross section data into a CSV.
 * @param cross_sections A set of continuous energy macroscopic cross sections.
 * @param filename The filename and path/relative path to save the data.
 * @return A boolean of the function's success.
 */
inline bool ExportCrossSectionsToCSV(const CrossSections *cross_sections, const std::string& filename) {
    constexpr unsigned int fidelity = 1000; //number of points on the plot

    std::vector<std::vector<float>> XSec(fidelity);

    const double delta_ln_e = std::log(ENERGY_MAX / ENERGY_MIN) / (fidelity - 1);
    double ln_e = std::log(ENERGY_MIN) - delta_ln_e;
    double energy{};
    std::vector<double> XSec_split_E(8);
    std::vector<float> XSec_E(5);

    for (unsigned int i = 0; i < fidelity; ++i) {
        ln_e += delta_ln_e;
        energy = std::exp(ln_e);

        XSec_split_E = cross_sections->getCrossSections(energy);
        XSec_E[0] = static_cast<float>(energy);
        XSec_E[1] = static_cast<float>(XSec_split_E[0] + XSec_split_E[3] + XSec_split_E[6]);//scatter
        XSec_E[2] = static_cast<float>(XSec_split_E[1] + XSec_split_E[4]);//capture
        XSec_E[3] = static_cast<float>(XSec_split_E[2] + XSec_split_E[5]);//fission
        XSec_E[4] = static_cast<float>(XSec_split_E[7]);//total

        XSec[i] = XSec_E;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Unable to open file " << filename << "\n";
        return false;
    }
    outfile.clear();

    for (const auto &row : XSec) {
        for (size_t i = 0; i < row.size(); ++i) {
            outfile << row[i];

            if (i < row.size() - 1)
                outfile << ",";

        }
        outfile << "\n";
    }

    outfile.close();

    return true;
}


#endif //NE684A2_CROSSSECTIONS_H
