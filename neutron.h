//
// Created by Trevor on 8/28/2026.
//

#pragma once
#import <iostream>
#import <cmath>
#import <vector>
#import "RandomManager.h"
#import "CrossSections.h"

#ifndef NE684A2_NEUTRON_H
#define NE684A2_NEUTRON_H
// k_eff optimal: 1e-3, 1
// general: 1e-5?, 10?
#define HYPER_WEIGHT_THRESH 1e-3
#define HYPER_ROUNDS INFINITY

#define FINE_FLUX_GROUPS 1000

/**
 * @brief Holds the tally data lumped from every thread.
 */
struct Global_Tallies {
    float N = 0.0f;
    float k_inf = 0.0f;
    float k_inf_2 = 0.0f;
    std::vector<double> group_flux = {0.0f, 0.0f, 0.0f};
    std::vector<double> group_flux_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr_2 = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> flux;
    std::vector<double> flux_2;
};

//Fission table will be a vector of pointers to Fission_Neutron(s)
//for asynchronicity, a neutron will compute its contribution to every tally then apply said contribution itself.
/**
 * @brief Holds the necessary information to store a neutron in the fission bank.
 */
struct Fission_Neutron {
    float weight = -1.0f;
    // This would also hold position, but in infinite homogenous, the transport problem simplifies to 0D
};

/**
 * @brief Holds the tally information collected per-thread.
 */
struct Semilocal_Results {
    float N = 0.0f;
    float k_inf = 0.0f;
    float k_inf_2 = 0.0f;
    std::vector<double> group_flux = {0.0f, 0.0f, 0.0f};
    std::vector<double> group_flux_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr_2 = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> flux;
    std::vector<double> flux_2;
};

/**
 * @brief Holds the tally contribution from a single neutron.
 */
struct Local_Results {
    float N = 0.0f;
    float k_inf = 0.0f;
    std::vector<double> group_flux = {0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> flux;
};

/**
 * @breif Neutron object used for Monte Carlo simulation.
 */
class neutron {
public:
    /**
     * @brief Creates a neutron object from a banked neutron.
     * @param n A neutron from the fission bank.
     */
    explicit neutron(const Fission_Neutron n) {
        weight = n.weight;
        energy = RandomManager::getRandomFrac() * (ENERGY_MAX - ENERGY_MIN) +ENERGY_MIN;
    }

    /**
     * @breif Simulates this neutron via Implicit Capture with Russian Roulette.
     * @param crossSections A set of continuous energy macroscopic cross sections.
     * @param results A mutable struct of this neutrons tally contributions.
     */
    void simulate(const CrossSections *crossSections, Local_Results *results) {
        bool alive = true;
        std::vector<double> XSec(8);
        std::vector<double> XSec_sums(3);
        const std::vector<float> atomic_weights = crossSections->targetMasses;
        double avg_nu = 0.0;
        double scatter;
        double energy_min = 0.0;
        unsigned int g{};
        while (alive) {
            XSec = crossSections->getCrossSections(energy);
            XSec_sums = {
                XSec[0] + XSec[3] + XSec[6],
                XSec[1] + XSec[4],
                XSec[2] + XSec[5],
            };
            avg_nu = (crossSections->neutronsFromFission[0] * XSec[2] +
                crossSections->neutronsFromFission[0] * XSec[5]) / XSec_sums[2];

            g = energy_to_group();

            //collision estimator of flux
            collision_flux_estimator(g, results, XSec[7]);

            //Fission Contribution
            results->k_inf += static_cast<float>(weight * avg_nu * XSec_sums[2] / XSec[7]);
            results->fission_rr[g] += weight * XSec_sums[2] / XSec[7];

            //Capture Contribution
            results->capture_rr[g] += weight * XSec_sums[1] / XSec[7];

            //Scatter Contribution
            scatter = RandomManager::getRandomFrac() * XSec_sums[0];
            for (int i = 0; i < 3; ++i) {
                if (scatter <= XSec[3 * i]) {
                    energy_min = std::pow((crossSections->targetMasses[i] - 1) / (crossSections->targetMasses[i] + 1), 2.0) * energy;;
                    break;
                }
                scatter -= XSec[3 * i];
            }

            new_energy = RandomManager::getRandomFrac() * (energy - energy_min) + energy_min;
            results->scatter_rr[scattering_index(g, results->scatter_rr.size())] += weight * XSec_sums[0] / XSec[7];

            weight *= XSec_sums[0] / XSec[7];
            energy = new_energy;

            //Russian Roulette
            if (weight <= HYPER_WEIGHT_THRESH) {
                if (1.0 - RandomManager::getRandomFrac() >= 1.0 / HYPER_ROUNDS) {
                    alive = false;
                }
                // neutron lives or is dead and weight doesn't matter
                weight *= HYPER_ROUNDS;
            }
        }
    }

private:
    double energy;

    double new_energy{};
    double weight;

    const std::vector<double> energy_groups = {1e2, 1.0, ENERGY_MIN};

    /**
     * @brief Returns the current energy group of this neutron.
     * @return The index of the current energy group
     */
    [[nodiscard]] unsigned int energy_to_group() const {
        unsigned int g{};
        for (g = 0; g < energy_groups.size(); ++g) {
            if (energy >= energy_groups[g]) break;
        }
        return g;
    }

    /**
     * @brief Converts the to and from energy into an index into the scattering matrix.
     * @return The index into the flattened scattering matrix.
     */
    [[nodiscard]] unsigned int scattering_index(const unsigned int from_group, const size_t scatter_size) const {
        int g{};
        for (g = 0; g < energy_groups.size(); ++g) {
            if (new_energy >= energy_groups[g]) break;
        }
        return from_group * (scatter_size - from_group + 1) / 2 + (g - from_group);
    }

    /**
     * @brief Estimates the flux contribution of this neutron via the collision estimator.
     * @param from_group The energy group the collision happens in.
     * @param results A mutable struct of this neutrons tally contributions.
     * @param total_XSec The total macroscopic cross section.
     */
    void collision_flux_estimator(const unsigned int from_group, Local_Results *results, const double total_XSec) const {
        const auto fine_energy = static_cast<unsigned int>(FINE_FLUX_GROUPS *
            std::log(ENERGY_MAX / energy) / std::log(ENERGY_MAX / ENERGY_MIN));
        results->group_flux[from_group] += weight / total_XSec;
        results->flux[fine_energy] += weight / total_XSec;
    }
};

/**
 * @brief Handles the entire lifetime of a single neutron. Simulates with Implicit Capture.
 * @param crossSections A set of continuous energy macroscopic cross sections.
 * @param fission_bank
 * @param results A mutable struct of this neutrons tally contributions.
 */
inline void runNeutron(const CrossSections *crossSections, const std::vector<Fission_Neutron> *fission_bank, Local_Results *results) {
    Fission_Neutron banked_neutron{.weight = 1.0f};
    if (!fission_bank->empty()) {
        const int index = static_cast<int>(RandomManager::getRandomFrac() * static_cast<double>(fission_bank->size()));
        banked_neutron = fission_bank->at(index);
    }

    results->N += banked_neutron.weight;

    neutron n(banked_neutron);
    n.simulate(crossSections, results);
}

/**
 * @breif Exports tally data into a CSV.
 * @param tallies The tally data to export.
 * @param filename The filename and path/relative path to save the data.
 * @return A boolean of the function's success.
 */
inline bool exportTallies(Global_Tallies *tallies, const std::string& filename) {
    std::cout <<std::defaultfloat;

    std::cout << "\nFission Cross Sections:\n";
    for (unsigned int g = 0; g < tallies->group_flux.size(); ++g) {
        const auto fiss_xsec = tallies->fission_rr[g] / tallies->group_flux[g];
        const auto fiss_var = static_cast<float>(std::abs(std::pow(tallies->fission_rr[g] / tallies->group_flux[g], 2.0) - tallies->fission_rr_2[g] / tallies->group_flux[g]) / (tallies->N - 1));
        std::cout << g + 1 << ":\t" << fiss_xsec << " +/- " << std::sqrt(fiss_var) << "\n";
    }
    std::cout << "\nCapture Cross Sections:\n";
    for (unsigned int g = 0; g < tallies->group_flux.size(); ++g) {
        const auto cap_xsec = tallies->capture_rr[g] / tallies->group_flux[g];
        const auto cap_var = static_cast<float>(std::abs(std::pow(tallies->capture_rr[g] / tallies->group_flux[g], 2.0) - tallies->capture_rr_2[g] / tallies->group_flux[g]) / (tallies->N - 1));
        std::cout << g + 1 << ":\t" << cap_xsec << " +/- " << std::sqrt(cap_var) << "\n";
    }
    std::cout << "\nScatter Cross Sections:\n";
    for (int g = 0; g < tallies->scatter_rr.size(); ++g) {
        int e_g = 0;
        if (g > 2)
            e_g = 1;
        if (g > 4)
            e_g = 2;

        const auto scat_xsec = tallies->scatter_rr[g] / tallies->group_flux[e_g];
        const auto scat_var = static_cast<float>(std::abs(std::pow(tallies->scatter_rr[g] / tallies->group_flux[e_g], 2.0) - tallies->scatter_rr_2[g] / tallies->group_flux[e_g]) / (tallies->N - 1));

        std::cout << e_g + 1 << "->" << g - e_g * (7 - e_g) / 2 + e_g + 1 << ":\t" << scat_xsec << " +/- " << std::sqrt(scat_var) << "\n";
    }
    std::cout << "\nFew-Group Flux:\n";
    const std::vector<double> energy_groups = {ENERGY_MAX, 1e2, 1.0, ENERGY_MIN};
    for (int g = 0; g < tallies->group_flux.size(); ++g) {
        const auto flux_var = static_cast<float>(std::abs(std::pow(tallies->group_flux[g] / tallies->N, 2.0) - tallies->group_flux_2[g] / tallies->N) / (tallies->N - 1));

        std::cout << g + 1 << ":\t" << tallies->group_flux[g] / tallies->N / (energy_groups[g] - energy_groups[g + 1]) << " +/- " << std::sqrt(flux_var) / (energy_groups[g] - energy_groups[g + 1]) << "\n";
    }

    std::vector<std::vector<float>> flux(FINE_FLUX_GROUPS);
    std::vector<float> single_group_flux(3);
    const double dlogE = std::log(ENERGY_MAX / ENERGY_MIN) / FINE_FLUX_GROUPS;
    for (int g = 0; g < FINE_FLUX_GROUPS; ++g) {
        const double logE = std::log(ENERGY_MAX) - g * dlogE;
        const double delta_energy = std::exp(logE) - std::exp(logE - dlogE);
        const double flux_var = std::abs(std::pow(tallies->flux[g] / tallies->N, 2.0) - tallies->flux_2[g] / tallies->N) / (tallies->N - 1);
        single_group_flux[0] = static_cast<float>(std::exp(logE));
        single_group_flux[1] = static_cast<float>(tallies->flux[g] / tallies->N / delta_energy);
        single_group_flux[2] = static_cast<float>(std::sqrt(flux_var) / delta_energy);

        flux[g] = single_group_flux;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Unable to open file " << filename << "\n";
        return false;
    }
    outfile.clear();

    for (const auto &row : flux) {
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

#endif //NE684A2_NEUTRON_H
