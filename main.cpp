#include <iostream>
#include <iomanip>
#include <vector>
#include <future>
#include <thread>
#include <cmath>
#include <chrono>
#include "Neutron.h"
#include "CrossSections.h"

int main() {
    constexpr unsigned int NUM_PARTICLES = 500000;
    constexpr double temperature = 1200.0; //K

    unsigned int num_cores = std::thread::hardware_concurrency();
    if (num_cores == 0) num_cores = 4; //fallback

    const unsigned int particles_per_core = NUM_PARTICLES / num_cores;

    DopplerCrossSections doppler_cross_sections;
    doppler_cross_sections.LoadCrossSections(temperature);
    {
        std::cout << "Running simulation at " << temperature << "K with " << NUM_PARTICLES << " particles in implicit mode...\n";

        std::vector<Fission_Neutron> fission_bank;
        std::vector<Fission_Neutron> fission_bank_new;
        Global_Tallies tallies;
        tallies.flux.resize(FINE_FLUX_GROUPS);
        tallies.flux_2.resize(FINE_FLUX_GROUPS);

        const auto start_time = std::chrono::high_resolution_clock::now();

        std::vector<std::future<Semilocal_Results>> futures;

        for (int i = 0; i < num_cores; ++i) {
            futures.push_back(std::async(std::launch::async, [&doppler_cross_sections, &fission_bank, particles_per_core]() {
                Semilocal_Results thread_results;
                thread_results.flux.resize(FINE_FLUX_GROUPS);
                thread_results.flux_2.resize(FINE_FLUX_GROUPS);
                for (unsigned int j = 0; j < particles_per_core; ++j) {
                    Local_Results results;
                    results.flux.resize(FINE_FLUX_GROUPS);
                    runNeutron(&doppler_cross_sections, &fission_bank, &results);

                    thread_results.N += results.N;
                    thread_results.k_inf += results.k_inf;
                    thread_results.k_inf_2 += std::pow(results.k_inf, 2.0f);
                    for (unsigned int k = 0; k < thread_results.group_flux.size(); ++k) {
                        thread_results.group_flux[k] += results.group_flux[k];
                        thread_results.group_flux_2[k] += std::pow(results.group_flux[k], 2.0);
                        thread_results.capture_rr[k] += results.capture_rr[k];
                        thread_results.capture_rr_2[k] += std::pow(results.capture_rr[k], 2.0);
                        thread_results.fission_rr[k] += results.fission_rr[k];
                        thread_results.fission_rr_2[k] += std::pow(results.fission_rr[k], 2.0);
                    }
                    for (unsigned int k = 0; k < thread_results.scatter_rr.size(); ++k) {
                        thread_results.scatter_rr[k] += results.scatter_rr[k];
                        thread_results.scatter_rr_2[k] += std::pow(results.scatter_rr[k], 2.0);
                    }
                    for (unsigned int k = 0; k < FINE_FLUX_GROUPS; ++k) {
                        thread_results.flux[k] += results.flux[k];
                        thread_results.flux_2[k] += std::pow(results.flux[k], 2.0);
                    }
                }
                return thread_results;
            }));
        }

        for (auto& f : futures) {
            const Semilocal_Results results = f.get();

            tallies.N += results.N;
            tallies.k_inf += results.k_inf;
            tallies.k_inf_2 += results.k_inf_2;
            for (unsigned int k = 0; k < tallies.group_flux.size(); ++k) {
                tallies.group_flux[k] += results.group_flux[k];
                tallies.group_flux_2[k] += results.group_flux_2[k];
                tallies.capture_rr[k] += results.capture_rr[k];
                tallies.capture_rr_2[k] += results.capture_rr_2[k];
                tallies.fission_rr[k] += results.fission_rr[k];
                tallies.fission_rr_2[k] += results.fission_rr_2[k];
            }
            for (unsigned int k = 0; k < tallies.scatter_rr.size(); ++k) {
                tallies.scatter_rr[k] += results.scatter_rr[k];
                tallies.scatter_rr_2[k] += results.scatter_rr_2[k];
            }
            for (unsigned int k = 0; k < FINE_FLUX_GROUPS; ++k) {
                tallies.flux[k] += results.flux[k];
                tallies.flux_2[k] += results.flux_2[k];
            }
        }

        fission_bank = fission_bank_new;
        fission_bank_new.clear();

        const auto end_time = std::chrono::high_resolution_clock::now();
        const auto delta_time = end_time - start_time;
        const auto dt = static_cast<double>(delta_time.count());

        const auto crit_var = static_cast<float>(std::abs(std::pow(tallies.k_inf / tallies.N, 2.0) - tallies.k_inf_2 / tallies.N) / (tallies.N - 1));
        std::cout << std::fixed << std::setprecision(5) << "k_inf = " << tallies.k_inf / tallies.N << " +/- " << std::sqrt(crit_var) << "\n";
        std::cout << "Elapsed time: " << dt / 1e6 << " ms \n";
        std::cout << "FOM: " << 1.0 / (dt / 1e9 * crit_var) << "\n";

        // if (!exportTallies(&tallies, std::format("finegroup_flux_{}K.csv", temperature))) {
        //     return 1;
        // }
    }

    return 0;
}
