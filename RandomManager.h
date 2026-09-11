//
// Created by Trevor on 8/28/2026.
//

#pragma once
#include <random>
#include <thread>

#ifndef NE684A2_RANDOMMANAGER_H
#define NE684A2_RANDOMMANAGER_H

/**
 * @brief Handles the generation of random floats.
 * @author Google Gemini 3.7 Flash
 */
class RandomManager {
public:
    /**
     * @brief Generates a random value from [0.0, 1.0). Thread local.
     * @return Random double in [0.0, 1.0).
     */
    static double getRandomFrac() {
        //thread_local is a static variable locked per thread, meaning it is only run the first time the function is
        //called per thread
        thread_local std::mt19937 engine(genUniqueSeed());
        thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(engine);
    }

private:
    /**
     * @brief I (Trevor) am not too sure what this does or why it works.
     * @return A thread-unique seed for a random number generator.
     */
    static unsigned int genUniqueSeed() {
        std::random_device rd;
        const size_t threadId = std::hash<std::thread::id>{}(std::this_thread::get_id());
        return rd() ^ static_cast<unsigned int>(threadId);
    }
};


#endif //NE684A2_RANDOMMANAGER_H
