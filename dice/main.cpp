#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "pcg/pcg_random.hpp"
#include "tclap/CmdLine.h"

using Mutation = std::pair<unsigned long, unsigned long>;
struct Mutation_hash {
    inline std::size_t operator()(const Mutation & m) const {
        return m.first*63+m.second;
    }
};

/* Collects results from threads for use back in sequential-land.
 * Data is protected by a std::mutex */
struct ResultHelper {
    ResultHelper(int size) {
        result.reserve(size);
    }
    void handle_result(std::vector<int> thread_result) {
        std::lock_guard<std::mutex> guard(mtx);
        result.insert(result.end(), thread_result.begin(), thread_result.end());
    }
    std::vector<int> result;
private:
    std::mutex mtx;
};

class mutation_generator {
public:
    mutation_generator(
            unsigned int seed,
            const std::vector<double> &_spectrum,
            const std::vector<unsigned long> &_opps) : engine { seed }
    {
        context_gen = std::discrete_distribution<>(_spectrum.begin(), _spectrum.end());
        for (unsigned long s : _opps) {
            position_gens.push_back(std::uniform_int_distribution<unsigned long>(0, s));
        }
    };

    mutation_generator(
            pcg_extras::seed_seq_from<std::random_device> &seed_source,
            const std::vector<double> &_spectrum,
            const std::vector<unsigned long> &_opps) : engine { seed_source }
    {
        context_gen = std::discrete_distribution<>(_spectrum.begin(), _spectrum.end());
        for (unsigned long s : _opps) {
            position_gens.push_back(std::uniform_int_distribution<unsigned long>(0, s));
        }
    };

    Mutation sample()
    {
        auto ctxt = context_gen(engine);
        auto pos = position_gens[ctxt](engine);
        return Mutation{ctxt, pos};
    }

    int run_single_rep(int n_mut) {
        std::unordered_set<Mutation, Mutation_hash> batch;
        for (int i = 0; i < n_mut; i++) {
            Mutation m = sample();
            if (batch.find(m) != batch.end()) {
                return 1;
            }
            else {
                batch.insert(m);
            }
        }
        return 0;
    }

    void run_reps(int n_reps, int n_mut, std::vector<int> &output) {
        output.reserve(n_reps);
        for (int i = 0; i < n_reps; i++) {
            output.push_back(run_single_rep(n_mut));
        }
    }

private:
    pcg32 engine;
    std::discrete_distribution<> context_gen;
    std::vector<std::uniform_int_distribution<unsigned long>> position_gens;
};


/* Read file line-by-line into vector of strings */
std::vector<std::string> read_file(const std::string &filename) {
    std::ifstream reader(filename);
    if (!reader.is_open()) {
        std::stringstream ss;
        ss << "failed to open " << filename << '\n';
        throw std::runtime_error(ss.str());
    }
    else {
        std::vector<std::string> output;
        std::string line;
        while(getline(reader, line)) {
            output.push_back(line);
        }
        return output;
    }
}


/* Read file and convert vector of strings to vector of doubles */
std::vector<double> read_spectrum(const std::string &filename) {
    std::vector<std::string> input = read_file(filename);
    std::vector<double> output;
    std::transform(input.begin(), input.end(),
                   std::back_inserter(output),
                   [](std::string s) { return std::stod(s); });
    return output;
}


/* Read file and convert vector of strings to vector of unsigned longs */
std::vector<unsigned long> read_opps(const std::string &filename) {
    std::vector<std::string> input = read_file(filename);
    std::vector<unsigned long> output;
    std::transform(input.begin(), input.end(),
                   std::back_inserter(output),
                   [](std::string s) { return (unsigned long) (std::stod(s)); });
    return output;
}


int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();
    try {
        TCLAP::CmdLine cmd("Expected mutational recurrence simulation", ' ', "0.9");

        TCLAP::UnlabeledValueArg<std::string> spectrumArg(
                "spectrum",
                "File containing spectrum (mutation probabilities), 1 per line",
                true, "", "string"
        );

        TCLAP::UnlabeledValueArg<std::string> oppsArg(
                "opps",
                "File containing opps (mutation opportunities), 1 per line",
                true, "", "string"
        );

        TCLAP::ValueArg<int> repsArg(
                "r", "reps",
                "Number of reps to simulate", true,
                1000, "integer"
        );

        TCLAP::ValueArg<unsigned int> threadsArg(
                "t", "threads",
                "Number of threads - default 1", false,
                1, "integer"
        );

        TCLAP::ValueArg<unsigned int> mutationArg(
                "m", "mutations",
                "Number of mutations to simulate in each rep", true,
                1000, "integer"
        );

        TCLAP::ValueArg<int> seedArg(
                "s", "seed",
                "Random number seed", false,
                0, "integer"
        );

        cmd.add( spectrumArg );
        cmd.add( oppsArg );
        cmd.add( repsArg );
        cmd.add( threadsArg );
        cmd.add( mutationArg );
        cmd.add( seedArg );

        // Parse the argv array.
        cmd.parse( argc, argv );

        std::random_device rd;

        // Read files
        std::vector<double> spec;
        std::vector<unsigned long> op;
        try {
            spec = read_spectrum(spectrumArg.getValue());
            op = read_opps(oppsArg.getValue());
        }
        catch (std::exception &ex) {
            std::cerr << "IOError: " << ex.what() << std::endl;
            exit(1);
        }

        auto nthreads = threadsArg.getValue();

        if (nthreads < 1) nthreads = 1;
        if (nthreads > std::thread::hardware_concurrency()) nthreads = std::thread::hardware_concurrency();

        auto reps = repsArg.getValue();
        auto n_mut = mutationArg.getValue();

        std::vector<int> result;
        result.reserve(reps);

        if (nthreads == 1) {
            std::cout << "Using 1 thread" << std::endl;
            // Seed RNG and get result (single thread version)
            pcg_extras::seed_seq_from<std::random_device> seed;
//            unsigned int seed = rd();
            mutation_generator gen(seed, spec, op);
            gen.run_reps(reps, n_mut, result);
        }

        else {
            std::cout << "Using " << nthreads << " threads" << std::endl;

            std::vector<std::thread> threads;
            int reps_per_thread = reps / nthreads;
            int remainder = reps - (reps_per_thread * nthreads);
            ResultHelper helper(reps);  // thread-safe collector of output

            for (int i = 0; i < nthreads; i++) {
                // Make final thread execute remainder of reps, in addition to reps_per_thread
                if (i == nthreads - 1) reps_per_thread = reps_per_thread + remainder;

                // Set up seed_source for thread-local RNG
                // pcg_extras::seed_seq_from<std::random_device> seed;
                unsigned int seed = rd();

                threads.push_back(std::thread( [&helper, seed, &spec, &op, reps_per_thread, n_mut]() {
                    mutation_generator gen(seed, spec, op);
                    std::vector<int> res_t;
                    res_t.reserve(reps_per_thread);
                    gen.run_reps(reps_per_thread, n_mut, res_t);
                    helper.handle_result(res_t);
                }));
            }

            for (auto &thread : threads) {
                thread.join();
            }

            result = helper.result;
        }

        double sum = std::accumulate(result.begin(), result.end(), 0);
        std::cout << "Estimated probability of at least one recurrent mutation = "
                  << sum / result.size() << std::endl;

    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Time taken: " << duration << "ms" << std::endl;

    return 0;
}