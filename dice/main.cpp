#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <thread>
#include <unordered_set>

#include "pcg/pcg_random.hpp"
#include "tclap/CmdLine.h"

using Mutation = std::pair<unsigned long, unsigned long>;
struct Mutation_hash {
    inline std::size_t operator()(const Mutation & m) const {
        return m.first*63+m.second;
    }
};


struct RepResult {
    RepResult(std::vector<int> &totals, std::vector<int> &duplicates) : tot_counts(totals), dup_counts(duplicates) {}
    std::vector<int> tot_counts;
    std::vector<int> dup_counts;
};


/* Collects results from threads for use back in sequential-land.
 * Data is protected by a std::mutex */
struct ResultHelper {
    ResultHelper(int size) {
        result.reserve(size);
    }
    void handle_result(std::vector<RepResult> thread_result) {
        std::lock_guard<std::mutex> guard(mtx);
        result.insert(result.end(), thread_result.begin(), thread_result.end());
    }
    std::vector<RepResult> result;
private:
    std::mutex mtx;
};


class ProgressBar {
public:
    ProgressBar(size_t total) : total(total)  {

    }
    void update(size_t done) {
        std::lock_guard<std::mutex> guard(mtx);
        total_done += done;
        std::cerr << "\33[2K\r[";
        int prop_done = ((linewidth - 2) * total_done) / total;
        for (int i = 0; i < prop_done; i++) {
            std::cerr << "#";
        }
        for (int i = 0; i < (linewidth - 2 - prop_done); i++) {
            std::cerr << " ";
        }
        std::cerr << "] " << total_done << " / " << total << std::flush;
        //if (total_done >= total) std::cerr << std::endl;
    }
private:
    size_t total;
    size_t total_done = 0;
    size_t linewidth = 80;
    std::mutex mtx;
};


class mutation_generator {
public:
    mutation_generator(
            unsigned int seed,
            const std::vector<double> &_spectrum,
            const std::vector<unsigned long> &_opps) : rng { seed }, size { _spectrum.size() }
    {
        if (_spectrum.size() != _opps.size()) {
            throw std::runtime_error("Spectrum and Opportunities must be same size");
        }
        context_generator = std::discrete_distribution<>(_spectrum.begin(), _spectrum.end());
        for (unsigned long s : _opps) {
            position_generators.push_back(std::uniform_int_distribution<unsigned long>(0, s));
        }
    };

    mutation_generator(
            pcg_extras::seed_seq_from<std::random_device> &seed_source,
            const std::vector<double> &_spectrum,
            const std::vector<unsigned long> &_opps) : rng { seed_source }, size { _spectrum.size() }
    {
        context_generator = std::discrete_distribution<>(_spectrum.begin(), _spectrum.end());
        for (unsigned long s : _opps) {
            position_generators.push_back(std::uniform_int_distribution<unsigned long>(0, s));
        }
    };

    Mutation sample()
    {
        auto ctxt = context_generator(rng);
        auto pos = position_generators[ctxt](rng);
        return Mutation{ctxt, pos};
    }

    RepResult run_single_rep(int n_mut) {
        std::unordered_set<Mutation, Mutation_hash> seen;
        std::vector<int> dup_counts(size, 0);
        std::vector<int> tot_counts(size, 0);

        for (int i = 0; i < n_mut; i++) {
            Mutation m = sample();
            tot_counts[m.first]++;

            if (seen.find(m) != seen.end()) {
                // Duplicate mutation!
                dup_counts[m.first]++;
            }

            else {
                seen.insert(m);
            }
        }

        return RepResult(tot_counts, dup_counts);
    }

    void run_reps(int n_reps, int n_mut, std::vector<RepResult> &output, ProgressBar &pbar) {
        output.reserve(n_reps);
        for (int i = 0; i < n_reps; i++) {
            output.push_back(run_single_rep(n_mut));
            pbar.update(1);
        }
    }

private:
    pcg32 rng;
    size_t size;
    std::discrete_distribution<> context_generator;
    std::vector<std::uniform_int_distribution<unsigned long>> position_generators;
};


void write_results(std::ofstream &reps_out, const RepResult &rep);

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


void write_csv(std::ofstream &outstream, const std::vector<int> &counts) {
    for (int i = 0; i < counts.size() - 1; i++) {
        auto val = counts[i];
        outstream << val << ",";
    }
    outstream << counts.back() << std::endl;
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

        TCLAP::UnlabeledValueArg<std::string> outpathArg(
                "outpath",
                "Directory to write result files (must already exist). Result files will be <path>/repeat_mutations.csv and <path>/total_mutations.csv",
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
        cmd.add( outpathArg );

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

        // Set up out files
        std::stringstream reps_filename;
        reps_filename << outpathArg.getValue() << "/repeat_mutations.csv";
        std::ofstream reps_out(reps_filename.str());
        if (!reps_out.is_open()) {
            std::stringstream ss;
            ss << "failed to open " << reps_filename.str() << '\n';
            throw std::runtime_error(ss.str());
        }

        std::stringstream totals_filename;
        totals_filename << outpathArg.getValue() << "/total_mutations.csv";
        std::ofstream totals_out(totals_filename.str());
        if (!totals_out.is_open()) {
            std::stringstream ss;
            ss << "failed to open " << totals_filename.str() << '\n';
            throw std::runtime_error(ss.str());
        }

        int ncat = spec.size();

        auto nthreads = threadsArg.getValue();

        if (nthreads < 1) nthreads = 1;
        if (nthreads > std::thread::hardware_concurrency()) nthreads = std::thread::hardware_concurrency();

        auto reps = repsArg.getValue();
        auto n_mut = mutationArg.getValue();

        std::vector<RepResult> result;

        ProgressBar pbar(reps);

        if (nthreads == 1) {
            std::cout << "Using 1 thread" << std::endl;
            // Seed RNG and get result (single thread version)
            pcg_extras::seed_seq_from<std::random_device> seed;
            mutation_generator gen(seed, spec, op);
            result.reserve(reps);
            gen.run_reps(reps, n_mut, result, pbar);
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

                threads.push_back(std::thread( [&helper, seed, &spec, &op, &pbar, reps_per_thread, n_mut]() {
                    mutation_generator gen(seed, spec, op);
                    std::vector<RepResult> res_t;
                    res_t.reserve(reps_per_thread);
                    gen.run_reps(reps_per_thread, n_mut, res_t, pbar);
                    helper.handle_result(res_t);
                }));
            }

            for (auto &thread : threads) {
                thread.join();
            }

            result = helper.result;
        }

        // Print some quick stats
        std::cout << std::endl;
        std::cout << "Expected counts of duplicate mutations in each category:" << std::endl;
        std::cout << "(R equivalent ` colMeans(dups) `)" << std::endl;
        std::vector<double> counts(ncat, 0.0);
        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < ncat; j++) {
                counts[j] += (double) result[i].dup_counts[j];
            }
        }
        for (auto &val : counts) {
            val /= reps;
            std::cout << val << '\t';
        }
        std::cout << std::endl;

        std::cout << "Estimated conditional probability of seeing duplicate mutations in each category (conditioned on category):" << std::endl;
        std::cout << "(R equivalent ` colMeans(dups / totals) `)" << std::endl;
        std::vector<double> cond_probs(ncat, 0.0);
        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < ncat; j++) {
                if (result[i].tot_counts[j] > 0) {
                    cond_probs[j] += (double) result[i].dup_counts[j] / (double) result[i].tot_counts[j];
                }
            }
        }
        for (auto &val : cond_probs) {
            val /= reps;
            std::cout << val << '\t';
        }
        std::cout << std::endl;

        std::cout << "Estimated overall probability of seeing duplicate mutations in each category:" << std::endl;
        std::cout << "(R equivalent ` colMeans(dups / sum(totals[1, ])) `)" << std::endl;
        std::vector<double> probs(ncat, 0.0);

        for (int i = 0; i < reps; i++) {
            for (int j = 0; j < result[i].dup_counts.size(); j++) {
                probs[j] += (double) result[i].dup_counts[j] / (double) n_mut;
            }
        }

        for (auto &val : probs) {
            val /= reps;
            std::cout << val << '\t';
        }
        std::cout << std::endl;

        std::cout << "Estimated probability of at least one recurrent mutation:" << std::endl;
        std::cout << "(R equivalent ` sum(rowSums(dups) > 0) / nrow(dups) `)" << std::endl;
        int seen = 0;
        for (auto &rep_result : result) {
            if (std::accumulate(rep_result.dup_counts.begin(), rep_result.dup_counts.end(), 0) > 0) {
                seen++;
            }
        }
        std::cout << (double) seen / (double) reps << std::endl;

        for (auto &rep : result) {
            write_csv(reps_out, rep.dup_counts);
            write_csv(totals_out, rep.tot_counts);
        }
        reps_out.close();
        totals_out.close();

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
