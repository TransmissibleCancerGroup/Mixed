import std.algorithm;
import std.array;
import std.conv: to;
import std.functional;
import std.parallelism;
import std.random: dice, uniform, Random, unpredictableSeed;
import std.range;
import std.stdio;
import std.traits;
import std.typecons;

import darg;


struct Options
{
    @Option("help", "h")
    @Help("Prints this help.")
    OptionFlag help;

    @Option("reps", "r")
    @Help("Number of simulation repetitions to run")
    size_t reps = 1000;

    @Option("threads", "t")
    @Help("Number of threads - default 1")
    size_t threads = 1;

    @Option("mutations", "m")
    @Help("Number of mutations per rep")
    size_t mutations = 100;

    @Option("seed", "s")
    @Help("Random number seed")
    int seed;
}

alias Mutation = Tuple!(ulong, "context", ulong, "position");


int mutate(const ulong _mut,
           const ref double[] _spectrum,
           const ref ulong[] _opps)
{
    bool[Mutation] results;
    for (int j; j < _mut; j++) {
        int context = to!int(dice(_spectrum));
        ulong opp = _opps[context];
        auto position = uniform(0, opp);
        Mutation mutation;
        mutation[0] = context;
        mutation[1] = position;
        if (mutation in results) {
            return 1;
        }
        else {
            results[mutation] = true;
        }
    }
    return 0;
}


immutable usage = usageString!Options("example");
immutable help = helpString!Options;

int main(string[] args) {
    Options options;

    try
    {
        options = parseArgs!Options(args[1 .. $]);
    }
    catch (ArgParseError e)
    {
        writeln(e.msg);
        writeln(usage);
        return 1;
    }
    catch (ArgParseHelp e)
    {
        // Help was requested
        writeln(usage);
        write(help);
        return 0;
    }

    if (options.seed > 0) Random(unpredictableSeed);
    else Random(options.seed);

    if (options.threads < 1) options.threads = 1;
    if (options.threads > totalCPUs) options.threads = totalCPUs;

    auto pool = new TaskPool(options.threads);
    scope(exit) pool.stop();

    ulong nrep = to!ulong(options.reps);
    ulong n_mutations = to!ulong(options.mutations);
    // auto spectrum = array([0.4, 0.6]);
    // auto opps = array([750_000_000,250_000_000]);
    double[] spectrum = array([1.109833e-02, 9.149341e-03, 1.490070e-03, 6.233885e-03, 6.595870e-03, 7.342368e-03, 8.928404e-04, 7.186582e-03, 8.232604e-03, 5.758021e-03,
              6.163352e-04, 4.459080e-03, 1.225006e-02, 1.116223e-02, 2.275496e-03, 1.525910e-02, 1.801068e-03, 2.580909e-03, 5.925480e-04, 2.963986e-03,
              1.284983e-03, 7.021348e-04, 5.062896e-04, 1.381543e-03, 6.021227e-04, 2.393352e-03, 2.485340e-07, 8.900807e-04, 1.874853e-03, 2.067419e-03,
              3.048970e-04, 3.151574e-03, 2.951453e-02, 1.432275e-02, 1.716469e-01, 1.262376e-02, 2.089645e-02, 1.850170e-02, 9.557722e-02, 1.711331e-02,
              2.494381e-02, 2.716149e-02, 1.035708e-01, 1.768985e-02, 1.449210e-02, 1.768078e-02, 7.600222e-02, 1.376170e-02, 4.021520e-03, 2.371144e-03,
              2.810910e-03, 8.360909e-03, 1.182587e-03, 1.903167e-03, 1.487961e-03, 2.179344e-03, 6.892894e-04, 5.524095e-04, 1.200229e-03, 2.107137e-03,
              5.600155e-03, 1.999079e-03, 1.090066e-03, 3.981023e-03, 1.391577e-02, 6.274961e-03, 1.013764e-02, 9.256316e-03, 4.176675e-03, 5.252593e-03,
              7.013225e-03, 6.713813e-03, 1.124784e-02, 6.999724e-03, 4.977593e-03, 1.066741e-02, 8.073616e-03, 4.857381e-03, 8.325454e-03, 6.257106e-03,
              1.587636e-03, 1.784091e-03, 1.385831e-03, 3.158539e-03, 3.026912e-04, 2.098502e-03, 1.599549e-03, 2.758538e-03, 9.904500e-05, 2.023656e-04,
              1.188353e-03, 8.007233e-04, 1.397554e-03, 1.291737e-03, 2.031077e-03, 4.030128e-03]);
    ulong[] opps = array([1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, 1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, 8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, 1.11e+08, 8.75e+07, 1.25e+07,
                   1.25e+08, 1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, 1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, 8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, 1.11e+08, 8.75e+07,
                   1.25e+07, 1.25e+08, 1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, 1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, 8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, 1.11e+08,
                   8.75e+07, 1.25e+07, 1.25e+08, 1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, 7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, 6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07,
                   1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08, 1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, 7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, 6.43e+07, 5.36e+07, 8.52e+07,
                   8.27e+07, 1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08, 1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, 7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, 6.43e+07, 5.36e+07,
                   8.52e+07, 8.27e+07, 1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08]).map!(to!ulong).array();

    auto results = new int[nrep];

    if (options.threads > 1) {
        foreach(ref elem; pool.parallel(results)) {
            elem = mutate(n_mutations, spectrum, opps);
        }
    }

    else {
        foreach(ref elem; results) {
            elem = mutate(n_mutations, spectrum, opps);
        }
    }

    double prob = to!double(results.sum()) / nrep;
    writefln("Estimated probability of at least one recurrent mutation = %f", prob);
    return 0;
}

