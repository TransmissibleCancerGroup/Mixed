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
import core.time: MonoTime;

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

    @Argument("spectrum")
    @Help("File containing spectrum (mutation probabilities, 1 per line)")
    string spectrum;

    @Argument("opps")
    @Help("File containing opps (mutation opportunities, 1 per line)")
    string opps;
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


string[] read_file(string filename) {
    auto file = File(filename, "r");
    return file.byLine().map!(to!string).array();
}


immutable usage = usageString!Options("example");
immutable help = helpString!Options;

int main(string[] args) {
    auto start = MonoTime.currTime;
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
    double[] spectrum = read_file(options.spectrum).map!(to!double).array();
    ulong[] opps = read_file(options.opps).map!(to!double).map!(to!ulong).array();

    auto results = new int[nrep];

    writefln("Using %d threads", options.threads);
    
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
    auto end = MonoTime.currTime;
    auto duration = end - start;
    writeln("Time taken: ", duration.total!"msecs", "ms");
    return 0;
}

