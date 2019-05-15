
#include <CLI11.hpp>
#include "fasttree/Options.h"
#include "fasttree/FastTree.h"
#include "fasttree/Utils.h"
#include <cmath>

std::string isNotEmpty(const std::string &input) {
    if (input.size() == 0) {
        return "Value is empty";
    }
    return std::string();
}

struct Min : CLI::Validator {

    template<typename T>
    Min(T min, bool equals = true) {
        func = [min, equals](std::string input) {
            T val;
            CLI::detail::lexical_cast(input, val);
            if (val < min)
                return "Min value " + std::to_string(min);

            if (!equals && val == min)
                return "Min value greater than " + std::to_string(min);

            return std::string();
        };
    }

};

void setDeprecated(CLI::Option *op) {
    op->check([&op](const std::string &in) {
        std::cerr << "Warning: " << op->get_name() << " is a deprecated option and it has no effect" << std::endl;
        return "";
    });
}

void cli(CLI::App &app, std::string &name, std::string &version, std::string &flags, fasttree::Options &options,
         std::vector<std::string> &args) {
    std::stringstream description;
    std::string padd = std::string(name.size() + 1, ' ');
    description << name << " " << version << " " << flags << std::endl;
    description << name << " [-nt] [-n 100] [-quote] [-pseudo | -pseudo 1.0]" << std::endl;
    description << padd << "[-boot 1000 | -nosupport]" << std::endl;
    description << padd << "[-intree starting_trees_file | -intree1 starting_tree_file]" << std::endl;
    description << padd << "[-quiet | -nopr]" << std::endl;
    description << padd << "[-nni 10] [-spr 2] [-noml | -mllen | -mlnni 10]" << std::endl;
    description << padd << "[-mlacc 2] [-cat 20 | -nocat] [-gamma]" << std::endl;
    description << padd << "[-slow | -fastest] [-2nd | -no2nd] [-slownni] [-seed 1253]" << std::endl;
    description << padd << "[-top | -notop] [-topm 1.0 [-close 0.75] [-refresh 0.8]]" << std::endl;
    description << padd << "[-matrix Matrix | -nomatrix] [-nj | -bionj]" << std::endl;
    description << padd << "[-lg] [-wag] [-nt] [-gtr] [-gtrrates ac ag at cg ct gt] [-gtrfreq A C G T]" << std::endl;
    description << padd << "[ -constraints constraintAlignment [ -constraintWeight 100.0 ] ]" << std::endl;
    description << padd << "[-log logfile]" << std::endl;
    description << padd << "[ alignment_file ]" << std::endl;
    description << padd << "[ -out output_newick_file | > newick_tree]" << std::endl;
    description << std::endl;
    description << name << " [-nt] [-matrix Matrix | -nomatrix] [-rawdist] -makematrix [alignment]" << std::endl;
    description << padd << "[-n 100] > phylip_distance_matrix" << std::endl;
    description << std::endl;
    description << "  FastTree supports fasta or phylip interleaved alignments" << std::endl;
    description << "  By default FastTree expects protein alignments,  use -nt for nucleotides" << std::endl;
    description << "  FastTree reads standard input if no alignment file is given" << std::endl;

    app.name(name);
    app.description(description.str());
    app.set_help_flag("-h,-help,--help", "Print this help message and exit")->group("");
    app.add_option("protein_alignment", options.inFileName, "input file instead of stdin")->
            type_name("")->check(isNotEmpty);


    auto io = "Input/output options";

    app.add_option("-out", options.outFileName,
                   "print tree in output file instead of stdout")->type_name("file")->
            check(isNotEmpty)->group(io);

    app.add_option("-n", options.nAlign,
                   "read in multiple alignments in. This only works with phylip interleaved format. For example,"
                   " you can use it with the output from phylip's seqboot. If you use -n, FastTree will write 1"
                   " tree per line to standard output.")->type_name("<number>")->check(Min(1))->group(io);

    app.add_flag_function("-nt", [&options](size_t) {
                              options.nCodes = 4;
                          },
                          "nucleotides instead of protein alignments")->group(io);

    app.add_option("-intree", options.intreeFile,
                   "read the starting tree in from newickfile. Any branch lengths in the starting trees are ignored."
                   " -intree with -n will read a separate starting tree for each alignment.")->
            type_name("newick_file")->check(isNotEmpty)->group(io);

    app.add_option("-intree1", options.intreeFile, "read the same starting tree for each alignment")->
            type_name("newick_file")->check(isNotEmpty)->check([&options](const std::string &) {
        options.intree1 = true;
        return "";
    })->group(io);

    app.add_option("-verbose", options.verbose, "level of details during normal operation")->type_name("lvl")->
            group(io);

    app.add_flag_function("-quiet", [&options](size_t) {
                              options.verbose = 0;
                              options.showProgress = false;
                          },
                          "do not write to standard error during normal operation (no progress indicator, no options"
                          " summary, no likelihood values, etc.)")->group(io);

    app.add_flag_function("-nopr", [&options](size_t) {
        options.showProgress = false;
    }, "do not write the progress indicator to stderr")->group(io);

    app.add_option("-log", options.logFileName,
                   "save intermediate trees so you can extract the trees and restart long-running jobs if they crash"
                   " -log also reports the per-site rates (1 means slowest category)")->type_name("logfile")->
            check(isNotEmpty)->group(io);

    app.add_flag("-safelog", options.safeLog, "to write log directly without using buffer, useful when the application"
                                              " instantly crash")->group(io);

    app.add_flag_function("-quote", [&options](size_t) {
                              options.bQuote = true;
                          },
                          "quote sequence names in the output and allow spaces, commas, parentheses, and colons in them"
                          " but not ' characters (fasta files only)")->group(io);


    std::stringstream dist_description;
    dist_description << "Distances:" << std::endl;
    dist_description << "  Default: For protein sequences, log-corrected distances and an" << std::endl;
    dist_description << "    or for nucleotide sequences, Jukes-Cantor distances" << std::endl;

    auto dist = dist_description.str();

    app.add_flag_function("-makematrix", [&options](size_t) { options.makeMatrix = true; },
                          "print distance matrix")->group(dist);

    app.add_flag_function("-rawdist", [&options](size_t) { options.logdist = false; }, "to turn the log-correction off")
            ->group(dist);

    app.add_option("-matrix", options.matrixPrefix, "to specify a different matrix")->type_name("file")->
            check(isNotEmpty)->group(dist);

    app.add_flag_function("-nomatrix", [&options](size_t) { options.useMatrix = false; },
                          "use %different as the uncorrected distance")->type_name("file")->check(
            isNotEmpty)->group(dist);


    app.add_option("-pseudo", [&args, &options](CLI::results_t in) {
                       if (std::find(args.begin(), args.end(), "-pseudo") != args.end()) {
                           if (!CLI::detail::lexical_cast(in[0], options.pseudoWeight)) {
                               throw CLI::ConversionError(in[0], "-pseudo");
                           }
                       }
                       return true;
                   },
                   "Use pseudocounts to estimate distances between sequences with little or no overlap. "
                   "(Off by default.) Recommended if analyzing the alignment has sequences with little or no overlap."
                   " If the weight is not specified, it is 1.0")->add_result("1.0")->
            multi_option_policy(CLI::MultiOptionPolicy::TakeLast)->check(Min(0.0))->group(dist);


    std::stringstream topo_description;
    topo_description << "Topology refinement:" << std::endl;
    topo_description << "  By default, FastTree tries to improve the tree with up to 4*log2(N)" << std::endl;
    topo_description << "  rounds of minimum-evolution nearest-neighbor interchanges (NNI)," << std::endl;
    topo_description << "  where N is the number of unique sequences, 2 rounds of" << std::endl;
    topo_description << "  subtree-prune-regraft (SPR) moves (also min. evo.), and" << std::endl;
    topo_description << "  up to 2*log(N) rounds of maximum-likelihood NNIs." << std::endl;

    auto topo = topo_description.str();

    app.add_option("-nni", [&options](CLI::results_t in) {
        if (!CLI::detail::lexical_cast(in[0], options.nni)) {
            throw CLI::ConversionError(in[0], "-nni");
        }
        if (options.nni == 0) {
            options.spr = 0;
        }
        return true;
    }, "to set the number of rounds of min. evo. NNIs")->group(topo);


    app.add_option("-spr", options.spr, "to set the rounds of SPRs")->group(topo);

    app.add_flag_function("-noml", [&options](size_t) {
                              options.MLnni = 0;
                          },
                          "to turn off both min-evo NNIs and SPRs (useful if refining an approximately "
                          "maximum-likelihood tree with further NNIs)")->group(topo);

    app.add_option("-sprlength", options.maxSPRLength, "set the maximum length of a SPR move (default 10)")->
            group(topo);

    app.add_option("-mlnni", options.MLnni, "to set the number of rounds of maximum-likelihood NNIs")->
            group(topo);

    app.add_option("-mlacc", options.mlAccuracy, "Use -mlacc 2 or -mlacc 3 to always optimize all 5 branches at each"
                                                 " NNI, and to optimize all 5 branches in 2 or 3 rounds")->
            check(Min(1))->group(topo);

    app.add_flag_function("-mllen", [&options](size_t) {
                              options.MLnni = 0;
                              options.MLlen = true;
                          },
                          "to optimize branch lengths without ML NNIs. Use -mllen -nome with -intree to optimize"
                          " branch lengths on a fixed topology ")->group(topo);

    app.add_flag_function("-approxml,-mlapprox", [&options](size_t) {
                              options.exactML = false;
                          },
                          "approximate posterior distributions for a.a.s")->group(topo);

    app.add_flag_function("-slownni", [&options](size_t) {
                              options.fastNNI = false;
                          },
                          "to optimize branch lengths without ML NNIs. Use -mllen -nome with -intree to optimize"
                          " branch lengths on a fixed topology ")->group(topo);


    auto model = "Maximum likelihood model options";

    app.add_flag_function("-lg", [&options](size_t) {
        options.bUseLg = true;
    }, "Le-Gascuel 2008 model instead of (default) Jones-Taylor-Thorton 1992 model (a.a. only)")->group(model);

    app.add_flag_function("-wag", [&options](size_t) {
        options.bUseWag = true;
    }, "Whelan-And-Goldman 2001 model instead of (default) Jones-Taylor-Thorton 1992 model (a.a. only)")->group(model);

    app.add_flag_function("-gtr", [&options](size_t) {
        options.bUseGtr = true;
    }, "generalized time-reversible instead of (default) Jukes-Cantor (nt only)")->group(model);

    app.add_option("-gtrrates", [&args, &options](CLI::results_t in) {
        options.bUseGtr = true;
        options.bUseGtrRates = true;
        for (size_t i = 0; i < in.size(); i++) {
            if (!CLI::detail::lexical_cast(in[i], options.gtrrates[i])) {
                throw CLI::ConversionError(in[i], "-gtrrates");
            }
        }
        return true;
    }, " set the gtr rates")->check(Min(1e-5))->type_size(-1)->expected(6)->group(model);

    app.add_option("-gtrfreq", [&args, &options](CLI::results_t in) {
        options.bUseGtr = true;
        options.bUseGtrFreq = true;
        double sum = 0;
        for (size_t i = 0; i < in.size(); i++) {
            if (!CLI::detail::lexical_cast(in[i], options.gtrfreq[i])) {
                throw CLI::ConversionError(in[i], "-gtrfreq");
            }
            sum += options.gtrrates[i];
        }
        if (std::fabs(1.0 - sum) > 0.01) {
            throw CLI::ValidationError("-gtrfreq", "values do not sum to 1");
        }
        for (size_t i = 0; i < in.size(); i++) {
            options.gtrfreq[i] /= sum;
        }
        return true;
    }, " set the gtr frequences")->check(Min(1e-5))->type_size(-1)->expected(4)->group(model);

    app.add_option("-cat", options.nRateCats, "specify the number of rate categories of sites (default 20)")->
            type_name("n")->check(Min(1))->group(model);

    app.add_flag_function("-nocat", [&options](size_t) {
        options.nRateCats = 1;
    }, "no CAT model (just 1 category)")->group(model);

    app.add_flag_function("-gamma", [&options](size_t) {
                              options.gammaLogLk = true;
                          },
                          "after the final round of optimizing branch lengths with the CAT model, report the likelihood"
                          " under the discrete gamma model with the same number of categories. FastTree uses the same"
                          " branch lengths but optimizes the gamma shape parameter and the scale of the lengths. "
                          "The final tree will have rescaled lengths. Used with -log, this also generates per-site"
                          " likelihoods for use with CONSEL, see GammaLogToPaup.pl and documentation on the FastTree"
                          " web site.")->group(model);

    std::stringstream support_description;
    support_description << "Support value options:" << std::endl;
    support_description << "  By default, FastTree computes local support values by resampling the site" << std::endl;
    support_description << "  likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome,"
                        << std::endl;
    support_description << "  it will compute minimum-evolution bootstrap supports instead" << std::endl;
    support_description << "  In either case, the support values are proportions ranging from 0 to 1" << std::endl;

    auto support = support_description.str();

    app.add_flag_function("-nome", [&options](size_t) {
                              options.spr = 0;
                              options.nni = 0;
                          },
                          "to compute minimum-evolution bootstrap supports")->group(support);

    app.add_flag_function("-nosupport", [&options](size_t) {
        options.nBootstrap = 0;
    }, "to turn off support values")->group(support);

    app.add_option("-boot", options.nBootstrap, "to use just n resamples")->type_name("n")->group(support);

    app.add_flag_function("-noboot", [&options](size_t) {
                              options.nBootstrap = 0;
                          },
                          "to no use resamples")->group(support);

    app.add_option("-seed", options.seed, "to initialize the random number generator")->type_name("n")->group(support);

    std::stringstream search_description;
    search_description << "Searching for the best join:" << std::endl;
    search_description << "  By default, FastTree combines the 'visible set' of fast neighbor-joining " << std::endl;
    search_description << "  with local hill-climbing as in relaxed neighbor-joining" << std::endl;

    auto search = search_description.str();

    app.add_flag_function("-slow", [&options](size_t) {
                              options.slow = true;
                          },
                          "exhaustive search (like NJ or BIONJ, but different gap handling) -slow takes half an hour"
                          " instead of 8 seconds for 1,250 proteins")->group(search);

    app.add_flag_function("-fastest", [&options](size_t) {
                              options.fastest = true;
                              options.tophitsRefresh = 0.5;
                              options.useTopHits2nd = true;
                          },
                          "search the visible set (the top hit for each node) only Unlike the original fast"
                          " neighbor-joining, -fastest updates visible(C) after joining A and B if join(AB,C)"
                          " is better than join(C,visible(C)) -fastest also updates out-distances in a very lazy way,"
                          " -fastest sets -2nd on as well, use -fastest -no2nd to avoid this")->
            group(search)->excludes("-slow");

    std::stringstream heuristics_description;
    heuristics_description << "Top-hit heuristics:" << std::endl;
    heuristics_description << "  By default, FastTree uses a top-hit list to speed up search" << std::endl;

    auto heuristics = heuristics_description.str();

    app.add_flag_function("-top", [&options](size_t) {
        if (options.tophitsMult < 0.01) {
            options.tophitsMult = 1.0;
        }
        return true;
    }, "set the top-hit list size to 1.0 if it is less than 0.01")->group(heuristics);

    app.add_flag_function("-notop", [&options](size_t) {
                              options.tophitsMult = 0.0;
                          },
                          "(or -slow) to turn this feature off and compare all leaves to each other, and all new "
                          "joined nodes to each other")->group(heuristics);

    app.add_option("-topm", options.tophitsMult, "set the top-hit list size to parameter*sqrt(N) FastTree estimates the"
                                                 " top m hits of a leaf from the top 2*m hits of a 'close' neighbor,"
                                                 " where close is defined as d(seed,close) < 0.75 * d(seed, hit of rank"
                                                 " 2*m), and updates the top-hits as joins proceed")->type_name("1.0")->
            group(heuristics);

    app.add_option("-close", options.tophitsMult, "modify the close heuristic, lower is more conservative")->
            type_name("0.75")->check([&options](const std::string &in) {
        if (options.tophitsMult <= 0) {
            return "Cannot use -close unless -top is set above 0";
        }
        if (options.tophitsClose <= 0 || options.tophitsClose >= 1) {
            return "-close argument must be between 0 and 1";
        }
        return "";
    })->group(heuristics);

    app.add_option("-refresh", options.tophitsMult, "compare a joined node to all other nodes if its top-hit list is "
                                                    "less than 80% of the desired length, or if the age of the top-hit"
                                                    " list is log2(m) or greater")->type_name("0.8")->
            check([&options](const std::string &in) {
        if (options.tophitsMult <= 0) {
            return "Cannot use -refresh unless -top is set above 0";
        }
        if (options.tophitsRefresh <= 0 || options.tophitsRefresh >= 1) {
            return "-refresh argument must be between 0 and 1";
        }
        return "";
    })->group(heuristics);


    app.add_flag_function("-2nd", [&options](size_t) { options.useTopHits2nd = true; },
                          "to turn 2nd-level top hits heuristic on. This reduces memory usage and running time but may"
                          " lead to marginal reductions in tree quality. (By default, -fastest turns on -2nd.)")->
            group(heuristics);

    app.add_flag_function("-no2nd", [&options](size_t) { options.useTopHits2nd = false; },
                          "to turn 2nd-level top hits heuristic off.")->group(heuristics);


    auto join = "Join options";

    app.add_flag_function("-nj", [&options](size_t) { options.bionj = false; },
                          "regular (unweighted) neighbor-joining (default)")->
            group(join);

    app.add_flag_function("-bionj", [&options](size_t) { options.bionj = true; },
                          "weighted joins as in BIONJ FastTree will also weight joins during NNIs")->
            group(join);


    auto constrains = "Constrained topology search options";

    app.add_option("-constraints", options.constraintsFile,
                   "an alignment with values of 0, 1, and Not all sequences need be present. A column of 0s and 1s"
                   " defines a constrained split. Some constraints may be violated (see 'violating constraints:' in"
                   " standard error).")->type_name("alignmentfile")->check(isNotEmpty)->group(constrains);

    app.add_option("-constraintWeight", options.constraintWeight,
                   "how strongly to weight the constraints. A value of 1 means a penalty of 1 in tree length for"
                   " violating a constraint Default: 100.0")->type_name("w")->check(Min(0.0, false))->
            group(constrains);

    auto optimizations = "Optimizations";

    app.add_option("-threads", options.threads, "to use a parallel version with multiple threads")->
            type_name("n")->check(Min(1))->envname("OMP_NUM_THREADS")->group(constrains);

    app.add_flag("-double-precision", options.doublePrecision, "use double precision instead of single")->
            group(optimizations);

    app.add_set_ignore_case ("-ext", options.extension, {"NONE","SSE","AVX128","AVX256","AVX512"},
                "compute multiple processing elements with one operation. Available: none, SSE3, AVX, AVX2 or AVX512")
                ->type_name("name")->default_val("NONE")->group(optimizations);
    for(auto& c:options.extension) {c = (char)std::toupper(c);}


    auto deprecated = "Deprecated";

    setDeprecated(app.add_flag("-logdist", "use logarithmic distances, now on by default and obsolete")->
            group(deprecated));

    setDeprecated(app.add_flag("-exactml", "Exact  posterior distributions, now on by default and obsolete")->
            group(deprecated));

    setDeprecated(app.add_flag("-mlexact", "Exact posterior distributions, now on by default and obsolete")->
            group(deprecated));

    app.footer(
            "For more information, see http://www.microbesonline.org/fasttree/ or the comments in the source code");
    app.add_flag("-expert", options.expert)->group("");
}

void basicCli(CLI::App &app, std::string &name, std::string &version, std::string &flags) {
    for (auto *op: app.get_options()) {
        ((CLI::Option *) op)->group("");
    }
    app.name(name);
    std::stringstream description;
    description << name << " " << version << " " << flags << std::endl;
    description << "  " << name << " protein_alignment > tree" << std::endl;
    description << "  " << name << " < protein_alignment > tree" << std::endl;
    description << "  " << name << " -out tree protein_alignment" << std::endl;
    description << "  " << name << " -nt nucleotide_alignment > tree" << std::endl;
    description << "  " << name << " -nt -gtr < nucleotide_alignment > tree" << std::endl;
    description << "  " << name << " < nucleotide_alignment > tree" << std::endl;
    description << name << " accepts alignments in fasta or phylip interleaved formats" << std::endl;
    app.description(description.str());

    auto common = "Common options";
    app.get_option("-n")->description("to analyze multiple alignments (phylip format only) (use for global"
                                      " bootstrap, with seqboot and CompareToBootstrap.pl)")->group(common);
    app.get_option("-intree")->description("to set the starting tree(s)")->group(common);
    app.get_option("-intree1")->description("to use this starting tree for all the alignments "
                                            "(for faster global bootstrap on huge alignments)")->group(common);
    app.get_option("-quiet")->description("to suppress reporting information")->group(common);
    app.get_option("-nopr")->description("to suppress progress indicator")->group(common);
    app.get_option("-log")->description("save intermediate trees, settings, and model details")->group(common);
    app.get_option("-quote")->description(
            "allow spaces and other restricted characters (but not ' ) in sequence names and quote names in the output "
            "tree (fasta input only; FastTree will not be able to read these trees back in)")->group(common);
    app.get_option("-pseudo")->description("to use pseudocounts (recommended for highly gapped sequences)")->
            group(common);
    app.get_option("-fastest")->description(
            "speed up the neighbor joining phase & reduce memory usage (recommended for >50,000 sequences)")->
            group(common);
    app.get_option("-nosupport")->description("to not compute support values")->group(common);
    app.get_option("-gtr")->description("generalized time-reversible model (nucleotide alignments only)")->
            group(common);
    app.get_option("-lg")->description("Le-Gascuel 2008 model (amino acid alignments only)")->group(common);
    app.get_option("-wag")->description("Whelan-And-Goldman 2001 model (amino acid alignments only)")->group(common);
    app.get_option("-noml")->description("to turn off maximum-likelihood")->group(common);
    app.get_option("-nome")->description(
            "to turn off minimum-evolution NNIs and SPRs (recommended if running additional ML NNIs with -intree),"
            " -nome -mllen with -intree to optimize branch lengths for a fixed topology")->
            group(common);
    app.get_option("-cat")->description("to specify the number of rate categories of sites (default 20)"
                                        " or -nocat to use constant rates")->
            type_name("n")->group(common);
    app.get_option("-gamma")->description("after optimizing the tree under the CAT approximation, rescale the lengths"
                                          " to optimize the Gamma20 likelihood")->group(common);
    app.get_option("-constraints")->description(
            "to constrain the topology search constraintAlignment should have 1s or 0s to indicates splits")->
            type_name("constraintAlignment")->group(common);
}


int main(int argc, char *argv[]) {
    fasttree::Options options;
    std::string name = "FastTree";
    std::string version = fasttree::Constants::version;
    std::string flags = fasttree::Constants::compileFlags;
    std::vector<std::string> args(argv + 1, argv + argc);
    CLI::App app;

    cli(app, name, version, flags, options, args);
    if (fasttree::isattyIn() && argc == 1) {
        std::cout << app.help() << std::endl;
        if (fasttree::isWindows()) {
            std::cerr << "Windows users: Please remember to run this inside a command shell" << std::endl;
            std::cerr << "Hit return to continue" << std::endl;
            std::cin.ignore(1);
        }
        return EXIT_SUCCESS;
    }

    CLI11_PARSE(app, argc, argv);
    if (options.expert) {
        CLI::App app2;
        cli(app2, name, version, flags, options, args);
        CLI11_PARSE(app2, "-h", false);
    }

    std::ifstream input;
    std::ofstream output;
    std::ofstream log;
    std::ostream &clog = options.safeLog ? std::cerr : std::clog;
    fasttree::TeeStream tee(log, clog);
    std::ostream teelog(&tee);

    if (!options.inFileName.empty()) {
        input.open(options.inFileName);
        if (input.fail()) {
            std::cerr << "Couldn't open the input file! " << options.inFileName << std::endl;
            return EXIT_FAILURE;
        }
    }

    if (!options.outFileName.empty()) {
        output.open(options.outFileName);
        if (output.fail()) {
            std::cerr << "Couldn't open the output file! " << options.outFileName << std::endl;
            return EXIT_FAILURE;
        }
    }

    if (!options.logFileName.empty()) {
        log.open(options.logFileName, std::ios::out | std::ios::app);
        if (log.fail()) {
            std::cerr << "Couldn't open the log file! " << options.logFileName << std::endl;
            return EXIT_FAILURE;
        }
    }

    fasttree::FastTree fastTree(options);
    std::ostream &applog = log ? teelog : clog;
    applog << "Command: ";
    std::copy(argv, argv + argc - 1, std::ostream_iterator<char *>(applog, ", "));
    applog << argv[argc - 1] << std::endl;

    try {
        fastTree.run(
                input ? input : std::cin,
                output ? output : std::cout,
                applog);

    } catch (std::invalid_argument &ex) {
        std::cerr << ex.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}