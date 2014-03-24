/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <csignal>
#include <iostream>
#include <exception>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "process.hpp"
#include "meta_pipe.hpp"
#include "po.hpp"
#include "string_arguments.hpp"
#include "name_to_stream.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace bloomrepeats {

void print_processor_tree(Processor* processor, int indent = 0) {
    const int SPACES_IN_TAB = 4;
    std::string tab(SPACES_IN_TAB * indent, ' ');
    *name_to_ostream(":cout") << tab << processor->key();
    *name_to_ostream(":cout") << ": " << processor->name() << "\n";
    BOOST_FOREACH (Processor* child, processor->children()) {
        print_processor_tree(child, indent + 1);
    }
}

void print_help(const std::string& output, const Processor* processor,
                const std::string& app, const std::string& positional) {
    boost::shared_ptr<std::ostream> output_ptr = name_to_ostream(output);
    std::ostream& out = *output_ptr;
    out << "Usage:" << std::endl;
    out << app << " [options]";
    if (!positional.empty()) {
        using namespace boost::algorithm;
        std::string::const_iterator it = std::find_if(positional.begin(),
                                         positional.end(), !is_any_of("-"));
        if (it != positional.end()) {
            out << ' ' << std::string(it, positional.end());
        }
    }
    po::options_description desc(processor->name());
    add_general_options(desc);
    processor->add_options(desc);
    out << std::endl << std::endl;
    out << desc << std::endl;
}

typedef void (*SignalHandler)(int);

static SignalHandler prev_handler_ = 0;
static Processor* signal_processor_ = 0;

static void process_handler(int) {
    BOOST_ASSERT(signal_processor_ != 0);
    signal_processor_->interrupt();
}

class SignalManager {
public:
    SignalManager(Processor* processor) {
        BOOST_ASSERT(prev_handler_ == 0);
        BOOST_ASSERT(signal_processor_ == 0);
        BOOST_ASSERT(processor != 0);
        signal_processor_ = processor;
        prev_handler_ = signal(SIGINT, process_handler);
    }

    ~SignalManager() {
        BOOST_ASSERT(signal_processor_ != 0);
        SignalHandler prev_handler = signal(SIGINT, prev_handler_);
        BOOST_ASSERT(prev_handler == process_handler);
        signal_processor_ = 0;
        prev_handler_ = 0;
    }
};

int process(int argc, char** argv,
            Processor* processor,
            const std::string& name,
            const std::string& positional,
            bool catch_sigint) {
    boost::shared_ptr<SignalManager> sm;
    if (catch_sigint) {
        sm.reset(new SignalManager(processor));
    }
    if (has_arg(argc, argv, "--tree")) {
        print_processor_tree(processor);
        return 0;
    }
    po::options_description desc(name);
    add_general_options(desc);
    processor->add_options(desc);
    po::positional_options_description pod;
    if (!positional.empty()) {
        pod.add(positional.c_str(), -1);
    }
    po::variables_map vm;
    int error = read_options(argc, argv, vm, desc, pod);
    if (error) {
        return error;
    }
    if (vm.count("debug")) {
        processor->apply_options(vm);
        if (vm.count("help")) {
            print_help(":cout", processor, argv[0], positional);
            return 1;
        }
    } else {
        try {
            processor->apply_options(vm);
            if (vm.count("help")) {
                print_help(":cout", processor, argv[0], positional);
                return 1;
            }
        } catch (std::exception& e) {
            std::cerr << argv[0];
            std::cerr << ": error while applying options" << std::endl;
            std::cerr << "  " << e.what() << std::endl;
            return 255;
        } catch (...) {
            std::cerr << argv[0];
            std::cerr << ": error while applying options" << std::endl;
            std::cerr << "  Unknown error" << std::endl;
            return 255;
        }
    }
    Strings errors = processor->options_errors();
    if (!errors.empty()) {
        std::cerr << argv[0];
        std::cerr << ": error while validating options" << std::endl;
        BOOST_FOREACH (const std::string& message, errors) {
            std::cerr << message << std::endl;
        }
        return 255;
    }
    Strings warnings = processor->options_warnings();
    if (!warnings.empty()) {
        std::cerr << argv[0];
        std::cerr << ": warnings while validating options" << std::endl;
        BOOST_FOREACH (const std::string& message, warnings) {
            std::cerr << message << std::endl;
        }
    }
    if (!processor->block_set()) {
        processor->set_empty_block_set();
    }
    int workers = processor->workers();
    uint32_t hash_1 = blockset_hash(*processor->block_set(), workers);
    if (vm.count("debug")) {
        processor->run();
    } else {
        try {
            processor->run();
        } catch (std::exception& e) {
            std::cerr << argv[0] << ": algorithm error" << std::endl;
            std::cerr << "  " << e.what() << std::endl;
            return 255;
        } catch (...) {
            std::cerr << argv[0];
            std::cerr << ": error while applying options" << std::endl;
            std::cerr << "  Unknown error" << std::endl;
            return 255;
        }
    }
    uint32_t hash_2 = blockset_hash(*processor->block_set(), workers);
    bool changed = (hash_1 != hash_2);
    std::cerr << processor->key() << ": ";
    std::cerr << (changed ? "changed" : "unchanged") << std::endl;
    return 0;
}

int process_and_delete(int argc, char** argv,
                       const std::vector<Processor*>& processors,
                       const std::string& positional) {
    int result = 0;
    std::vector<SharedProcessor> ps;
    BOOST_FOREACH (Processor* p, processors) {
        ps.push_back(SharedProcessor(p));
    }
    BOOST_FOREACH (SharedProcessor p, ps) {
        int r = process(argc, argv, p.get(),
                        /* name */ "", positional);
        if (r) {
            result = r;
        }
    }
    return result;
}

int execute_script(const std::string& script,
                   const std::string& output,
                   int argc, char** argv, Meta* meta, bool debug,
                   const std::string& positional) {
    int result = 0;
    boost::shared_ptr<std::ostream> output_ptr = name_to_ostream(output);
    std::ostream& output_stream = *output_ptr;
    std::vector<Processor*> raw_ps;
    if (debug) {
        raw_ps = parse_script_to_processors(script, meta);
        result |= process_and_delete(argc, argv, raw_ps, positional);
    } else {
        try {
            raw_ps = parse_script_to_processors(script, meta);
            result |= process_and_delete(argc, argv, raw_ps,
                                         positional);
        } catch (std::exception& e) {
            output_stream << e.what() << std::endl;
            result = 15;
        } catch (...) {
            output_stream << "Unknown error" << std::endl;
            result = 15;
        }
    }
    return result;
}

static bool has_opt(std::string buffer, const std::string& opt) {
    // FIXME
    using namespace boost::algorithm;
    trim_right(buffer);
    buffer.resize(buffer.size() - 1); // remove ';'
    buffer += " "; // append ' '
    return buffer.find(" " + opt + " ") != std::string::npos;
}

static SignalHandler prev_handler_2_ = 0;
static std::ostream* signal_ostream_ = 0;

static void process_handler_2(int) {
    BOOST_ASSERT(signal_ostream_ != 0);
    (*signal_ostream_) << "SIGINT catched. Enter quit;" << "\n";
}

class SignalManager2 {
public:
    SignalManager2(std::ostream* signal_ostream) {
        BOOST_ASSERT(prev_handler_2_ == 0);
        BOOST_ASSERT(signal_ostream_ == 0);
        BOOST_ASSERT(signal_ostream != 0);
        signal_ostream_ = signal_ostream;
        prev_handler_2_ = signal(SIGINT, process_handler_2);
    }

    ~SignalManager2() {
        BOOST_ASSERT(signal_ostream_ != 0);
        SignalHandler prev_handler = signal(SIGINT, prev_handler_2_);
        BOOST_ASSERT(prev_handler == process_handler_2);
        signal_ostream_ = 0;
        prev_handler_2_ = 0;
    }
};

void get_line(std::istream& input_stream, std::ostream& output_stream,
              std::string& line) {
    SignalManager2 sm2(&output_stream);
    std::getline(input_stream, line);
}

int interactive_loop(const std::string& input, const std::string& output,
                     int argc, char** argv, Meta* meta) {
    int result = 0;
    boost::shared_ptr<std::istream> input_ptr = name_to_istream(input);
    std::istream& input_stream = *input_ptr;
    boost::shared_ptr<std::ostream> output_ptr = name_to_ostream(output);
    std::ostream& output_stream = *output_ptr;
    StringToArgv args0(argc, argv);
    args0.remove_argument("--help"); // TODO DRY see po.cpp add_general_options
    args0.remove_argument("-h");
    args0.remove_argument("--debug");
    args0.remove_argument("--tree");
    args0.remove_argument("-i");
    bool debug0 = has_arg(argc, argv, "--debug");
    std::string buffer;
    std::string line;
    while (input_stream) {
        output_stream << (buffer.empty() ? "% " : ". ");
        line.clear();
        get_line(input_stream, output_stream, line);
        using namespace boost::algorithm;
        trim_right(line);
        buffer += line;
        if (ends_with(buffer, ";")) {
            if (buffer == "quit;") {
                break;
            }
            bool debug = debug0;
            StringToArgv args(args0);
            // TODO DRY
            if (has_opt(buffer, "--help")) {
                args.add_argument("--help");
            }
            if (has_opt(buffer, "-h")) {
                args.add_argument("-h");
            }
            if (has_opt(buffer, "--debug")) {
                args.add_argument("--debug");
                debug = true;
            }
            if (has_opt(buffer, "--tree")) {
                args.add_argument("--tree");
            }
            int r = execute_script(buffer, output,
                                   args.argc(), args.argv(),
                                   meta, debug);
            buffer.clear();
            if (r) {
                result = r;
                continue;
            }
        }
    }
    output_stream << std::endl;
    return result;
}

void copy_processor_options(Processor& dest, const Processor& source) {
    StringToArgv dummy;
    po::options_description desc;
    source.add_options(desc);
    po::variables_map vm;
    po::store(po::command_line_parser(dummy.argc(), dummy.argv()).options(desc)
              .allow_unregistered().run(), vm);
    // po::notify(vm); // to pass required options check
    dest.apply_options(vm);
    dest.set_timing(source.timing());
}

}

