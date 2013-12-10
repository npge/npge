/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <exception>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "process.hpp"
#include "meta_pipe.hpp"
#include "po.hpp"
#include "string_arguments.hpp"
#include "name_to_stream.hpp"

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

int process(int argc, char** argv,
            Processor* processor,
            const std::string& name,
            const std::string& positional) {
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
    } else {
        try {
            processor->apply_options(vm);
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
    std::vector<std::string> errors = processor->options_errors();
    if (!errors.empty()) {
        std::cerr << argv[0];
        std::cerr << ": error while validating options" << std::endl;
        BOOST_FOREACH (const std::string& message, errors) {
            std::cerr << message << std::endl;
        }
        return 255;
    }
    std::vector<std::string> warnings = processor->options_warnings();
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
    return 0;
}

int process_and_delete(int argc, char** argv,
                       const std::vector<Processor*>& processors) {
    int result = 0;
    std::vector<SharedProcessor> ps;
    BOOST_FOREACH (Processor* p, processors) {
        ps.push_back(SharedProcessor(p));
    }
    BOOST_FOREACH (SharedProcessor p, ps) {
        int r = process(argc, argv, p.get());
        if (r) {
            result = r;
        }
    }
    return result;
}

int execute_script(const std::string& script, const std::string& output,
                   int argc, char** argv, Meta* meta, bool debug) {
    int result = 0;
    boost::shared_ptr<std::ostream> output_ptr = name_to_ostream(output);
    std::ostream& output_stream = *output_ptr;
    std::vector<Processor*> raw_ps;
    if (debug) {
        raw_ps = parse_script_to_processors(script, meta);
        result |= process_and_delete(argc, argv, raw_ps);
    } else {
        try {
            raw_ps = parse_script_to_processors(script, meta);
            result |= process_and_delete(argc, argv, raw_ps);
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
        std::getline(input_stream, line);
        using namespace boost::algorithm;
        trim_right(line);
        buffer += line;
        if (ends_with(buffer, ";")) {
            bool debug = debug0;
            StringToArgv args(args0);
            // TODO DRY
            if (buffer.find(" --help") != std::string::npos) {
                args.add_argument("--help");
            }
            if (buffer.find(" -h") != std::string::npos) {
                args.add_argument("-h");
            }
            if (buffer.find(" --debug") != std::string::npos) {
                args.add_argument("--debug");
                debug = true;
            }
            if (buffer.find(" --tree") != std::string::npos) {
                args.add_argument("--tree");
            }
            int r = execute_script(buffer, output, args.argc(), args.argv(),
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

