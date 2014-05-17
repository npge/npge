/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <cstring>
#include <boost/foreach.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/fusion/include/at_c.hpp>

#include "meta_pipe.hpp"
#include "Pipe.hpp"
#include "Meta.hpp"
#include "tss_meta.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

typedef std::vector<char> String;

std::string to_string(const String& k) {
    return std::string(k.begin(), k.end());
}

void set_k(Pipe* pipe, const String& k) {
    pipe->set_key(to_string(k));
}

void set_n(Pipe* pipe, const String& k) {
    pipe->set_name(to_string(k));
}

typedef boost::fusion::vector<String, String> TwoStrings;

void set_bs(Pipe* pipe, const TwoStrings& n_d) {
    const String& name = boost::fusion::at_c<0>(n_d);
    const String& desc = boost::fusion::at_c<1>(n_d);
    pipe->declare_bs(to_string(name), to_string(desc));
}

void add_p(Pipe* pipe, const Meta* meta, const TwoStrings& processor) {
    std::string key = to_string(boost::fusion::at_c<0>(processor));
    std::string options = to_string(boost::fusion::at_c<1>(processor));
    BOOST_ASSERT_MSG(meta->has(key), ("No such processor: " + key).c_str());
    Processor* p = meta->get_plain(key);
    pipe->add(p, options);
}

template <typename Iterator>
bool parse_pipe(Iterator& first, Iterator last,
                Pipe* pipe, const Meta* meta) {
    namespace qi = boost::spirit::qi;
    using qi::bool_;
    using qi::char_;
    using qi::int_;
    using qi::lexeme;
    using qi::lit;
    using qi::phrase_parse;
    using boost::spirit::iso8859_1::space;
    using boost::spirit::eol;
    bool r = phrase_parse(
                 first, last,
                 //  Begin grammar
                 (
                     lit("pipe") >> lexeme[+char_("a-zA-Z0-9")]
                     [boost::bind(set_k, pipe, _1)]
                     >> '{' >> *(
                         lit("name") >> lexeme['"' >> +(char_ - '"') >> '"']
                         [boost::bind(set_n, pipe, _1)] >> ';'
                         || (lit("bs") >> +(char_ - '"' - ' ') >>
                             lexeme['"' >> +(char_ - '"') >> '"'])
                         [boost::bind(set_bs, pipe, _1)] >> ';'
                         || lit("max_loops") >> int_
                         [boost::bind(&Pipe::set_max_loops, pipe, _1)] >> ';'
                         || lit("workers") >> int_
                         [boost::bind(&Processor::set_workers, pipe, _1)] >> ';'
                         || lit("no_options") >> bool_
                         [boost::bind(&Processor::set_no_options, pipe, _1)] >>
                         ';'
                         || lit("timing") >> bool_
                         [boost::bind(&Processor::set_timing, pipe, _1)] >> ';'
                         || lit("add") >> (lexeme[+char_("a-zA-Z0-9")] >>
                                           lexeme[*(char_ - ';')])
                         [boost::bind(add_p, pipe, meta, _1)] >> ';'
                     ) >> '}' >> ';'
                 )
                 ,
                 //  End grammar
                 space | '#' >> *(char_ - eol) >> eol // comment skipper
             );
    return r;
}

Pipe* create_pipe_c(const char*& begin, const char* end,
                    const Meta* meta) {
    if (meta == 0) {
        meta = tss_meta();
    }
    Pipe* result = new Pipe;
    bool ok = parse_pipe(begin, end, result, meta);
    if (!ok) {
        delete result;
        throw Exception("Can't parse pipe from " + std::string(begin, end));
    }
    if (begin > end) {
        delete result;
        throw Exception("Pipe creation error: begin > end");
    }
    return result;
}

Pipe* create_pipe(const std::string& script,
                  const Meta* meta, std::string* tail) {
    const char* begin = script.c_str();
    const char* end = &(*script.end());
    Pipe* result = create_pipe_c(begin, end, meta);
    if (tail) {
        tail->assign(begin, end);
    }
    return result;
}

Pipe* create_pipe_c(const char* script, const Meta* meta,
                    const char** tail) {
    const char* begin = script;
    const char* end = begin + strlen(begin);
    Pipe* result = create_pipe_c(begin, end, meta);
    if (tail) {
        *tail = begin;
    }
    return result;
}

static void trim_begin(const char*& begin) {
    while (isspace(*begin)) {
        begin++;
    }
}

static void trim_end(const char* begin, const char*& end) {
    end--;
    for (; end >= begin && isspace(*end); end--) {
    }
    end++;
}

static std::string remove_comments(const std::string& script) {
    using namespace boost::algorithm;
    Strings lines;
    split(lines, script, is_any_of("\n"));
    BOOST_FOREACH (std::string& line, lines) {
        size_t sharp_pos = line.find('#');
        if (sharp_pos != std::string::npos) {
            line.resize(sharp_pos);
        }
    }
    return join(lines, "\n");
}

std::vector<Processor*> parse_script_to_processors(const std::string& script0,
        Meta* meta) {
    std::vector<Processor*> result;
    std::string script = remove_comments(script0);
    const char* begin = script.c_str();
    const char* end = &(*script.end());
    trim_end(begin, end);
    while (begin != end) {
        using namespace boost::algorithm;
        trim_begin(begin);
        if (strncmp(begin, "run", 3) == 0) {
            const char* run_end = strchr(begin, ';');
            BOOST_ASSERT_MSG(run_end, "No ';' found after 'run' command");
            std::string processor_text(begin + 4, run_end);
            trim(processor_text);
            size_t space_pos = processor_text.find(' ');
            std::string processor_name = processor_text.substr(0, space_pos);
            BOOST_ASSERT_MSG(meta->has(processor_name),
                             ("No such processor: " + processor_name).c_str());
            Processor* p = meta->get_plain(processor_name);
            p->set_options("", meta->placeholder_processor());
            if (space_pos != std::string::npos) {
                std::string processor_opts = processor_text.substr(space_pos);
                p->set_options(processor_opts, meta->placeholder_processor());
            }
            result.push_back(p);
            begin = run_end + 1;
        } else {
            const char* script_begin = begin;
            Processor* new_pipe = create_pipe_c(begin, end, meta);
            delete new_pipe;
            const char* script_end = begin;
            std::string beginning(script_begin, script_end);
            std::string* tail = 0;
            meta->set_returner(boost::bind(create_pipe, beginning, meta, tail));
        }
    }
    return result;
}

Processor* parse_script(const std::string& script, Meta* meta) {
    std::vector<Processor*> ps = parse_script_to_processors(script, meta);
    if (ps.empty()) {
        return 0;
    } else if (ps.size() == 1) {
        return ps[0];
    } else {
        Pipe* main_pipe = new Pipe;
        main_pipe->set_name("Main pipe");
        BOOST_FOREACH (Processor* p, ps) {
            main_pipe->add(p);
        }
        return main_pipe;
    }
}

}

