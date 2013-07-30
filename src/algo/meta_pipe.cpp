/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <cstring>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "meta_pipe.hpp"
#include "Pipe.hpp"
#include "Meta.hpp"
#include "tss_meta.hpp"
#include "throw_assert.hpp"

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
    namespace ascii = boost::spirit::ascii;
    using qi::bool_;
    using qi::char_;
    using qi::int_;
    using qi::lexeme;
    using qi::lit;
    using qi::phrase_parse;
    using ascii::space;
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
    BOOST_ASSERT_MSG(ok, ("Can't parse pipe description: " +
                          std::string(begin, end)).c_str());
    BOOST_ASSERT(begin <= end);
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

static void trim_begin(const char*& begin) {
    while (isspace(*begin)) {
        begin++;
    }
}

static void trim_end(const char* begin, const char*& end) {
    end--;
    for (; end >= begin && isspace(*end); end--)
    { }
    end++;
}

Processor* parse_script(const std::string& script, Meta* meta) {
    const char* begin = script.c_str();
    const char* end = &(*script.end());
    trim_end(begin, end);
    while (begin != end) {
        using namespace boost::algorithm;
        trim_begin(begin);
        if (strncmp(begin, "run", 3) == 0) {
            std::string processor_name(begin + 4, end - 1);
            trim(processor_name);
            BOOST_ASSERT_MSG(meta->has(processor_name),
                             ("No such processor: " + processor_name).c_str());
            return meta->get_plain(processor_name);
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
    return 0;
}

}

