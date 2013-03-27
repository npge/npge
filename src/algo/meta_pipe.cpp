/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "meta_pipe.hpp"
#include "Pipe.hpp"
#include "Meta.hpp"
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
    ProcessorPtr p = meta->get(key);
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
    bool r = phrase_parse(first, last,
                          //  Begin grammar
                          (
                              lit("pipe") >> lexeme[+char_("a-zA-Z0-9")][boost::bind(set_k, pipe, _1)]
                              >> '{' >> *(
                                  lit("name") >> lexeme['"' >> +(char_ - '"') >> '"']
                                  [boost::bind(set_n, pipe, _1)] >> ';'
                                  || lit("max_loops") >> int_
                                  [boost::bind(&Pipe::set_max_loops, pipe, _1)] >> ';'
                                  || lit("workers") >> int_
                                  [boost::bind(&Processor::set_workers, pipe, _1)] >> ';'
                                  || lit("no_options") >> bool_
                                  [boost::bind(&Processor::set_no_options, pipe, _1)] >> ';'
                                  || lit("timing") >> bool_
                                  [boost::bind(&Processor::set_timing, pipe, _1)] >> ';'
                                  || lit("add") >> (lexeme[+char_("a-zA-Z0-9")] >> *(char_ - ';'))
                                  [boost::bind(add_p, pipe, meta, _1)] >> ';'
                              ) >> '}' >> ';'
                          )
                          ,
                          //  End grammar
                          space | '#' >> *(char_ - eol) >> eol // comment skipper
                         );
    return r;
}

static Meta default_meta;

boost::shared_ptr<Pipe> create_pipe(const std::string& script,
                                    const Meta* meta, std::string* tail) {
    if (meta == 0) {
        meta = &default_meta;
    }
    boost::shared_ptr<Pipe> result(new Pipe);
    typedef std::string::const_iterator It;
    It first = script.begin();
    bool ok = parse_pipe(first, script.end(), result.get(), meta);
    BOOST_ASSERT_MSG(ok, ("Can't parse pipe description: " + script).c_str());
    BOOST_ASSERT(first <= script.end());
    if (tail) {
        tail->assign(first, script.end());
    }
    return result;
}

ProcessorPtr parse_script(const std::string& script0, Meta* meta) {
    std::string script = script0;
    while (!script.empty()) {
        // TODO N**2, where N = len(script0)
        using namespace boost::algorithm;
        trim(script);
        if (starts_with(script, "run")) {
            std::string processor_name = script.substr(4, script.size() - 5);
            trim(processor_name);
            BOOST_ASSERT_MSG(meta->has(processor_name),
                             ("No such processor: " + processor_name).c_str());
            return meta->get(processor_name);
        } else {
            std::string tail;
            ProcessorPtr new_pipe = create_pipe(script, meta, &tail);
            std::string beginning(script, 0, script.size() - tail.size());
            meta->set_returner(boost::bind(create_pipe, beginning, meta,
                                           /* tail */ static_cast<std::string*>(0)));
            script.swap(tail); // script = tail
        }
    }
    return ProcessorPtr();
}

}

