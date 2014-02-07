/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include "boost-xtime.hpp"
#include <boost/thread/mutex.hpp>

#include "AbstractOutput.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

static bool file_and_mask_check(AbstractOutput* p,
        std::string& message) {
    if (p->opt_value("file").as<std::string>() != "" &&
            p->opt_value("mask").as<std::string>() != "") {
        message = "both '" + p->opt_prefixed("file") +
                  "' and '" + p->opt_prefixed("mask") +
                  "' were specified";
        return false;
    } else {
        return true;
    }
}

static bool mask_check(AbstractOutput* p, std::string& message) {
    std::string mask = p->opt_value("mask").as<std::string>();
    if (mask != "" && mask.find("${block}") == std::string::npos) {
        message = "'" + p->opt_prefixed("mask") +
            "' must contain '${block}'";
        return false;
    } else {
        return true;
    }
}

typedef boost::shared_ptr<std::ostringstream> SstreamPtr;
typedef std::map<Block*, SstreamPtr> Block2SP;
typedef std::vector<Block*> Blocks;

struct AbstractOutput::Impl {
    boost::mutex mutex_;
    Block2SP block2op_;
    Blocks blocks_;
    Blocks::const_iterator blocks_it_;
    boost::shared_ptr<std::ostream> out_;
    bool main_thread_;

    void move_text() {
        // do not call this from concurent threads
        std::vector<SstreamPtr> sstreams;
        {
            boost::mutex::scoped_lock lock(mutex_);
            while (blocks_it_ != blocks_.end()) {
                Block* block = *blocks_it_;
                Block2SP::iterator it = block2op_.find(block);
                if (it != block2op_.end()) {
                    sstreams.push_back(it->second);
                    block2op_.erase(it);
                    blocks_it_++;
                } else {
                    break;
                }
            }
        }
        BOOST_ASSERT(out_);
        std::ostream& out = *out_;
        BOOST_FOREACH (SstreamPtr& sstr, sstreams) {
            out << sstr->str();
        }
    }

    void add_text(Block* block, const SstreamPtr& sstr) {
        boost::mutex::scoped_lock lock(mutex_);
        block2op_[block] = sstr;
    }
};

AbstractOutput::AbstractOutput():
        impl_(new Impl) {
    add_opt("file", "output file with all blocks", std::string());
    add_opt("mask", "mask of output files (${block} is "
            "replaced with block name)", std::string());
    add_opt_check(boost::bind(file_and_mask_check, this, _1));
    add_opt_check(boost::bind(mask_check, this, _1));
}

AbstractOutput::~AbstractOutput() {
    delete impl_;
}

bool AbstractOutput::one_file() const {
    return opt_value("mask").as<std::string>().empty();
}

static struct BlockCompareName2 {
    bool operator()(const Block* b1, const Block* b2) const {
        typedef boost::tuple<int, const std::string&> Tie;
        return Tie(-b1->size(), b1->name()) < Tie(-b2->size(), b2->name());
    }
} bcn2;

bool AbstractOutput::change_blocks_impl(std::vector<Block*>& blocks) const {
    if (one_file()) {
        std::sort(blocks.begin(), blocks.end(), bcn2);
        if (workers() >= 2) {
            impl_->blocks_ = blocks;
            impl_->blocks_it_ = impl_->blocks_.begin();
        }
    }
    return false;
}

bool AbstractOutput::initialize_work_impl() const {
    impl_->main_thread_ = false;
    if (one_file()) {
        std::string file = opt_value("file").as<std::string>();
        impl_->out_ = name_to_ostream(file);
        print_header(*impl_->out_);
    }
    prepare();
    return false;
}

ThreadData* AbstractOutput::before_thread_impl() const {
    if (!impl_->main_thread_) {
        impl_->main_thread_ = true;
        return new ThreadData;
    }
    return 0;
}

bool AbstractOutput::process_block_impl(Block* block, ThreadData* data) const {
    if (one_file()) {
        if (workers() >= 2) {
            SstreamPtr sstr(new std::ostringstream);
            print_block(*sstr, block);
            impl_->add_text(block, sstr);
            if (data) {
                // existance of data mark this thread as moving
                // buffers to the output file
                impl_->move_text();
            }
        } else {
            print_block(*impl_->out_, block);
        }
    } else {
        using namespace boost::algorithm;
        std::string mask = opt_value("mask").as<std::string>();
        std::string path = replace_all_copy(mask,
                "${block}", block->name());
        boost::shared_ptr<std::ostream> o = name_to_ostream(path);
        BOOST_ASSERT(o);
        print_header(*o);
        print_block(*o, block);
        print_footer(*o);
    }
    return false;
}

bool AbstractOutput::finish_work_impl() const {
    if (one_file()) {
        if (workers() >= 2) {
            impl_->move_text();
        }
        print_footer(*impl_->out_);
        impl_->out_.reset(); // close file
    }
}

void AbstractOutput::prepare() const
{ }

void AbstractOutput::print_header(std::ostream& o) const
{ }

void AbstractOutput::print_footer(std::ostream& o) const
{ }

}

