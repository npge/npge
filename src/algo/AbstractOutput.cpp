/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>
#include "boost-xtime.hpp"
#include <boost/thread/mutex.hpp>

#include "AbstractOutput.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

typedef boost::shared_ptr<std::ostringstream> SstreamPtr;
typedef std::map<Block*, SstreamPtr> Block2SP;

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
        ASSERT_TRUE(out_);
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
    add_opt("file", "output file with all blocks",
            std::string());
}

AbstractOutput::~AbstractOutput() {
    delete impl_;
}

void AbstractOutput::change_blocks_impl(Blocks& blocks) const {
    sort_blocks(blocks);
    if (workers() >= 2) {
        impl_->blocks_ = blocks;
        impl_->blocks_it_ = impl_->blocks_.begin();
    }
}

void AbstractOutput::initialize_work_impl() const {
    impl_->main_thread_ = false;
    std::string file = opt_value("file").as<std::string>();
    impl_->out_ = name_to_ostream(file);
    print_header(*impl_->out_);
    prepare();
}

ThreadData* AbstractOutput::before_thread_impl() const {
    if (!impl_->main_thread_) {
        impl_->main_thread_ = true;
        return new ThreadData;
    }
    return 0;
}

void AbstractOutput::process_block_impl(Block* block, ThreadData* data) const {
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
}

void AbstractOutput::finish_work_impl() const {
    if (workers() >= 2) {
        impl_->move_text();
    }
    print_footer(*impl_->out_);
    impl_->out_.reset(); // close file
}

void AbstractOutput::prepare() const {
}

void AbstractOutput::print_header(std::ostream& o) const {
}

void AbstractOutput::print_footer(std::ostream& o) const {
}

}

