/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PRINT_OVERLAPS_HPP_
#define BR_PRINT_OVERLAPS_HPP_

#include "AbstractOutput.hpp"

namespace bloomrepeats {

/** Print ASCII diagram with all fragments overlapping with a block.
Fragments must be \ref Connector "connected"

It is recommended to use this processor if blocks have alignment.
*/
class PrintOverlaps : public AbstractOutput {
public:
    /** Constructor */
    PrintOverlaps();

    /** Get if block name is printed */
    bool print_block() const {
        return print_block_;
    }

    /** Set if block name is printed.
    Defaults to true;
    */
    void set_print_block(bool print_block) {
        print_block_ = print_block;
    }

    /** Get if fragment name is printed */
    bool print_fragment() const {
        return print_fragment_;
    }

    /** Set if fragment name is printed.
    Defaults to false;
    */
    void set_print_fragment(bool print_fragment) {
        print_fragment_ = print_fragment;
    }

    /** Get if top scale is printed */
    bool top_scale() const {
        return top_scale_;
    }

    /** Set if top scale is printed.
    Defaults to true.
    */
    void set_top_scale(bool top_scale) {
        top_scale_ = top_scale;
    }

    /** Get if bottom scale is printed */
    bool bottom_scale() const {
        return bottom_scale_;
    }

    /** Set if bottom scale is printed.
    Defaults to true.
    */
    void set_bottom_scale(bool bottom_scale) {
        bottom_scale_ = bottom_scale;
    }

    /** Get max allowed line width of output */
    int width() const {
        return width_;
    }

    /** Set max allowed line width of output.
    Defaults to 76.
    */
    void set_width(int width) {
        width_ = width;
    }

    /** Get char used to mark fragment */
    char marker() const {
        return marker_;
    }

    /** Set char used to mark fragment.
    Defaults to '*'.
    */
    void set_marker(char marker) {
        marker_ = marker;
    }

    /** Print ASCII diagram with all fragments overlapping with a block */
    void print_block(std::ostream& o, Block* block) const;

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    const char* name_impl() const;

private:
    bool print_block_;
    bool print_fragment_;
    bool top_scale_;
    bool bottom_scale_;
    int width_;
    char marker_;
};

}

#endif

