/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ROW_STORAGE_HPP_
#define BR_ROW_STORAGE_HPP_

#include "OptionsPrefix.hpp"

namespace bloomrepeats {

/** Utility class storing type of alignment row storage */
class RowStorage : public OptionsPrefix {
public:
    /** Constructor */
    RowStorage(bool keep_alignment = false, RowType row_type = COMPACT_ROW);

    /** Get if alignments is extracted too (not only blocks) */
    bool keep_alignment() const {
        return keep_alignment_;
    }

    /** Set if alignments is extracted too (not only blocks) */
    void set_keep_alignment(bool keep_alignment) {
        keep_alignment_ = keep_alignment;
    }

    /** Get alignment row type */
    const RowType row_type() const {
        return row_type_;
    }

    /** Set alignment row type */
    void set_row_type(RowType row_type) {
        row_type_ = row_type;
    }

    /** Create new alignment row.
    \warning Asserts keep_alignment().
    */
    AlignmentRow* create_row() const;

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

private:
    bool keep_alignment_;
    RowType row_type_;
};

}

#endif

