/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_IMPORT_BLAST_HITS_HPP_
#define BR_IMPORT_BLAST_HITS_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "OtherBlockSet.hpp"
#include "FileReader.hpp"

namespace bloomrepeats {

/** Add blocks from blast hits (blast output -m 8).
\note This processor depends on AddBlocks(keep_alignment = true).
*/
class ImportBlastHits : public Processor, public OtherBlockSet,
    public FileReader {
public:
    /** Constructor.
    \param block_set The block set, passed to blast.
    \param min_length Min accepted length of blast hit
    \param min_ident Min accepted identity of blast hit
    \param max_evalue Max accepted e-value of blast hit
    */
    ImportBlastHits(const BlockSetPtr& block_set = BlockSetPtr(),
                    int min_length = 100,
                    float min_ident = 0.95, float max_evalue = 1e-6);

    /** Get min accepted length of blast hit */
    int min_length() const {
        return min_length_;
    }

    /** Set min accepted length of blast hit */
    void set_min_length(int min_length) {
        min_length_ = min_length;
    }

    /** Get min accepted identity of blast hit */
    float min_ident() const {
        return min_ident_;
    }

    /** Set min accepted identity of blast hit */
    void set_min_ident(float min_ident) {
        min_ident_ = min_ident;
    }

    /** Get max accepted e-value of blast hit */
    float max_evalue() const {
        return max_evalue_;
    }

    /** Set max accepted e-value of blast hit */
    void set_max_evalue(float max_evalue) {
        max_evalue_ = max_evalue;
    }

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;

private:
    int min_length_;
    float min_ident_;
    float max_evalue_;
};

}

#endif

