/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Meta.hpp"
#include "AddSequences.hpp"
#include "SequencesFromOther.hpp"
#include "AddBlocks.hpp"
#include "BlastRunner.hpp"
#include "ImportBlastHits.hpp"
#include "AddBlastBlocks.hpp"
#include "ResolveBlast.hpp"
#include "ResolveAnchors.hpp"
#include "AnchorFinder.hpp"
#include "Filter.hpp"
#include "Connector.hpp"
#include "OriByMajority.hpp"
#include "StickBoundaries.hpp"
#include "OverlapsResolver.hpp"
#include "OverlapsResolver2.hpp"
#include "CheckNoOverlaps.hpp"
#include "Joiner.hpp"
#include "Union.hpp"
#include "OverlaplessUnion.hpp"
#include "Subtract.hpp"
#include "BlocksExpander.hpp"
#include "FragmentsExpander.hpp"
#include "Pipe.hpp"
#include "CleanUp.hpp"
#include "Output.hpp"
#include "OutputPipe.hpp"
#include "PrintOverlaps.hpp"
#include "BlockInfo.hpp"
#include "Stats.hpp"
#include "IsPangenome.hpp"
#include "Consensus.hpp"
#include "UniqueNames.hpp"
#include "Rest.hpp"
#include "ExternalAligner.hpp"
#include "RemoveAlignment.hpp"
#include "ConSeq.hpp"
#include "DeConSeq.hpp"
#include "MoveGaps.hpp"
#include "CutGaps.hpp"
#include "MetaProcessor.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

Meta::Meta() {
    set_processor<AddSequences>();
    set_processor<SequencesFromOther>();
    set_processor<AddBlocks>();
    set_processor<BlastRunner>();
    set_processor<ImportBlastHits>();
    set_processor<AddBlastBlocks>();
    set_processor<ResolveBlast>();
    set_processor<ResolveAnchors>();
    set_processor<AnchorFinder>();
    set_processor<Filter>();
    set_processor<Connector>();
    set_processor<OriByMajority>();
    set_processor<StickBoundaries>();
    set_processor<OverlapsResolver>();
    set_processor<OverlapsResolver2>();
    set_processor<CheckNoOverlaps>();
    set_processor<Joiner>();
    set_processor<Union>();
    set_processor<OverlaplessUnion>();
    set_processor<Subtract>();
    set_processor<BlocksExpander>();
    set_processor<FragmentsExpander>();
    set_processor<Pipe>();
    set_processor<CleanUp>();
    set_processor<Output>();
    set_processor<OutputPipe>();
    set_processor<PrintOverlaps>();
    set_processor<BlockInfo>();
    set_processor<Stats>();
    set_processor<IsPangenome>();
    set_processor<Consensus>();
    set_processor<UniqueNames>();
    set_processor<Rest>();
    set_processor<ExternalAligner>();
    set_processor<RemoveAlignment>();
    set_processor<ConSeq>();
    set_processor<DeConSeq>();
    set_processor<MoveGaps>();
    set_processor<CutGaps>();
    set_processor<MetaProcessor>();
}

bool Meta::has(const std::string& key) const {
    return map_.find(key) != map_.end();
}

ProcessorPtr Meta::get(const std::string& key) const {
    ReturnerMap::const_iterator it = map_.find(key);
    if (it == map_.end()) {
        throw Exception("No such proessor: " + key);
    }
    const ProcessorReturner& returner = it->second;
    ProcessorPtr processor = returner();
    BOOST_ASSERT(processor->key() == key);
    processor->set_meta(const_cast<Meta*>(this));
    return processor;
}

std::vector<std::string> Meta::keys() const {
    std::vector<std::string> result;
    BOOST_FOREACH (const ReturnerMap::value_type& key_and_func, map_) {
        result.push_back(key_and_func.first);
    }
    return result;
}

bool Meta::empty() const {
    return map_.empty();
}

void Meta::clear() {
    map_.clear();
}

}

