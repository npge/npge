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
#include "AddGenes.hpp"
#include "BlastFinder.hpp"
#include "BlastRunner.hpp"
#include "ImportBlastHits.hpp"
#include "AddBlastBlocks.hpp"
#include "ResolveBlast.hpp"
#include "ResolveAnchors.hpp"
#include "AnchorFinder.hpp"
#include "Filter.hpp"
#include "Stem.hpp"
#include "SameChr.hpp"
#include "Connector.hpp"
#include "OriByMajority.hpp"
#include "StickBoundaries.hpp"
#include "OverlapsResolver.hpp"
#include "OverlapsResolver2.hpp"
#include "CheckNoOverlaps.hpp"
#include "Joiner.hpp"
#include "Union.hpp"
#include "Clear.hpp"
#include "OverlaplessUnion.hpp"
#include "OneByOne.hpp"
#include "Partition.hpp"
#include "FindGeneGroups.hpp"
#include "FindGeneConversion.hpp"
#include "Subtract.hpp"
#include "BlocksExpander.hpp"
#include "FragmentsExpander.hpp"
#include "Hash.hpp"
#include "FileRemover.hpp"
#include "Pipe.hpp"
#include "CleanUp.hpp"
#include "Output.hpp"
#include "OutputPipe.hpp"
#include "FragmentDistance.hpp"
#include "PrintTree.hpp"
#include "ConsensusTree.hpp"
#include "PrintOverlaps.hpp"
#include "PrintPartition.hpp"
#include "PrintGeneGroups.hpp"
#include "PrintMutations.hpp"
#include "BlockInfo.hpp"
#include "Stats.hpp"
#include "Info.hpp"
#include "MakePrePangenome.hpp"
#include "MakePangenome.hpp"
#include "IsPangenome.hpp"
#include "Consensus.hpp"
#include "UniqueNames.hpp"
#include "Rest.hpp"
#include "ExternalAligner.hpp"
#include "RemoveAlignment.hpp"
#include "RemoveNames.hpp"
#include "MarkNonWeak.hpp"
#include "LinkEqualFragments.hpp"
#include "ConSeq.hpp"
#include "DeConSeq.hpp"
#include "MoveGaps.hpp"
#include "CutGaps.hpp"
#include "Align.hpp"
#include "MetaProcessor.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

Meta::Meta() {
    placeholder_processor_ = new Processor;
    set_processor<AddSequences>();
    set_processor<SequencesFromOther>();
    set_processor<AddBlocks>();
    set_processor<AddGenes>();
    set_processor<BlastFinder>();
    set_processor<BlastRunner>();
    set_processor<ImportBlastHits>();
    set_processor<AddBlastBlocks>();
    set_processor<ResolveBlast>();
    set_processor<ResolveAnchors>();
    set_processor<AnchorFinder>();
    set_processor<Filter>();
    set_processor<Stem>();
    set_processor<SameChr>();
    set_processor<Connector>();
    set_processor<OriByMajority>();
    set_processor<StickBoundaries>();
    set_processor<OverlapsResolver>();
    set_processor<OverlapsResolver2>();
    set_processor<CheckNoOverlaps>();
    set_processor<Joiner>();
    set_processor<Union>();
    set_processor<Clear>();
    set_processor<OverlaplessUnion>();
    set_processor<OneByOne>();
    set_processor<Partition>();
    set_processor<FindGeneGroups>();
    set_processor<FindGeneConversion>();
    set_processor<Subtract>();
    set_processor<BlocksExpander>();
    set_processor<FragmentsExpander>();
    set_processor<Hash>();
    set_processor<FileRemover>();
    set_processor<Pipe>();
    set_processor<CleanUp>();
    set_processor<Output>();
    set_processor<OutputPipe>();
    set_processor<FragmentDistance>();
    set_processor<PrintTree>();
    set_processor<ConsensusTree>();
    set_processor<PrintOverlaps>();
    set_processor<PrintPartition>();
    set_processor<PrintGeneGroups>();
    set_processor<PrintMutations>();
    set_processor<BlockInfo>();
    set_processor<Stats>();
    set_processor<Info>();
    set_processor<MakePrePangenome>();
    set_processor<MakePangenome>();
    set_processor<IsPangenome>();
    set_processor<Consensus>();
    set_processor<UniqueNames>();
    set_processor<Rest>();
    set_processor<ExternalAligner>();
    set_processor<RemoveAlignment>();
    set_processor<RemoveNames>();
    set_processor<MarkNonWeak>();
    set_processor<LinkEqualFragments>();
    set_processor<ConSeq>();
    set_processor<DeConSeq>();
    set_processor<MoveGaps>();
    set_processor<CutGaps>();
    set_processor<Align>();
    set_processor<MetaProcessor>();
}

Meta::~Meta() {
    delete placeholder_processor_;
}

bool Meta::has(const std::string& key) const {
    return map_.find(key) != map_.end();
}

Processor* Meta::get_plain(const std::string& key) const {
    ReturnerMap::const_iterator it = map_.find(key);
    if (it == map_.end()) {
        throw Exception("No such proessor: " + key);
    }
    const ProcessorReturner& returner = it->second;
    Processor* processor = returner();
    BOOST_ASSERT(processor->key() == key);
    processor->set_meta(const_cast<Meta*>(this));
    return processor;
}

SharedProcessor Meta::get(const std::string& key) const {
    return SharedProcessor(get_plain(key));
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

std::string Meta::get_key_and_delete(const Processor* p) {
    std::string key = p->key();
    delete p;
    return key;
}

}

