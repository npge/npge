/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "meta_lib.hpp"
#include "Meta.hpp"
#include "Processor.hpp"
#include "SequencesFromOther.hpp"
#include "In.hpp"
#include "AddGenes.hpp"
#include "BlastFinder.hpp"
#include "BlastRunner.hpp"
#include "ImportBlastHits.hpp"
#include "AddBlastBlocks.hpp"
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
#include "SelfOverlapsResolver.hpp"
#include "Joiner.hpp"
#include "MergeUnique.hpp"
#include "SplitRepeats.hpp"
#include "Union.hpp"
#include "Move.hpp"
#include "MoveUnchanged.hpp"
#include "Clear.hpp"
#include "OverlaplessUnion.hpp"
#include "OneByOne.hpp"
#include "Partition.hpp"
#include "FindGeneGroups.hpp"
#include "FindGeneConversion.hpp"
#include "Subtract.hpp"
#include "RemoveWithSameName.hpp"
#include "BlocksExpander.hpp"
#include "FragmentsExpander.hpp"
#include "FragmentsExtender.hpp"
#include "SplitExtendable.hpp"
#include "Hash.hpp"
#include "FileRemover.hpp"
#include "MkDir.hpp"
#include "FileCopy.hpp"
#include "Pipe.hpp"
#include "CleanUp.hpp"
#include "Output.hpp"
#include "OutputPipe.hpp"
#include "FragmentDistance.hpp"
#include "FragmentFinder.hpp"
#include "OverlapFinder.hpp"
#include "BlockFinder.hpp"
#include "PrintTree.hpp"
#include "ConsensusTree.hpp"
#include "FindBSA.hpp"
#include "ChrBSA.hpp"
#include "PrintBSA.hpp"
#include "InputBSA.hpp"
#include "FastaBSA.hpp"
#include "ExactStemBSA.hpp"
#include "PrintOverlaps.hpp"
#include "PrintPartition.hpp"
#include "PrintGeneGroups.hpp"
#include "PrintMutations.hpp"
#include "MutationsSequences.hpp"
#include "FindLowSimilar.hpp"
#include "BlockInfo.hpp"
#include "Stats.hpp"
#include "Info.hpp"
#include "OldMakePrePangenome.hpp"
#include "OldMakePangenome.hpp"
#include "AreBlocksGood.hpp"
#include "IsPangenome.hpp"
#include "Consensus.hpp"
#include "UniqueNames.hpp"
#include "Rest.hpp"
#include "ExternalAligner.hpp"
#include "MultipleAligner.hpp"
#include "SimilarAligner.hpp"
#include "DummyAligner.hpp"
#include "MetaAligner.hpp"
#include "RemoveAlignment.hpp"
#include "RemoveNames.hpp"
#include "MarkNonWeak.hpp"
#include "RemoveWeak.hpp"
#include "LinkEqualFragments.hpp"
#include "ConSeq.hpp"
#include "DeConSeq.hpp"
#include "MoveGaps.hpp"
#include "CutGaps.hpp"
#include "FixEnds.hpp"
#include "Align.hpp"
#include "ReAlign.hpp"
#include "MetaProcessor.hpp"
#include "TrySmth.hpp"
#include "AllProcessors.hpp"
#include "GetData.hpp"
#include "ReplaceNames.hpp"

namespace npge {

void add_meta_lib(Meta* meta) {
    meta->set_processor<Processor>();
    meta->set_processor<SequencesFromOther>();
    meta->set_processor<In>();
    meta->set_processor<AddGenes>();
    meta->set_processor<BlastFinder>();
    meta->set_processor<BlastRunner>();
    meta->set_processor<ImportBlastHits>();
    meta->set_processor<AddBlastBlocks>();
    meta->set_processor<AnchorFinder>();
    meta->set_processor<LiteFilter>();
    meta->set_processor<Filter>();
    meta->set_processor<Stem>();
    meta->set_processor<SameChr>();
    meta->set_processor<Connector>();
    meta->set_processor<OriByMajority>();
    meta->set_processor<StickBoundaries>();
    meta->set_processor<OverlapsResolver>();
    meta->set_processor<OverlapsResolver2>();
    meta->set_processor<CheckNoOverlaps>();
    meta->set_processor<SelfOverlapsResolver>();
    meta->set_processor<Joiner>();
    meta->set_processor<MergeUnique>();
    meta->set_processor<SplitRepeats>();
    meta->set_processor<Union>();
    meta->set_processor<Move>();
    meta->set_processor<MoveUnchanged>();
    meta->set_processor<Clear>();
    meta->set_processor<OverlaplessUnion>();
    meta->set_processor<OneByOne>();
    meta->set_processor<Partition>();
    meta->set_processor<FindGeneGroups>();
    meta->set_processor<FindGeneConversion>();
    meta->set_processor<Subtract>();
    meta->set_processor<RemoveWithSameName>();
    meta->set_processor<BlocksExpander>();
    meta->set_processor<FragmentsExpander>();
    meta->set_processor<FragmentsExtender>();
    meta->set_processor<SplitExtendable>();
    meta->set_processor<Hash>();
    meta->set_processor<FileRemover>();
    meta->set_processor<MkDir>();
    meta->set_processor<FileCopy>();
    meta->set_processor<Pipe>();
    meta->set_processor<CleanUp>();
    meta->set_processor<Output>();
    meta->set_processor<OutputPipe>();
    meta->set_processor<FragmentDistance>();
    meta->set_processor<FragmentFinder>();
    meta->set_processor<OverlapFinder>();
    meta->set_processor<BlockFinder>();
    meta->set_processor<PrintTree>();
    meta->set_processor<ConsensusTree>();
    meta->set_processor<FindBSA>();
    meta->set_processor<ChrBSA>();
    meta->set_processor<PrintBSA>();
    meta->set_processor<InputBSA>();
    meta->set_processor<FastaBSA>();
    meta->set_processor<ExactStemBSA>();
    meta->set_processor<PrintOverlaps>();
    meta->set_processor<PrintPartition>();
    meta->set_processor<PrintGeneGroups>();
    meta->set_processor<PrintMutations>();
    meta->set_processor<MutationsSequences>();
    meta->set_processor<FindLowSimilar>();
    meta->set_processor<BlockInfo>();
    meta->set_processor<Stats>();
    meta->set_processor<Info>();
    meta->set_processor<OldMakePrePangenome>();
    meta->set_processor<OldMakePangenome>();
    meta->set_processor<AreBlocksGood>();
    meta->set_processor<IsPangenome>();
    meta->set_processor<Consensus>();
    meta->set_processor<UniqueNames>();
    meta->set_processor<Rest>();
    meta->set_processor<ExternalAligner>();
    meta->set_processor<MafftAligner>();
    meta->set_processor<MuscleAligner>();
    meta->set_processor<MultipleAligner>();
    meta->set_processor<SimilarAligner>();
    meta->set_processor<DummyAligner>();
    meta->set_processor<MetaAligner>();
    meta->set_processor<RemoveAlignment>();
    meta->set_processor<RemoveNames>();
    meta->set_processor<MarkNonWeak>();
    meta->set_processor<RemoveWeak>();
    meta->set_processor<LinkEqualFragments>();
    meta->set_processor<ConSeq>();
    meta->set_processor<DeConSeq>();
    meta->set_processor<MoveGaps>();
    meta->set_processor<CutGaps>();
    meta->set_processor<FixEnds>();
    meta->set_processor<LiteAlign>();
    meta->set_processor<Align>();
    meta->set_processor<ReAlign>();
    meta->set_processor<MetaProcessor>();
    meta->set_processor<TrySmth>();
    meta->set_processor<AddingLoopBySize>();
    meta->set_processor<AllProcessors>();
    meta->set_processor<AllOptions>();
    meta->set_processor<GetData>();
    meta->set_processor<ReplaceNames>();
}

}

