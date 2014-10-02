/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#define NPGE_SCRIPT(...) #__VA_ARGS__

const char* meta_lua = NPGE_SCRIPT(

function simple_terminal()
    while true do
        io.write("> ")
        io.flush()
        local line = io.read()
        if not line then
            break
        end
        local chunk, load_error = loadstring(line)
        if chunk then
            local success, pcall_result = pcall(chunk)
            print(pcall_result)
        else
            print(load_error)
        end
    end
end

function get(key)
    return meta:get_opt(key)
end

function set(key, value, description)
    if meta:has_opt(key) and
        type(value) ~= type(meta:get_opt(key)) then
        print("Type mismatch, key: " .. key)
        return
    end
    if description then
        meta:set_opt(key, value, description)
    else
        meta:set_opt(key, value)
    end
end

function new_p(name)
    return meta:get_plain(name)
end

function run(name, opts)
    local p = new_p(name)
    opts = opts or ""
    p:set_options(opts, meta:placeholder_processor())
    p:run()
    Processor.delete(p)
end

function arg_has(a, pattern)
    for k, v in next, a do
        if v == pattern then
            return true
        end
    end
    return false
end

function arg_value(a, pattern)
    for i, v in ipairs(a) do
        if v == pattern and i < #a then
            return a[i + 1]
        end
    end
end

function file_exists(name)
    local f = io.open(name, "r")
    if f ~= nil then
        io.close(f)
        return true
    else
        return false
    end
end

function is_first_upper(name)
    local first_letter = name:sub(1, 1)
    return string.upper(first_letter) == first_letter
end

// based on http://lua-users.org/wiki/SplitJoin
function string:split(sep, nMax, plain)
    if not sep then
        sep = '%s+'
    end
    assert(sep ~= '')
    assert(nMax == nil or nMax >= 1)
    local aRecord = {}
    if self:len() > 0 then
        nMax = nMax or -1
        local nField=1 nStart=1
        local nFirst,nLast = self:find(sep, nStart, plain)
        while nFirst and nMax ~= 0 do
            aRecord[nField] = self:sub(nStart, nFirst-1)
            nField = nField+1
            nStart = nLast+1
            nFirst,nLast = self:find(sep, nStart, plain)
            nMax = nMax-1
        end
        aRecord[nField] = self:sub(nStart)
    end
    return aRecord
end

function string:starts_with(prefix)
   return self:sub(1, prefix:len()) == prefix
end

function string:ends_with(suffix)
   return suffix == '' or self:sub(-suffix:len()) == suffix
end

function read_config(fname0)
    local fname = resolve_home_dir(fname0)
    if #fname == 0 or not file_exists(fname) then
        return
    end
    local keys_before = {}
    for key, value in next, _G do
        keys_before[key] = true
    end
    dofile(fname)
    local new_keys = {}
    for key, value in next, _G do
        if keys_before[key] == nil and meta:has_opt(key) then
            set(key, value)
            table.insert(new_keys, key)
        end
    end
    for i, key in next, new_keys do
        _G[key] = nil
    end
end

function run_main(name, opts)
    local p = new_p(name)
    opts = opts or ""
    p:set_options(opts, meta:placeholder_processor())
    p:apply_vector_options(arg)
    if arg_has(arg, '-h') or arg_has(arg, '--help') then
        p:print_help()
    elseif arg_has(arg, '--tree') then
        p:print_tree()
    else
        p:run()
    end
    Processor.delete(p)
end

function main()
    assert(arg)
    if arg_has(arg, '-g') then
        local conf = arg_value(arg, '-g')
        if conf then
            meta:print_config(conf)
        else
            meta:print_config()
        end
        return
    end
    // remove argv[0] from argument list
    table.remove(arg, 1)
    local fname = arg[1]
    if fname and fname:sub(1, 1) ~= '-' then
        // remove fname from arguments list
        table.remove(arg, 1)
        if file_exists(fname) then
            io.stderr:write("Running script " .. fname .. "\n")
            io.stderr:flush()
            local f, message = loadfile(fname)
            if f then
                f()
            else
                error(message)
            end
        elseif is_first_upper(fname) then
            io.stderr:write("Running " .. fname .. "..\n")
            io.stderr:flush()
            local pp = {}
            for p in string.gmatch(fname, '([^,]+)') do
                table.insert(pp, p)
            end
            if #pp == 1 then
                run_main(pp[1])
            else
                register_p('MainPipe', function()
                    local main_pipe = Pipe.new()
                    for i, p in next, pp do
                        main_pipe:add(p)
                    end
                    return main_pipe
                end)
                run_main('MainPipe')
            end
        else
            error('No such file: ' .. fname)
        end
        io.stderr:write(".. done\n")
        io.stderr:flush()
        if not arg_has(arg, '-i') then
            return
        end
    end
    run_main('Processor')
    if arg_has(arg, '-h') or arg_has(arg, '--help') or
            arg_has(arg, '--tree') then
        return
    end
    if terminal then
        terminal()
    else
        simple_terminal()
    end
end

function block_set()
    return meta:placeholder_processor():block_set()
end

function other()
    return meta:placeholder_processor():other()
end

function register_p(name, returner)
    local returner1 = function()
        local p = returner()
        p:set_key(name)
        return p
    end
    meta:set_returner(returner1, name)
end

// pipes

register_p('RemoveMinorBlocks', function()
    local p = LuaProcessor.new()
    p:declare_bs('target', 'Target blockset')
    p:set_name('Remove minor blocks (name starts with "m")')
    p:set_action(function()
        local bs = p:block_set()
        for _, block in next, bs:blocks() do
            local name = block:name()
            if name ~= "" and name:sub(1, 1) == 'm' then
                bs:erase(block)
            end
        end
    end)
    return p
end)

register_p('JoinerP', function()
    local p = Pipe.new()
    p:add('LiteFilter')
    p:add('OriByMajority')
    p:add('Rest', 'target=target other=target')
    p:add('MergeUnique')
    p:add('MetaAligner')
    p:add('Joiner')
    p:add('Rest', 'target=target other=target')
    return p
end)

register_p('MergeAndJoin', function()
    local p = Pipe.new()
    p:add('Rest', 'target=target other=target')
    p:add('MergeUnique')
    p:add('MetaAligner')
    p:add('Joiner')
    return p
end)

register_p('LiteJoinerP', function()
    local p = Pipe.new()
    p:add('LiteFilter')
    p:add('OriByMajority')
    p:add('Joiner')
    return p
end)

register_p('ExtendAndAlign', function()
    local p = Pipe.new()
    p:add('FragmentsExtender', '--extend-length-portion:=0.5')
    p:add('Align')
    return p
end)

register_p('ExtendLoop', function()
    local p = Pipe.new()
    p:set_max_loops(-1)
    p:add('MoveUnchanged', 'target=unchanged other=target')
    p:add('ExtendAndAlign')
    p:add('Move', 'target=target other=unchanged')
    p:add('AddingLoopBySize', 'target=ol other=target')
    p:add('Clear', 'target=target --clear-seqs:=1 no_options')
    p:add('Move', 'target=target other=ol')
    p:add('Clear', 'target=ol --clear-seqs:=1 no_options')
    return p
end)

register_p('ExtendAndFix', function()
    local p = Pipe.new()
    p:add('FragmentsExtender', '--extend-length-portion:=0.5')
    p:add('FixEnds')
    return p
end)

register_p('ExtendLoopFast', function()
    local p = Pipe.new()
    p:set_max_loops(-1)
    p:add('MoveUnchanged', 'target=unchanged other=target')
    p:add('ExtendAndFix')
    p:add('Move', 'target=target other=unchanged')
    p:add('OverlaplessUnion',
        'target=ol other=target --ou-move:=1')
    p:add('Clear', 'target=target --clear-seqs:=0 no_options')
    p:add('Move', 'target=target other=ol')
    p:add('Clear', 'target=ol --clear-seqs:=1 no_options')
    return p
end)

register_p('AnchorLoop', function()
    local p = Pipe.new()
    p:add('Filter')
    p:add('Rest', 'target=target other=target')
    p:add('ConSeq', 'target=cons other=target')
    p:add('AnchorFinder', 'target=cons')
    p:add('MoveUnchanged', 'target=null other=cons')
    p:add('Clear', 'target=null')
    p:add('DummyAligner', 'target=cons')
    p:add('UniqueNames', 'target=cons no_options')
    p:add('Union', 'target=anchors other=cons')
    p:add('ExtendAndAlign', 'target=cons')
    p:add('RemoveWithSameName', 'target=anchors other=cons')
    p:add('SplitExtendable', 'other=anchors target=cons')
    p:add('RemoveNames', 'target=cons '..
        '--remove-seqs-names:=0 no_options')
    p:add('DeConSeq', 'target=deconseq other=cons')
    p:add('ExtendLoop', 'target=cons')
    p:add('ExtendLoop', 'target=deconseq')
    p:add('DeConSeq', 'target=target other=cons')
    p:add('Align')
    p:add('Move', 'target=target other=deconseq')
    p:add('Clear', 'target=cons --clear-seqs:=1 no_options')
    p:add('Clear', 'target=anchors --clear-seqs:=1 no_options')
    p:add('Clear', 'target=deconseq --clear-seqs:=1 no_options')
    return p
end)

register_p('AnchorLoopFast', function()
    local p = Pipe.new()
    p:add('Filter')
    p:add('Rest', 'target=target other=target')
    p:add('ConSeq', 'target=cons other=target')
    p:add('AnchorFinder', 'target=cons')
    p:add('MoveUnchanged', 'target=null other=cons')
    p:add('Clear', 'target=null')
    p:add('DummyAligner', 'target=cons')
    p:add('ExtendAndAlign', 'target=cons')
    p:add('ExtendLoopFast', 'target=cons')
    p:add('DeConSeq', 'target=target other=cons')
    p:add('Align')
    p:add('Clear', 'target=cons --clear-seqs:=1 no_options')
    return p
end)

register_p('AddBlastBlocksToSelf', function()
    local p = Pipe.new()
    p:add('Filter')
    p:add('Rest', 'target=target other=target')
    p:add('AddBlastBlocks', 'target=target other=target')
    return p
end)

register_p('AnchorJoinerFast', function()
    local p = Pipe.new()
    p:set_max_loops(4)
    p:add('TrySmth', '--smth-processor:=AnchorLoopFast')
    p:add('TrySmth', '--smth-processor:=JoinerP')
    p:add('TrySmth', '--smth-processor:=FragmentsExtender')
    p:add('Info', '--short-stats:=true')
    return p
end)

register_p('AnchorJoiner', function()
    local p = Pipe.new()
    p:set_max_loops(3)
    p:add('TrySmth', '--smth-processor:=AnchorLoop')
    p:add('TrySmth', '--smth-processor:=JoinerP')
    p:add('TrySmth', '--smth-processor:=FragmentsExtender')
    p:add('Info', '--short-stats:=true')
    return p
end)

register_p('FragmentsExtender3', function()
    local p = new_p('FragmentsExtender')
    p:set_opt_value('extend-length', 3)
    return p
end)

register_p('ShortExtender', function()
    local p = Pipe.new()
    p:set_max_loops(-1)
    p:add('TrySmth', '--smth-processor:=FragmentsExtender3')
    return p
end)

register_p('AnchorBlastJoiner', function()
    local p = Pipe.new()
    p:set_max_loops(-1)
    p:add('TrySmth', '--smth-processor:=AddBlastBlocksToSelf')
    p:add('TrySmth', '--smth-processor:=JoinerP')
    p:add('TrySmth', '--smth-processor:=MergeAndJoin')
    p:add('TrySmth', '--smth-processor:=LiteJoinerP')
    p:add('TrySmth', '--smth-processor:=FragmentsExtender')
    p:add('TrySmth', '--smth-processor:=ShortExtender')
    p:add('Info', '--short-stats:=true')
    return p
end)

register_p('MakePangenome', function()
    local p = Pipe.new()
    p:add('AnchorJoinerFast')
    p:add('AnchorJoiner')
    p:add('AnchorBlastJoiner')
    p:add('RemoveNames', '--remove-seqs-names:=0 '..
        '--remove-blocks-names:=1')
    p:add('UniqueNames')
    p:add('Rest', 'target=target other=target')
    p:add('MergeUnique')
    p:add('MetaAligner')
    return p
end)

register_p('MergeStemJoin', function()
    local p = Pipe.new()
    p:add('Rest', 'target=target other=target')
    p:add('MergeUnique')
    p:add('MetaAligner')
    p:add('Stem', '--exact:=1')
    p:add('Joiner')
    return p
end)

function un_out_info(p, target, base_name)
    p:add('UniqueNames', 'target=' .. target)
    p:add('Output', 'target=' .. target ..
        ' --out-file:=' .. base_name .. '.bs')
    p:add('BlockInfo', 'target=' .. target ..
        ' --info-count-seqs:=1 '..
        '--info-file:=' .. base_name .. '.bi')
end

register_p('PostProcessing', function()
    local p = Pipe.new()
    p:set_name("Postprocess pangenome")

    p:add('In', '--in-blocks:=pangenome.bs')
    p:add('SequencesFromOther', 'target=features other=target')
    p:add('In', 'target=features --in-blocks:=features.bs')

    p:add('BlockInfo', '--info-count-seqs:=1 '..
        '--info-file:=pangenome.bi')
    p:add('Info', '--out-stats=pangenome.info')
    p:add('Hash', '--hash-file=pangenome.hash')

    p:add('Union', 'target=stem other=target')
    p:add('Stem', 'target=stem --exact:=1')

    // check

    p:add('MkDir', '--dirname:=check')

    p:add('IsPangenome',
        '--out-is-pangenome=check/isgood '..
        '--blast-cons-dst=check/consensuses.fasta '..
        '--blast-hits-dst=check/hits.blast '..
        'blast-hits=blast-hits joined=joined '..
        'all-blast-hits=all-blast-hits '..
        'non-internal-hits=non-internal-hits')

    un_out_info(p, 'blast-hits', 'check/good-blast-hits')
    un_out_info(p, 'all-blast-hits', 'check/all-blast-hits')
    un_out_info(p, 'non-internal-hits',
                'check/non-internal-hits')
    un_out_info(p, 'joined', 'check/joined')

    // mutations

    p:add('MkDir', '--dirname:=mutations')

    p:add('PrintMutations', '--file:=mutations/mut.tsv')
    p:add('MutationsSequences', '--mutation-distance=1 '..
        'target=mut other=stem')
    p:add('Output', 'target=mut '..
        '--out-dump-seq:=1 --out-dump-block:=0 '..
        '--out-file:=mutations/mutseq.fasta')
    p:add('Output', 'target=mut '..
        '--out-dump-seq:=1 --out-dump-block:=1 '..
        '--out-file:=mutations/mutseq-with-blocks.bs')

    // trees

    p:add('MkDir', '--dirname:=trees')

    p:add('PrintTree', '--tree-file=trees/all_trees.tsv')
    p:add('ConsensusTree', 'prefix|nj- --nj-tree-method:=nj '..
        '--nj-out-consensus-tree=trees/nj-constree.tre '..
        '--nj-out-branch=trees/nj.branch '..
        '--nj-tree-pseudo-leafs=1 '..
        '--nj-bootstrap-print=before-length')
    p:add('ConsensusTree', 'prefix|upgma- '..
        '--upgma-tree-method:=upgma '..
        '--upgma-out-consensus-tree='..
        'trees/upgma-constree.tre '..
        '--upgma-out-branch=trees/upgma.branch '..
        '--upgma-tree-pseudo-leafs=1 '..
        '--upgma-bootstrap-print=before-length')
    p:add('FragmentDistance',
        '--distance-file=trees/distances.tsv')

    // genes

    p:add('MkDir', '--dirname:=genes')

    p:add('FindGoodGeneGroups',
          'target=ggg pangenome=target features=features')
    p:add('UniqueNames', 'target=ggg')
    p:add('Output', 'target=ggg --out-file:=genes/good.bs '..
          '--out-export-contents:=0')

    // BSA

    p:add('ChrBSA')
    p:add('PrintBSA', '--out-bsa:=pangenome.bsa')
    p:add('ExactStemBSA')
    p:add('PrintBSA', '--out-bsa:=pangenome-stem.bsa')

    // split

    p:add('SplitRepeats', 'target=split other=target')
    un_out_info(p, 'split', 'split')

    // low

    p:add('FindLowSimilar', 'target=low other=target')
    un_out_info(p, 'low', 'low')

    return p
end)

// shortcuts

register_p('GetFasta', function()
    local p = Pipe.new()
    p:add('GetData', '--type:=fasta --table=genomes.tsv '..
        '--data=genomes-raw.fasta')
    return p
end)

register_p('GetGenes', function()
    local p = Pipe.new()
    p:add('GetData', '--type:=features --table=genomes.tsv '..
        '--data=features.embl')
    return p
end)

register_p('Rename', function()
    local p = Pipe.new()
    p:add('In', '--in-blocks=genomes-raw.fasta')
    p:add('ReplaceNames', '--table=genomes.tsv')
    p:add('Output',
        '--out-dump-seq:=1 --out-file=genomes-renamed.fasta')
    return p
end)

register_p('ExtractGenes', function()
    local p = Pipe.new()
    p:add('In', '--in-blocks=genomes-renamed.fasta')
    p:add('AddGenes', '--in-genes=features.embl')
    p:add('Output', '--out-file=features.bs '..
        '--out-export-contents:=0')
    return p
end)

register_p('Prepare', function()
    local p = Pipe.new()
    p:add('GetFasta', '--data:=genomes-raw.fasta')
    p:add('Rename')
    p:add('GetGenes', '--data:=features.embl')
    p:add('AddGenes', '--in-genes=features.embl')
    p:add('Output', '--out-file:=features.bs '..
        '--out-export-contents:=0')
    return p
end)

register_p('Pangenome', function()
    local p = Pipe.new()
    p:add('In', '--in-blocks=genomes-renamed.fasta')
    p:add('MakePangenome')
    p:add('OutputPipe', '--out-file=pangenome.bs')
    return p
end)

register_p('CheckPangenome', function()
    local p = Pipe.new()
    p:add('In', '--in-blocks=pangenome.bs')
    p:add('IsPangenome')
    return p
end)

register_p('MakeSubPangenome', function()
    local p = LuaProcessor.new()
    p:declare_bs('target', 'Where to write sub-pangenome')
    p:declare_bs('other', 'Input pangenome')
    p:add_opt('genomes', 'Remaining genomes', {}, true)
    p:set_name('Make sub-pangenome on set of genomes')
    p:set_action(function()
        local genomes = {}
        for _, genome in pairs(p:opt_value('genomes')) do
            genomes[genome] = 1
        end
        for _, seq in pairs(p:other():seqs()) do
            if genomes[seq:genome()] == 1 then
                p:block_set():add_sequence(seq)
            end
        end
        for _, block in pairs(p:other():blocks()) do
            local new_block = Block.new()
            for _, fragment in pairs(block:fragments()) do
                if genomes[fragment:seq():genome()] == 1 then
                    new_block:insert(fragment:clone())
                end
            end
            if new_block:size() == 0 then
                Block.delete(new_block)
            else
                if new_block:size() == 1 then
                    new_block:front():set_row(nil)
                else
                    remove_pure_gap_columns(new_block)
                end
                p:block_set():insert(new_block)
            end
        end
    end)
    return p
end)

register_p('SubPangenome', function()
    local p = Pipe.new()
    p:add('In', 'target=other --in-blocks=pangenome.bs')
    p:add('Filter', 'target=other --find-subblocks:=0')
    p:add('MakeSubPangenome')
    p:add('Joiner')
    p:add('Rest', 'target=target other=target')
    p:add('MergeUnique')
    p:add('MetaAligner')
    p:add('OutputPipe')
    return p
end)

register_p('FindGoodGeneGroups', function()
    local p = LuaProcessor.new()
    p:set_name("Find groups of CDS's located equaly in one " ..
               "stem block")
    p:declare_bs('target', 'Where to write good gene groups')
    p:declare_bs('pangenome', 'pangenome')
    p:declare_bs('features', 'features')
    p:set_action(function()
        local fc = VectorFc()
        fc:add_bs(p:get_bs('pangenome'))
        fc:prepare()
        // key: block_pos1_pos2
        // value: array of gene fragments
        local groups = {}
        for _, gene in pairs(p:get_bs('features'):blocks()) do
            if gene:size() == 1
                    and gene:name():starts_with('CDS') then
                local gene_f = gene:front()
                local p_ff = fc:find_overlap_fragments(gene_f)
                if #p_ff == 1 then
                    local p_f = p_ff[1]
                    local p_b = p_f:block()
                    if p_b:name():starts_with('s') then
                        local block_pos = function(s_pos)
                            local f_pos = p_f:seq_to_frag(s_pos)
                            return p_f:block_pos(f_pos)
                        end
                        local bp1 = block_pos(gene_f:min_pos())
                        local bp2 = block_pos(gene_f:max_pos())
                        local bpa = math.min(bp1, bp2)
                        local bpb = math.max(bp1, bp2)
                        local key = p_b:name() .. bpa .. bpb
                        if groups[key] == nil then
                            groups[key] = {}
                        end
                        table.insert(groups[key], gene_f)
                    end
                end
            end
        end
        local genomes_n = p:get_bs('pangenome'):genomes_number()
        for _, group in pairs(groups) do
            if #group == genomes_n then
                local new_block = Block.new()
                new_block:set_weak(true)
                for _, gene_f in pairs(group) do
                    new_block:insert(gene_f)
                end
                p:get_bs('target'):insert(new_block)
            end
        end
    end)
    return p
end)

register_p('Upstreams', function()
    local p = LuaProcessor.new()
    p:set_name("Replace each block in 'other' with a block " ..
               "in 'target' such that its fragments are " ..
               "upstreams of fragments of original blocks")
    p:declare_bs('target', 'Where to write upstream blocks')
    p:declare_bs('other', 'Source blocks')
    p:add_gopt('upstream-length', 'Length of upstream',
               'UPSTREAM_LENGTH')
    p:set_action(function()
        local l = p:opt_value('upstream-length')
        for _, block in pairs(p:other():blocks()) do
            local new_block = Block.new()
            new_block:set_name('upstream' .. block:name())
            for _, f in pairs(block:fragments()) do
                local new_f = nil
                if f:ori() == 1 and f:min_pos() - l > 0 then
                    new_f = Fragment.new(f:seq())
                    new_f:set_min_pos(f:min_pos() - l)
                    new_f:set_max_pos(f:min_pos() - 1)
                elseif f:ori() == -1 and
                        f:max_pos() + l < f:seq():size() then
                    new_f = Fragment.new(f:seq())
                    new_f:set_ori(-1)
                    new_f:set_min_pos(f:max_pos() + 1)
                    new_f:set_max_pos(f:max_pos() + l)
                end
                if new_f then
                    new_block:insert(new_f)
                end
            end
            if not new_block:empty() then
                p:get_bs('target'):insert(new_block)
            else
                Block.delete(new_block)
            end
        end
    end)
    return p
end)

);

