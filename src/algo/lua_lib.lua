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
        end;
    end
end

function get(key)
    return meta:get_opt(key)
end

function set(key, value, description)
    if meta:has_opt(key) and
        type(value) ~= type(meta:get_opt(key)) then
        print("Type mismatch, key: " .. key)
        return;
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

function arg_value(a, opt)
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
            local f, message = loadfile(fname)
            if f then
                f()
            else
                error(message)
            end
        elseif is_first_upper(fname) then
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
    end;
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
    p:declare_bs('target', 'Target block set')
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
    p:add('Filter')
    p:add('OriByMajority')
    p:add('MergeUnique')
    p:add('MetaAligner')
    p:add('Joiner')
    p:add('Rest', 'target=target other=target')
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
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=AnchorLoopFast')
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=JoinerP')
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=FragmentsExtender')
    p:add('Hash')
    p:add('Info')
    return p
end)

register_p('AnchorJoiner', function()
    local p = Pipe.new()
    p:set_max_loops(3)
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=AnchorLoop')
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=JoinerP')
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=FragmentsExtender')
    p:add('Hash')
    p:add('Info')
    return p
end)

register_p('AnchorBlastJoiner', function()
    local p = Pipe.new()
    p:set_max_loops(-1)
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=AddBlastBlocksToSelf')
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=JoinerP')
    p:add('Hash')
    p:add('TrySmth', '--smth-processor:=FragmentsExtender')
    p:add('Hash')
    p:add('Info')
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
    p:set_name("Postprocess pangenome");

    p:add('In', '--in-blocks=pangenome.bs')

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
    p:add('GetData', '--type:=genes --table=genomes.tsv '..
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

);

