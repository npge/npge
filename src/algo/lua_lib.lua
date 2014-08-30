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
    p:run()
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
                print(message)
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
            print('No such file: ' .. fname)
        end
        if not arg_has(arg, '-i') then
            return
        end
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
    p:set_max_loops(-1)
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
    p:set_max_loops(-1)
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

register_p('PostProcessing', function()
    local p = Pipe.new()
    p:set_name("Postprocess pangenome");
    p:add('AddBlocks', '--in-blocks=pangenome0.fasta')
    p:add('RemoveNames', '--remove-seqs-names:=0 '..
        '--remove-blocks-names:=1')
    p:add('UniqueNames')
    p:add('OutputPipe', '--out-file=pangenome.fasta')

    p:add('IsPangenome',
        '--out-is-pangenome=pangenome.isgood '..
        '--blast-cons-dst=pangenome-cons.fasta '..
        '--blast-hits-dst=pangenome-hits.blast '..
        'blast-hits=blast-hits joined=joined '..
        'all-blast-hits=all-blast-hits '..
        'non-internal-hits=non-internal-hits')

    // TODO DRY

    p:add('UniqueNames', 'target=blast-hits')
    p:add('Output', 'prefix|blast-good-fasta- '..
        'target=blast-hits '..
        '--blast-good-fasta-file=pangenome-blast-good.fasta')
    p:add('BlockInfo', 'prefix|blast-good-bi- '..
        'target=blast-hits --blast-good-bi-count-seqs=1 '..
        '--blast-good-bi-file=pangenome-blast-good.bi')

    p:add('UniqueNames', 'target=all-blast-hits')
    p:add('Output', 'prefix|blast-all-fasta- '..
        'target=all-blast-hits '..
        '--blast-all-fasta-file=pangenome-blast-all.fasta')
    p:add('BlockInfo', 'prefix|blast-all-bi- '..
        'target=all-blast-hits '..
        '--blast-all-bi-count-seqs=1 '..
        '--blast-all-bi-file=pangenome-blast-all.bi')

    p:add('UniqueNames', 'target=non-internal-hits')
    p:add('Output', 'prefix|blast-non-internal-fasta- '..
        'target=non-internal-hits '..
        '--blast-non-internal-fasta-file='..
        'pangenome-blast-non-internal.fasta')
    p:add('BlockInfo', 'prefix|blast-non-internal-bi- '..
        '--blast-non-internal-bi-count-seqs=1 '..
        'target=non-internal-hits '..
        '--blast-non-internal-bi-file='..
        'pangenome-blast-non-internal.bi')

    p:add('UniqueNames', 'target=joined')
    p:add('Output', 'prefix|joined-fasta- target=joined '..
        '--joined-fasta-file=pangenome-joined.fasta')
    p:add('BlockInfo', 'prefix|joined-bi- target=joined '..
        '--joined-bi-count-seqs=1 '..
        '--joined-bi-file=pangenome-joined.bi')

    p:add('BlockInfo', 'prefix|bi- --bi-count-seqs=1 '..
        '--bi-file=pangenome.bi')
    p:add('Info', '--out-stats=pangenome.info')
    p:add('PrintMutations', 'prefix|mut- '..
        '--mut-file=pangenome.mutations')
    p:add('Union', 'target=stem other=target')
    p:add('Stem', 'target=stem --exact:=1')
    p:add('MutationsSequences', '--mutation-distance=1 '..
        'target=mut other=stem')
    p:add('Output', 'prefix|mut-fasta- target=mut '..
        '--mut-fasta-dump-seq:=1 --mut-fasta-dump-block:=0 '..
        '--mut-fasta-file=pangenome-mutations.fasta')
    p:add('Output', 'prefix|mut-blocks- target=mut '..
        '--mut-fasta-dump-seq:=1 --mut-fasta-dump-block:=1 '..
        '--mut-fasta-file='..
        'pangenome-mutations-with-blocks.fasta')
    p:add('PrintTree', '--tree-file=pangenome.trees')
    p:add('ConsensusTree', 'prefix|nj- --nj-tree-method:=nj '..
        '--nj-out-consensus-tree=pangenome.constree-nj '..
        '--nj-tree-pseudo-leafs=1 '..
        '--nj-bootstrap-print=before-length')
    p:add('ConsensusTree', 'prefix|upgma- '..
        '--upgma-tree-method:=upgma '..
        '--upgma-out-consensus-tree='..
        'pangenome.constree-upgma '..
        '--upgma-tree-pseudo-leafs=1 '..
        '--upgma-bootstrap-print=before-length')
    p:add('FragmentDistance',
        '--distance-file=pangenome.distances')
    p:add('Hash', '--hash-file=pangenome.hash')

    p:add('MergeUnique')
    p:add('MetaAligner')
    p:add('OutputPipe', 'prefix|merged- '..
        '--merged-out-file=pangenome-merged.fasta')
    p:add('BlockInfo', 'prefix|merged-bi- '..
        '--merged-bi-count-seqs=1 '..
        '--merged-bi-file=pangenome-merged.bi')
    p:add('Info', 'prefix|merged- '..
        '--merged-out-stats=pangenome-merged.info')
    p:add('Hash', 'prefix|merged- '..
        '--merged-hash-stats=pangenome-merged.hash')
    p:add('PrintMutations', 'prefix|merged-mut- '..
        '--merged-mut-file=pangenome-merged.mutations')
    p:add('Union', 'target=merged-stem other=target')
    p:add('Stem', 'target=merged-stem --exact:=1')
    p:add('MutationsSequences', '--mutation-distance=1 '..
        'target=merged-mut other=merged-stem')
    p:add('Output', 'prefix|merged-mut-fasta- '..
        'target=merged-mut '..
        '--merged-mut-fasta-dump-seq:=1 '..
        '--merged-mut-fasta-dump-block:=0 '..
        '--merged-mut-fasta-file='..
        'pangenome-merged-mutations.fasta')
    p:add('Output', 'prefix|merged-mut-blocks- '..
        'target=merged-mut '..
        '--merged-mut-blocks-dump-seq:=1 '..
        '--merged-mut-blocks-dump-block:=1 '..
        '--merged-mut-blocks-file='..
        'pangenome-merged-mutations-with-blocks.fasta')

    p:add('ChrBSA')
    p:add('PrintBSA', '--out-bsa=pangenome-merged.bsa')
    p:add('ExactStemBSA')
    p:add('PrintBSA', 'prefix|exact-stem- '..
        '--exact-stem-out-bsa=pangenome-merged-exact-stem.bsa')

    p:add('SplitRepeats', 'target=split other=target')
    p:add('UniqueNames', 'target=split')
    p:add('Output', 'prefix|split-fasta- target=split '..
        '--split-fasta-file=pangenome-merged-split.fasta')
    p:add('BlockInfo', 'prefix|split-bi- target=split '..
        '--split-bi-count-seqs=1 '..
        '--split-bi-file=pangenome-merged-split.bi')

    p:add('FindLowSimilar', 'target=low other=target')
    p:add('UniqueNames', 'target=low')
    p:add('Output', 'prefix|low-fasta- target=low '..
        '--low-fasta-file=pangenome-merged-low.fasta')
    p:add('BlockInfo', 'prefix|low-bi- target=low '..
        '--low-bi-count-seqs=1 '..
        '--low-bi-file=pangenome-merged-low.bi')

    return p
end)

);

