--
-- NPG-explorer, Nucleotide PanGenome explorer
-- Copyright (C) 2014 Boris Nagaev
--
-- See the LICENSE file for terms of use.
--

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

string.extract_value = extract_value
extract_value = nil

function string.trim(self)
    local text = self:gsub("%s+$", "")
    text = text:gsub("^%s+", "")
    return text
end

if not unpack then
    unpack = table.unpack
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

-- based on http://lua-users.org/wiki/SplitJoin
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
    local fname = file.resolve_home_dir(fname0)
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
    -- remove argv[0] from argument list
    table.remove(arg, 1)
    local fname = arg[1]
    if fname and fname:sub(1, 1) ~= '-' then
        -- remove fname from arguments list
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
                _G.pp = pp
                register_p('MainPipe', function()
                    local main_pipe = Pipe.new()
                    for i, p in next, _G.pp do
                        main_pipe:add(p)
                    end
                    return main_pipe
                end)
                run_main('MainPipe')
                _G.pp = nil
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
    meta:set_returner(returner, name)
end

-- pipes

register_p('RemoveMinorBlocks', function()
    local p = LuaProcessor.new()
    p:declare_bs('target', 'Target blockset')
    p:set_name('Remove minor blocks (name starts with "m")')
    p:set_action(function(p)
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
    p:add('OutputPipe', '--out-file:=pre-pangenome.bs')
    return p
end)

register_p('AnchorJoiner', function()
    local p = Pipe.new()
    p:set_max_loops(3)
    p:add('TrySmth', '--smth-processor:=AnchorLoop')
    p:add('TrySmth', '--smth-processor:=JoinerP')
    p:add('TrySmth', '--smth-processor:=FragmentsExtender')
    p:add('Info', '--short-stats:=true')
    p:add('OutputPipe', '--out-file:=pre-pangenome.bs')
    return p
end)

register_p('FragmentsExtender3', function()
    local p = new_p('FragmentsExtender')
    p:set_opt_value('extend-length', 3)
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
    p:add('TrySmth', '--smth-processor:=FragmentsExtender3')
    p:add('Info', '--short-stats:=true')
    p:add('OutputPipe', '--out-file:=pre-pangenome.bs')
    return p
end)

register_p('ShortUniqueToMinor', function()
    local p = LuaProcessor.new()
    p:declare_bs('target', 'Target blockset')
    p:add_gopt('min-length', 'Min length of non-minor block',
               'MIN_LENGTH')
    p:set_name('Rename short unique blocks to minor')
    p:set_action(function(p)
        local min_length = p:opt_value('min-length')
        local bs = p:block_set()
        for _, block in next, bs:blocks() do
            if block:size() == 1 and
                    block:front():length() < min_length then
                block:set_name('m' .. block:id())
            end
        end
    end)
    return p
end)

register_p('Pangenome', function()
    local p = Pipe.new()
    p:add('AnchorJoinerFast')
    p:add('AnchorJoiner')
    p:add('AnchorBlastJoiner')
    p:add('RemoveNames', '--remove-seqs-names:=0 '..
        '--remove-blocks-names:=1')
    p:add('UniqueNames')
    p:add('Rest', 'target=target other=target')
    p:add('MergeUnique')
    p:add('ShortUniqueToMinor')
    p:add('MetaAligner')
    return p
end)

register_p('RenameBlocksToGlobal', function()
    local p = LuaProcessor.new()
    p:set_action(function(p)
        local bs = p:block_set()
        for _, block in pairs(bs:blocks()) do
            local name = block:name()
            if name:starts_with('s') then
                name = name:gsub('^s', 'g')
                block:set_name(name)
            end
        end
    end)
    return p
end)

register_p('MergeStemJoin', function()
    local p = Pipe.new()
    p:add('LiteFilter')
    p:add('RemoveNames', '--remove-seqs-names:=0')
    p:add('Stem', '--exact:=true')
    p:add('RemoveAlignment')
    p:add('Joiner')
    p:add('Rest', 'target=target other=target')
    p:add('MergeUnique')
    p:add('Joiner')
    p:add('Stem', '--exact:=true')
    p:add('Rest', 'target=target other=target')
    p:add('MergeUnique', "--both-neighbours:=false");
    p:add('Joiner')
    p:add('UniqueNames')
    p:add('RenameBlocksToGlobal')
    return p
end)

register_p('CheckNoRest', function()
    local p = LuaProcessor.new()
    p:set_action(function(p)
        local target_bs = p:block_set()
        local rest_bs = BlockSet.new()
        local rest_p = new_p('Rest')
        rest_p:set_bs('other', target_bs)
        rest_p:set_bs('target', rest_bs)
        rest_p:run()
        assert(rest_bs:empty())
    end)
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

    -- dirs

    p:add('MkDir', '--dirname:=check')
    p:add('MkDir', '--dirname:=mutations')
    p:add('MkDir', '--dirname:=trees')
    p:add('MkDir', '--dirname:=genes')
    p:add('MkDir', '--dirname:=global-blocks')
    p:add('MkDir', '--dirname:=extra-blocks')

    -- check

    p:add('IsPangenome',
        '--out-is-pangenome=check/isgood '..
        '--blast-hits-dst=check/hits.blast '..
        'blast-hits=blast-hits joined=joined '..
        'all-blast-hits=all-blast-hits '..
        'non-internal-hits=non-internal-hits')

    un_out_info(p, 'blast-hits', 'check/good-blast-hits')
    un_out_info(p, 'all-blast-hits', 'check/all-blast-hits')
    un_out_info(p, 'non-internal-hits',
                'check/non-internal-hits')
    un_out_info(p, 'joined', 'check/joined')

    -- mutations

    p:add('Consensus',
        '--cons-file:=mutations/consensuses.fasta')
    p:add('PrintMutations', '--file:=mutations/mut.tsv')
    p:add('MutationsSequences', '--mutation-distance=1 '..
        'target=mut other=stem')
    p:add('Output', 'target=mut '..
        '--out-dump-seq:=1 --out-dump-block:=0 '..
        '--out-file:=mutations/mutseq.fasta')
    p:add('Output', 'target=mut '..
        '--out-dump-seq:=1 --out-dump-block:=1 '..
        '--out-file:=mutations/mutseq-with-blocks.bs')

    -- trees

    p:add('GlobalTree',
        '--out-global-tree=trees/nj-global-tree.tre')

    -- genes

    p:add('FindGoodGeneGroups',
          'target=ggg pangenome=target features=features')
    p:add('UniqueNames', 'target=ggg')
    p:add('Output', 'target=ggg --out-file:=genes/good.bs '..
          '--out-export-contents:=0')
    p:add('Upstreams', 'target=ggg-upstreams other=ggg')
    p:add('Output', 'target=ggg-upstreams '..
          '--out-file:=genes/good-upstreams.bs')
    p:add('UniqueFragments', 'target=ggg-upstreams')
    p:add('Output', 'target=ggg-upstreams '..
          '--out-file:=genes/good-upstreams-unique.bs')

    -- BSA

    p:add('Union', 'target=global-blocks other=target')
    p:add('MergeStemJoin', 'target=global-blocks')
    p:add('CheckNoRest', 'target=global-blocks')
    p:add('ChrBSA', 'target=global-blocks')
    p:add('PrintBSA', [[target=global-blocks
        --out-bsa:=global-blocks/pangenome.ba]])
    p:add('Output', [[target=global-blocks
        --out-export-contents:=0
        --out-file:=global-blocks/pangenome.bs]])
    p:add('LocalBSA', 'target=target other=global-blocks')
    p:add('PrintBSA', "--out-bsa:=pangenome.ba")

    -- split

    p:add('SplitRepeats', 'target=split other=target')
    un_out_info(p, 'split', 'extra-blocks/split')

    -- low

    p:add('FindLowSimilar', 'target=low other=target')
    un_out_info(p, 'low', 'extra-blocks/low')

    return p
end)

-- shortcuts

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

register_p('Examine', function()
    local p = Pipe.new()
    p:add('In', '--in-blocks=genomes-renamed.fasta')
    p:add('MkDir', '--dirname:=examine')
    p:add('GenomeLengths')
    p:add('DraftAndRecommend')
    return p
end)

register_p('MakePangenome', function()
    local p = Pipe.new()
    p:add('In', '--in-blocks=genomes-renamed.fasta')
    p:add('Pangenome')
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
    p:set_action(function(p)
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
    p:set_action(function(p)
        local fc = VectorFc()
        fc:add_bs(p:get_bs('pangenome'))
        fc:prepare()
        -- key: block_pos1_pos2
        -- value: array of gene fragments
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
    p:set_action(function(p)
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

register_p('UniqueFragments', function()
    local p = BlocksJobs.new()
    p:set_name("Remove fragments from blocks so that all " ..
               "fragments in a block are of different " ..
               "sequence")
    p:declare_bs('target', 'Target blockset')
    p:set_process_block(function(block)
        local fragments = block:fragments()
        local seen_seqs = {}
        for _, f in pairs(fragments) do
            if not seen_seqs[f:str()] then
                seen_seqs[f:str()] = 1
            else
                block:erase(f)
            end
        end
    end)
    return p
end)

register_p('ReadMutations', function()
    local p = LuaProcessor.new()
    p:set_name("Read table file with mutations (output of " ..
               "PrintMutations), takes consensuses from " ..
               "and constructs blockset in target blockset")
    p:declare_bs('target', 'Target blockset')
    p:declare_bs('other', 'Consensus sequences')
    p:add_opt('mutations', 'File with mutations ' ..
              '(output of PrintMutations)', '', true)
    p:set_action(function(p)
        local name2cons = {}
        local mut = {} --mut[block_name][fr_name][pos] = change
        for _, seq in pairs(p:other():seqs()) do
            name2cons[seq:name()] = seq
            local mut1 = {}
            mut[seq:name()] = mut1
            -- look for "fragments=f1,f2,f3..."
            local descr = seq:description()
            local fragments = descr:extract_value('fragments')
            for _, f in pairs(fragments:split(',')) do
                mut1[f] = {}
            end
        end
        --
        local fname = p:opt_value('mutations')
        local mut_file = file.name_to_istream(fname)
        local prev_b, prev_f
        while mut_file:good() do
            local line = mut_file:readline():trim()
            local b, f, pos, c = unpack(line:split())
            if b == '.' then
                b = prev_b
            end
            if f == '.' then
                f = prev_f
            end
            prev_b = b
            prev_f = f
            if b and b ~= 'block' then -- header
                local mut1 = assert(mut[b])
                local mut2 = assert(mut1[f])
                pos = tonumber(pos)
                if tonumber(c) then
                    -- long gaps
                    local stop = tonumber(c)
                    c = '-'
                    for i = pos, stop do
                        mut2[i] = c
                    end
                else
                    mut2[pos] = c
                end
            end
        end
        -- write fragments to temp fasta file
        local tmp_fname = ':read_mut_' .. rand_name(10)
        file.set_sstream(tmp_fname, '')
        local fasta = file.name_to_ostream(tmp_fname)
        for b, v in pairs(mut) do
            local cons = assert(name2cons[b])
            local length = cons:size()
            for f, m in pairs(v) do
                local chars = {}
                for pos = 0, length - 1 do
                    local c = m[pos] or cons:char_at(pos)
                    table.insert(chars, c)
                end
                local str = table.concat(chars)
                fasta:write_fasta(f, 'block=' .. b, str)
            end
        end
        fasta:flush()
        -- read fasta file
        local in_p = new_p('In')
        in_p:set_parent(p)
        in_p:fix_opt_value('in-blocks', tmp_fname)
        in_p:set_bs('target', p:block_set())
        in_p:run()
        -- remove fasta file
        file.remove_stream(tmp_fname)
    end)
    return p
end)

register_p('DownloadGenomesTables', function()
    local p = LuaProcessor.new()
    p:set_name("Download genomes tables from EBI server " ..
               "and writes files like bacteria/Escherichia.tsv")
    p:add_opt('no-drafts', 'Skip draft genomes', true)
    p:add_opt('min-genomes',
              'Minimum number of genomes in genus', 10)
    p:add_opt('kingdoms', 'Kingdoms downloaded',
            {'archaea', 'archaealvirus', 'bacteria', 'virus',
            'eukaryota', 'organelle', 'phage', 'plasmid'})
    p:set_action(function(p)
        -- download genomes from EBI
        local no_drafts = p:opt_value('no-drafts')
        local genomes = {}
        for _, k in ipairs(p:opt_value('kingdoms')) do
            genomes[k] = {}
            local fname = k .. '.details.txt'
            local base_url = 'http://www.ebi.ac.uk/genomes/'
            local url = base_url .. fname
            file.download_file(url, fname)
            local input = file.name_to_istream(fname)
            while input:good() do
                local line = input:readline():trim()
                local draft = line:lower():find('draft')
                local header = line:starts_with('#')
                if not header and
                        not (draft and no_drafts) then
                    local id, ver, date, tax, descr =
                        unpack(line:split('\t'))
                    if descr then
                        local genus = descr:split()[1]
                        if genomes[k][genus] == nil then
                            genomes[k][genus] = {}
                        end
                        table.insert(genomes[k][genus],
                            {id, descr, tax})
                    end
                end
            end
        end
        -- download taxons from Uniprot
        local taxonomy = 'taxonomy-all.tab'
        if not file_exists(taxonomy) then
            local url = 'http://www.uniprot.org/' ..
                'taxonomy/?format=tab&force=yes'
            file.download_file(url, taxonomy)
        end
        -- read taxons
        local taxons = {}
        local taxon_parent = {}
        local input = file.name_to_istream(taxonomy)
        while input:good() do
            local line = input:readline()
            local fields = line:split('\t')
            local tax, name = unpack(fields)
            local parent = fields[#fields - 1]
            if tax then
                taxons[tax] = name
                taxon_parent[tax] = parent
            end
        end
        -- write .tsv files in subdirs
        function guess_chromosome(descr)
            local d = descr:lower()
            d = d:gsub('chromosome: ', 'chromosome ')
            for i = 1, 50 do
                if d:ends_with('chromosome ' .. i) then
                    return 'chr' .. i
                end
            end
            if d:ends_with('chromosome i') then
                return 'chr1'
            end
            if d:ends_with('chromosome ii') then
                return 'chr2'
            end
            if d:ends_with('chromosome iii') then
                return 'chr3'
            end
            if d:find('plasmid') then
                return 'plasmid'
            end
            return 'chr'
        end
        local circularity = {
            archaea = 'c', archaealvirus = 'l',
            bacteria = 'c', virus = 'l',
            eukaryota = 'l', organelle = 'c',
            phage = 'c', plasmid = 'c'
        }
        function guess_circular(descr, k)
            local d = descr:lower()
            if d:find('circular') then
                return 'c'
            end
            if d:find('linear') then
                return 'l'
            end
            return circularity[k]
        end
        function taxon_of(tax)
            if taxons[tax] and taxons[tax] ~= '' then
                return taxons[tax]
            end
            local parent = taxon_parent[tax]
            if parent and taxons[parent] and
                    taxons[parent] ~= '' then
                return taxons[parent]
            end
            return tax
        end
        function write_genus(k, genus, g2)
            collectgarbage() -- close opened files
            local names = {}
            genus = genus:gsub(' ', '_'):gsub('/', '_')
            local tsv = file.cat_paths(k, genus .. '.tsv')
            local output = file.name_to_ostream(tsv)
            for _, genome in ipairs(g2) do
                local id, descr, tax = unpack(genome)
                -- CP000380.1 -> CP000380
                id = id:split('%.')[1]
                local r = {}
                local mnem = taxon_of(tax)
                local chr = guess_chromosome(descr)
                local name = mnem .. '&' .. chr
                if names[name] ~= nil then
                    table.insert(r, '#')
                end
                names[name] = 1
                table.insert(r, 'all:embl:' .. id)
                table.insert(r, mnem)
                table.insert(r, chr)
                table.insert(r, guess_circular(descr, k))
                table.insert(r, descr)
                output:write(table.concat(r, ' ') .. '\n')
            end
        end
        function count_genomes(g2)
            local cc = 0
            local set = {}
            for _, genome in ipairs(g2) do
                local id, descr, tax = unpack(genome)
                local mnem = taxon_of(tax)
                if set[mnem] == nil then
                    cc = cc + 1
                end
                set[mnem] = 1
            end
            return cc
        end
        local min_genomes = p:opt_value('min-genomes')
        for k, g1 in pairs(genomes) do
            file.make_dir(k)
            for genus, g2 in pairs(g1) do
                if count_genomes(g2) >= min_genomes then
                    write_genus(k, genus, g2)
                end
            end
        end
    end)
    return p
end)

function genome_seqs(bs, genome)
    local seqs = {}
    for _, seq in ipairs(bs:seqs()) do
        if seq:genome() == genome then
            table.insert(seqs, seq)
        end
    end
    return seqs
end

register_p('DraftPangenome', function()
    local p = LuaProcessor.new()
    p:set_name("Build draft pangenome on random subset " ..
               "of genomes (stem blocks, one iteration)")
    p:add_opt('ngenomes', 'Number of genomes (-1 is all)', -1)
    p:declare_bs('other', 'Input genomes')
    p:declare_bs('target', 'Where draft is written')
    p:set_action(function(p)
        -- select 4 genomes
        local ngenomes = p:opt_value('ngenomes')
        local other = p:other()
        local bs = p:block_set()
        local genomes = other:genomes_list()
        if ngenomes == -1 then
            ngenomes = #genomes
        end
        ngenomes = math.min(ngenomes, #genomes)
        math.randomseed(os.time())
        for i = 1, ngenomes do
            local r = math.random(1, #genomes)
            local genome = table.remove(genomes, r)
            bs:add_sequences(genome_seqs(other, genome))
        end
        --
        local finder = new_p('AnchorFinder')
        finder:set_parent(p)
        finder:apply(bs)
        Processor.delete(finder)
        --
        local stem = new_p('Stem')
        stem:set_parent(p)
        stem:set_opt_value('exact', true)
        stem:apply(bs)
        Processor.delete(stem)
        --
        local aligner = new_p('DummyAligner')
        aligner:set_parent(p)
        aligner:apply(bs)
        Processor.delete(aligner)
        --
        local extender = new_p('ExtendLoopFast')
        extender:set_parent(p)
        extender:set_max_loops(10)
        extender:apply(bs)
        Processor.delete(extender)
        --
        local filter = new_p('Filter')
        filter:set_parent(p)
        filter:apply(bs)
        Processor.delete(filter)
    end)
    return p
end)

register_p('RecommendIdentity', function()
    local p = LuaProcessor.new()
    p:set_name("Recommend value of MIN_IDENTITY of " ..
               "pangenome to be built")
    p:add_opt('recommendation',
        'Output file of recommended MIN_IDENTITY ' ..
        'of pangenome to be built', '', true)
    p:declare_bs('target', 'Draft pangenome')
    p:set_action(function(p)
        -- find average identity
        local bs = p:block_set()
        local sum_identity = 0
        local sum_weights = 0
        for _, block in pairs(bs:blocks()) do
            local weight = block:alignment_length()
            local identity = block:identity():to_d()
            sum_identity = sum_identity + weight * identity
            sum_weights = sum_weights + weight
        end
        local avg_identity = sum_identity / sum_weights
        local recommended = avg_identity - 0.1
        local fname = p:opt_value('recommendation')
        local out = file.name_to_ostream(fname)
        out:write(string.format(
[[Estimation of average identity in core blocks: %0.3f
Recommended value of MIN_IDENTITY: %0.3f
]], avg_identity, recommended))
        out:flush()
    end)
    return p
end)

register_p('DraftAndRecommend', function()
    local p = Pipe.new()
    p:add('DraftPangenome', 'target=draft other=target')
    p:add('RecommendIdentity', 'target=draft ' ..
        '--recommendation=examine/identity_recommended.txt')
    p:add('OutputPipe', 'target=draft ' ..
        '--out-file=examine/draft.bs --skip-rest:=1')
    return p
end)

register_p('MakeDraftPangenome', function()
    local p = Pipe.new()
    p:add('In', [[target=other
        --in-blocks=genomes-renamed.fasta]])
    p:add('DraftAndRecommend')
    return p
end)

register_p('FindBlock', function()
    local p = LuaProcessor.new()
    p:set_name('Find block by block name')
    p:declare_bs('other', 'Input blocks')
    p:declare_bs('target', 'Where found block is copied')
    p:add_opt('block-name', 'Block name pattern ' ..
        '(See Lua String Patterns)', 'u1x1')
    p:set_action(function(p)
        local pattern = p:opt_value('block-name')
        for _, block in ipairs(p:other():blocks()) do
            if block:name():match(pattern) then
                p:block_set():insert(block:clone())
            end
        end
    end)
    return p
end)

register_p('MakeFindBlock', function()
    local p = Pipe.new()
    p:add('In', 'target=other --in-blocks=pangenome.bs')
    p:add('FindBlock')
    p:add('Output')
    return p
end)

register_p('GenomeLengths', function()
    local p = LuaProcessor.new()
    p:set_name('Print lengths of all genomes')
    p:declare_bs('target', 'Target blockset')
    p:add_opt('genomes-info', 'Output file',
              'examine/genomes-info.tsv')
    p:set_action(function(p)
        local fname = p:opt_value('genomes-info')
        local out = file.name_to_ostream(fname)
        local bs = p:block_set()
        for _, genome in ipairs(bs:genomes_list()) do
            local length = 0
            for _, seq in ipairs(genome_seqs(bs, genome)) do
                length = length + seq:size()
            end
            out:write(genome .. '\t' .. length .. '\n')
        end
    end)
    return p
end)

register_p('AllProcessors', function()
    local p = LuaProcessor.new()
    p:set_name('Print table of all processors')
    p:add_opt('out', 'Output file', ':stdout')
    p:set_action(function(p)
        local fname = p:opt_value('out')
        local out = file.name_to_ostream(fname)
        out:write("<table border='1'>")
        out:write("<tr><td>Summary</td><td>Help</td></tr>")
        local key2pr = {}
        local deleters = {}
        for _, key in ipairs(meta:keys()) do
            local pr = new_p(key)
            table.insert(deleters, Processor.deleter(pr))
            if pr:name() ~= "" and pr:name() ~= key then
                key2pr[key] = pr
            end
        end
        local tr = [[
        <tr valign="top">
            <td>
                <a name="%s"></a>
                <b>%s</b>
                <br/>
                %s

                %s
            </td>
            <td>
                %s

                %s
                &nbsp;
            </td>
        </tr>
        ]]
        local get_bss = function(pr)
            local bss = pr:get_block_sets()
            local out = {}
            local ok = false
            table.insert(out, "<ul>")
            for _, bs_name in ipairs(bss) do
                local descr = pr:bs_description(bs_name)
                if descr ~= "" then
                    ok = true
                    table.insert(out, "<li>")
                    table.insert(out, "<u>")
                    table.insert(out, bs_name)
                    table.insert(out, "</u>")
                    table.insert(out, ": ")
                    table.insert(out, descr)
                    table.insert(out, "</li>")
                end
            end
            table.insert(out, "</ul>")
            if not ok then
                return ''
            end
            local header = "<i>Blocksets</i>:"
            return header .. table.concat(out)
        end
        local get_opts = function(pr)
            local opts = pr:opts()
            if #opts == 2 then
                -- 2 because 'timing', 'workers'
                return ''
            end
            local out = {}
            table.insert(out, "<table>")
            for _, opt in ipairs(opts) do
                if opt ~= 'timing' and opt ~= 'workers' then
                    local full_opt = '--' ..
                        pr:opt_prefixed(opt)
                    local descr = pr:opt_description(opt)
                    table.insert(out, "<tr valign='top'>")
                    table.insert(out, "<td width='20px'></td>")
                    table.insert(out, "<td><nobr>")
                    table.insert(out, full_opt)
                    table.insert(out, "</nobr></td>")
                    table.insert(out, "<td>")
                    table.insert(out, descr)
                    table.insert(out, "</td>")
                    table.insert(out, "</tr>")
                end
            end
            table.insert(out, "</table>")
            local header = "<i>Options</i>:"
            return header .. table.concat(out)
        end
        local get_childs = function(pr)
            local childs = pr:children()
            if #childs == 0 then
                return ''
            end
            local childs_set = {}
            for _, child in ipairs(childs) do
                childs_set[child:key()] = child
            end
            local out = {}
            for ch, child in pairs(childs_set) do
                if key2pr[ch] then
                    local t = '<a href="#%s">%s</a>'
                    table.insert(out, t:format(ch, ch))
                else
                    table.insert(out, ch)
                end
            end
            local header =
                "<br/><br/><i>Child processors</i>:<br/>"
            return header .. table.concat(out, ', ')
        end
        for key, pr in pairs(key2pr) do
            local name = pr:name()
            local bss = get_bss(pr)
            local opts = get_opts(pr)
            local ch = get_childs(pr)
            out:write(tr:format(key, key, name, ch, bss, opts))
        end
        out:write("</table>")
        out:flush()
    end)
    return p
end)

