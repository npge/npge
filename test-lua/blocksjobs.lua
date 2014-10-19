
register_p('RemoveOverlapsWithOther', function()
    local p = BlocksJobs.new()

    p:set_before_work(function()
        local work_data = {}
        work_data.other_fc = VectorFc()
        work_data.other_fc:add_bs(p:other())
        work_data.other_fc:prepare()
        return work_data
    end)

    p:set_initialize_thread(function(work_data)
        assert(work_data)
        assert(work_data.processor)
        assert(work_data.processor:key() ==
            'RemoveOverlapsWithOther')
        assert(work_data.other_fc)
        local thread_data = {}
        thread_data.bad_blocks = {}
        return thread_data
    end)

    p:set_process_block(function(block, thread_data, work_data)
        assert(work_data.other_fc)
        assert(thread_data.bad_blocks)
        if work_data.other_fc:block_has_overlap(block) then
            table.insert(thread_data.bad_blocks, block)
        end
    end)

    p:set_after_thread(function(thread_data, work_data)
        assert(work_data.other_fc)
        assert(thread_data.bad_blocks)
        local bs = work_data.processor:block_set()
        for _, block in ipairs(thread_data.bad_blocks) do
            bs:erase(block)
        end
    end)

    return p
end)

local SEQ_LENGTH = 1000000

local seq = DummySequence.new('A', SEQ_LENGTH)

local test_p = new_p('RemoveOverlapsWithOther')

math.randomseed(os.time())

local fillRandom = function(bs, number, length)
    for i = 1, number do
        local start = math.random(0, SEQ_LENGTH - length * 2)
        local stop = start + length
        local f = Fragment.new(seq, start, stop)
        local block = Block.new()
        block:insert(f)
        bs:insert(block)
    end
end

fillRandom(test_p:other(), 20, 10)
fillRandom(test_p:block_set(), 1000, 1000)

test_p:set_workers(8)
test_p:run()

local other_fc = SetFc()
other_fc:add_bs(test_p:other())
other_fc:prepare()

assert(not other_fc:bs_has_overlap(test_p:block_set()))

Processor.delete(test_p)

