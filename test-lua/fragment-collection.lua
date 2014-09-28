local Fcs = {
    VectorFc, SetFc
}

for _, Fc in pairs(Fcs) do
    local s = Sequence.new('ATGCATGC')
    local f1 = Fragment.new(s, 1, 2)
    local f1_d = Fragment.deleter(f1)
    local f2 = Fragment.new(s, 2, 3, -1)
    local f2_d = Fragment.deleter(f2)
    local f3 = Fragment.new(s, 3, 3)
    local f3_d = Fragment.deleter(f3)
    local f4 = Fragment.new(s, 0, 3)
    local f4_d = Fragment.deleter(f4)

    local fc = Fc()
    fc:add_fragment(f1)
    fc:prepare()
    assert(fc:has_overlap(f1))
    assert(fc:has_overlap(f2))
    assert(not fc:has_overlap(f3))
    assert(fc:has_overlap(f4))

    -- circular, 1 fragment
    s:set_name('A&cI&c')
    assert(fc:next(f1) == f1)
    assert(fc:prev(f1) == f1)
    -- linear, 1 fragment
    s:set_name('A&cI&l')
    assert(fc:next(f1) == nil)
    assert(fc:prev(f1) == nil)

    -- 2 fragments
    fc:add_fragment(f2)
    fc:prepare()
    -- circular, 2 fragments
    s:set_name('A&cI&c')
    assert(fc:next(f1) == f2)
    assert(fc:prev(f1) == f2)
    assert(fc:next(f2) == f1)
    assert(fc:prev(f2) == f1)
    -- linear, 2 fragments
    s:set_name('A&cI&l')
    assert(fc:next(f1) == f2)
    assert(fc:prev(f1) == nil)
    assert(fc:next(f2) == nil)
    assert(fc:prev(f2) == f1)

    -- 3 fragments
    fc:add_fragment(f3)
    fc:prepare()
    -- linear, 3 fragments
    s:set_name('A&cI&l')
    assert(fc:next(f2) == f3)
    assert(fc:prev(f2) == f1)
    assert(fc:next(f1) == f2)
    assert(fc:prev(f3) == f2)
    assert(fc:neighbor(f2, 1) == f3)
    assert(fc:neighbor(f2, -1) == f1)
    assert(fc:logical_neighbor(f2, 1) == f1)
    assert(fc:logical_neighbor(f2, -1) == f3)
    assert(fc:are_neighbors(f1, f3) == 0)
    assert(fc:are_neighbors(f3, f1) == 0)
    assert(fc:are_neighbors(f1, f2) == 1)
    assert(fc:are_neighbors(f2, f1) == -1)
    assert(fc:are_neighbors(f2, f3) == 1)
    assert(fc:are_neighbors(f3, f2) == -1)
    assert(fc:another_neighbor(f1, f2) == nil)
    assert(fc:another_neighbor(f1, f3) == nil)
    assert(fc:another_neighbor(f2, f3) == f1)
    assert(fc:another_neighbor(f2, f1) == f3)
    -- TODO
    -- local seqs = fc:seqs()
    -- assert(#seqs == 1)
    -- assert(#seqs[1] == s)
end

