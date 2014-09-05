local Fcs = {
    VectorFc, SetFc
}

for _, Fc in pairs(Fcs) do
    local s = Sequence.new('ATGCATGC')
    local f1 = Fragment.new(s, 1, 2)
    local f1_d = Fragment.deleter(f1)
    local f2 = Fragment.new(s, 2, 3)
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
end

