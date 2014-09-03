s = Sequence.new()
s:set_name('abc')
s:push_back("ATTC")
f1 = Fragment.new(s, 1,1)
f2 = Fragment.new(s, 2,2)
b = Block.new()
b:insert(f1)
b:insert(f2)

for f in b:iter_fragments() do
    assert(f:id() == 'abc_1_1' or f:id() == 'abc_2_2')
end

for _, f in next, b:fragments() do
    assert(f:id() == 'abc_1_1' or f:id() == 'abc_2_2')
end

assert(b:identity() == Decimal(1))
Block.delete(b)

