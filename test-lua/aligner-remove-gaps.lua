
local s = Sequence.new("ATACT")
local f1 = Fragment.new(s, 0, 1)
f1:set_row(AlignmentRow.new("A--T"))
local f2 = Fragment.new(s, 2, 4)
f2:set_row(AlignmentRow.new("A-CT"))
local block = Block.new()
block:insert(f1)
block:insert(f2)

remove_pure_gap_columns(block)

assert(f1:contents() == 'A-T')
assert(f2:contents() == 'ACT')

Block.delete(block)

