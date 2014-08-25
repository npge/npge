assert(proportion(100, 200, 2000) == 1000)

assert(size_to_char(char_to_size('A')) == 'A')
assert(size_to_char(char_to_size('T')) == 'T')
assert(size_to_char(char_to_size('G')) == 'G')
assert(size_to_char(char_to_size('C')) == 'C')

assert(complement_char('A') == 'T')
assert(complement_char('T') == 'A')
assert(complement_char('G') == 'C')
assert(complement_char('C') == 'G')

assert(complement_str('ATGC') == 'GCAT')
assert(complement_str('A') == 'T')
assert(complement_str('AN-G') == 'C-NT')

local ATGC_hash = make_hash('ATGC')
local ATGC_rev_hash = make_hash('ATGC', -1)
local ATGC_hash_rev = complement_hash(ATGC_hash, 4)
local GCAT_hash = make_hash('GCAT')
assert(GCAT_hash == ATGC_rev_hash)
assert(GCAT_hash == ATGC_hash_rev)

local TGCC_hash = reuse_hash(ATGC_hash, 4, 'A', 'C', true)
local GGCA_hash = reuse_hash(GCAT_hash, 4, 'T', 'G', false)
assert(complement_hash(TGCC_hash, 4) == GGCA_hash)

-- long sequences

local seq = 'CATAGTTAGTCCAAGCAGAGGTCGGCTTAGGCGCCTAATGTAGTT'
local seq_rev = complement_str(seq)
assert(complement_str(seq_rev) == seq)

local seq_hash = make_hash(seq)
local seq_rev_hash = make_hash(seq, -1)
local seq_hash_rev = complement_hash(seq_hash, #seq)
local qes_hash = make_hash(seq_rev)
assert(qes_hash == seq_rev_hash)
assert(qes_hash == seq_hash_rev)

local remove = 'C'
local add = 'T'

local eqT_hash = reuse_hash(seq_hash, #seq, remove, add, true)
local Aqe_hash = reuse_hash(qes_hash, #seq,
    complement_char(remove),
    complement_char(add),
    false)
assert(complement_hash(Aqe_hash, #seq) == eqT_hash)

-- random numbers

assert(type(rand_name(5)) == 'string')
assert(#rand_name(5) == 5)
assert(type(make_seed()) == 'number')
