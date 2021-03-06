one = Decimal(1)
assert(one:to_s() == '1')
assert(Decimal(one) == one)
assert(Decimal(1) == Decimal('1'))
assert(Decimal(1) ~= Decimal('1.7'))
assert(Decimal(10) * Decimal('0.1') == Decimal(1))
assert(Decimal('20.1') > Decimal('20'))
assert(Decimal('20.1') >= Decimal('20'))
assert(Decimal('-20.1') < Decimal('20'))
assert(Decimal('-20.1') < Decimal('-20'))
assert(Decimal('-20.0') == Decimal('-20'))
assert(Decimal('0') == Decimal(0))
assert(Decimal('100') * Decimal(0) == Decimal(0))
assert(Decimal('1.2'):round() == 1)
assert(Decimal('1.2'):to_d() >= 1.199)
assert(Decimal('1.2'):to_d() <= 1.201)
assert(-Decimal('1.2') == Decimal('-1.2'))
assert(-Decimal() == Decimal())
assert(Decimal(0) == Decimal('0'))
assert(Decimal(0) == Decimal())
assert(Decimal('1.2') + Decimal('3.4') == Decimal('4.6'))
assert(Decimal('1.2') - Decimal('3.4') == Decimal('-2.2'))
assert(Decimal('1.2') * Decimal('3.4') == Decimal('4.08'))
assert(Decimal('1.7') / Decimal('3.4') == Decimal('0.5'))
assert(Decimal('1.7'):to_s() == '1.7')

