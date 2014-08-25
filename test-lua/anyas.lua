assert(type(get('MIN_LENGTH')) == 'number')
assert(type(get('MIN_IDENTITY')) == 'userdata')
d = get('MIN_IDENTITY')
d2 = Decimal(d)
assert(d2 == d)
set('MIN_IDENTITY', d2)

set('test', {'1', '2'}, "Test descr")
v = get('test')
assert(#v == 2)
assert(v[1] == '1')
assert(v[2] == '2')
assert(meta:get_description('test') == "Test descr")

set('test2', true)
v = get('test2')
assert(v == true)

