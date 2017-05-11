#!/usr/bin/env python

import csdrpy

csdrpy_syms = filter(lambda x: x.startswith('_') == False, dir(csdrpy))
csdrpy_syms = set(csdrpy_syms)

assert( csdrpy.is_nan(1.0) == False)
assert( csdrpy.is_nan(float('NaN')) == True)
assert( csdrpy.log2n(1) == 0 )
assert( csdrpy.log2n(256) == 8)

print 'Loaded {0} symbols'.format(len(csdrpy_syms))
print 'csdr python binding works! :)'
