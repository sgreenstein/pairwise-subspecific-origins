"""
File: subspecies.pyx
Author: Seth Greenstein
Purpose: Provide methods for working with subspecies and combinations of subspecies

>>> for prox in _subspecies_ints:
...    for dist in _subspecies_ints:
...        if not is_distal(dist, combine(prox, dist)):
...            print 'Failed'
...        if not is_proximal(prox, combine(prox, dist)):
...            print 'Failed'
...        for test in _subspecies_ints:
...            if test != prox:
...                if is_proximal(test, combine(prox, dist)):
...                    print 'Failed'
...            if test != dist:
...                if is_distal(test, combine(prox, dist)):
...                    print 'Failed'
"""

cdef int DOM = 0b001
cdef int MUS = 0b010
cdef int CAS = 0b100
cdef int UNK = 0b1000
NUM_SUBSPECIES = 3

UNKNOWN = UNK

cdef int _SHIFT = 4
cdef int _PROXIMAL_MASK = 0b1111 << _SHIFT
cdef int _DISTAL_MASK = 0b1111


cpdef int combine(int proximal_origin, int distal_origin):
    """ Converts two int species to an int combination
    :param proximal_origin: int representing subspecific origin of proximal locus
    :param distal_origin: int representing subspecific origin of distal locus
    :return: int representing the combination
    """
    return (proximal_origin << _SHIFT) | distal_origin

_subspecies_ints = [DOM, MUS, CAS, UNK]
_subspecies_names = ['dom', 'mus', 'cas', 'unk']
_int_to_str = {i: s for i, s in zip(_subspecies_ints, _subspecies_names)}
for prox_int, prox_name in zip(_subspecies_ints, _subspecies_names):
    for dist_int, dist_name in zip(_subspecies_ints, _subspecies_names):
        _int_to_str[combine(prox_int, dist_int)] = (prox_name, dist_name)
_str_to_int = {v: k for k, v in _int_to_str.iteritems()}
_combos = []


def iter_subspecies(include_unknown=False):
    """ an iterator for all subspecies
    :param include_unknown: Include combos with UNK as one of the sources, default False
    """
    return _subspecies_ints[:len(_subspecies_ints) - (not include_unknown)]

for prox in iter_subspecies(False):
    for dist in iter_subspecies(False):
        _combos.append(combine(prox, dist))
for known in iter_subspecies(False):
    _combos.append(combine(known, UNK))
    _combos.append(combine(UNK, known))
_combos.append(combine(UNK, UNK))


def iter_combos(include_unknown=False):
    """ an iterator for all combinations of subspecies
    :param include_unknown: Include combos with UNK as one of the sources, default False
    """
    if include_unknown:
        return _combos
    else:
        return _combos[:-7]

combo_nums = {}
combo_names = []
for n, c in enumerate(iter_combos(True)):
    combo_nums[c] = n
    combo_names.append(_int_to_str[c])

_COLORS = {DOM: 0x0000ff, MUS: 0xff0000, CAS: 0x00ff00}
_ORDINAL_COLORS = []

for combo in iter_combos():
    prox_color = _COLORS[proximal(combo)]
    dist_color = _COLORS[distal(combo)]
    color = ((prox_color >> 1) & prox_color) + ((dist_color >> 1) & dist_color)
    _ORDINAL_COLORS.append(color)
    _COLORS[combo] = color


cpdef int to_color(integer, ordinal=False):
    """
    :param integer: int representation of subspecies or subspecies combo
    :param ordinal: if the integer is the ordinal representation, default False
    :return: integer representing the rgb color
    >>> to_color(combine(DOM, MUS)) == to_color(to_ordinal(combine(DOM, MUS)), True)
    """
    if ordinal:
        return _ORDINAL_COLORS[integer]
    return _COLORS[integer]



cpdef int to_ordinal(int integer):
    """
    :param integer: int representation of subspecies combo
    :return: integer in the range 0...num combos/subspecies
    >>> to_ordinal(MUS) == iter_subspecies().index(MUS)
    True
    >>> to_ordinal(combine(DOM, CAS)) == iter_combos().index(combine(DOM, CAS))
    True
    """
    if proximal(integer):
        return combo_nums[integer]
    return iter_subspecies(True).index(integer)


def to_string(integer, ordinal=False):
    """
    :param integer: int representation of subspecies
    :param ordinal: if the integer is the ordinal representation, default False
    :return: subspecies name
    >>> to_string(DOM)
    'dom'
    >>> to_string(to_ordinal(combine(DOM, MUS)), True)
    ('dom', 'mus')
    """
    if ordinal:
        return combo_names[integer]
    return _int_to_str[integer]


cpdef int to_int(string):
    """
    :param string: subspecies name
    :return: int representation of subspecies
    >>> to_int('dom') == DOM
    True
    """
    # TODO: remove stuff about "bad" (residual heterozygosity)
    if len(string) == 2:
        return -999
    return _subspecies_ints[_subspecies_names.index(string.lower())]


cpdef int is_distal(int subspecies, int combo):
    """
    :param subspecies: int representation of single subspecies
    :param combo: int representing combination of two subspecies
    :return: true if subspecies is the distal of the combination
    >>> is_distal(CAS, combine(DOM, CAS))
    True
    >>> is_distal(DOM, combine(DOM, CAS))
    False
    """
    return (subspecies & _DISTAL_MASK) == (combo & _DISTAL_MASK)


cpdef int is_proximal(int subspecies, int combo):
    """
    :param subspecies: int representation of single subspecies
    :param combo: int representing combination of two subspecies
    :return: true if subspecies is the proximal of the combination
    >>> is_proximal(DOM, combine(DOM, CAS))
    True
    >>> is_proximal(CAS, combine(DOM, CAS))
    False
    """
    return ((subspecies << _SHIFT) & _PROXIMAL_MASK) == (combo & _PROXIMAL_MASK)


cpdef int is_homo(int combo):
    """
    :param combo: integer representation of a combination
    :return: True if both origins are the same
    >>> is_homo(combine(DOM, MUS))
    False
    >>> is_homo(combine(MUS, MUS))
    True
    """
    return ((combo >> _SHIFT) & _DISTAL_MASK) == (combo & _DISTAL_MASK)

cdef int _NONE_NONE = combine(UNK, UNK)


cpdef int is_known(int combo):
    """ Checks to see if either origin of combo is unknown
    :param combo: int representation of origin combination
    :return: True if neither origin is unknown
    >>> is_known(combine(DOM, CAS))
    True
    >>> is_known(combine(UNK, CAS))
    False
    """
    return not _NONE_NONE & combo


cpdef int proximal(int combo):
    """ Returns the proximal subspecies from a combo
    :param combo: int representation of origin combination
    :return: int representation of the proximal origin
    >>> proximal(combine(CAS, DOM)) == CAS
    True
    """
    return (combo & _PROXIMAL_MASK) >> _SHIFT


cpdef int distal(int combo):
    """ Returns the distal subspecies from a combo
    :param combo: int representation of origin combination
    :return: int representation of the distal origin
    >>> distal(combine(CAS, DOM)) == DOM
    True
    """
    return combo & _DISTAL_MASK


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    print to_string(DOM)
    print to_string(combine(DOM, UNK))
