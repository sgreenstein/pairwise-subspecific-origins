"""
File: subspecies.py
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

DOM = 0b001
MUS = 0b010
CAS = 0b100
NONE = 0b1000
NUM_SUBSPECIES = 3

_SHIFT = 4
_PROXIMAL_MASK = 0b111 << _SHIFT
_DISTAL_MASK = 0b111


def combine(proximal_origin, distal_origin):
    """ Converts two int species to an int combination
    :param proximal_origin: int representing subspecific origin of proximal locus
    :param distal_origin: int representing subspecific origin of distal locus
    :return: int representing the combination
    """
    return (proximal_origin << _SHIFT) | distal_origin

_subspecies_ints = [DOM, MUS, CAS, NONE]
_subspecies_names = ['dom', 'mus', 'cas', '???']
_int_to_str = {i: s for i, s in zip(_subspecies_ints, _subspecies_names)}
for prox_int, prox_name in zip(_subspecies_ints, _subspecies_names):
    for dist_int, dist_name in zip(_subspecies_ints, _subspecies_names):
        _int_to_str[combine(prox_int, dist_int)] = prox_name + ' :: ' + dist_name
_str_to_int = {v: k for k, v in _int_to_str.iteritems()}
_combos = []
_combos_with_none = []


def iter_subspecies(include_unknown=False):
    """ an iterator for all subspecies (does not include unknown)
    """
    return _subspecies_ints[:len(_subspecies_ints) - (not include_unknown)]

for proximal in iter_subspecies(False):
    for distal in iter_subspecies(False):
        _combos.append(combine(proximal, distal))

for proximal in iter_subspecies(True):
    for distal in iter_subspecies(True):
        _combos_with_none.append(combine(proximal, distal))


def iter_combos(include_unknown=False):
    """ an iterator for all combinations of subspecies (does not include unknown)
    """
    if include_unknown:
        return _combos_with_none
    else:
        return _combos


def to_string(integer):
    """
    :param integer: int representation of subspecies
    :return: subspecies name
    >>> to_string(DOM) == 'dom'
    True
    """
    return _int_to_str[integer]


def to_int(string):
    """
    :param string: subspecies name
    :return: int representation of subspecies
    >>> to_int('dom') == DOM
    True
    """
    return _subspecies_ints[_subspecies_names.index(string.lower())]


def is_distal(subspecies, combo):
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


def is_proximal(subspecies, combo):
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


def is_homo(combo):
    """
    :param combo: integer representation of a combination
    :return: True if both origins are the same
    >>> is_homo(combine(DOM, MUS))
    False
    >>> is_homo(combine(MUS, MUS))
    True
    """
    return ((combo >> _SHIFT) & _DISTAL_MASK) == (combo & _DISTAL_MASK)


def is_known(combo):
    """ Checks to see if either origin of combo is unknown
    :param combo: int representation of origin combination
    :return: True if neither origin is unknown
    >>> is_known(combine(DOM, CAS))
    True
    >>> is_known(combine(NONE, CAS))
    False
    """
    return not combine(NONE, NONE) & combo


def proximal(combo):
    """ Returns the proximal subspecies from a combo
    :param combo: int representation of origin combination
    :return: int representation of the proximal origin
    >>> proximal(combine(CAS, DOM)) == CAS
    True
    """
    return (combo & _PROXIMAL_MASK) >> _SHIFT


def distal(combo):
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
    print to_string(combine(DOM, NONE))
