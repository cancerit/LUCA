from luca.codec import Codec


def test_codec():
    c = Codec(start=10)
    assert c.add('A') == 10
    assert c.add('B') == 11
    assert c.add('B') == 11

    assert c.decode(10) == 'A'
    assert c.decode(11) == 'B'

    assert c.encode('A') == 10
    assert c.encode('B') == 11

    assert c.items == ['A', 'B']

    assert len(c) == 2
    assert 'A' in c
    assert 'B' in c
    assert 'C' not in c
