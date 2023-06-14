from myutils.jupyter import hide_code
from IPython.display import HTML


def test_hide_code():
    hc = hide_code()
    assert isinstance(hc, HTML)
