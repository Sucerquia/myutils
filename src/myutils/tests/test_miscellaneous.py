from myutils.miscellaneous import (optimized_e,
                                   output_terminal,
                                   time_g09)


from myutils.tests.variables4tests import gpa_endo_log


def test_output_terminal():
    out = output_terminal("echo 'working'")
    h_msg = output_terminal("myutils -h").split('\n')
    assert out == 'working\n'
    assert isinstance(h_msg, list)
    assert isinstance(h_msg[0], str)
    assert len(h_msg) > 1


def test_time_g09():
    time = time_g09(gpa_endo_log)
    assert time == 103.75


def test_optimized_e():
    e = optimized_e(gpa_endo_log)
    assert e == -27965.4226678293
