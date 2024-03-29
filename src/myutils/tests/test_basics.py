from myutils.tests.variables4tests import basics_bash
from myutils.miscellaneous import output_terminal
from pytest import raises
from os.path import exists


def test_adjust():
    outcome = output_terminal(f"{basics_bash} "
                              "adjust simple trial")
    outcome = outcome.split("\n")

    assert outcome[1][:21] == '++++++++ test-basics:'
    assert 'simple trial' in outcome[1]
    assert len(outcome[1]) == 79
    assert len(outcome[-1]) == 0


def test_verbose():
    outcome = output_terminal(f"{basics_bash} "
                              "verbose simple trial")
    outcome = outcome.split("\n")
    assert outcome[1][:29] == '++++++++ test-basics: VERBOSE'
    assert 'simple trial' in outcome[1]
    assert len(outcome[0]) == 79
    assert len(outcome[-1]) == 0


def test_warning():
    outcome = output_terminal(f"{basics_bash} "
                              "warning simple trial")
    outcome = outcome.split("\n")
    assert outcome[1][:29] == '++++++++ test-basics: WARNING'
    assert 'simple trial' in outcome[1]
    assert len(outcome[0]) == 79
    assert len(outcome[-1]) == 0


def test_finish():
    outcome = output_terminal(f"{basics_bash} "
                              "finish simple trial")
    outcome = outcome.split("\n")
    assert outcome[1][:21] == '++++++++ test-basics:'
    assert 'simple trial' in outcome[1]
    assert len(outcome[0]) == 79
    assert len(outcome[-1]) == 0


def test_fail():
    with raises(Exception) as e_info:
        output_terminal(f"{basics_bash} "
                        "fail simple trial")
    outcome = str(e_info.value).split("\n")
    assert "ERROR executing the command" in outcome[0]
    assert outcome[1][:27] == '++++++++ test-basics: ERROR'
    assert 'simple trial' in outcome[1]
    assert len(outcome[1]) == 79


def test_create_bck():
    # file
    output_terminal("touch remove-test.log remove-test.com; "
                    f"{basics_bash} create_bck remove-test.log "
                    "remove-test.com")
    output_terminal("touch remove-test.com; "
                    f"{basics_bash} create_bck remove-test.com")
    assert exists('remove-test-bck_1.com')
    assert exists('remove-test-bck_2.com')
    assert exists('remove-test-bck_1.log')
    assert not exists('remove-test.com')
    assert not exists('remove-test.log')
