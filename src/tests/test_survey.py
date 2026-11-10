#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2026  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later).
#
# The full copyright notice, including terms governing use, modification,
# and redistribution, is contained in the files LICENSE and COPYRIGHT,
# which can be found at the root of the source code distribution tree:
# - LICENSE <https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE>
# - COPYRIGHT <https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT>
#
"""
Tests of Survey, which turns an Input into the parameters and the commands that run Galaxia, and runs them. The backend is not
run: execute is replaced by a recording one that also makes the files that the backend would (the work of the backend is for
the integration tests), and the photometric systems are small ones that do not need the isochrones of the installation. What is
written goes under temporary directories, as the output directory of a survey is the current one unless it is told otherwise.
"""
import os
import re
import sys
import types

import numpy as np
import pytest

from ..galaxia_ananke import Input, Output, Survey, photometry
from ..galaxia_ananke._constants import GALAXIA
from ..galaxia_ananke._defaults import DEFAULTS_FOR_PARFILE, DEFAULT_CMD
from ..galaxia_ananke._templates import FTTAGS, PARFILE_TEMPLATE
from ..galaxia_ananke.photometry.Photometry import nested_dict
from ..galaxia_ananke.photometry.PhotoSystem import PhotoSystem
from ..galaxia_ananke.utils import hash_iterable
from .utils import minimal_kernels, minimal_particles, patch_nbody1

SURVEY_MODULE = sys.modules[Survey.__module__]


class FakeSystem:
    """A photometric system as Survey uses one, without the isochrones of the installation: it has the real check of the colors."""
    _required_cmd_magnames_dictkeys = PhotoSystem._required_cmd_magnames_dictkeys
    check_cmd_magnames = PhotoSystem.check_cmd_magnames

    def __init__(self, category, name, mag_names):
        self.category, self.name, self.mag_names = category, name, mag_names

    @property
    def categ_and_name(self):
        return [self.category, self.name]

    @property
    def key(self):
        return f"{self.category}/{self.name}"


GAIA = FakeSystem('padova', 'GAIA__DR2', ['G', 'Gbp', 'Grp'])
LSST = FakeSystem('padova', 'LSST', ['u', 'g', 'r', 'i', 'z', 'Y'])
MINE = FakeSystem('py_custom', 'Mine', ['A', 'B', 'C'])


@pytest.fixture(autouse=True)
def systems(monkeypatch):
    monkeypatch.setattr(photometry, 'available_photo_systems',
                        nested_dict({'padova': {'GAIA__DR2': GAIA, 'LSST': LSST}, 'py_custom': {'Mine': MINE}}))


@pytest.fixture(autouse=True)
def nbody1(tmp_path, monkeypatch):
    return patch_nbody1(monkeypatch, tmp_path)


@pytest.fixture
def directories(tmp_path):
    (tmp_path / 'in').mkdir()
    (tmp_path / 'out').mkdir()
    return tmp_path / 'in', tmp_path / 'out'


def make_input(directories, **kwargs):
    inp = Input(minimal_particles(), minimal_kernels(), **{'input_dir': directories[0], 'name': 'sim', 'ngb': 64, **kwargs})
    inp.input_sorter = None                           # as a survey does when it runs, before anything asks the input for its hash
    return inp


@pytest.fixture
def an_input(directories):
    return make_input(directories)


@pytest.fixture
def backend(monkeypatch, directories):
    """execute, recording what it is given and making the files that the backend makes: one for each generation that it runs."""
    calls = []

    def execute(cmds, max_workers=None, verbose=None, **kwargs):
        calls.append(types.SimpleNamespace(cmds=list(cmds), max_workers=max_workers, verbose=verbose))
        for cmd in cmds:
            name = re.search(r"/([^/ ]*)_params$", cmd)
            made = re.search(r"--nfile=(\S+) --ngen=(\d+)", cmd)
            if " -r " in cmd and name and made:
                (directories[1] / f"{name.group(1)}.{made.group(1)}.{made.group(2)}.ebf").write_text('')
    monkeypatch.setattr(SURVEY_MODULE, 'execute', execute)
    return calls


def prepare(survey, directories, cmd=DEFAULT_CMD, n_gens=(0,), **kwargs):
    survey._prepare_survey_parameters_and_output(cmd, n_gens, **{'output_dir': directories[1], 'fsample': 1, **kwargs})
    return survey


# -- how a survey is made -----------------------------------------------------------

def test_the_defaults_of_a_survey(an_input):
    survey = Survey(an_input)
    assert survey.surveyname == 'survey' and survey.verbose is True and survey.input is an_input
    assert survey.photosystems == [GAIA] and survey.photo_sys == {'padova/GAIA__DR2'}
    assert repr(survey) == "Survey(surveyname=survey, photo_sys={'padova/GAIA__DR2'})"


def test_what_is_given_is_what_the_survey_says(an_input):
    survey = Survey(an_input, photo_sys='py_custom/Mine', surveyname='mine', verbose=False)
    assert (survey.surveyname, survey.verbose, survey.photosystems) == ('mine', False, [MINE])


def test_photometric_systems_can_be_given_as_a_list_and_keep_their_order(an_input):
    survey = Survey(an_input, photo_sys=['padova/LSST', 'padova/GAIA__DR2', 'py_custom/Mine'])
    assert survey.photosystems == [LSST, GAIA, MINE]
    assert survey.photo_sys == {'padova/LSST', 'padova/GAIA__DR2', 'py_custom/Mine'}


def test_a_photometric_system_that_does_not_exist_is_a_key_error(an_input):
    with pytest.raises(KeyError, match='Nothing'):
        Survey(an_input, photo_sys='padova/Nothing')


def test_the_default_photometric_system_is_gaia_dr2_of_padova():
    from ..galaxia_ananke._defaults import DEFAULT_PSYS
    assert DEFAULT_PSYS == ['padova/GAIA__DR2']


def test_the_old_names_of_the_photometric_systems_still_work_and_say_they_are_old(an_input):
    survey = Survey(an_input)
    with pytest.warns(DeprecationWarning, match='photosystems'):
        assert survey.isochrones == [GAIA]
    with pytest.warns(DeprecationWarning, match='prepare_photosystems'):
        assert Survey.set_isochrones_from_photosys('padova/LSST') == [LSST]


def test_the_verbosity_can_be_changed_but_not_unset(an_input):
    survey = Survey(an_input, verbose=False)
    survey.verbose = True
    assert survey.verbose is True
    survey.verbose = None
    assert survey.verbose is True


def test_a_survey_says_what_it_takes_from_its_input(an_input, directories):
    survey = Survey(an_input)
    assert (survey.hdim, survey.ngb, survey.caching, survey.append_hash) == (6, 64, False, True)
    assert survey.inputdir == directories[0] and survey.inputname_hash == an_input.name_hash


# -- before it has been made --------------------------------------------------------

@pytest.mark.parametrize("attribute", ["fileparam", "parameters", "hash", "surveyname_hash", "output"])
def test_a_survey_that_was_not_made_says_to_make_it_first(attribute, an_input):
    with pytest.raises(RuntimeError, match="hasn't been made yet, run method `make_survey` first"):
        getattr(Survey(an_input), attribute)


@pytest.mark.parametrize("attribute", ["n_gens", "fsample"])
def test_a_survey_that_was_not_made_says_to_make_it_first_for_its_generations_too(attribute, an_input):
    with pytest.raises(RuntimeError, match="hasn't been made yet"):
        getattr(Survey(an_input), attribute)


def test_a_survey_that_was_not_made_has_no_parameters_and_no_output(an_input):
    survey = Survey(an_input)
    assert survey.has_no_fileparam and survey.has_no_output


def test_the_parameter_file_of_a_survey_that_was_not_made_cannot_be_written(an_input):
    with pytest.raises(RuntimeError, match="hasn't been made yet"):
        Survey(an_input)._write_parameter_file()


# -- the parameters -----------------------------------------------------------------

def test_preparing_a_survey_gives_it_parameters_and_an_output(an_input, directories):
    survey = prepare(Survey(an_input), directories)
    assert not survey.has_no_fileparam and not survey.has_no_output
    assert isinstance(survey.output, Output) and survey.output.survey is survey


def test_the_parameters_are_the_defaults_with_the_photometry_the_neighbours_and_what_is_given(an_input, directories):
    survey = prepare(Survey(an_input), directories, fsample=0.25, rand_seed=5)
    given = {FTTAGS.photo_categ: 'padova', FTTAGS.photo_sys: 'GAIA__DR2', FTTAGS.mag_color_names: 'G,Gbp-Grp', FTTAGS.nres: 64,
             FTTAGS.output_dir: directories[1], FTTAGS.fsample: 0.25, FTTAGS.rand_seed: 5}
    assert survey.fileparam == {**DEFAULTS_FOR_PARFILE, **given}


def test_preparing_a_survey_leaves_the_defaults_as_they_were(an_input, directories):
    before = dict(DEFAULTS_FOR_PARFILE)
    prepare(Survey(an_input), directories, fsample=0.25, rand_seed=5)
    assert DEFAULTS_FOR_PARFILE == before


def test_the_number_of_neighbours_is_that_of_the_input(directories):
    survey = prepare(Survey(make_input(directories, ngb=32)), directories)
    assert survey.fileparam[FTTAGS.nres] == 32


def test_what_the_survey_gives_cannot_be_used_to_change_it(an_input, directories):
    survey = prepare(Survey(an_input), directories)
    survey.fileparam[FTTAGS.fsample] = 99
    survey.parameters[FTTAGS.rand_seed] = 99
    assert survey.fileparam[FTTAGS.fsample] == 1 and survey.fileparam[FTTAGS.rand_seed] == DEFAULTS_FOR_PARFILE[FTTAGS.rand_seed]


def test_the_colors_of_the_diagram_can_be_given_as_a_dictionary(an_input, directories):
    cmd = {'magnitude': 'Grp', 'color_minuend': 'Gbp', 'color_subtrahend': 'G'}
    assert prepare(Survey(an_input), directories, cmd=cmd).fileparam[FTTAGS.mag_color_names] == 'Grp,Gbp-G'


@pytest.mark.parametrize("cmd, error", [
    ('G,Gbp-Nothing', "CMD magnitude names 'G,Gbp-Nothing' don't match isochrone names"),
    ({'magnitude': 'G'}, "Given cmd_magnames dict covers wrong set of keys"),
    (5, "cmd_magnames should be either a string or a dict"),
])
def test_colors_that_the_system_does_not_have_are_refused(cmd, error, an_input, directories):
    with pytest.raises(ValueError, match=error):
        prepare(Survey(an_input), directories, cmd=cmd)


def test_the_colors_are_those_of_the_first_photometric_system(an_input, directories):
    survey = Survey(an_input, photo_sys=['padova/GAIA__DR2', 'padova/LSST'])
    with pytest.raises(ValueError, match="don't match isochrone names"):
        prepare(survey, directories, cmd='g,r-i')                        # those of LSST, which is the second one


@pytest.mark.xfail(raises=pytest.fail.Exception, reason="the parameters that are given are put in the parameters whatever their "
                   "name, and the parameter file ignores those that it has no place for: a misspelled one (rand_seeed) is "
                   "silently not used, and changes the hash of the survey")
def test_a_parameter_that_does_not_exist_is_refused(an_input, directories):
    with pytest.raises((TypeError, ValueError, KeyError)):
        prepare(Survey(an_input), directories, rand_seeed=5)


# -- the generations and the other photometric systems ------------------------------

def test_each_generation_has_its_number_in_the_parameters(an_input, directories):
    survey = prepare(Survey(an_input), directories, n_gens=(3, 5))
    assert survey._extraparam == {'n_gen_0': 3, 'n_gen_1': 5} and survey.n_gens == [3, 5]


def test_the_tags_of_the_generations_are_as_wide_as_the_largest_number(an_input, directories):
    survey = prepare(Survey(an_input), directories, n_gens=range(12))
    assert list(survey._extraparam) == [f"n_gen_{i:02d}" for i in range(12)] and survey.n_gens == list(range(12))
    assert list(prepare(Survey(an_input), directories, n_gens=(0, 10))._extraparam) == ['n_gen_00', 'n_gen_01']


def test_the_generations_come_back_in_the_order_of_their_tags(an_input, directories):
    assert prepare(Survey(an_input), directories, n_gens=(7, 2, 9)).n_gens == [7, 2, 9]


@pytest.mark.parametrize("n_gens", [(), [], 0])
def test_no_generation_is_refused_with_a_message_about_the_generations(n_gens, an_input, directories):
    with pytest.raises(ValueError, match='n_gens must hold at least one generation'):
        Survey(an_input)._vanilla_survey(n_gens=n_gens, output_dir=directories[1])


def test_the_photometric_systems_after_the_first_are_in_the_parameters_by_number(an_input, directories):
    survey = prepare(Survey(an_input, photo_sys=['padova/GAIA__DR2', 'padova/LSST', 'py_custom/Mine']), directories)
    assert survey._extraparam == {'n_gen_0': 0, 'photo_categ_1': 'padova', 'photo_sys_1': 'LSST',
                                  'photo_categ_2': 'py_custom', 'photo_sys_2': 'Mine'}
    assert survey.fileparam[FTTAGS.photo_sys] == 'GAIA__DR2'            # the first one is the one that the diagram is of


def test_all_the_parameters_are_those_of_the_file_and_the_other_ones(an_input, directories):
    survey = prepare(Survey(an_input, photo_sys=['padova/GAIA__DR2', 'padova/LSST']), directories, n_gens=(0, 1))
    assert survey.parameters == {**survey.fileparam, **survey._extraparam}


def test_the_sampling_is_that_of_each_generation_times_their_number(an_input, directories):
    survey = prepare(Survey(an_input), directories, n_gens=range(3), fsample=0.01)
    assert survey.fsample == pytest.approx(0.03)


# -- the name of a survey -----------------------------------------------------------

def test_the_hash_of_a_survey_is_that_of_its_parameters_in_the_order_of_their_names_as_text(an_input, directories):
    survey = prepare(Survey(an_input, photo_sys=['padova/GAIA__DR2', 'padova/LSST']), directories, n_gens=(0, 1))
    expected = hash_iterable(str(value).encode('ascii') for _, value in sorted(survey.parameters.items()))
    assert survey.hash == expected.decode()
    assert re.fullmatch(r"[0-9a-f]{64}", survey.hash)


def test_the_name_of_a_survey_has_the_start_of_its_hash(an_input, directories):
    survey = prepare(Survey(an_input, surveyname='mine'), directories)
    assert survey.surveyname_hash == f"mine_{survey.hash[:7]}"


def test_the_hash_can_be_left_out_of_the_name_of_a_survey(directories):
    survey = prepare(Survey(make_input(directories, append_hash=False), surveyname='mine'), directories)
    assert survey.surveyname_hash == 'mine' and survey.append_hash is False


@pytest.mark.parametrize("change", ["fsample", "seed", "generations", "diagram", "second system", "output directory"])
def test_the_hash_of_a_survey_changes_with_its_parameters(change, an_input, directories):
    base = prepare(Survey(an_input), directories).hash
    other = {'fsample': lambda: prepare(Survey(an_input), directories, fsample=0.5),
             'seed': lambda: prepare(Survey(an_input), directories, rand_seed=1),
             'generations': lambda: prepare(Survey(an_input), directories, n_gens=(0, 1)),
             'diagram': lambda: prepare(Survey(an_input), directories, cmd='Grp,Gbp-G'),
             'second system': lambda: prepare(Survey(an_input, photo_sys=['padova/GAIA__DR2', 'padova/LSST']), directories),
             'output directory': lambda: prepare(Survey(an_input), directories, output_dir=directories[0])}[change]()
    assert other.hash != base


def test_surveys_with_the_same_parameters_have_the_same_hash(an_input, directories):
    assert prepare(Survey(an_input), directories).hash == prepare(Survey(an_input), directories).hash


def test_a_survey_that_is_made_again_with_other_parameters_has_another_name(an_input, directories):
    survey = prepare(Survey(an_input), directories, fsample=1)
    first = survey.surveyname_hash
    prepare(survey, directories, fsample=0.5)
    assert survey.surveyname_hash != first
    assert survey.output.output_name == f"{survey.surveyname_hash}.{an_input.name_hash}"


def test_a_survey_that_is_made_again_with_the_same_parameters_keeps_its_name(an_input, directories):
    survey = prepare(Survey(an_input), directories, fsample=1)
    first = survey.surveyname_hash
    prepare(survey, directories, fsample=1)
    assert survey.surveyname_hash == first


def test_the_metadata_of_a_survey_is_its_name_and_its_hash(an_input, directories):
    survey = prepare(Survey(an_input, surveyname='mine'), directories)
    assert survey.metadata == {'surveyname': 'mine', 'hash': survey.hash}


def test_the_output_is_named_after_the_survey_and_the_input(an_input, directories):
    survey = prepare(Survey(an_input), directories)
    assert survey.output.output_name == f"{survey.surveyname_hash}.{an_input.name_hash}"


# -- the parameter file -------------------------------------------------------------

def test_the_parameter_file_is_written_in_the_input_directory_as_the_backend_reads_it(an_input, directories):
    survey = prepare(Survey(an_input), directories)
    parfile, parameters = survey._write_parameter_file()
    assert parfile == directories[0] / f"{survey.surveyname_hash}_params"
    assert parfile.read_text() == PARFILE_TEMPLATE.substitute({FTTAGS.output_file: survey.surveyname_hash, **survey.fileparam})
    assert parameters == survey.fileparam


def test_the_parameter_file_is_written_again_without_caching(an_input, directories):
    survey = prepare(Survey(an_input), directories)
    parfile, _ = survey._write_parameter_file()
    os.utime(parfile, (1_000_000_000, 1_000_000_000))
    survey._write_parameter_file()
    assert parfile.stat().st_mtime > 1_000_001_000


def test_with_caching_a_parameter_file_that_is_right_is_left_alone_and_one_that_is_not_is_written_over(directories):
    with pytest.warns(DeprecationWarning, match="caching mode"):
        inp = make_input(directories, caching=True)
    survey = prepare(Survey(inp), directories)
    parfile, _ = survey._write_parameter_file()
    os.utime(parfile, (1_000_000_000, 1_000_000_000))
    survey._write_parameter_file()
    assert parfile.stat().st_mtime == 1_000_000_000
    parfile.write_text('something else')
    survey._write_parameter_file()
    assert parfile.read_text() == PARFILE_TEMPLATE.substitute({FTTAGS.output_file: survey.surveyname_hash, **survey.fileparam})


# -- the commands -------------------------------------------------------------------

def test_the_commands_to_run_the_backend_have_the_input_the_generation_and_the_parameter_file(an_input, directories, backend):
    survey = prepare(Survey(an_input, verbose=False), directories)
    survey._run_survey(directories[0] / 'params', n_gens=(0, 2), max_gen_workers=2)
    assert backend[0].cmds == [f"{GALAXIA} -r --hdim=6 --nfile={an_input.name_hash} --ngen={n} {directories[0] / 'params'}" for n in (0, 2)]
    assert (backend[0].max_workers, backend[0].verbose) == (2, False)


def test_the_commands_to_run_the_backend_for_kernels_in_position_space_say_three_dimensions(directories, backend):
    inp = Input(minimal_particles(), np.array([1.0, 2.0]), input_dir=directories[0], name='sim', ngb=64)
    inp.input_sorter = None
    prepare(Survey(inp), directories)._run_survey(directories[0] / 'params', n_gens=(0,), max_gen_workers=1)
    assert '--hdim=3' in backend[0].cmds[0]


def test_without_a_dimension_the_commands_have_no_such_option(an_input, directories, backend, monkeypatch):
    monkeypatch.setattr(Survey, 'hdim', property(lambda self: None))
    prepare(Survey(an_input), directories)._run_survey(directories[0] / 'params', n_gens=(0,), max_gen_workers=1)
    assert backend[0].cmds == [f"{GALAXIA} -r --nfile={an_input.name_hash} --ngen=0 {directories[0] / 'params'}"]


def test_the_commands_to_append_a_photometric_system_are_for_each_file_of_the_output(an_input, directories, backend):
    survey = prepare(Survey(an_input, photo_sys=['padova/GAIA__DR2', 'padova/LSST']), directories)
    files = [directories[1] / f"{survey.output.output_name}.{n}.ebf" for n in (0, 1)]
    for path in files:
        path.write_text('')                                                   # what the backend makes
    survey._append_survey(LSST, max_gen_workers=None)
    assert sorted(backend[0].cmds) == sorted(f"{GALAXIA} -a --pcat=padova --psys=LSST {path}" for path in files)
    assert backend[0].max_workers == 2                                         # as many as there are files by default
    survey._append_survey(LSST, max_gen_workers=1)
    assert backend[1].max_workers == 1


# -- running a survey ---------------------------------------------------------------

def test_a_survey_prepares_its_input_and_runs_the_backend_once_by_default(an_input, directories, backend):
    survey = Survey(an_input, verbose=False)
    survey._vanilla_survey(output_dir=directories[1])
    parfile = directories[0] / f"{survey.surveyname_hash}_params"
    assert [call.cmds for call in backend] == [[f"{GALAXIA} -r --hdim=6 --nfile={an_input.name_hash} --ngen=0 {parfile}"]]
    assert parfile.is_file() and an_input.pname.is_file() and an_input.kname.is_file()       # all prepared before it ran
    assert (backend[0].max_workers, backend[0].verbose) == (1, False)


def test_a_survey_runs_each_generation_it_is_asked_for(an_input, directories, backend):
    Survey(an_input)._vanilla_survey(n_gens=3, output_dir=directories[1])
    assert [cmd.split('--ngen=')[1].split(' ')[0] for cmd in backend[0].cmds] == ['0', '1', '2']
    assert backend[0].max_workers == 3


def test_a_survey_runs_the_generations_that_are_listed(an_input, directories, backend):
    Survey(an_input)._vanilla_survey(n_gens=(2, 5), output_dir=directories[1])
    assert [cmd.split('--ngen=')[1].split(' ')[0] for cmd in backend[0].cmds] == ['2', '5']


def test_n_jobs_is_the_old_name_of_the_number_of_generations(an_input, directories, backend):
    with pytest.warns(DeprecationWarning, match='n_jobs will be deprecated'):
        Survey(an_input)._vanilla_survey(n_jobs=2, output_dir=directories[1])
    assert len(backend[0].cmds) == 2


def test_the_number_of_workers_of_the_generations_is_said_to_be_not_implemented_and_given_to_execute_as_it_is(an_input, directories, backend):
    with pytest.warns(UserWarning, match='max_gen_workers is currently not implemented'):
        Survey(an_input)._vanilla_survey(n_gens=3, max_gen_workers=2, output_dir=directories[1])
    assert backend[0].max_workers == 2


def test_the_sorter_of_the_input_is_set_by_the_survey(directories, backend):
    inp = Input({**minimal_particles(), 'partitionid': np.array([1, 0])}, minimal_kernels(), input_dir=directories[0], ngb=64)
    Survey(inp)._vanilla_survey(output_dir=directories[1])
    assert inp.input_sorter.tolist() == [1, 0]
    given = np.array([0, 1])
    Survey(inp)._vanilla_survey(input_sorter=given, output_dir=directories[1])
    assert inp.input_sorter is given


def test_the_parameters_that_are_given_to_the_run_are_those_of_the_survey(an_input, directories, backend):
    survey = Survey(an_input)
    survey._vanilla_survey(fsample=0.1, rand_seed=3, output_dir=directories[1])
    assert survey.fileparam[FTTAGS.fsample] == 0.1 and survey.fileparam[FTTAGS.rand_seed] == 3


def test_the_other_photometric_systems_are_appended_to_each_file_after_the_run(an_input, directories, backend):
    survey = Survey(an_input, photo_sys=['padova/GAIA__DR2', 'padova/LSST', 'py_custom/Mine'])
    survey._vanilla_survey(n_gens=2, output_dir=directories[1])
    run, *appended = backend
    assert len(run.cmds) == 2 and ' -r ' in run.cmds[0]
    assert [sorted(call.cmds) for call in appended] == [
        sorted(f"{GALAXIA} -a --pcat={system.category} --psys={system.name} {directories[1] / (survey.output.output_name + f'.{n}.ebf')}"
               for n in (0, 1)) for system in (LSST, MINE)]


def test_without_caching_a_survey_runs_again(an_input, directories, backend):
    survey = Survey(an_input)
    survey._vanilla_survey(output_dir=directories[1])
    survey._vanilla_survey(output_dir=directories[1])
    assert len(backend) == 2


def test_with_caching_a_survey_that_was_run_is_not_run_again(directories, backend):
    with pytest.warns(DeprecationWarning, match="caching mode"):
        inp = make_input(directories, caching=True)
    survey = Survey(inp, photo_sys=['padova/GAIA__DR2', 'padova/LSST'])
    survey._vanilla_survey(output_dir=directories[1])
    first = len(backend)
    assert first == 2                                                          # the run and the append of the second system
    survey._vanilla_survey(output_dir=directories[1])
    assert len(backend) == first


def test_making_a_survey_runs_it_and_then_reads_and_post_processes_its_output(an_input, directories, backend, monkeypatch):
    seen = []
    monkeypatch.setattr(Output, 'read_galaxia_output', lambda self, *arguments: seen.append(('read', arguments)))
    monkeypatch.setattr(Output, 'post_process_output', lambda self: seen.append(('post-process',)))
    rule = lambda df: 0
    survey = Survey(an_input, verbose=True)
    output = survey.make_survey(verbose=False, partitioning_rule=rule, max_pp_workers=2, pp_auto_flush=False,
                                output_dir=directories[1], rand_seed=4)
    assert output is survey.output and survey.verbose is False
    assert len(backend) == 1 and survey.fileparam[FTTAGS.rand_seed] == 4
    assert seen == [('read', (rule, 2, False)), ('post-process',)]


def test_making_a_survey_reads_its_output_the_default_way(an_input, directories, backend, monkeypatch):
    seen = []
    monkeypatch.setattr(Output, 'read_galaxia_output', lambda self, *arguments: seen.append(arguments))
    monkeypatch.setattr(Output, 'post_process_output', lambda self: None)
    Survey(an_input).make_survey(output_dir=directories[1])
    assert seen == [(None, 1, True)]


def test_the_documentation_of_making_a_survey_has_what_it_is_made_of_filled_in():
    documentation = Survey.make_survey.__doc__
    assert "TTAGS_" not in documentation and "{DEFAULT" not in documentation
    assert f"``'{DEFAULT_CMD}'``" in documentation and "Default to 17052" in documentation
    assert "'abs_mag'" in documentation                                         # the box of the diagram, as a dictionary


if __name__ == '__main__':
    pass
