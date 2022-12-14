#!/usr/bin/env python
"""
Docstring
"""
from .constants import *
from .Input import Input
from .Survey import Survey

__all__ = []


def make_survey_from_particles(particles, rho_pos, rho_vel, photo_sys=DEFAULT_PSYS, cmd_magnames=DEFAULT_CMD, simname='sim', surveyname='survey', fsample=1, ngb=64, knorm=0.596831, **kwargs):
    input = Input(particles, rho_pos, rho_vel, simname, knorm, ngb)
    survey = Survey(input, photo_sys=photo_sys, surveyname=surveyname)
    output = survey.make_survey(cmd_magnames=cmd_magnames, fsample=fsample, **kwargs)
    return output


if __name__ == '__main__':
    pass
