#!/usr/bin/python
# -*- coding: utf-8 -*-

from time import localtime, strftime
import subprocess


def log(string):
    """Prints a log"""
    print("\n{}: {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), string))


def run(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    proc.communicate()
    return proc.returncode
