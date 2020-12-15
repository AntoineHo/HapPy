#!/usr/bin/python
# -*- coding: utf-8 -*-

# General
from time import localtime, strftime
import os, subprocess

# Stats & representation
import numpy as np
import matplotlib.pyplot as plt

def log(string: str):
    """Prints a log"""
    print("\n{}: {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), string))


def run(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    proc.communicate()
    return proc.returncode


def check_files(file):
    """Returns absolute file paths and raise exception if file does not exist"""
    if not os.path.isfile(file):
        raise Exception("ERROR: {0} is not found!".format(file))
    return file


def size_from_string(string: str) -> int:
    size_multiplier = {"K": 1000, "M": 1000000, "G": 1000000000}
    if string[-1].upper() in ["K", "M", "G"] and isfloat(string[:-1]):
        return int(float(string[:-1]) * size_multiplier[string[-1]])
    elif isfloat(string):
        if int(float(string)) < 1000 :
            print("WARNING: Small genome size, did you forget a multiplier?")
        return int(float(string))
    else:  # ERROR
        raise Exception("Size argument is not a valid number.")

def gauss(x, mu, sigma, A) :
    """Returns a gaussian curve"""
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return False

def isfloat(string) :
    try :
        float(string)
        return True
    except :
        return False
