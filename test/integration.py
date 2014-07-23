from __future__ import print_function
import subprocess as sp
import shlex
import unittest
from sys import argv
from os import path

PWD = path.dirname(__file__)
R1_PLAIN = path.join(PWD, "r1.fastq")

def run_fdb(exe, *args, **kwargs):
    cmd_parts = [exe,]
    for arg, val in kwargs:
        cmd_parts.append("-" + arg)
        cmd_parts.append(val)
    cmd_parts.extend(args)
    cmd = " ".join(cmd_parts)
    print("Running:\n{}".format(cmd))
    cmd = shlex.split(cmd)
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    return proc


class TestFDBOutput(unittest.TestCase):
    def setUp(self):
        exe = argv[1]
        if path.exists(exe):
            self.exe = exe

    def test_noargs(self):
        run_fdb()

if __name__ == "__main__":
    unittest.main()
