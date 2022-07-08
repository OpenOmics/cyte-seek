#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from subprocess import CalledProcessError
import os, subprocess, time

# Local imports
from utils import fatal, err


def retry(times=5, exceptions=(Exception)):
    """
    Decorator to retry running a function. Retries the wrapped function 
    N times with an exponential backoff stragety. A tuple of Exceptions
    can be passed that trigger a retry attempt. When times is equal to 
    4 the back-off strategy will be {4, 16, 64, 256} seconds. Calls fatal
    if the function cannot be run within the defined number of times.
    @param times <int>: 
        The number of times to repeat the wrapped function,
        default: 5
    @param exceptions tuple(<Exception>):
        Tuple of Python Exceptions that will trigger a retry attempt,
        default: (Exception)
    @return <object func>:
        Calls fatal when func cannot be run after N times.
    """
    def decorator(func):
        def do(*args, **kwargs):
            # Begin the attempt loop and return if successful 
            attempt = 0
            delay = 1
            while attempt < times:
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    # Default expection, Exception
                    err('Function failed: {0}'.format(func.__name__))
                    err('\t@args: {0}'.format(args))
                    err('\t@kwargs: {0}'.format(kwargs))
                    err('\t@reason: {0}'.format(e))
                    # Increase backoff: 4, 16, 64, 256, 1024, 4096...
                    attempt += 1
                    delay = 4**attempt
                    err("Attempt: {0}/{1}".format(attempt, times))
                    err("Trying again in {0} seconds!\n".format(delay))
                    time.sleep(delay)
            # Could not succesfully run the given function 
            # in its alloted number of tries, exit with error 
            fatal('Fatal: failed after max retry attempts!')
        return do
    return decorator


def set_options(strict):
    """
    Changes behavior of default shell and get overrides options 
    to run bash in a strict mode. 
    @param strict <bool>:
        Overrides default shell options and runs shell in strict or 
        less permissive mode.
    @return prefix <int>:
        Returns overrides options to run bash in a strict mode
    """
    prefix = ''  # permissive shell option
    if strict: 
        # Changes behavior of default shell
        # set -e: exit immediately upon error
        # set -u: treats unset variables as an error
        # set -o pipefail: exits if a error occurs in any point of a pipeline
        prefix = 'set -euo pipefail; '

    return prefix


@retry(times=5, exceptions=(CalledProcessError))
def bash(cmd, interpreter='/bin/bash', strict=set_options(True), cwd=os.getcwd(), **kwargs):
    """
    Interface to run a process or bash command. Using subprocess.call_check()
    due to portability across most python versions. It was introduced in python 2.5
    and it is also interoperabie across all python 3 versions. 
    @param cmd <str>:
        Shell command to run
    @param interpreter <str>:
        Interpreter for command to run [default: bash]
    @pararm strict <bool>:
        Prefixes any command with 'set -euo pipefail' to ensure process fail with
        the expected exit-code  
    @params kwargs <check_call()>:
        Keyword arguments to modify subprocess.check_call() behavior
    @return exitcode <int>:
        Returns the exit code of the run command, failures return non-zero exit codes
    """
    exitcode = subprocess.check_call(strict + cmd, 
        shell=True, 
        executable=interpreter, 
        cwd=cwd, 
        **kwargs
    )
    
    return exitcode


if __name__ == '__main__':
    # Tests
    bash('ls -la /home/')
    bash('ls -la /fake/dne/path')

