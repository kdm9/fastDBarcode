import platform
from os import path, mkdir

build_platform = platform.system().lower()
build_dir = "-".join(["build", build_platform])
if not path.exists(build_dir):
    mkdir(build_dir)


SConscript('main.scons', variant_dir=build_dir)

Clean('src/main.c', build_dir)

