import setuptools
import os
import re
import sys
import platform
import subprocess
import warnings

from distutils.version import LooseVersion
from setuptools.command.build_ext import build_ext


class CMakeExtension(setuptools.Extension):
    def __init__(self, name, sourcedir='', cmake_args=(), exclude_arch=False):
        setuptools.Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        self.cmake_args = cmake_args
        self.exclude_arch = exclude_arch


class CMakeBuild(build_ext):
    def run(self):
        if os.path.exists('.git'):
            subprocess.check_call(['git', 'submodule', 'update', '--init', '--recursive'])

        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.2.0':
                raise RuntimeError("CMake >= 3.2.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        extdir = os.path.join(extdir, "point_cloud_utils")
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir, '-DPYTHON_EXECUTABLE=' + sys.executable]
        cmake_args.extend(ext.cmake_args)

        cfg = 'Debug' if self.debug or os.environ.get("PCU_DEBUG") else 'Release'
        build_args = ['--config', cfg]

        if cfg == 'Debug':
            warnings.warn("Building extension %s in debug mode" % ext.name)

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if os.environ.get('CMAKE_GENERATOR') != "NMake Makefiles":
                if sys.maxsize > 2 ** 32 and not ext.exclude_arch:
                    cmake_args += ['-A', 'x64']
                    print("64 bit!")
                else:
                    print("32 bit!")
                    cmake_args += ['-A', 'Win32']
                build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''), self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake'] + cmake_args + [ext.sourcedir], cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
        print()  # Add an empty line for cleaner output


def main():
    with open("README.md", "r") as fh:
        long_description = fh.read()

    cmake_args = []
    exclude_arch = False
    if 'USE_MKL' in os.environ or '--use-mkl' in sys.argv:
        cmake_args.append('-DEIGEN_WITH_MKL=ON')
        sys.argv.remove('--use-mkl')
    if 'EXCLUDE_ARCH' in os.environ or '--exclude-arch' in sys.argv:
        exclude_arch = True
        sys.argv.remove('--exclude-arch')

    setuptools.setup(
        name="point-cloud-utils",
        version="0.24.1",
        author="Francis Williams",
        author_email="francis@fwilliams.info",
        description="A Python library for common tasks on 3D point clouds and meshes",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/fwilliams/point-cloud-utils",
        packages=setuptools.find_packages(exclude=["tests"]),
        classifiers=[
            "Programming Language :: C++",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        ],
        ext_modules=[CMakeExtension('point_cloud_utils', cmake_args=cmake_args, exclude_arch=exclude_arch)],
        cmdclass=dict(build_ext=CMakeBuild),
        zip_safe=False,
        install_requires=[
            'numpy',
            'scipy'
        ],
        test_suite="tests"
    )


if __name__ == "__main__":
    main()

