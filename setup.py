long_description =  """\
WRITE ME
"""

import os, re, sys, sysconfig, shutil, subprocess, site
from setuptools import setup, Command, Extension
from distutils.util import get_platform
from glob import glob


# Get version number from module
version = re.search("__version__ = '(.*)'",
                    open('python_src/__init__.py').read()).group(1)

extra_link_args = []
extra_compile_args=['-msse2', '-O3', '-funroll-loops']

F4base = Extension(
    name = 'snaffour.F4base',
    sources = ['cython_src/F4base.c', 'c_src/F4term.c', 'c_src/F4poly.c'],
    include_dirs = ['cython_src', 'c_src'], 
    extra_compile_args = extra_compile_args,
    extra_link_args = extra_link_args
)

class SnaffourClean(Command):
    """
    Clean *all* the things!
    """
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        for dir in ['build', 'dist']:
            shutil.rmtree(dir, ignore_errors=True)
        for file in glob('*.pyc') + glob('cython_src/*.c') + glob('*.egg-info'):
            if os.path.exists(file):
                os.remove(file)

class SnaffourTest(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        build_lib_dir = os.path.join(
            'build',
            'lib.{platform}-{version_info[0]}.{version_info[1]}'.format(
                platform=sysconfig.get_platform(),
                version_info=sys.version_info)
        )
        sys.path.insert(0, build_lib_dir)
        import snaffour.test

def check_call(args):
    try:
        subprocess.check_call(args)
    except subprocess.CalledProcessError:
        executable = args[0]
        command = [a for a in args if not a.startswith('-')][-1]
        raise RuntimeError(command + ' failed for ' + executable)

# For manylinux1 wheels, need to set the platform name manually
try:
    from wheel.bdist_wheel import bdist_wheel
    class SnaffourBuildWheel(bdist_wheel):
        def initialize_options(self):
            bdist_wheel.initialize_options(self)
            if sys.platform.startswith('linux'):
                plat = get_platform().replace('linux', 'manylinux1')
                plat = plat.replace('-', '_')
                self.plat_name = plat
        
except ImportError:
    SnaffourBuildWheel = None
    
        
class SnaffourRelease(Command):
    user_options = [('install', 'i', 'install the release into each Python')]
    def initialize_options(self):
        self.install = False
    def finalize_options(self):
        pass
    def run(self):
        if os.path.exists('build'):
            shutil.rmtree('build')
        if os.path.exists('dist'):
            shutil.rmtree('dist')

        pythons = os.environ.get('RELEASE_PYTHONS', sys.executable).split(',')
        for python in pythons:
            check_call([python, 'setup.py', 'build'])
            check_call([python, 'setup.py', 'test'])
            if sys.platform.startswith('linux'):
                check_call([python, 'setup.py', 'bdist_egg'])
            if self.install:
                check_call([python, 'setup.py', 'pip_install'])
            else:
                check_call([python, 'setup.py', 'bdist_wheel'])

        # Build sdist using the *first* specified Python
        check_call([pythons[0], 'setup.py', 'sdist'])

        # Double-check the Linux wheels
        if sys.platform.startswith('linux'):
            for name in os.listdir('dist'):
                if name.endswith('.whl'):
                    subprocess.check_call(['auditwheel', 'repair', os.path.join('dist', name)])

class SnaffourPipInstall(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        python = sys.executable
        check_call([python, 'setup.py', 'bdist_wheel'])
        egginfo = 'Snaffour.egg-info'
        if os.path.exists(egginfo):
            shutil.rmtree(egginfo)
        check_call([python, '-m', 'pip', 'install', '--upgrade',
                    '--no-index', '--no-cache-dir', '--find-links=dist',
                    'snaffour'])

# Check that Cython's .c files are up to date:

try:
    from Cython.Build import cythonize
    if 'clean' not in sys.argv:
        file = 'cython_src/F4base.pyx'
        if os.path.exists(file):
            cythonize([file])
except ImportError:
    pass 


setup(
    name = 'Snaffour',
    version = version,
    description = 'Finds a grevlex Groebner basis using the F4 algorithm.',
    long_description = long_description,
    url = 'http://t3m.computop.org/Snaffour',
    author = 'Marc Culler and Nathan M. Dunfield and Matthias Goerner',
    author_email = 'culler@uic.edu, nathan@dunfield.info',
    license='GPLv2+',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: C',
        'Programming Language :: Cython',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics',
        ],

    packages = ['snaffour'],
    package_dir = {'snaffour':'python_src'}, 
    ext_modules = [F4base],
    cmdclass = {'clean':SnaffourClean,
                'test':SnaffourTest,
                'release':SnaffourRelease,
                'bdist_wheel':SnaffourBuildWheel,
                'pip_install':SnaffourPipInstall,
    },
    zip_safe=False, 
)
