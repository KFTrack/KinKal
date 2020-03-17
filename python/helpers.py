#!/usr/bin/env python
#
# Define some helper functions and one helper class, all used to build the shared libraries.
#
# The functions are:
# def validateEnvironment():
#     - Check that some important environment variables are set.  Exit if they are not.
#
# def validateDebugLevel( ):
#     - Check that the debug level has one of the supported values
#
# def defineMergeFlags( debug_level):
#     - Define the compiler and linker flags; some flags depend on the debug leve.
#
# def defineExportedOSEnvironment():
#     - Prepare a set a of environment variables, plus their values, that will
#       be added to the scons environment.
#
# def locateSConscriptFiles(sourceRoot):
#     - Walk the directory tree rooted at sourceRoot and return a list of the paths to
#       all SConscript files under that root.
#
# def dispatchSConscriptFiles( ss, sourceRoot, build_base ):
#     - Given a list of SConscript files, ss, tell scons to build the targets specified
#       in this files and to put the outputs into the defined spot in the build heirarchy.
#
# The class build_helper has four user callable entry points:
#
#  1) make_mainlib( self, userlibs, cppf=[], pf=[], addfortran=False ):
#     Find all non-plugin .cc files in the current directory, compile
#     them and add them to a library named after the path to the current diretory.
#     The user must provide a link list (userlibs), which may be empty.
#     Optionally provide:
#        - additional compiler flags (cppf)
#        - additional flags to the parse flags feature
#     Optionally compile all fortran files and add them to the library.
#
#  2) make_plugin( self, cc, userlibs, cppf = [], pf = []):
#     Make single plugin library, specified by the named .cc file.
#     The user must provide a link list (userlibs), which may be empty.
#     The library is named after that path to the source file.
#     The optional arguments are the same as for 1)
#
#  3) make_plugins( self, userlibs, exclude_cc = [], cppf = [], pf = [] ):
#     Find all plugin files in the current directory and build them against the
#     same link list; the link list may be empty.  Exclude any plugins specified
#     on the exclude_cc list.  For each plugin this makes a shared library
#     named after the path to the source file.
#     This does not make the _dict and _map plugins.
#
#
#  4) make_dict_and_map( self, userlibs ):
#     Make the _dict and _map plugins for this directory. This looks for the files
#     classes_def.xml and classes.h.
#

from glob import glob
import os, re, string

import sys
import subprocess

# Check that some of the required environment variables have been set.
def validateEnvironment():
    if not 'PACKAGE_SOURCE' in os.environ:
        sys.exit('You have not specified PACKAGE_SOURCE for this build.\nExiting.')
    if not 'DEBUG_LEVEL' in os.environ:
        sys.exit('You have not specified DEBUG_LEVEL for this build.\nExiting.')
    if not 'BUILD_BASE' in os.environ:
        sys.exit('You have not specified BUILD_BASE for this build.\nExiting.')

# Check that the debug level is one of the known values
def validateDebugLevel():
    level = os.environ['DEBUG_LEVEL']
    allowedLevels = ['prof', 'debug' ]
    if not level in allowedLevels:
        print ('Unrecognized value for --debuglevel ', level)
        print ('   The value must be one of the allowed levels: ', str(allowedLevels) )
        raise Exception('Illegal value for --debuglevel')
        pass
    print ("Recognized debug level is: ", level)
    return level

# Define the compiler and linker options.  Scons knows how to sort out which is which.
# These are given to scons using its Evironment.MergeFlags call.
def defineMergeFlags( debug_level):
#   commonFlags = [ '-std=c++17', '-rdynamic', '-Wall', '-Wno-unused-local-typedefs', '-g', '-Werror',
   commonFlags = [ '-std=c++17', '-Wall', '-Wno-unused-local-typedefs', '-g', '-Werror',
                   '-gdwarf-2', '-Werror=return-type', '-Winit-self', '-Woverloaded-virtual', '-fPIC', '-ftrapping-math' ]
                   # Linker flags
                   #'-Wl,--no-undefined',
#                   '-Wl,--as-needed' ]

   # Hopefully the need for this will go away in later compiler versions.
   gcc_version = os.environ['GCC_VERSION']
   if gcc_version == 'v7_3_0':
       commonFlags = [ commonFlags, '-Wno-strict-overflow' ]
       pass

   extraFlags = []
   if debug_level == 'prof':
       extraFlags = [ '-O3', '-fno-omit-frame-pointer', '-DNDEBUG' ]
       pass
   if debug_level == 'debug':
       extraFlags = [ '-O0' ]
       pass
   return [ commonFlags, extraFlags ]

# Prepare some shell environment in a form to be pushed into the scons environment.
def defineExportedOSEnvironment():
    osenv = {}
    for var in [ 'LD_LIBRARY_PATH',  'GCC_FQ_DIR',  'PATH', 'PYTHONPATH',  'ROOTSYS', 'PYTHON_ROOT', 'PYTHON_DIR' ]:
        if var in os.environ.keys():
            osenv[var] = os.environ[var]
            pass
        pass
    return osenv

# Walk the directory tree to locate all SConscript files.
def locateSConscriptFiles(sourceRoot):
    ss=[]
    for root,dirs,files in os.walk(sourceRoot):
        for file in files:
            if file == 'SConscript': ss.append('%s/%s'%(root,file))
            pass
        pass
    return ss

# Tell scons to do the work found in the set of SConscript files.
def dispatchSConscriptFiles( env, ss, sourceRoot, build_base ):
    for sourceFile in ss:
        relpath = os.path.relpath(sourceFile, sourceRoot )
        tokens = relpath.split('/')
        tokens.pop()
        sep = '/'
        objPath = build_base + '/' + sep.join(tokens)
        env.SConscript ( sourceFile, variant_dir=objPath)

# An instance of this class available to the SConscript files via the scons environment.
class build_helper:
    """build_helper: class to build libraries from source"""
#   This behaves like c++ static data member and is initialized at construction time.
    def __init__(self,env):
        self.env       = env
        self.buildBase = os.environ['BUILD_BASE']
#
#   Accesor
#
    def base(self):
        return self.buildBase
#
#   Build the name of the shared library into which non-plugin compiled code will be inserted.
#   Two versions: with and without the '#/lib' path prefix.
#
    def libname(self):
        relpath = os.path.relpath('.',self.buildBase)
        tokens = relpath.split('/')
        if len(tokens) > 1:
            if tokens[len(tokens)-1] == 'src':
                tokens.pop()
                pass
            pass
        sep = '_'
        return sep.join(tokens)

    def prefixed_libname(self):
        return '#/lib/' + self.libname()
#
#   Build the name of the shared library into which plugin code will be inserted.
#   Two versions: with and without the '#/lib' path prefix.
#
    def plugin_libname(self,sourcename):
        return self.libname() + '_' + sourcename[:sourcename.find('.cc')]
    def prefixed_plugin_libname(self,sourcename):
        return '#/lib/' + self.plugin_libname(sourcename)
#
#   Build a list of plugins to be built.
#
    def plugin_cc(self):
        return ( self.env.Glob('*_module.cc',  strings=True) +
                 self.env.Glob('*_service.cc', strings=True) +
                 self.env.Glob('*_source.cc',  strings=True) )

#
#   Build a list of .cc files that are not plugings; these go into the
#   library named after the directory.
#
    def non_plugin_cc(self):
        tmp = non_plugin_cc = self.env.Glob('*.cc', strings=True)
        for cc in self.plugin_cc(): tmp.remove(cc)
        return tmp
#
#   Names need to build the _dict and _map libraries.
#
    def dict_tmp_name(self):
        relpath = os.path.relpath('.',self.buildBase)
        return '#/tmp/src/' + relpath + '/' + self.libname() + '_dict.cpp'

    def map_tmp_name(self):
        relpath = os.path.relpath('.',self.buildBase)
        return '#/tmp/src/' + relpath + '/' + self.libname() + '_map.cpp'

    def dict_libname(self):
        relpath = os.path.relpath('.',self.buildBase)
        return self.libname() + '_dict'

    def map_libname(self):
        relpath = os.path.relpath('.',self.buildBase)
        return self.libname() + '_map'

    def prefixed_dict_libname(self):
        return '#/lib/' + self.dict_libname()

    def prefixed_map_libname(self):
        return '#/lib/' + self.map_libname()
#
#   Make the main library.
#
    def make_mainlib( self, userlibs, cppf=[], pf=[], addfortran=False ):
        non_plugin_cc = self.non_plugin_cc()
        if addfortran:
            fortran = self.env.Glob('*.f')
            non_plugin_cc = [ non_plugin_cc, fortran]
            pass
        libs = []
        if non_plugin_cc:
            self.env.SharedLibrary( self.prefixed_libname(),
                                    non_plugin_cc,
                                    LIBS=[ userlibs ],
                                    CPPFLAGS=cppf,
                                    parse_flags=pf
                                )
            libs = [ self.libname() ]
            pass
        return libs
#
#   Make one plugin library ( but does not work for _dict and _map plugins )
#
    def make_plugin( self, cc, userlibs, cppf = [], pf = []):
        self.env.SharedLibrary( self.prefixed_plugin_libname(cc),
                                cc,
                                LIBS=[ userlibs, ],
                                CPPFLAGS=cppf,
                                parse_flags=pf
                            )
#
#   Make all plugin libraries, excluding _dict and _map; this works if
#   all libraries need the same link list.
#
    def make_plugins( self, userlibs, exclude_cc = [], cppf = [], pf = [] ):
        plugin_cc = self.plugin_cc()
        for cc in exclude_cc: plugin_cc.remove(cc)
        for cc in plugin_cc:
            self.env.SharedLibrary( self.prefixed_plugin_libname(cc),
                                    cc,
                                    LIBS=[ userlibs ],
                                    CPPFLAGS=cppf,
                                    parse_flags=pf
            )

#
#   Make the dictionary and map plugins.
#
    def make_dict_and_map( self, userlibs ):
        if os.path.exists('classes.h'):
            if os.path.exists('classes_def.xml'):
                self.env.DictionarySource([ self.dict_tmp_name(),
                                            self.map_tmp_name() ],
                                          [ 'classes.h', 'classes_def.xml'] )
                self.env.SharedLibrary( self.prefixed_dict_libname(),
                                        self.dict_tmp_name(),
                                        LIBS=[ userlibs ]
                                    )
                self.env.SharedLibrary( self.prefixed_map_libname(),
                                        self.map_tmp_name()
                                    )

#
#   Build a list of the source filenames for unit tests to be built and run
#
    def unittest_cc(self):
        return ( self.env.Glob('*_unit.cc',  strings=True) )


#
#   filename of executable to be made from the unit test source.
#   Just the filename without other path elements.
#
    def executable_name(self, name ):
       return name[:name.find('_unit.cc')]

#
#   Run one unit test - Fixme: still under development
#
    def run_unit_test( self, executable ):
        cmd = "../bin/"+executable   # Fixme: use proper path not hard coded ../
        print ("\n\nRunning unit test ...: ", cmd )
        p=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()
        print ( "   Status code: ", p_status)  # Fixme: modify to write to file
        print ( "   cout: ", output)
        print ( "   cerr: ", err)
        return p_status

#
#   Build all unit test
#
    def make_unit_tests( self, linkLists ):
        unitTests = self.unittest_cc()
        for u in unitTests:
            exe = self.executable_name(u)
            libs  = linkLists.get(exe,[])
            self.env.Program(
                target = "#/bin/"+exe,
                source = [ u ],
                LIBS   = libs
            )

#
#   Run all unit tests.
#
    def run_unit_tests( self ):
        unitTests = self.unittest_cc()
        for u in unitTests:
            exe = self.executable_name(u)
            status = self.run_unit_test( exe )
            print ( "   Done: ", exe, "  status: ", status)
